#!/bin/bash

#Run script from f5c directory. f5c, results, and data folders should all be located in the same directory
#The resultdata directory will have all its contents erased

#default arguments
exedir=../f5c
#where is the raw data for methylation calling located
testdir=../data
#what is the name of the folder containing the data and where the results will be stored
dataset=CHR22
#where is the directory containing the dataset directory
resultdir=../results
#path to the directory where results are stored
resultdata="${resultdir}/${dataset}"
#number of times to repeat test for each set of parameters
num_runs=3
#current individual run number initialize to -1
run_number=-1
#current group of runs
test_number=0
#number of times B has been modified
modified_B=0
#number of times K has been modified
modified_K=0
#factor to reduce parameter by each time it crashes
crash_sf=0.7

#whether the current run produced no warnings
no_warnings=0
#whether there has been a crash so far
crashed_before=0
#flag for whether the most recent run of param_test crashed or not
recently_crashed=0
#parameter that we are currently trying to optimise
NUM_PARAMETERS=7

#flag for whether to only run the divide and conquer algorithm or not
binary_only=0
#flag for whether to skip the warning step or not
skip_warnings=0
#flag for whether to run the crash test or not
run_crash_test=0
#fix max_bases
max_bases=0
#batchsize
batchsize=0

# terminate script
die() {
	echo "$1" >&2
	exit 1
}

#sets the run and test numbers to the correct values by looking in the appropriate files
set_run_and_test_number(){
    run_number=$(($(cat "${resultdata}/run.config")))
    if [ -f "${resultdata}/averages.txt" ]; then
        test_number=$(wc -l "${resultdata}/averages.txt" | awk '{ print $1 }')
    elif [ ${run_number} -eq -1  ]; then
        test_number=0
    else
        test_number=$((${run_number} / ${num_runs}))
    fi
    echo
    echo "Run number adjusted: ${run_number}, test number adjusted: ${test_number}"
    echo
}

#modify a parameter value when a crash occurs
#1st arg: path to file to be modified, 2nd arg: index of parameter to modify (3 or 4), 3rd arg: optional arguments to pass
crash_modify() {
    cd ${resultdir}
    python3 process_results.py modify "$1" "$2" "$3"
    if [ ${2} -eq 3 ]; then
        modified_K=$((${modified_K} + 1))
    elif [ ${2} -eq 4 ]; then
        modified_B=$((${modified_B} + 1))
    fi
    # debug
    # echo
    # echo "K: ${modified_K}, B: ${modified_B}"
    # echo
}

#1st arg: string containing the warning to be parsed
handle_gpu_warning() {
    cd ${resultdir}
    #get the recommended and current values for B
    recommended_B=$(echo $1 | sed -e 's/.* \(.*\)M bases.*/\1/')
    current_B=$(echo $1 | sed -e 's/.*currently \(.*\)M.*/\1/')
    # echo "Rec: ${recommended_B}"
    # echo "Cur: ${current_B}"
    #compare recommended and current value of B
    python3 process_results.py compare ${recommended_B} ${current_B} > "${dataset}/comparison.txt"
    comparison=$(head -n 1 "${dataset}/comparison.txt")
    if [ "${comparison}" = '>2.0' ]; then
        #adjust the B to be 2 * the recommended value
        new_B=$(tail -n 1 "${dataset}/comparison.txt")
        crash_modify "${dataset}/test_${dataset}.profile" 4 "-v ${new_B}"
    elif [ "${comparison}" = '>0.75' ]; then
        crash_modify "${dataset}/test_${dataset}.profile" 4 "-s ${crash_sf}"
    elif [ "${comparison}" = '<=0.75' ]; then
        #ratio is already below the recommended value so alternate between K and B
        if [ ${modified_B} -ge ${modified_K} ]; then
            #modify K
            crash_modify "${dataset}/test_${dataset}.profile" 3 "-s ${crash_sf}"
        else
            #modify B
            crash_modify "${dataset}/test_${dataset}.profile" 4 "-s ${crash_sf}"
        fi
    else
        echo "Comparison failed"
    fi
}

#handle failure of a given test if there was one. assumes that there was a crash.
#1st arg: path to the file containing output of param_test.sh where there was a crash.
handle_crash() {
    cd ${exedir}

    #look at what the warnings are saying
    gpu_warning=$(cat $1 | grep "Your GPU can accommodate")

    #modify value of either K or B depending on whether there is a warning for the B parameter.
    cd ${resultdir}

    if [ -z "${gpu_warning}" ]; then
        #no warning has been seen so far so simply modify whichever parameter has been modified less
        if [ ${modified_B} -ge ${modified_K} ]; then
            #modify K
            crash_modify "${dataset}/test_${dataset}.profile" 3 "-s ${crash_sf}"
        else
            #modify B
            crash_modify "${dataset}/test_${dataset}.profile" 4 "-s ${crash_sf}"
        fi

    else
        #gpu warning present so handle the crash
        handle_gpu_warning "${gpu_warning}"
    fi
    crashed_before=1

    #reset the run and test numbers to the appropriate values
    set_run_and_test_number

    # remove the crash file and the failed test
    cd ${exedir}
    rm ${1}
}

#all arguments have to be passed in to allow flexibility
#1st arg: dataset, 2nd arg: number of runs, 3rd arg: profile file (optional)
param_test() {
    cd ${exedir}
    #set up the command to run first
    if [ -z "${3}" ]; then
        command="scripts/param_test.sh -D ${1} -n ${2}"
    else
        command="scripts/param_test.sh -D ${1} -n ${2} -x ${3}"
        # echo "Parameter file:"
        # cat $3
        # echo
    fi

    #need to add the flag in case of override
    if [ $batchsize -ne 0 ]; then
        command="${command} -K ${batchsize} -0"
    fi

    if [ $max_bases -ne 0 ]; then
        command="${command} -B ${max_bases} -1"
    fi

    #run the command
    eval ${command} || echo "f5c crashed. Readjusting parameters..."

    #handle any crashes that may have occurred
    if [ -s "${resultdir}/${dataset}/crash_output.txt" ]; then
        handle_crash "${resultdir}/${dataset}/crash_output.txt"
        recently_crashed=1
    else
        recently_crashed=0
    fi
}

#run before starting the actual program to see if it crashes
crash_test() {
    cd ${exedir}
    recently_crashed=1
    first_run=1
    while [ ${recently_crashed} -eq 1 ]; do
        #either use default parameters or profile
        if [ ${first_run} -eq 1 ]; then
            param_test ${dataset} 1
            first_run=0
        else
            param_test ${dataset} 1 "${resultdir}/${dataset}/test_${dataset}.profile"
        fi
    done
    #remove the crash tests to start the real testing
    cd ${exedir}
    ./scripts/remove_test.sh -d ${dataset} -a
}

run_test() {
    cd ${exedir}
    recently_crashed=1
    #keep running until no crash
    while [ ${recently_crashed} -eq 1 ]; do
        param_test ${dataset} ${num_runs} ${1}
        if [ ${recently_crashed} -eq 0 ]; then
            run_number=$((${run_number} + ${num_runs}))
            test_number=$((${test_number} + 1))
        fi
    done
}

#first run: default settings and no profile
first_run() {
    cd ${exedir}
    # need to manually run because there is no profile file on the first time
    recently_crashed=1
    first_run=1
    #keep running until no crash
    while [ ${recently_crashed} -eq 1 ]; do
        #don't pass profile on first go
        if [ ${first_run} -eq 1 ]; then
            param_test ${dataset} ${num_runs}
            first_run=0
        else
            param_test ${dataset} ${num_runs} "${resultdir}/${dataset}/test_${dataset}.profile"
        fi
        if [ ${recently_crashed} -eq 0 ]; then
            run_number=$((${run_number} + ${num_runs}))
            test_number=$((${test_number} + 1))
        fi
    done

    cd ${resultdir}
    python3 process_results.py warning "${dataset}/${run_number}/warning_${run_number}.txt" "${dataset}/${run_number}/useful_${run_number}.txt" -x -n | tee evaluate_result.txt
    #if first run has no warnings then skip the warning stage by setting the flag
    evaluate_result=$(head -n 1 evaluate_result.txt)
    if [ "${evaluate_result}" = "No warnings found." ] || [ "${evaluate_result}" = "No parameters modified." ]; then
        no_warnings=1
    fi
    cd ${exedir}
}

#warning run: keep running until no warnings produced
warning_run() {
    cd ${exedir}
    run_test "${resultdata}/test_${dataset}.profile"
    cd ${resultdir}
    python3 process_results.py average "${dataset}/parameters.txt" -n $num_runs
    python3 process_results.py warning "${dataset}/${run_number}/warning_${run_number}.txt" "${dataset}/test_${dataset}.profile" | tee evaluate_result.txt
    evaluate_result=$(head -n 1 evaluate_result.txt)
    if [ "${evaluate_result}" = "No warnings found." ] || [ "${evaluate_result}" = "No parameters modified." ]; then
        no_warnings=1
    fi
    cd ${exedir}
}

#the binary search algorithm in one direction
#1st arg: which parameter, 2nd arg: initial multiplier value (either 0.5 or 2)
#3rd arg: initial test number for test 1, 4th arg: initial test number for test 2
binary_search() {
    # Modify the parameter and run
    cd ${resultdir}
    python3 process_results.py modify "${dataset}/test_${dataset}.profile" $1 -s $2
    echo "False" > "${dataset}/stop.test"
    #initialise values
    stop_test=$(cat "${dataset}/stop.test")
    streak=0
    num_worse=0
    evaluate_args="${dataset}/averages.txt -p $1 $3 $4 ${streak} ${num_worse}"
    #loop until the algorithm tells us to stop
    while [ ${stop_test} = 'False' ]; do
        run_test "${resultdata}/test_${dataset}.profile"
        cd ${resultdir}
        #evaluate the result
        python3 process_results.py average "${dataset}/parameters.txt" -n ${num_runs}
        if [ -z "${evaluate_args}" ];then
            die "Evaluate args empty: ${evaluate_args}. Terminating"
        else
            python3 process_results.py evaluate ${evaluate_args}
        fi
        stop_test=$(cat "${dataset}/stop.test")
        evaluate_args=$(cat "${dataset}/evaluate.parameters")
        #move the files
        echo "Profile: "
        cat "${dataset}/test_${dataset}.profile"
        cat "${dataset}/test_${dataset}.profile" > "${dataset}/test_${dataset}.profile.old"
        if [ -f "${dataset}/new.profile" ]; then
            cat "${dataset}/new.profile" > "${dataset}/test_${dataset}.profile"
        else
            die "Something went wrong with process_results.py evaluate function. Terminating."
        fi
        cd ${exedir}
    done
}

help_msg() {
	echo "Parameter tuning script for f5c."
	echo "Usage: f5c_dir/scripts/tune_parameters.sh args"
	echo
	echo "-t [test data dir]   Directory where test data is located. Default is ../data"
	echo "-d [dataset]   	   Name of dataset to be used. Default is CHR22"
    echo "-r [resultdir]       Directory where test results are located. Default is ../results"
    echo "-n [num_runs]        How many runs to do for each set of parameters. Default is 3"
    echo "-B [max_bases]       Allows user to specify a set value for max_bases. This will not allow the parameter to be tuned in case of crashes."
    echo "-K [batchsize]       Allows user to specify a set value for batchsize. This will not allow the parameter to be tuned in case of crashes."
    echo "-b                   Flag to indicate that only the binary search should be performed. Requires previous results to exist."
    echo "-w                   Flag to indicate that we want to skip the warning stage."
    echo "-c                   Flag to indicate whether we should test for crashes before starting the parameter tuning."
    echo "-h                   Show this help message."
}

# parse options
while getopts t:d:r:n:B:K:bwhc opt
do
	case $opt in
		t) testdir="$OPTARG";;
		d) dataset="$OPTARG";;
        r) resultdir="$OPTARG";;
        n) num_runs=$((${OPTARG}));;
        b) binary_only=1;;
        w) skip_warnings=1;;
        B) max_bases=$((${OPTARG}));;
        K) batchsize=$((${OPTARG}));;
        c) run_crash_test=1;;
		h) help_msg
		   exit 0;;
		?) printf "Usage: %s args\n" "$0"
		   exit 2;;
	esac
done
#shift $(($OPTIND - 1))

resultdata="${resultdir}/${dataset}"

if [ ${run_crash_test} -eq 1 ]; then
    echo
    echo 'Running crash test'
    crash_test
fi

#set the run number to the correct value
if [ -d ${resultdata} ] && [ -f "${resultdata}/run.config" ] && [ -f "${resultdata}/averages.txt" ]; then
    set_run_and_test_number
fi

if [ ${binary_only} -eq 0 ]; then
    #if tests have already been done, check to see if the latest one has a warning
    run_warnings=1
    if [ -f "${resultdata}/parameters.txt" ]; then
        latest_run=$((${run_number}-1))
        warnings_exist=$(cat "${resultdata}/${latest_run}/warning_${latest_run}.txt")
        if [ "${warnings_exist}" = "No warnings." ]; then
            echo 'No warnings. Skipping first and warning runs.'
            run_warnings=0
        else
            echo 'Warnings found. Continuing optimisation.'
        fi
    else
        echo
        echo 'Doing first run'
        first_run
    fi

    #user flag overrides automatic behaviour. If any crashes have occurred so far, skip the warnings
    if [ ${skip_warnings} -eq 0 ] && [ "${no_warnings}" -eq 0 ] && [ ${modified_K} -eq 0 ] && [ ${modified_B} -eq 0 ] ; then
        echo
        echo 'Doing warnings'

        #only do warning rounds if we need to
        if [ "$run_warnings" -eq 1 ]; then
            #Keep running and following warnings until no more warnings
            while [ "${no_warnings}" -eq 0 ]
            do
                warning_run
            done
        fi
    fi
else
    echo
    echo "Skipping first and warning runs"
fi

echo
echo 'Doing divide and conquer'
echo

#list of indices for parameters to optimize
optimize_indices=(0 1 2)

#if overriding, skip optimization of this parameter
if [ $batchsize -eq 0 ] && [ ${modified_K} -eq 0 ]; then
    optimize_indices+=(3)
fi

if [ $max_bases -eq 0 ] && [ ${modified_B} -eq 0 ]; then
    optimize_indices+=(4)
fi

optimize_indices+=(5 6)

#tune each parameter by applying the binary search approach
for i in ${optimize_indices[@]}; do
    #skip optimization of B or K if it was adjusted due to a crash
    if ( [ ${i} -eq 3 ] && [ ${modified_K} -ne 0 ] ) ||
    ( [ ${i} -eq 4 ] && [ ${modified_B} -ne 0 ] ); then
        continue
    fi

    #no need to do run_test on 1st round
    if [ ${i} -gt 0 ]; then
        run_test "${resultdata}/test_${dataset}.profile"
    fi
    cd ${resultdir}

    #update averages
    python3 process_results.py average "${dataset}/parameters.txt" -n $num_runs

    #save the current parameter values, test time and the line number
    tail -n 1 "${dataset}/averages.txt" > "${dataset}/averages${i}.txt"
    cat "${dataset}/test_${dataset}.profile" > "${dataset}/test_${dataset}.profile.orig${i}"
    original_test_number=$((${test_number} - 1))

    #binary search to the left
    binary_search "$i" 0.5 ${test_number} ${original_test_number}

    #save the new optimal, restore the original values
    cd ${resultdir}
    tail -n 1 "${dataset}/averages.txt" >> "${dataset}/averages${i}.txt"
    cat "${dataset}/test_${dataset}.profile.orig${i}" > "${dataset}/test_${dataset}.profile"

    #binary search to the right
    binary_search "$i" 2.0 ${test_number} ${original_test_number}

    #save the results
    cd ${resultdir}
    tail -n 1 "${dataset}/averages.txt" >> "${dataset}/averages${i}.txt"
    echo >> "${dataset}/averages${i}.txt"

    python3 process_results.py choose "${dataset}/averages${i}.txt"
    cat "${dataset}/test_${dataset}.profile" > "${dataset}/test_${dataset}.profile.old"
    cat "${dataset}/best.profile" > "${dataset}/test_${dataset}.profile"
    cat "${dataset}/best.profile" >> "${dataset}/all_best.profile"

    #need new line
    echo >> "${dataset}/all_best.profile"
done

#compare the best profiles from each set of parameter optimizations and select the best performer
cd ${resultdir}

#remove duplicate profiles from the all_best file
python3 process_results.py remove_duplicates "${dataset}/all_best.profile"

#for each of the best parameters, re-run the test again
while read p; do
    echo "$p" > "${dataset}/best.profile"
    cd ${exedir}
    run_test "${resultdir}/${dataset}/best.profile"
    cd ${resultdir}
    python3 process_results.py average "${dataset}/parameters.txt" -n ${num_runs}
    tail -n 1 "${dataset}/averages.txt" >> "${dataset}/final_averages.txt"
done < "${dataset}/all_best.profile"

python3 process_results.py choose "${dataset}/final_averages.txt"
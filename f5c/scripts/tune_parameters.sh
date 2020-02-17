#!/bin/bash

#Run script from f5c directory. f5c, results, and data folders should all be located in the same directory
#All functions assume we are starting in the f5c directory.
#The resultdata directory will have all its contents erased

# exit when command fails
set -e

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
#whether the current run produced no warnings
no_warnings=0
#parameter that we are currently trying to optimise
NUM_PARAMETERS=7
#flag for whether to only run the divide and conquer algorithm or not
binary_only=0
#flag for whether to skip the warning step or not
skip_warnings=0
#fix max_bases
max_bases=0

run_test() {
    cd ${exedir}
    if [ $max_bases -eq 0 ]; then
        scripts/param_test.sh -D ${dataset} -n ${num_runs} -x ${1}
    else
        scripts/param_test.sh -D ${dataset} -n ${num_runs} -x ${1} -B ${max_bases}
    fi
    run_number=$((${run_number} + ${num_runs}))
    test_number=$((${test_number} + 1))
}

#first run: default settings and no profile
first_run() {
    cd ${exedir}
    if [ $max_bases -eq 0 ]; then
        scripts/param_test.sh -D ${dataset} -n ${num_runs}
    else
        scripts/param_test.sh -D ${dataset} -n ${num_runs} -B ${max_bases}
    fi
    run_number=$((${run_number} + ${num_runs}))
    test_number=$((${test_number} + 1))
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
#3rd arg: initial value for test 1, 4th arg: initial value for test 2
binary_search() {
    # Modify the parameter and run
    cd ${resultdir}
    python3 process_results.py modify "${dataset}/test_${dataset}.profile" $1 -s $2
    echo "False" > "${dataset}/stop.test"
    #initialise values
    stop_test=$(cat "${dataset}/stop.test")
    streak=0
    num_worse=0
    # w_limit=2
    # s_limit=5
    # midpoint=0
    evaluate_args="${dataset}/averages.txt -p $1 $3 $4 ${streak} ${num_worse}"
    #loop until the algorithm tells us to stop
    while [ ${stop_test} = 'False' ]; do
        run_test "${resultdata}/test_${dataset}.profile"
        cd ${resultdir}
        #evaluate the result
        python3 process_results.py average "${dataset}/parameters.txt" -n ${num_runs}
        python3 process_results.py evaluate ${evaluate_args}
        stop_test=$(cat "${dataset}/stop.test")
        evaluate_args=$(cat "${dataset}/evaluate.parameters")
        #move the files
        cat "${dataset}/test_${dataset}.profile" > "${dataset}/test_${dataset}.profile.old"
        cat "${dataset}/new.profile" > "${dataset}/test_${dataset}.profile"
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
    echo "-b                   Flag to indicate that only the binary search should be performed. Requires previous results to exist."
    echo "-w                   Flag to indicate that we want to skip the warning stage."
    echo "-h                   Show this help message."
}

# parse options
while getopts t:d:r:n:B:bwh opt
do
	case $opt in
		t) testdir="$OPTARG";;
		d) dataset="$OPTARG";;
        r) resultdir="$OPTARG";;
        n) num_runs=$((${OPTARG}));;
        b) binary_only=1;;
        w) skip_warnings=1;;
        B) max_bases=$((${OPTARG}));;
		h) help_msg
		   exit 0;;
		?) printf "Usage: %s args\n" "$0"
		   exit 2;;
	esac
done
#shift $(($OPTIND - 1))

resultdata="${resultdir}/${dataset}"

#set the run number to the correct value
if [ -d ${resultdata} ] && [ -f "${resultdata}/run.config" ] && [ -f "${resultdata}/parameters.txt" ]; then
    run_number=$(($(cat "${resultdata}/run.config") - 1))
    test_number=$(wc -l "${resultdata}/parameters.txt" | awk '{ print $1 }')
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

    #user flag overrides automatic behaviour
    if [ ${skip_warnings} -eq 0 ] && [ "${no_warnings}" -eq 0 ]; then
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

#tune each parameter by applying the divide and conquer approach
# for i in $(seq 0 $((${NUM_PARAMETERS} - 1))); do
for i in 0 1 2 3 5 6; do
    #run with the original parameter value before changing
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

    python3 process_results.py choose "${dataset}/averages${i}.txt"
    cat "${dataset}/test_${dataset}.profile" > "${dataset}/test_${dataset}.profile.old"
    cat "${dataset}/best.profile" > "${dataset}/test_${dataset}.profile"
done

#test
# divide_and_conquer 0 0.5

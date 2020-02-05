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

run_test() {
    scripts/param_test.sh -D ${dataset} -n ${num_runs} -x "${resultdata}/test_${dataset}.profile"
    run_number=$((${run_number} + ${num_runs}))
    test_number=$((${test_number} + 1))
}

#first run: default settings and no profile
first_run() {
    scripts/param_test.sh -D ${dataset} -n ${num_runs}
    run_number=$((${run_number} + ${num_runs}))
    cd ${resultdir}
    python3 process_results.py warning "${dataset}/${run_number}/warning_${run_number}.txt" "${dataset}/${run_number}/useful_${run_number}.txt" -x
    cd ${exedir}
}

#warning run: keep running until no warnings produced
warning_run() {
    run_test;
    cd ${resultdir}
    python3 process_results.py average "${dataset}/parameters.txt"
    python3 process_results.py warning "${dataset}/${run_number}/warning_${run_number}.txt" "${dataset}/test_${dataset}.profile" | tee evaluate_result.txt
    evaluate_result=$(cat evaluate_result.txt)
    if [ "${evaluate_result}" = "No warnings found." ] || [ "${evaluate_result}" = "No parameters modified." ]; then
        no_warnings=1
    fi
    cd ${exedir}
}

divide_and_conquer() {
    # Modify the parameter and run
    cd ${resuldir}
    python3 process_results.py modify "${dataset}/test_${dataset}.profile" $1 -s $2
    echo "False" > "${resultdata}/stop.test"
    #initialise values
    stop_test=$(cat "${resultdata}/stop.test")
    streak=0
    num_worse=0
    w_limit=2
    s_limit=5
    midpoint=0
    #loop until the algorithm tells us to stop
    while [ ${stop_test} = 'False']; do
        run_test;
        #evaluate the result
        python3 process_results.py average "${dataset}/parameters.txt"
        python3 process_results.py evaluate "${dataset}/averages.txt" $test_number $((${test_number} - 1)) ${streak} ${num_worse}
        cd ${exedir}
    done
    #If improvement, repeat
    #If no change or worsens, try the midpoint of the old and new value
        # i. If midpoint shows no improvement or worse performance, go to step c.

        # ii. If midpoint is an improvement over both, take the midpoint of the current value and whichever one had better performance. Go to step d.
}

help_msg() {
	echo "Parameter tuning script for f5c."
	echo "Usage: f5c_dir/scripts/tune_parameters.sh args"
	echo
	echo "-t [test data dir]   Directory where test data is located."
	echo "-d [dataset]   	   Name of dataset to be used."
    echo "-r [resultdir]       Directory where test results are located."
    echo "-n [num_runs]        How many runs to do for each set of parameters"
	echo "-h                   Show this help message."
}

# parse options
while getopts t:d:r:n:h opt
do
	case $opt in
		t) testdir="$OPTARG";;
		d) dataset="$OPTARG";;
        r) resultdir="$OPTARG";;
        n) num_runs=$((${OPTARG}));;
		h) help_msg
		   exit 0;;
		?) printf "Usage: %s args\n" "$0"
		   exit 2;;
	esac
done
#shift $(($OPTIND - 1))

resultdata="${resultdir}/${dataset}"

first_run

#Keep running and following warnings until no more warnings
while [ ${no_warnings} -eq 0 ]
do
    warning_run
done

#tune each parameter by applying the divide and conquer approach
# for i in $(seq 0 $((${NUM_PARAMETERS} - 1))); do
#     divide_and_conquer "$i" 0.5
# done
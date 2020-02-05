#!/bin/sh

# exit when command fails
set -e

dataset=LIGATIONFAB42804
resultdir=../results

help_msg() {
	echo "Parameter tuning script for f5c."
	echo "Usage: f5c_dir/scripts/remove_test.sh args"
	echo
	echo "-d [dataset]         Name of dataset to be used."
    echo "-r [result dir]      Directory containing all the results."
	echo "-h                   Show this help message."
}

# parse options
while getopts d:r:h opt
do
	case $opt in
		d) dataset="$OPTARG";;
        r) resultdir="$OPTARG";;
		h) help_msg
		   exit 0;;
		?) printf "Usage: %s args\n" "$0"
		   exit 2;;
	esac
done
shift $(($OPTIND - 1))

#initialize variables
result_path="${resultdir}/${dataset}"
line_number=$(find "${result_path}" -maxdepth 1 -type d -print| wc -l | xargs -n1 expr -1 +)
test_number=$((${line_number} - 1))

#remove test from relevant files
rm -r "${result_path}/${test_number}"
sed -e "${test_number}d" "${result_path}/parameters.txt" > "${result_path}/parameters.txt"

#read and update the run number
run_number=$( cat "${result_path}/run.config" )
run_number=$((${run_number} - 1))
echo $run_number >  "${result_path}/run.config"
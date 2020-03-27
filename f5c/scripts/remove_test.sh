#!/bin/bash

#removes the most recent test

# exit when command fails
set -e

dataset=LIGATIONFAB42804
resultdir=../results
remove_all=0

help_msg() {
	echo "Parameter tuning script for f5c."
	echo "Usage: f5c_dir/scripts/remove_test.sh args"
	echo
	echo "-d [dataset]         Name of dataset to be used."
    echo "-r [result dir]      Directory containing all the results."
	echo "-a				   Remove all the data in the result folder for the dataset"
	echo "-h                   Show this help message."
}

# parse options
while getopts d:r:ha opt
do
	case $opt in
		d) dataset="$OPTARG";;
        r) resultdir="$OPTARG";;
		a) remove_all=1;;
		h) help_msg
		   exit 0;;
		?) printf "Usage: %s args\n" "$0"
		   exit 2;;
	esac
done
shift $(($OPTIND - 1))

#initialize variables
result_path="${resultdir}/${dataset}"

if [ $remove_all -eq 1 ]; then
	rm -r "${result_path}" || echo "Directory does not exist. Creating now."
	mkdir $result_path
else
	#remove test from relevant files
	if [ -f "${result_path}/parameters.txt" ]; then
		head -n -1 "${result_path}/parameters.txt" > temp.txt
		cat temp.txt > "${result_path}/parameters.txt"
	fi
fi

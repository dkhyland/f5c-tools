#!/bin/bash

#Run script from f5c directory. f5c, results, and data folders should all be located in the same directory

# defaults
exepath=./f5c
testdir=../data
dataset=LIGATIONFAB42804
bamfile=${testdir}/${dataset}/${dataset}.bam
ref=${testdir}/humangenome.fa
reads=${testdir}/${dataset}/${dataset}.fastq
resultdir=../results
numruns=1
override_B=0
override_K=0

#GPU parameters
batchsize=256
max_bases=2000000
max_lf=3.0
avg_epk=2.0
max_epk=5.0
ultra_thresh=100000
if command -v nproc > /dev/null; then
	threads=$(nproc --all)
else
	threads=8
fi

#Whether or not to index the dataset before running
index=0

# execution mode (valgrind/gdb/cpu/cuda/echo)
mode=
#profile file name
profile=
#testset_url="http://genome.cse.unsw.edu.au/tmp/f5c_ecoli_2kb_region_test.tgz"
#fallback_url="https://ndownloader.figshare.com/files/13784075?private_link=b04e3976eaed2225b848"
# empty by default
# testset_url=
# fallback_url=

# terminate script
die() {
	echo "$1" >&2
	exit 1
}

clear_fscache() {
	sync
	echo 3 | tee /proc/sys/vm/drop_caches
}

#checks the output of f5c and the return code from calling f5c and removes all the runs related to this test if there was an error
check_failed() {
	#grep for the last line of f5c if it successfully runs
	check_output=$(cat ${1} | grep "Real time:")
	error_message=$(tail -n 1 ${1})
	echo "check: ${check_output}, err: ${error_message}"
	echo
	#check for common error messages
	if [[ -z ${check_output} || ! -z $(echo ${error_message} | grep "Killed") || ! -z $(echo ${error_message} | grep "Segmentation fault") || ! -z $(echo ${error_message} | grep "meth_main::ERROR") ]]; then
		# output the raw data to a crash_output.txt file
		cat ${1} > "${resultdir}/${dataset}/crash_output.txt"

		#debug
		# echo "Removing ${2} rounds."
		# echo "Run.config and parameters.txt before:"
		# cat "${resultdir}/${dataset}/run.config"
		# cat "${resultdir}/${dataset}/parameters.txt"

		# remove all the previous successful tests for this round
		for j in $( seq 1 $((${2} - 1)) ); do
			scripts/remove_test.sh -d ${dataset}
		done

		#read and update the run number
		run_number=$( cat "${resultdir}/${dataset}/run.config" )
		run_number=$((${run_number} - ${2}))
		#don't output values lower than 0
		config_number=${run_number}
		if [ ${run_number} -lt -1 ]; then
			config_number=-1
		fi
		echo $config_number >  "${resultdir}/${dataset}/run.config"

		#debug
		# echo "Run.config and parameters.txt after:"
		# cat "${resultdir}/${dataset}/run.config"
		# cat "${resultdir}/${dataset}/parameters.txt"
		# echo
		exit $2
	fi
}

# download test set given url
#
# $1 - URL
# $2 - Fallback URL (optional)
# download_test_set() {
# 	#no data to download
# 	if [ -z $1 ] && [ -z $2 ]; then
# 		return
# 	fi
#
# 	# data set exists
# 	if [ -d ${testdir} ]; then
# 		return
# 	fi
#
# 	mkdir -p test
# 	tar_path=test/data.tgz
# 	wget -O $tar_path "$1" || wget -O $tar_path "$2" || rm -rf $tar_path ${testdir}
# 	echo "Extracting. Please wait."
# 	tar -xf $tar_path -C test || rm -rf $tar_path ${testdir}
# 	rm -f $tar_path
# }

#run in a specified mode
mode_test() {
	if [ -z ${profile} ]; then
		cmd="${exepath} call-methylation -b ${bamfile} -g ${ref} -r ${reads} --cuda-max-lf ${max_lf} --cuda-avg-epk ${avg_epk} --cuda-max-epk ${max_epk} -K ${batchsize} -B ${max_bases} -t ${threads} --ultra-thresh ${ultra_thresh}"
	else
		cmd="${exepath} call-methylation -b ${bamfile} -g ${ref} -r ${reads} -x ${profile}"
	fi

	case $1 in
		valgrind) valgrind $cmd > /dev/null;;
		gdb) gdb --args "$cmd";;
		cpu) $cmd --disable-cuda=yes;;
		cuda) $cmd --disable-cuda=no;;
		# echo) echo "$cmd -t $threads > ${testdir}/result.txt";;
		# nvprof) nvprof  -f --analysis-metrics -o profile.nvprof "$cmd" --disable-cuda=no --debug-break=5 > /dev/null;;
		# custom) shift; $cmd "$@" > ${testdir}/result.txt; execute_test;;
		*) die "Unknown mode: $1";;
	esac
}

help_msg() {
	echo "Test script for f5c."
	echo "Usage: f5c_dir/scripts/param_test.sh [-c] [-b bam file] [-g reference genome] [-r fastq/fasta read] mode"
	echo
	echo "mode                 one of:valgrind/gdb/cpu/cuda"
	echo
	echo "-d [test data dir]   Directory where test data is located. Default is ../data"
	echo "-D [dataset]   	   Name of dataset to be used. Default is LIGATIONFAB42804"
	echo "-b [bam file]        Same as f5c -b."
	echo "-r [read file]       Same as f5c -r."
	echo "-g [ref genome]      Same as f5c -g."
	echo "-x [config file]	   Same as f5c -x (specifies config file for GPU parameters)."
	echo "-l [n]			   Same as f5c --cuda-max-lf (sets boundary for CPU/GPU read allocation)."
	echo "-a [n]	   		   Same as f5c --cuda-avg-epk (average number of events per kmer)."
	echo "-m [n]	   		   Same as f5c --cuda-max-epk (max number of events per kmer)."
	echo "-K [n]               Same as f5c -K (number of reads loaded per batch)."
	echo "-B [n]               Same as f5c -B (number of bases loaded per batch)."
	echo "-t [n]               Same as f5c -t (number of threads)."
	echo "-u [n]			   Same as f5c --ultra-thresh (threshold for skipping ultra long reads)."
	echo "-n [n]			   Number of times to run f5c with the given parameters"
	echo "-0				   Override flag for K"
	echo "-1				   Override flag for B"
	# echo "-d                   Download chr22_meth_example data set and exit."
	echo "-i 		   		   Index the dataset."
	echo "-h                   Show this help message."
}

# parse options
while getopts d:D:b:g:r:x:l:a:m:K:B:t:u:n:ih01 opt
do
	case $opt in
		d) testdir="$OPTARG";;
		D) dataset=$OPTARG;;
		b) bamfile="$OPTARG";;
		g) ref="$OPTARG";;
		r) reads="$OPTARG";;
		x) profile="$OPTARG";;
		l) max_lf="$OPTARG";;
		a) avg_epk="$OPTARG";;
		m) max_epk="$OPTARG";;
		K) batchsize="$OPTARG";;
		B) max_bases="$OPTARG";;
		t) threads="$OPTARG";;
		u) ultra_thresh="$OPTARG";;
		n) numruns=$((${OPTARG}));;
		# c) testdir=test/chr22_meth_example
		#    bamfile=${testdir}/reads.sorted.bam
		#    ref=${testdir}/humangenome.fa
		#    reads=${testdir}/reads.fastq
		#    testset_url="http://genome.cse.unsw.edu.au/tmp/f5c_na12878_test.tgz"
		#    fallback_url="https://ndownloader.figshare.com/files/13784792?private_link=5dd2077f1041412a9518";;
		# d) download_test_set "http://genome.cse.unsw.edu.au/tmp/f5c_na12878_test.tgz" "https://ndownloader.figshare.com/files/13784792?private_link=5dd2077f1041412a9518"
		#    exit 0;;
		0) override_K=1;;
		1) override_B=1;;
		h) help_msg
		   exit 0;;
		i) index=1;;
		?) printf "Usage: %s [-c] [-b bam file] [-g reference genome] [-r fastq/fasta read] args\n" "$0"
		   exit 2;;
	esac
done
shift $(($OPTIND - 1))
mode=$1

# download_test_set $testset_url $fallback_url

bamfile=${testdir}/${dataset}/${dataset}.bam
ref=${testdir}/humangenome.fa
reads=${testdir}/${dataset}/${dataset}.fastq

# validate files and folders
[ -d $testdir ] || die "${testdir}: Directory does not exist"
[ -d "${testdir}/${dataset}" ] || die "${testdir}/${dataset}: Directory does not exist"
[ -d "${testdir}/${dataset}/fast5" ] || die "${testdir}/${dataset}/fast5: fast5 directory does not exist"
for file in ${bamfile} ${ref} ${reads}; do
	[ -f ${file} ] || die "${file}: File does not exist"
done

#result directory configuration
if [ ! -d "${resultdir}/${dataset}" ]; then
	echo "Results directory does not exist. Creating one now..."
	[ -d ${resultdir} ] || mkdir ${resultdir}
	mkdir "${resultdir}/${dataset}"
	echo -1 > "${resultdir}/${dataset}/run.config"
fi

#file to keep track of current run number
[ -f "${resultdir}/${dataset}/run.config" ] || echo -1 > "${resultdir}/${dataset}/run.config"

for i in $( seq 1 ${numruns} );
do
	#read and update the run number
	run_number=$( cat "${resultdir}/${dataset}/run.config" )
	run_number=$((${run_number} + 1))
	echo $run_number >  "${resultdir}/${dataset}/run.config"

	#make directory for the current test
	if [ ! -d "${resultdir}/${dataset}/${run_number}" ]; then
		mkdir "${resultdir}/${dataset}/${run_number}"
	fi
	run_folder="${resultdir}/${dataset}/${run_number}"

	#runs if user specified -i flag
	if [ ${index} -eq 1 ] && [ ${i} -eq 1 ]; then
		echo
		echo "Indexing..."
		echo
		${exepath} index -d ${testdir}/${dataset}/fast5 ${reads}
	fi

	echo
	echo "Methylation calling..."
	echo

	#variable for storing the execution result of f5c in case it crashes
	exe_result=0

	#set the command to be run
	if [ -z "$mode" ]; then
		if [ -z "${profile}" ]; then
			#create .profile file if one is not specified
			echo "${max_lf} ${avg_epk} ${max_epk} ${batchsize} ${max_bases} ${threads} ${ultra_thresh}" > "${resultdir}/${dataset}/test_${dataset}.profile"
			cat "${resultdir}/${dataset}/test_${dataset}.profile"
			command="${exepath} call-methylation -b ${bamfile} -g ${ref} -r ${reads} --cuda-max-lf ${max_lf} --cuda-avg-epk ${avg_epk} --cuda-max-epk ${max_epk} -K ${batchsize} -B ${max_bases} -t ${threads} --ultra-thresh ${ultra_thresh}"
		else
			command="${exepath} call-methylation -b ${bamfile} -g ${ref} -r ${reads} -x ${profile}"
			if [ ${override_B} -eq 1 ]; then
				command="${command} -B ${max_bases}"
			fi

			if [ ${override_K} -eq 1 ]; then
				command="${command} -K ${batchsize}"
			fi
		fi
	else
		command="mode_test "$@""
	fi

	#run f5c
	eval ${command} >"${testdir}/${dataset}/result.txt" 2> "${run_folder}/raw_${run_number}.txt"

	#check whether the execution was successful
	check_failed "${run_folder}/raw_${run_number}.txt" ${i}

	#extract data
	echo
	echo "Extracting data..."
	echo

	#useful data
	cat "${run_folder}/raw_${run_number}.txt" | grep 'init_cuda\|align_cuda\|load_balance\|memory_balance\|total entries\|total bases\|sec\|max-lf\|Real time:\|advisor::INFO\|cuda::info' | tee "${run_folder}/useful_${run_number}.txt"

	#parameters used for the test and key metrics
	printf "${run_number}: " >> "${resultdir}/${dataset}/parameters.txt"
	cat "${run_folder}/raw_${run_number}.txt" | grep 'max-lf:' | tr '\n' ' ' >> "${resultdir}/${dataset}/parameters.txt"
	printf ",Wall: " >> "${resultdir}/${dataset}/parameters.txt"
	cat "${run_folder}/raw_${run_number}.txt" | grep 'Real time:' | grep -o -P '(?<=Real time: ).*(?= sec; CPU)' >> "${resultdir}/${dataset}/parameters.txt"

	#warnings for parameter tuning
	cat "${run_folder}/useful_${run_number}.txt" | grep '::INFO' | grep -v 'Performance bounded by file I/O' > "${run_folder}/warning_${run_number}.txt" || echo "No warnings." > "${run_folder}/warning_${run_number}.txt"

	#summary of test results
	tail -n 24 "${run_folder}/raw_${run_number}.txt" > "${run_folder}/summary_${run_number}.txt" || echo "Failed to create summary" > "${run_folder}/summary_${run_number}.txt"

	echo
	echo "Data extraction complete. Useful data can be found in ${run_folder}"
	echo

	# clear_fscache
done
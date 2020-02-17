#!/bin/bash

#Run script from directory containing f5c and results folders

# terminate script
die() {
	echo "$2" >&2
	echo
	exit 1
}

XAVIER=xavier
NANOJET=nanojet
JETSON=jetson
DAV=/storage/david
F5CTOOLS=/storage/f5c-tools

source=''
target=''
all=0

#1st arg: which machine. Choices are: XAVIER, JETSON, NANOJET
#2nd arg: upload ('u') or download ('d')
#3rd arg: first word of filename to sync: param_test.sh ('param') or process_results.py ('process') or tune_parameters.sh ('tune') or 'all' to upload all files


if [[ $1 != $XAVIER ]] && [[ $1 != $NANOJET ]] && [[ $1 != $JETSON ]]; then
	die "Please use 'xavier','nanojet', or 'jetson' as the first argument"
fi

#set the server folder
if [[ $1 = $XAVIER ]]; then
	server_folder="${XAVIER}:${DAV}"
elif [[ $1 = $NANOJET ]]; then
	server_folder="${NANOJET}:${F5CTOOLS}"
elif [[ $1 = $JETSON ]]; then
	server_folder="${JETSON}:${F5CTOOLS}"
fi

if [[ $2 = 'd' ]]; then
	if [[ -z $3 ]] || [[ $3 = 'param' ]]; then
		source="${server_folder}/f5c/scripts/param_test.sh"
		target="f5c/scripts"
	elif [[ $3 = 'process' ]]; then
		source="${server_folder}/results/process_results.py"
		target="results"
	elif [[ $3 = 'tune' ]]; then
		source="${server_folder}/f5c/scripts/tune_parameters.sh"
		target="f5c/scripts"
	elif [[ $3 = 'results' ]]; then
		source="${server_folder}/results/${3}"
		target="results"
	else
		die "Failed. Please use 'param' for param_test.py or 'process' for process_results.py or 'tune' for tune_parameters.sh"
	fi
elif [[ $2 = 'u' ]]; then
	if [[ -z $3 ]] || [[ $3 = 'param' ]]; then
		source="f5c/scripts/param_test.sh"
		target="${server_folder}/f5c/scripts"
	elif [[ $3 = 'process' ]]; then
		source="results/process_results.py"
		target="${server_folder}/results"
	elif [[ $3 = 'tune' ]]; then
		source="f5c/scripts/tune_parameters.sh"
		target="${server_folder}/f5c/scripts"
	elif [[ $3 = 'all' ]]; then
		all=1
		source="f5c/scripts/param_test.sh"
		target="${server_folder}/f5c/scripts"
	else
		die "Failed. Please use 'param' for param_test.py or 'process' for process_results.py or 'tune' for tune_parameters.sh or 'all' to upload all"
	fi
	dos2unix.exe ${source}
else
	die "Failed. Please use 'u' for upload and 'd' for download."
fi

scp -r "${source}" "${target}"

if [[ $all -eq 1 ]]; then
	source="results/process_results.py"
	target="${server_folder}/results"
	dos2unix.exe ${source}
	scp "${source}" "${target}"
	source="f5c/scripts/tune_parameters.sh"
	target="${server_folder}/f5c/scripts"
	dos2unix.exe ${source}
	scp "${source}" "${target}"
fi
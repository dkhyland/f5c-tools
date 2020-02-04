#!/bin/bash

#Run script from directory containing f5c and results folders

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

source=''
target=''
all=0

#1st arg: upload ('u') or download ('d'), 2nd arg: first word of filename to sync: param_test.sh ('param') or process_results.py ('process') or tune_parameters.sh ('tune') or 'all' to upload all files
if [[ $1 = 'd' ]]; then
	if [[ -z $2 ]] || [[ $2 = 'param' ]]; then
		source="${XAVIER}:${DAV}/f5c/scripts/param_test.sh"
		target="/f5c/scripts"
	elif [[ $2 = 'process' ]]; then
		source="${XAVIER}:${DAV}/results/process_results.py"
		target="tools/process_results"
	elif [[ $2 = 'tune' ]]; then
		source="${XAVIER}:${DAV}/f5c/scripts/tune_parameters.sh"
		target="f5c/scripts"
	else
		die "Failed. Please use 'param' for param_test.py or 'process' for process_results.py or 'tune' for tune_parameters.sh"
	fi
elif [[ $1 = 'u' ]]; then
	if [[ -z $2 ]] || [[ $2 = 'param' ]]; then
		source="f5c/scripts/param_test.sh"
		target="${XAVIER}:${DAV}/f5c/scripts"
	elif [[ $2 = 'process' ]]; then
		source="tools/process_results/process_results.py"
		target="${XAVIER}:${DAV}/results"
	elif [[ $2 = 'tune' ]]; then
		source="f5c/scripts/tune_parameters.sh"
		target="${XAVIER}:${DAV}/f5c/scripts"
	elif [[ $2 = 'all' ]]; then
		all=1
		source="f5c/scripts/param_test.sh"
		target="${XAVIER}:${DAV}/f5c/scripts"
	else
		die "Failed. Please use 'param' for param_test.py or 'process' for process_results.py or 'tune' for tune_parameters.sh or 'all' to upload all"
	fi
	dos2unix.exe ${source}
else
	die "Failed. Please use 'u' for upload and 'd' for download."
fi

scp "${source}" "${target}"

if [[ $all -eq 1 ]]; then
	source="tools/process_results/process_results.py"
	target="${XAVIER}:${DAV}/results"
	dos2unix.exe ${source}
	scp "${source}" "${target}"
	source="f5c/scripts/tune_parameters.sh"
	target="${XAVIER}:${DAV}/f5c/scripts"
	dos2unix.exe ${source}
	scp "${source}" "${target}"
fi
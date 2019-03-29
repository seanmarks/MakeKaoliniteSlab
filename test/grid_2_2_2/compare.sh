#!/bin/bash

# Read list of files to check
read -d '' -ra output_files < "output_files_to_check.input"

for file in ${output_files[@]}; do
	diff_out=$( diff $file ref/$file )
	if [[ ! -z "$diff_out" ]]; then
		echo "ERROR: File $file does not match!"
		echo "$diff_out"
	fi
done

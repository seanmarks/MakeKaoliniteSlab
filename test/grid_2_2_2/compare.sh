#!/bin/bash

declare -a files=( \
	"freeze_groups_monolayer.ndx" \
	"freeze_groups_bilayer.ndx" \
	"kaolinite_slab.gro" \
	"kaolinite_reflected_slab.gro" \
	"kaolinite_slab.itp" \
)

for file in ${files[@]}; do
	diff_out=$( diff $file ref/$file )
	if [[ ! -z "$diff_out" ]]; then
		echo "ERROR: File $file does not match!"
		echo "$diff_out"
	fi
done

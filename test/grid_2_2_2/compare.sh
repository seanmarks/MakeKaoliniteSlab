#!/bin/bash

declare -a files=( \
	"freeze_groups_monolayer.ndx" \
	"freeze_groups_bilayer.ndx" \
	"kaolinite_slab.gro" \
	"kaolinite_reflected_slab.gro" \
	"kaolinite_slab.itp" \
)

for file in ${files[@]}; do
	diff $file ref/$file
done

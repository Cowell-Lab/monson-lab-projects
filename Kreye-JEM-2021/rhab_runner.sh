#!/bin/bash

# Define your function, takes one argument
my_function() {
	local arg=$1
	echo "Function called with argument: $arg"
	
	python3 /Users/s236922/code/cowell-lab/monson-lab-projects/RHAB/rhab_summary_v5.py \
		/Users/s236922/code/projects/Kreye-JEM-2021/jobs/repcalc/45497c38-dbf6-435d-b895-aea4a2c6ceec-007/45497c38-dbf6-435d-b895-aea4a2c6ceec-007/${arg}.igblast.makedb.gene.clone.mutations.airr.tsv \
		/Users/s236922/code/projects/Kreye-JEM-2021/jobs/igblast/2e15daec-f10c-4a21-8cc4-de55bbed4000-007/2e15daec-f10c-4a21-8cc4-de55bbed4000-007/${arg}.igblast.airr.tsv \
		/Users/s236922/code/projects/Kreye-JEM-2021/jobs/repcalc/45497c38-dbf6-435d-b895-aea4a2c6ceec-007/45497c38-dbf6-435d-b895-aea4a2c6ceec-007/${arg}.igblast.makedb.gene.clone.aa_properties.airr.tsv \
		/Users/s236922/code/projects/Kreye-JEM-2021/summary/${arg}.summary.tsv
}

# Put all arguments into an array
args=("$@")

# Loop through each argument and call the function
for arg in "${args[@]}"; do
	my_function "$arg"
done

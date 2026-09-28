#!/bin/zsh
conda activate ordering
python3 /Users/s236922/code/cowell-lab/monson-lab-projects/ordering_script/ordering.py \
	--data /Users/s236922/code/cowell-lab/monson-lab-projects/Kreye-JEM-2021/jobs/igblast/2e15daec-f10c-4a21-8cc4-de55bbed4000-007/2e15daec-f10c-4a21-8cc4-de55bbed4000-007/a4a761a2-a8f2-447a-a184-6abe37d9c8db.igblast.makedb.airr.tsv \
	--v_call /Users/s236922/code/cowell-lab/monson-lab-projects/ordering_script/data/genes_v_call.csv
python3 ~/code/cowell-lab/monson-lab-projects/ordering_script/ordering.py \
	--data /Users/s236922/code/cowell-lab/monson-lab-projects/Kreye-JEM-2021/jobs/igblast/2e15daec-f10c-4a21-8cc4-de55bbed4000-007/2e15daec-f10c-4a21-8cc4-de55bbed4000-007/d3bbcbdd-92b4-4a8b-8c9a-57e263d81932.igblast.makedb.airr.tsv \
	--v_call ~/code/cowell-lab/monson-lab-projects/ordering_script/data/genes_v_call.csv
conda deactivate

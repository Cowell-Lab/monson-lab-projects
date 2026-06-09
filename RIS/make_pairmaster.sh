#!/bin/bash
# python3 ../ordering_script/ordering.py --data "./<FILE_NAME>.igblast.makedb.airr.tsv" --v_call "../ordering_script/data/genes_v_call.csv"
DATA_DIR=/mnt/md0/Projects/MonsonLab/TM-HC-RIS/vdjserver/analysis_v3/stats
python3 ../ordering_script/ordering.py --data "${DATA_DIR}/UTSW33_S42_L001_R1_001.fastq.merged.unique.gene.mutations.airr.tsv" --v_call "../ordering_script/data/genes_v_call.csv"
# python3 ../ordering_script/pairmaster.py --summary "UTSW33_S42_L001_R1_001.summary.tsv" --order "<FILE_NAME>.to_order.csv" --output "RIS_pairmaster.csv"
# python3 ./txt_2_xlsx_pairmaster.py
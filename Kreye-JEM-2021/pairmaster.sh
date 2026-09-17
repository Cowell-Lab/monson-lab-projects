#!/bin/bash
conda activate ordering
python3 /Users/s236922/code/cowell-lab/monson-lab-projects/ordering_script/pairmaster.py \
    --summary /Users/s236922/code/cowell-lab/monson-lab-projects/Kreye-JEM-2021/summary/a4a761a2-a8f2-447a-a184-6abe37d9c8db.positive.summary.tsv \
    --order /Users/s236922/code/cowell-lab/monson-lab-projects/Kreye-JEM-2021/order_v2/a4a761a2-a8f2-447a-a184-6abe37d9c8db.to_order.csv \
    --output ./a4a761a2-a8f2-447a-a184-6abe37d9c8db.positive.pairmaster.csv
python3 /Users/s236922/code/cowell-lab/monson-lab-projects/ordering_script/pairmaster.py \
    --summary /Users/s236922/code/cowell-lab/monson-lab-projects/Kreye-JEM-2021/summary/d3bbcbdd-92b4-4a8b-8c9a-57e263d81932.negative.summary.tsv \
    --order /Users/s236922/code/cowell-lab/monson-lab-projects/Kreye-JEM-2021/order_v2/d3bbcbdd-92b4-4a8b-8c9a-57e263d81932.to_order.csv \
    --output ./d3bbcbdd-92b4-4a8b-8c9a-57e263d81932.negative.pairmaster.csv
conda deactivate
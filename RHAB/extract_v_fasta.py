#
# extra V gene AA into fasta
#
# Author: Scott Christley
# Date: Sept 17, 2020
#

from __future__ import print_function
import json
import argparse
import airr
import csv
import sys

# heavy chain
data_dir = '/work/data/immune/MonsonLab/RHAB/'
airr_file = 'rhab_heavy_results.airr.tsv'

reader = airr.read_rearrangement(data_dir + airr_file)
for row in reader:
    print('>' + row['sequence_id'])
    print(row['sequence_alignment_aa'])

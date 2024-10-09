#
# pull FASTA out of spreadsheet
#
# Author: Scott Christley
# Date: May 17, 2021
#

from __future__ import print_function
import json
import argparse
import airr
import csv
import sys

data_file = '/work/data/immune/MonsonLab/Sanger/5.6.2021/5-6-21 YW.csv'
fasta_file = '/work/data/immune/MonsonLab/Sanger/5.6.2021/seqs-5-6-21.fasta'

writer = open(fasta_file, 'w')
reader = csv.DictReader(open(data_file, 'r'))
#seq = 0
for row in reader:
    #print(row)
    seq_id = row['seq_id']
    #seq_id = seq_id + '.seq.' + str(seq)
    #seq += 1
    writer.write('>' + seq_id + '\n')
    writer.write(row['sequence'].strip() + '\n')


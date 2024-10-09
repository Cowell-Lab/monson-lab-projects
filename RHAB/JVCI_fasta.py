#
# combine data files into one
#
# Author: Scott Christley
# Date: Sept 16, 2020
#

from __future__ import print_function
import json
import argparse
import airr
import csv
import sys

data_file = '/work/data/immune/MonsonLab/Sanger/JCVI/JCVI.csv'
fasta_file = '/work/data/immune/MonsonLab/Sanger/JCVI/JCVI.fasta'

writer = open(fasta_file, 'w')
reader = csv.DictReader(open(data_file, 'r'))
seq = 0
for row in reader:
    #print(row)
    seq_id = row['BR CODE'] + '_' + row['PLATE'] + '_' + row['WELL']
    seq_id = seq_id + '.seq.' + str(seq)
    seq += 1
    writer.write('>' + seq_id + '\n')
    writer.write(row['SEQUENCE'].strip() + '\n')


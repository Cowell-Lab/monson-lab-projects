#
# Make FASTA file from TSVready CSV
# Monson Lab
#
# Author: Scott Christley
# Date: Sept 3, 2021
#

from __future__ import print_function
import json
import argparse
import csv
import sys

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Make FASTA file from TSVready.')
    parser.add_argument('tsv_ready', type=str, help='TSVready CSV filename')
    parser.add_argument('out_name', type=str, help='Output filename')
    args = parser.parse_args()

    seqs = {}

    if args:
        reader = csv.DictReader(open(args.tsv_ready, 'r'))
        writer = open(args.out_name, 'w')
        for row in reader:
            #seq_id = row['PCR product name'] + '_' + row['Order ID']
            seq_id = row['PCR product name']
            seq_id = seq_id.replace(' ','')
            if len(seq_id) == 0:
                continue
            if seqs.get(seq_id) is not None:
                print('Duplicate sequence ID:', row['PCR product name'])
                sys.exit(1)
            seqs[seq_id] = row
            writer.write('>' + seq_id + '\n')
            writer.write(row['Sequence'] + '\n')

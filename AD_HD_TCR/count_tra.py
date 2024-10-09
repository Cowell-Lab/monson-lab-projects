#
# Check TRA locus in sequences
#
# Author: Scott Christley
# Date: Sept 12, 2019
#

from __future__ import print_function
import json
import argparse
import os
import sys
import airr

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Count TRA locus in AIRR.')
    parser.add_argument('input_tsv', nargs='*', type=str, help='Input AIRR TSV file names')
    args = parser.parse_args()

    if args:
        for fname in args.input_tsv:
            total_cnt = 0
            np_total_cnt = 0
            tra_cnt = 0
            trb_cnt = 0
            other_cnt = 0
            np_tra_cnt = 0
            np_trb_cnt = 0
            np_other_cnt = 0
            print('Processing: ' + fname)
            reader = airr.read_rearrangement(fname)
            for row in reader:
                total_cnt += 1
                if row['productive']:
                    if row['locus'] == 'TRA':
                        tra_cnt += 1
                    elif row['locus'] == 'TRB':
                        trb_cnt += 1
                    else:
                        other_cnt += 1
                else:
                    np_total_cnt += 1
                    if row['locus'] == 'TRA':
                        np_tra_cnt += 1
                    elif row['locus'] == 'TRB':
                        np_trb_cnt += 1
                    else:
                        np_other_cnt += 1
            print(fname)
            print('productive')
            print('TRA: ' + str(tra_cnt))
            print('TRB: ' + str(trb_cnt))
            print('other: ' + str(other_cnt))
            print('non-productive')
            print('TRA: ' + str(np_tra_cnt))
            print('TRB: ' + str(np_trb_cnt))
            print('other: ' + str(np_other_cnt))
            print('  total: ' + str(total_cnt))

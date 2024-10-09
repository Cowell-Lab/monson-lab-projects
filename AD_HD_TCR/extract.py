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
            print('Processing: ' + fname)
            reader = airr.read_rearrangement(fname)
            writer = airr.derive_rearrangement('output.airr.tsv', fname)
            for row in reader:
                if row['productive']:
                    if row['locus'] == 'TRA':
                        pass
                    elif row['locus'] == 'TRB':
                        writer.write(row)
                    else:
                        writer.write(row)

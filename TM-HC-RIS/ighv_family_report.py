#
# Summary gene family usage report
#
# Author: Scott Christley
# Date: Aug 21, 2021
#

from __future__ import print_function
import argparse
import os
import sys
import csv

names = ['file', 'IGHV1', 'IGHV2', 'IGHV3', 'IGHV4', 'IGHV5', 'IGHV6', 'IGHV7']
if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Generate gene family usage summary report.')
    parser.add_argument('input_files', nargs='*', type=str, help='Gene family usage TSVs')
    args = parser.parse_args()

    if args:
        summary_cnt = {}
        summary_rel = {}
        for f in args.input_files:
            print(f)
            fields = f.split('.v_subgroup_usage.tsv')
            sample_id = fields[0]
            summary_cnt[sample_id] = {}
            summary_rel[sample_id] = {}
            summary_cnt[sample_id]['file'] = f
            summary_rel[sample_id]['file'] = f

            reader = csv.DictReader(open(f, 'r'), dialect='excel-tab')
            for row in reader:
                if row['gene'] in names:
                    summary_cnt[sample_id][row['gene']] = row['copy_count']
                    summary_rel[sample_id][row['gene']] = row['copy_freq']

        writer = csv.DictWriter(open('gene_family_count_summary.csv', 'w'), fieldnames = names)
        writer.writeheader()
        for s in summary_cnt:
            writer.writerow(summary_cnt[s])

        writer = csv.DictWriter(open('gene_family_freq_summary.csv', 'w'), fieldnames = names)
        writer.writeheader()
        for s in summary_rel:
            writer.writerow(summary_rel[s])

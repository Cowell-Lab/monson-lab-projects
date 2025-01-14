#
# just a test
#

import argparse
import airr
import csv

filename = "2070_S19_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.summary.mutations.airr.tsv"
reader = airr.read_rearrangement(filename)
for row in reader:
    print(float(row['mu_count']) / float(row['mu_freq']))
    print(len(row['germline_alignment']))
    print(len(row['cdr1']))
    print(len(row['cdr2']))
    print(len(row['fwr1']))
    print(len(row['fwr2']))
    print(len(row['fwr3']))
    cnt = 0
    for c in row['cdr1']:
        if c != '.':
            cnt += 1
    print(cnt)
    cnt = 0
    for c in row['cdr2']:
        if c != '.':
            cnt += 1
    print(cnt)
    cnt = 0
    for c in row['fwr1']:
        if c != '.':
            cnt += 1
    print(cnt)
    cnt = 0
    for c in row['fwr2']:
        if c != '.':
            cnt += 1
    print(cnt)
    cnt = 0
    for c in row['fwr3']:
        if c != '.':
            cnt += 1
    print(cnt)
    break

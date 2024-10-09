#
# AD vs HD analysis
# accumulate clone counts for clonality graphs
#
# Author: Scott Christley
# Date: Apr 8, 2020
#

from __future__ import print_function
import json
import argparse
import os
import sys
import airr
import csv

data_dir = './library5/7c83e049-bb46-4fe1-8877-93fbdecd1af1-007/cdr3_sharing_data/'
data_dir2 = './library7/c17f5818-ae03-4efb-b6de-e31b2345ba54-007/cdr3_sharing_data/'

AD_files = [
    data_dir + 'file23_cdr3_vj_aa_sharing.tsv',
    data_dir + 'file22_cdr3_vj_aa_sharing.tsv',
    data_dir + 'file1_cdr3_vj_aa_sharing.tsv',
    data_dir + 'file20_cdr3_vj_aa_sharing.tsv',
    data_dir + 'file14_cdr3_vj_aa_sharing.tsv',
    data_dir + 'file15_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample1_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample44_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample32_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample4_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample15_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample34_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample35_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample47_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample5_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample37_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample51_cdr3_vj_aa_sharing.tsv'
    ]

AD_CD8_files = [
    data_dir2 + 'sample17_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample7_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample14_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample49_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample33_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample9_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample10_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample48_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample20_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample52_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample27_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample24_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample6_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample53_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample12_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample0_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample41_cdr3_vj_aa_sharing.tsv'
    ]

MCI_files = [
    data_dir + 'file0_cdr3_vj_aa_sharing.tsv',
    data_dir + 'file3_cdr3_vj_aa_sharing.tsv',
    data_dir + 'file2_cdr3_vj_aa_sharing.tsv',
    data_dir + 'file19_cdr3_vj_aa_sharing.tsv',
    data_dir + 'file13_cdr3_vj_aa_sharing.tsv',
    data_dir + 'file16_cdr3_vj_aa_sharing.tsv'
    ]

HD_files = [
    data_dir + 'file5_cdr3_vj_aa_sharing.tsv',
    data_dir + 'file4_cdr3_vj_aa_sharing.tsv',
    data_dir + 'file7_cdr3_vj_aa_sharing.tsv',
    data_dir + 'file6_cdr3_vj_aa_sharing.tsv',
    data_dir + 'file9_cdr3_vj_aa_sharing.tsv'
    ]

HD_CD8_files = [
    data_dir2 + 'sample38_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample29_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample18_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample43_cdr3_vj_aa_sharing.tsv',
    data_dir2 + 'sample28_cdr3_vj_aa_sharing.tsv',
    ]

def accumulate_counts(clone_dict, reader):
    for row in reader:
        if clone_dict.get(row['CDR3 AA (imgt)']) is None:
            clone_dict[row['CDR3 AA (imgt)']] = {}
            clone_dict[row['CDR3 AA (imgt)']][row['LEVEL']] = int(row['TOTAL_COUNT'])
        else:
            if clone_dict[row['CDR3 AA (imgt)']].get(row['LEVEL']) is None:
                clone_dict[row['CDR3 AA (imgt)']][row['LEVEL']] = int(row['TOTAL_COUNT'])
            else:
                clone_dict[row['CDR3 AA (imgt)']][row['LEVEL']] += int(row['TOTAL_COUNT'])

def accumulate_ranges(result_dict, clone_dict):
    file_cnt = 0
    for CDR3 in clone_dict:
        for level in clone_dict[CDR3]:
            file_cnt += clone_dict[CDR3][level]
    print('file total:', file_cnt)
    result_dict['total'] = file_cnt
    result_dict['>= 1%'] = 0
    result_dict['>= 0.01% and < 1%'] = 0
    result_dict['< 0.01%'] = 0
    for CDR3 in clone_dict:
        for level in clone_dict[CDR3]:
            pct = clone_dict[CDR3][level] / file_cnt
            if pct >= 0.01:
                result_dict['>= 1%'] += clone_dict[CDR3][level]
            elif pct >= 0.0001:
                result_dict['>= 0.01% and < 1%'] += clone_dict[CDR3][level]
            else:
                result_dict['< 0.01%'] += clone_dict[CDR3][level]

AD_clones = {}
for file in AD_files:
    print('Processing:',file)
    reader = csv.DictReader(open(file, 'r'), dialect='excel-tab')
    accumulate_counts(AD_clones, reader)

outfile = open('ACS_CD4_cdr3_vj_aa_sharing.tsv', 'w')
print(len(AD_clones))
for CDR3 in AD_clones:
    for level in AD_clones[CDR3]:
        outfile.write(CDR3 + '\t' + level + '\t' + str(AD_clones[CDR3][level]) + '\n')
outfile.close()

sample_clones = []
for file in AD_files:
    file_dict = {}
    sample_clones.append(file_dict)
    file_dict['sample'] = 'ACS CD4'
    file_dict['file'] = file
    file_clones = {}
    print('Processing:',file)
    reader = csv.DictReader(open(file, 'r'), dialect='excel-tab')
    accumulate_counts(file_clones, reader)
    print(len(file_clones))
    accumulate_ranges(file_dict, file_clones)

AD_CD8_clones = {}
for file in AD_CD8_files:
    print('Processing:',file)
    reader = csv.DictReader(open(file, 'r'), dialect='excel-tab')
    accumulate_counts(AD_CD8_clones, reader)

outfile = open('ACS_CD8_cdr3_vj_aa_sharing.tsv', 'w')
print(len(AD_CD8_clones))
for CDR3 in AD_CD8_clones:
    for level in AD_CD8_clones[CDR3]:
        outfile.write(CDR3 + '\t' + level + '\t' + str(AD_CD8_clones[CDR3][level]) + '\n')
outfile.close()

for file in AD_CD8_files:
    file_dict = {}
    sample_clones.append(file_dict)
    file_dict['sample'] = 'ACS CD8'
    file_dict['file'] = file
    file_clones = {}
    print('Processing:',file)
    reader = csv.DictReader(open(file, 'r'), dialect='excel-tab')
    accumulate_counts(file_clones, reader)
    print(len(file_clones))
    accumulate_ranges(file_dict, file_clones)

#MCI_clones = {}
#for file in MCI_files:
#    print('Processing:',file)
#    reader = csv.DictReader(open(file, 'r'), dialect='excel-tab')
#    for row in reader:
#        if MCI_clones.get(row['CDR3 AA (imgt)']) is None:
#            MCI_clones[row['CDR3 AA (imgt)']] = {}
#            MCI_clones[row['CDR3 AA (imgt)']][row['LEVEL']] = int(row['TOTAL_COUNT'])
#        else:
#            if MCI_clones[row['CDR3 AA (imgt)']].get(row['LEVEL']) is None:
#                MCI_clones[row['CDR3 AA (imgt)']][row['LEVEL']] = int(row['TOTAL_COUNT'])
#            else:
#                MCI_clones[row['CDR3 AA (imgt)']][row['LEVEL']] += int(row['TOTAL_COUNT'])

#outfile = open('MCI_only_cdr3_vj_aa_sharing.tsv', 'w')
#print(len(MCI_clones))
#for CDR3 in MCI_clones:
#    for level in MCI_clones[CDR3]:
#        outfile.write(CDR3 + '\t' + level + '\t' + str(MCI_clones[CDR3][level]) + '\n')
#outfile.close()


HD_clones = {}
for file in HD_files:
    print('Processing:',file)
    reader = csv.DictReader(open(file, 'r'), dialect='excel-tab')
    accumulate_counts(HD_clones, reader)

outfile = open('HC_cdr3_vj_aa_sharing.tsv', 'w')
print(len(HD_clones))
for CDR3 in HD_clones:
    for level in HD_clones[CDR3]:
        outfile.write(CDR3 + '\t' + level + '\t' + str(HD_clones[CDR3][level]) + '\n')
outfile.close()

for file in HD_files:
    file_dict = {}
    sample_clones.append(file_dict)
    file_dict['sample'] = 'HC CD4'
    file_dict['file'] = file
    file_clones = {}
    print('Processing:',file)
    reader = csv.DictReader(open(file, 'r'), dialect='excel-tab')
    accumulate_counts(file_clones, reader)
    print(len(file_clones))
    accumulate_ranges(file_dict, file_clones)

HD_CD8_clones = {}
for file in HD_CD8_files:
    print('Processing:',file)
    reader = csv.DictReader(open(file, 'r'), dialect='excel-tab')
    accumulate_counts(HD_CD8_clones, reader)

outfile = open('HC_CD8_cdr3_vj_aa_sharing.tsv', 'w')
print(len(HD_CD8_clones))
for CDR3 in HD_CD8_clones:
    for level in HD_CD8_clones[CDR3]:
        outfile.write(CDR3 + '\t' + level + '\t' + str(HD_CD8_clones[CDR3][level]) + '\n')
outfile.close()

for file in HD_CD8_files:
    file_dict = {}
    sample_clones.append(file_dict)
    file_dict['sample'] = 'HC CD8'
    file_dict['file'] = file
    file_clones = {}
    print('Processing:',file)
    reader = csv.DictReader(open(file, 'r'), dialect='excel-tab')
    accumulate_counts(file_clones, reader)
    print(len(file_clones))
    accumulate_ranges(file_dict, file_clones)
#print(sample_clones)

outfile = open('ACS_HC_sample_clones.csv', 'w')
writer = csv.DictWriter(outfile, fieldnames=['sample', 'file', 'total', '>= 1%', '>= 0.01% and < 1%', '< 0.01%'])
writer.writeheader()
for entry in sample_clones:
    writer.writerow(entry)


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

data_dir = '/work/data/immune/MonsonLab/RHAB/'

# heavy chain
#aa_file = '589bfc04-8129-44d4-9ce6-ab2692b8d281-007/rhab_heavy.igblast.aa_props.tsv'
#mut_file = '589bfc04-8129-44d4-9ce6-ab2692b8d281-007/rhab_heavy.igblast.mut.tsv'
#airr_file = 'rhab_heavy.igblast.airr.tsv'
#output_file = 'rhab_heavy_results.airr.tsv'

# light chain
aa_file = 'f1c3d5d0-e24f-4f71-9ecb-728bc69edaf2-007/rhab_light.igblast.aa_props.tsv'
mut_file = 'f1c3d5d0-e24f-4f71-9ecb-728bc69edaf2-007/rhab_light.igblast.mut.tsv'
airr_file = 'rhab_light.igblast.airr.tsv'
output_file = 'rhab_light_results.airr.tsv'

# read the AA properties
aa_fields = ['CDR3_AA_GRAVY', 'CDR3_AA_BULK', 'CDR3_AA_ALIPHATIC', 'CDR3_AA_POLARITY', 'CDR3_AA_CHARGE', 'CDR3_AA_BASIC', 'CDR3_AA_ACIDIC', 'CDR3_AA_AROMATIC']
aa_data = {}
reader = csv.DictReader(open(data_dir + aa_file, 'r'), dialect='excel-tab')
for row in reader:
    aa_data[row['SEQUENCE_ID']] = row

# read the mutations data
mut_fields = ['num_aa_mutations', 'num_nt_mutations']
mut_data = {}
reader = csv.DictReader(open(data_dir + mut_file, 'r'), dialect='excel-tab')
for row in reader:
    mut_data[row['SEQUENCE_ID']] = row
    nt_cnt = 0
    aa_cnt = 0
    for i in range(1,105):
        f = 'MU_COUNT_' + str(i) + '_R'
        if int(row[f]) > 0:
            nt_cnt += int(row[f])
            aa_cnt += 1
    row['num_nt_mutations'] = nt_cnt
    row['num_aa_mutations'] = aa_cnt

fields = aa_fields + mut_fields
reader = airr.read_rearrangement(data_dir + airr_file)
writer = airr.derive_rearrangement(output_file, data_dir + airr_file, fields=fields)
for row in reader:
    aa_row = aa_data.get(row['sequence_id'])
    if aa_row is None:
        print(row['sequence_id'], 'is missing from AA table')
        sys.exit(1)
    for f in aa_fields:
        row[f] = aa_row[f]

    mut_row = mut_data.get(row['sequence_id'])
    if mut_row is None:
        print(row['sequence_id'], 'is missing from mutations table, ignoring...')
    else:
        for f in mut_fields:
            row[f] = mut_row[f]
    writer.write(row)

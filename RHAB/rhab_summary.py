#
# Combine gene calls, CDR3 charge and mutation info
# into a single summary TSV.
#
# This is based upon VDJServer Change-O files as we cannot
# use IgBlast's AIRR TSV for mutational analysis yet. We make
# the output look like AIRR TSV though.
#
# This is MonsonLab specific for "single cell" Sanger
# sequencing data.
#
# Author: Scott Christley
# Date: Feb 3, 2021
#

from __future__ import print_function
import json
import argparse
import airr
import csv
import sys

base_fields = ['sequence_id', 'sequence', 'productive', 'v_call', 'd_call', 'j_call', 'cdr3_aa']
aa_fields = ['CDR3_AA_GRAVY', 'CDR3_AA_BULK', 'CDR3_AA_ALIPHATIC', 'CDR3_AA_POLARITY', 'CDR3_AA_CHARGE', 'CDR3_AA_BASIC', 'CDR3_AA_ACIDIC', 'CDR3_AA_AROMATIC']
mut_fields = ['num_aa_R_mutations', 'num_nt_R_mutations', 'num_aa_S_mutations', 'num_nt_S_mutations']
fieldnames = base_fields + aa_fields +  mut_fields

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='RHAB summary report.')
    parser.add_argument('mutation_file', type=str, help='VDJServer mutation analysis TSV file')
    parser.add_argument('airr_file', type=str, help='IgBlast AIRR TSV file')
    parser.add_argument('output_file', type=str, help='Output summary TSV file')
    args = parser.parse_args()

    if args:
        # base AIRR TSV
        reader = airr.read_rearrangement(args.airr_file)
        airr_seq = {}
        final_seq = {}
        for row in reader:
            airr_seq[row['sequence_id']] = row
            entry = {}
            for f in base_fields:
                entry[f] = row[f]
            final_seq[row['sequence_id']] = entry

        # the mutational analysis files should be subset of sequences that
        # are productive and valid for the analysis
        reader = csv.DictReader(open(args.mutation_file, 'r'), dialect='excel-tab')
        for row in reader:
            airr_entry = airr_seq[row['SEQUENCE_ID']]
            entry = final_seq[row['SEQUENCE_ID']]

            for f in aa_fields:
                entry[f] = row[f]

            # sum the mutations data
            nt_cnt = 0
            aa_cnt = 0
            for i in range(1,105):
                f = 'MU_COUNT_' + str(i) + '_R'
                if int(row[f]) > 0:
                    nt_cnt += int(row[f])
                    aa_cnt += 1
            entry['num_nt_R_mutations'] = nt_cnt
            entry['num_aa_R_mutations'] = aa_cnt
            nt_cnt = 0
            aa_cnt = 0
            for i in range(1,105):
                f = 'MU_COUNT_' + str(i) + '_S'
                if int(row[f]) > 0:
                    nt_cnt += int(row[f])
                    aa_cnt += 1
            entry['num_nt_S_mutations'] = nt_cnt
            entry['num_aa_S_mutations'] = aa_cnt

        writer = csv.DictWriter(open(args.output_file, 'w'), dialect='excel-tab', fieldnames=fieldnames)
        writer.writeheader()
        for row in final_seq:
            writer.writerow(final_seq[row])

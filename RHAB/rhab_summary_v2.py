#
# Combine gene calls, CDR3 charge and mutation info
# into a single summary TSV.
#
# This is newer version for the new makedb IgBlast output
# and the new Immcantation tools
#
# This is MonsonLab specific for "single cell" Sanger
# sequencing data.
#
# Author: Scott Christley
# Date: Jun 17, 2021
#

from __future__ import print_function
import json
import argparse
import airr
import csv
import sys

base_fields = ['sequence_id', 'sequence', 'productive', 'v_call', 'd_call', 'j_call', 'cdr3_aa']
aa_fields = ['cdr3_aa_gravy','cdr3_aa_bulk','cdr3_aa_aliphatic','cdr3_aa_polarity','cdr3_aa_charge','cdr3_aa_basic','cdr3_aa_acidic','cdr3_aa_aromatic']
mut_fields = ['mu_count_v_r','mu_count_v_s']
nt_mut_fields = ['mu_aa_count_v_r', 'mu_aa_count_v_s']
fieldnames = base_fields + aa_fields +  mut_fields + nt_mut_fields

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='RHAB V2 summary report.')
    parser.add_argument('mutation_file', type=str, help='Mutation analysis TSV file')
    parser.add_argument('airr_file', type=str, help='IgBlast AIRR TSV file with AA properties')
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
            for f in aa_fields:
                entry[f] = row[f]
            final_seq[row['sequence_id']] = entry

        # the mutational analysis files should be subset of sequences that
        # are productive and valid for the analysis
        reader = airr.read_rearrangement(args.mutation_file)
        for row in reader:
            airr_entry = airr_seq[row['sequence_id']]
            entry = final_seq[row['sequence_id']]

            for f in mut_fields:
                entry[f] = row[f]

            # generate AA sums
            aa_cnt = 0
            for i in range(1,105):
                f = 'mu_count_' + str(i) + '_r'
                if int(row[f]) > 0:
                    aa_cnt += 1
            entry['mu_aa_count_v_r'] = aa_cnt
            aa_cnt = 0
            for i in range(1,105):
                f = 'mu_count_' + str(i) + '_s'
                if int(row[f]) > 0:
                    aa_cnt += 1
            entry['mu_aa_count_v_s'] = aa_cnt

        writer = csv.DictWriter(open(args.output_file, 'w'), dialect='excel-tab', fieldnames=fieldnames)
        writer.writeheader()
        for row in final_seq:
            writer.writerow(final_seq[row])

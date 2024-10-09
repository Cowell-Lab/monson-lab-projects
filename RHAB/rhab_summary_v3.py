#
# Combine gene calls, CDR3 charge and mutation info
# into a summary and all TSVs.
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
from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna, Gapped

ags_pos = [ 36, 45, 64, 65, 90, 101 ]

base_fields = ['sequence_id', 'productive', 'v_call', 'd_call', 'j_call', 'fwr1', 'cdr1', 'fwr2', 'cdr2', 'fwr3', 'cdr3', 'v_identity', 'sequence_alignment']
mut_fields = ['cdr3_aa_charge']
ags_fields = ['mu_count_36_r_aa', 'mu_count_45_r_aa', 'mu_count_64_r_aa', 'mu_count_65_r_aa', 'mu_count_90_r_aa', 'mu_count_101_r_aa']
fieldnames = ['sequence_id', 'productive', 'v_call', 'd_call', 'j_call', 'fwr1', 'cdr1', 'fwr2', 'cdr2', 'fwr3', 'cdr3', 'cdr3_aa_length', 'cdr3_aa_charge', 'mu_aa_count_v_r', 'v_identity', 'mu_count_36_r_aa', 'ags_36_codon', 'mu_count_45_r_aa', 'ags_45_codon', 'mu_count_64_r_aa', 'ags_64_codon', 'mu_count_65_r_aa', 'ags_65_codon', 'mu_count_90_r_aa', 'ags_90_codon', 'mu_count_101_r_aa', 'ags_101_codon', 'sequence_alignment']

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='RHAB V3 summary report.')
    parser.add_argument('mutation_file', type=str, help='Mutation analysis TSV file')
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
            entry['cdr3_aa_length'] = len(row['cdr3_aa'])
            final_seq[row['sequence_id']] = entry

        # the mutational analysis files should be subset of sequences that
        # are productive and valid for the analysis
        reader = airr.read_rearrangement(args.mutation_file)
        for row in reader:
            airr_entry = airr_seq[row['sequence_id']]
            entry = final_seq[row['sequence_id']]

            for f in mut_fields:
                entry[f] = row[f]

            if row['locus'] == 'IGH':
                for f in ags_fields:
                    entry[f] = row[f]
                seq_aa = Seq(row['sequence_alignment'].replace('-','N').replace('.','N')).translate()
                germ_aa = Seq(row['germline_alignment'].replace('-','N').replace('.','N')).translate()
#                seq_aa = Seq(row['sequence_alignment'].replace('-','N').replace('.','N'), generic_dna).translate()
#                germ_aa = Seq(row['germline_alignment'].replace('-','N').replace('.','N'), generic_dna).translate()

                for pos in ags_pos:
                    if float(row['mu_count_'+str(pos)+'_r_aa']) > 0:
                        entry['ags_'+str(pos)+'_codon'] = germ_aa[pos-1] + ' > ' + seq_aa[pos-1]
                        #print(pos, seq_aa[pos-1], germ_aa[pos-1])
                        #print(seq_aa)
                        #print(germ_aa)

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
            #entry['mu_aa_count_v_s'] = aa_cnt

        writer = csv.DictWriter(open(args.output_file, 'w'), dialect='excel-tab', fieldnames=fieldnames)
        writer.writeheader()
        for row in final_seq:
            writer.writerow(final_seq[row])

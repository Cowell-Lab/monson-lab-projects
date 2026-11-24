#
# Extract sequence and sequence_id into FASTA
#
# VDJServer Analysis Portal
# VDJServer Tapis applications
# https://vdjserver.org
#
# Author: Scott Christley
# Copyright (C) 2020-2023 The University of Texas Southwestern Medical Center
# Date: Sept 27, 2020
#

from __future__ import print_function
import argparse
import os
import sys
import airr

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Extract FASTA from AIRR TSV.')
    parser.add_argument('input_airr', type=str, help='Input AIRR TSV filename')
    parser.add_argument('output_fasta', type=str, help='Output FASTA filename prefix')
    args = parser.parse_args()

    if args:
        input_data = airr.read_rearrangement(args.input_airr)

        # cellranger is run on a single locus
        filename = args.output_fasta + ".fasta"

        print('Extracting FASTA for ' + args.input_airr)
        output = open(filename, 'w')
        for row in input_data:
            output.write('>' + row['adc_annotation_sequence_id'] + '|CELL_ID=' + row['adc_annotation_cell_id'] + '\n')
            output.write(row['sequence'] + '\n')

        output.close()

#
# Attach metadata fields to file
#
# Author: Scott Christley
# Date: Mar 11, 2022
#

from __future__ import print_function
import json
import argparse
import os
import sys
import airr
import csv

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Output sample table.')
    parser.add_argument('airr_json', type=str, help='Repertoire AIRR JSON file')
    parser.add_argument('input_file', type=str, help='Input file')
    parser.add_argument('output_prefix', type=str, help='Output prefix')
    args = parser.parse_args()

    if args:
        data = airr.read_airr(args.airr_json)
        reps = { obj['repertoire_id'] : obj for obj in data['Repertoire'] }

        fields = args.input_file.split('.')
        if fields[-1] == 'tsv':
            reader = csv.DictReader(open(args.input_file, 'r'), dialect='excel-tab')
            f = reader.fieldnames.copy()
            f.remove('repertoire_id')
            fieldnames = ['repertoire_id', 'library', 'sample_id', 'subject_id', 'diagnosis', 'cell type', 'template']
            fieldnames.extend(f)
            output_file = args.output_prefix + args.input_file
            output_writer = csv.DictWriter(open(output_file, 'w'), fieldnames=fieldnames, dialect='excel-tab', lineterminator='\n')
            output_writer.writeheader()
        elif fields[-1] == 'csv':
            reader = csv.DictReader(open(args.input_file, 'r'))
            f = reader.fieldnames.copy()
            f.remove('repertoire_id')
            fieldnames = ['repertoire_id', 'library', 'sample_id', 'subject_id', 'diagnosis', 'cell type', 'template']
            fieldnames.extend(f)
            output_file = args.output_prefix + args.input_file
            output_writer = csv.DictWriter(open(output_file, 'w'), fieldnames=fieldnames, lineterminator='\n')
            output_writer.writeheader()
        else:
            print('Unknown input file type, not tsv or csv.')
            sys.exit(1)

        objs = []
        for row in reader:
            rep = reps[row['repertoire_id']]
            row['library'] = rep['sample'][0]['sequencing_run_id']
            row['sample_id'] = rep['sample'][0]['sample_id']
            row['subject_id'] = rep['subject']['subject_id']
            row['diagnosis'] = rep['sample'][0]['disease_state_sample']
            row['cell type'] = rep['sample'][0]['cell_subset']['label']
            row['template'] = rep['sample'][0]['template_class']
            objs.append(row)

        # sort by subject
        decorated = [(row['subject_id'], i, row) for i, row in enumerate(objs)]
        decorated.sort()
        output_objs = [row for subject_id, i, row in decorated]

        for row in output_objs:
            output_writer.writerow(row)

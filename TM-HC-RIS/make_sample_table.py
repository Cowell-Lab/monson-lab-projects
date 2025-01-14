#
# Output sample table with custom fields
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
    parser.add_argument('output_prefix', type=str, help='Output prefix')
    parser.add_argument('--group', type=str, nargs=2, help='Repertoire group file and id')
    args = parser.parse_args()

    if args:
        rep_groups = None
        rep_group_id = None
        data = airr.load_repertoire(args.airr_json)
        reps = { obj['repertoire_id'] : obj for obj in data['Repertoire'] }

        # processing a specific group instead of all the repertoires
        if args.group:
            with open(args.group[0], 'r', encoding='utf-8') as handle:
                data = json.load(handle)
                rep_groups = { obj['repertoire_group_id'] : obj for obj in data['RepertoireGroup'] }
            rep_group_id = args.group[1]
            if rep_groups.get(rep_group_id) is None:
                print('ERROR: Cannot find repertoire group id:', rep_group_id)
                sys.exit(1)
            else:
                data = {}
                for rep in rep_groups[rep_group_id]['repertoires']:
                    data[rep['repertoire_id']] = reps[rep['repertoire_id']]
                reps = data

        fieldnames = ['repertoire_id', 'library', 'sample_id', 'subject_id', 'diagnosis', 'cell type', 'template']
        output_file = args.output_prefix + '.sample_table.csv'
        output_writer = csv.DictWriter(open(output_file, 'w'), fieldnames=fieldnames, lineterminator='\n')
        output_writer.writeheader()
        for rep_id in reps:
            rep = reps[rep_id]
            entry = { 'repertoire_id':rep['repertoire_id'] }
            entry['library'] = rep['sample'][0]['sequencing_run_id']
            entry['sample_id'] = rep['sample'][0]['sample_id']
            entry['subject_id'] = rep['subject']['subject_id']
            entry['diagnosis'] = rep['sample'][0]['disease_state_sample']
            entry['cell type'] = rep['sample'][0]['cell_subset']['label']
            entry['template'] = rep['sample'][0]['template_class']
            output_writer.writerow(entry)

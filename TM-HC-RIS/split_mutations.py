#
# Monson Lab TM/HC/RIS
#
# split mutation repertoires into groups
#

import argparse
import airr
import csv
import sys
import json

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Split mutation files.')
    parser.add_argument('airr_metadata', type=str, help='AIRR repertoire metadata file name')
    parser.add_argument('airr_group', type=str, help='AIRR repertoire group file name')
    parser.add_argument('group_name', type=str, help='AIRR repertoire group name for split')
    parser.add_argument('--frequency', dest='frequency', default=False, action='store_true', help='split frequency file')
    args = parser.parse_args()

    if args:
        data = airr.load_repertoire(args.airr_metadata)
        all_reps = { obj['repertoire_id'] : obj for obj in data['Repertoire'] }
        print('Loaded', len(all_reps), 'repertoires.')

        with open(args.airr_group, 'r', encoding='utf-8') as handle:
            data = json.load(handle)
        groups = { obj['repertoire_group_id'] : obj for obj in data['RepertoireGroup'] }
        print('Loaded', len(groups), 'repertoire groups.')

        group = groups.get(args.group_name)
        if group is None:
            print('Cannot find group:', args.group_name)
            sys.exit(1)

        reps = [ obj['repertoire_id'] for obj in group['repertoires'] ]
        print('Group', args.group_name, 'contains', len(reps), 'repertoires.')

        if args.frequency:
            reader = csv.DictReader(open('mutational_report.frequency.repertoire.csv', 'r'))
            fieldnames = reader.fieldnames.copy()
            writer = csv.DictWriter(open(args.group_name + '.mutational_report.frequency.repertoire.csv', 'w'), fieldnames=fieldnames)
        else:
            reader = csv.DictReader(open('mutational_report.repertoire.csv', 'r'))
            fieldnames = reader.fieldnames.copy()
            writer = csv.DictWriter(open(args.group_name + '.mutational_report.repertoire.csv', 'w'), fieldnames=fieldnames)

        writer.writeheader()
        for row in reader:
            if row['repertoire_id'] in reps:
                writer.writerow(row)


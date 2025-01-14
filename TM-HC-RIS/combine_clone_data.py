#
# Combine mutation frequency and abundance for clones
#
from __future__ import print_function
import json
import yaml
import argparse
import os
import sys
import csv
import airr

def load_data_file(filename):
    ext = filename.split('.')[-1]
    if ext in ('yaml', 'yml'):
        with open(filename, 'r', encoding='utf-8') as handle:
            data = yaml.load(handle, Loader=yamlordereddictloader.Loader)
    elif ext == 'json':
        with open(filename, 'r', encoding='utf-8') as handle:
            data = json.load(handle)
    else:
        if debug:
            sys.stderr.write('Unknown file type: %s. Supported file extensions are "yaml", "yml" or "json"\n' % (ext))
        raise TypeError('Unknown file type: %s. Supported file extensions are "yaml", "yml" or "json"\n' % (ext))
    return data

clones = {}
fields = ['seq_count','copy_count','seq_freq','copy_freq']

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Combine clones.')
    parser.add_argument('airr_metadata', type=str, help='AIRR repertoire metadata file name')
    parser.add_argument('airr_group', type=str, help='AIRR repertoire group file name')
    args = parser.parse_args()

    if args:
        data = airr.load_repertoire(args.airr_metadata)
        repertoires = { obj['repertoire_id'] : obj for obj in data['Repertoire'] }

        data = load_data_file(args.airr_group)
        groups = { obj['repertoire_group_id'] : obj for obj in data['RepertoireGroup'] }

        # extract mutational frequencies
        first = True
        reader = csv.DictReader(open('mutational_report.frequency.clone.csv', 'r'))
        for row in reader:
            rep_id = row['repertoire_id']
            if first:
                first = False
                fieldnames = reader.fieldnames.copy()
                fieldnames.extend(fields)
                writer = csv.DictWriter(open('combine_report.clone.csv','w'), fieldnames=fieldnames)
                writer.writeheader()

            # load clone abundances
            if clones.get(rep_id) is None:
                clones[rep_id] = {}
                r = repertoires[rep_id]
                filename = rep_id + '.ighv4.ge3ags.mutations.count.tsv'
                #filename = r['data_processing'][0]['data_processing_files'][0]
                #filename = filename.replace('.airr.tsv.gz', '.gene.clone.count.tsv')
                print('Processing:', filename)
                reader2 = csv.DictReader(open(filename, 'r'), dialect='excel-tab')
                for row2 in reader2:
                    clones[rep_id][row2['clone_id']] = row2

            clone = clones[rep_id][row['clone_id']]
            row['seq_count'] = clone['seq_count']
            row['copy_count'] = clone['copy_count']
            row['seq_freq'] = clone['seq_freq']
            row['copy_freq'] = clone['copy_freq']
            writer.writerow(row)


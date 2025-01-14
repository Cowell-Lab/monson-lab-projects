#
# Extract mutation frequency for clones for bubbleplot
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

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Extract clones.')
    parser.add_argument('airr_metadata', type=str, help='AIRR repertoire metadata file name')
    parser.add_argument('airr_group', type=str, help='AIRR repertoire group file name')
    parser.add_argument('mf_threshold', type=str, help='Mutation frequency threshold')
    args = parser.parse_args()

    if args:
        data = airr.load_repertoire(args.airr_metadata)
        repertoires = { obj['repertoire_id'] : obj for obj in data['Repertoire'] }

        data = load_data_file(args.airr_group)
        groups = { obj['repertoire_group_id'] : obj for obj in data['RepertoireGroup'] }

        # collect clones by mutation frequency
        for group in groups:
            clones[group] = []
        reader = csv.DictReader(open('mutational_report.frequency.clone.csv', 'r'))
        for row in reader:
            if float(row['mu_freq_r']) < float(args.mf_threshold):
                continue
            for group in clones:
                for rep in groups[group]['repertoires']:
                    if rep['repertoire_id'] == row['repertoire_id']:
                        clone = {}
                        clone['repertoire_id'] = row['repertoire_id']
                        clone['clone_id'] = row['clone_id']
                        clone['mu_freq_r'] = row['mu_freq_r']
                        clone['mu_freq_s'] = row['mu_freq_s']
                        clone['mu_freq_r_aa'] = row['mu_freq_r_aa']
                        clone['mu_freq_s_aa'] = row['mu_freq_s_aa']
                        clones[group].append(clone)
                        #print(group, clone['repertoire_id'], clone['clone_id'])

        # extract clonal abundance
        for group in clones:
            for rep in groups[group]['repertoires']:
                rep_id = rep['repertoire_id']
                r = repertoires[rep_id]
                filename = r['data_processing'][0]['data_processing_files'][0]
                filename = filename.replace('.airr.tsv.gz', '.gene.clone.count.tsv')
                print('Processing:', filename)
                reader = csv.DictReader(open(filename, 'r'), dialect='excel-tab')
                for row in reader:
                    for clone in clones[group]:
                        if clone['repertoire_id'] == row['repertoire_id'] and clone['clone_id'] == row['clone_id']:
                            print(group, clone['repertoire_id'], clone['clone_id'])
                            clone['seq_count'] = row['seq_count']
                            clone['copy_count'] = row['copy_count']
                            clone['seq_freq'] = row['seq_freq']
                            clone['copy_freq'] = row['copy_freq']

        # write output
        for group in clones:
            filename = group + '.mf_clones_' + args.mf_threshold + '.csv'
            writer = csv.DictWriter(open(filename,'w'), fieldnames = clones[group][0].keys())
            writer.writeheader()
            for clone in clones[group]:
                writer.writerow(clone)

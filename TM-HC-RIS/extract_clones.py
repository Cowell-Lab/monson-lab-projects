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
    parser.add_argument('processing_stage', type=str, help='Processing stage')
    parser.add_argument('top_clones', type=str, help='Top clones to extract')
    args = parser.parse_args()

    if args:
        data = airr.read_airr(args.airr_metadata)
        repertoires = { obj['repertoire_id'] : obj for obj in data['Repertoire'] }

        data = airr.read_airr(args.airr_group)
        groups = { obj['repertoire_group_id'] : obj for obj in data['RepertoireGroup'] }

        # collect clones
        for group in groups:
            if clones.get(group) is None:
                clones[group] = []
            for rep in groups[group]['repertoires']:
                rep_id = rep['repertoire_id']
                filename = rep_id + '.' + args.processing_stage + '.count.tsv'
                #r = repertoires[rep_id]
                #filename = r['data_processing'][0]['data_processing_files'][0]
                #filename = filename.replace('.airr.tsv.gz', '.gene.clone.count.tsv')
                print('Processing:', filename)
                reader = csv.DictReader(open(filename, 'r'), dialect='excel-tab')
                i = 0
                for row in reader:
                    clones[group].append(row)
                    i = i + 1
                    if i >= int(args.top_clones):
                        break

                # extract V/J genes
                filename = rep_id + '.' + args.processing_stage + '.airr.tsv'
                reader = airr.read_rearrangement(filename)
                print('Getting V/J calls:', filename)
                for row in reader:
                    for clone in clones[group]:
                        if clone['repertoire_id'] == row['repertoire_id'] and clone['clone_id'] == row['clone_id']:
                            #print(group, clone['repertoire_id'], clone['clone_id'])
                            clone['v_call'] = row['v_call']
                            clone['j_call'] = row['j_call']
                            break

        # extract mutational frequencies
        reader = csv.DictReader(open(args.processing_stage + '.clone.frequency.mutational_report.csv', 'r'))
        for row in reader:
            for group in clones:
                for clone in clones[group]:
                    if clone['repertoire_id'] == row['repertoire_id'] and clone['clone_id'] == row['clone_id']:
                        print(group, clone['repertoire_id'], clone['clone_id'])
                        clone['mu_freq_r'] = row['mu_freq_r']
                        clone['mu_freq_s'] = row['mu_freq_s']
                        clone['mu_freq_r_aa'] = row['mu_freq_r_aa']
                        clone['mu_freq_s_aa'] = row['mu_freq_s_aa']
                        break

        # write output
        for group in clones:
            filename = group + '.top_clones_' + args.top_clones + '.csv'
            writer = csv.DictWriter(open(filename,'w'), fieldnames = clones[group][0].keys())
            writer.writeheader()
            for clone in clones[group]:
                writer.writerow(clone)

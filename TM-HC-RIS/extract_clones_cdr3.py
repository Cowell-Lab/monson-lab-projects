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

cdr3s = ['CAREGGEYGDNTALDVW',
    'CARMDCNSRTCKSMDVW',
    'CGIGYSAVAAGTVDYW',
    'CARWGALLGDYYYGLDVW',
    'CGGMGLGSGSDLEDYW',
    'CAKPNAFGVVSNFDYW',
    'CAREGLYFEKEAFDIW',
    'CARDLRIPIRYNWNYGFNVLDYW',
    'CATYYYDNKHYFDYW',
    'CARGRNWEGEFDPW',
    'CASDSNNSWFYYW',
    'CARSHYIVVVVAASTCFDCW',
    'CAREVSNDFWSGYFTRWFDPW',
    'CARQPLDYYDPGRYYSAYFDYW',
    'CARHSKGGYYDSSGYNFHFDIW',
    'CARAEYYYYGMDVW',
    'CAREEYTTSSVDYW',
    'CATDLIGGSFFMYW',
    'CARDVDIAMLSLTGTTHGGLDYW',
    'CARGPGDCSSTRCYYYHYW',
    'CASIRGQHNTFDYW'
    ]

clones = {}

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Extract clones.')
    parser.add_argument('airr_metadata', type=str, help='AIRR repertoire metadata file name')
    parser.add_argument('airr_group', type=str, help='AIRR repertoire group file name')
    parser.add_argument('top_clones', type=str, help='Top clones to extract')
    args = parser.parse_args()

    if args:
        data = airr.load_repertoire(args.airr_metadata)
        repertoires = { obj['repertoire_id'] : obj for obj in data['Repertoire'] }

        data = load_data_file(args.airr_group)
        groups = { obj['repertoire_group_id'] : obj for obj in data['RepertoireGroup'] }

        # collect clones
        for group in groups:
            if clones.get(group) is None:
                clones[group] = []
            for rep in groups[group]['repertoires']:
                rep_id = rep['repertoire_id']
                r = repertoires[rep_id]
                filename = r['data_processing'][0]['data_processing_files'][0]
                filename = filename.replace('.airr.tsv.gz', '.gene.clone.count.tsv')
                print('Processing:', filename)
                reader = csv.DictReader(open(filename, 'r'), dialect='excel-tab')
                i = 0
                for row in reader:
                    clones[group].append(row)
                    i = i + 1
                    if i >= int(args.top_clones):
                        break

        new_clones = clones.copy()
        for group in clones:
            for rep in groups[group]['repertoires']:
                rep_id = rep['repertoire_id']
                r = repertoires[rep_id]
                filename = r['data_processing'][0]['data_processing_files'][0]
                filename = filename.replace('.airr.tsv.gz', '.gene.clone.airr.tsv')
                print('Processing:', filename)
                reader = airr.read_rearrangement(filename)
                i = 0
                for row in reader:
                    for clone in clones[group]:
                        if clone['repertoire_id'] == row['repertoire_id'] and clone['clone_id'] == row['clone_id']:
                            if row['junction_aa'] in cdr3s:
                                new_clones[group].remove(clone)
                                print(group, row['repertoire_id'], row['clone_id'], row['junction_aa'])
        clones = new_clones

        # extract mutational frequencies
        reader = csv.DictReader(open('mutational_report.frequency.clone.csv', 'r'))
        for row in reader:
            for group in clones:
                for clone in clones[group]:
                    if clone['repertoire_id'] == row['repertoire_id'] and clone['clone_id'] == row['clone_id']:
                        print(group, clone['repertoire_id'], clone['clone_id'])
                        clone['mu_freq_r'] = row['mu_freq_r']
                        clone['mu_freq_s'] = row['mu_freq_s']
                        clone['mu_freq_r_aa'] = row['mu_freq_r_aa']
                        clone['mu_freq_s_aa'] = row['mu_freq_s_aa']

        # write output
        for group in clones:
            filename = group + '.top_clones_' + args.top_clones + '_cdr3.csv'
            if len(clones[group]) == 0:
                continue
            writer = csv.DictWriter(open(filename,'w'), fieldnames = clones[group][0].keys())
            writer.writeheader()
            for clone in clones[group]:
                writer.writerow(clone)


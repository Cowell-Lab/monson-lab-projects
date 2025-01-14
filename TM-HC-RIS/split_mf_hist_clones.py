#
# filter AIRR TSV records based upon clones
#
from __future__ import print_function
import json
import yaml
import argparse
import os
import sys
import csv
import airr

dir = '../stats/'
clones = {}

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Extract AIRR TSV for clones.')
    parser.add_argument('airr_metadata', type=str, help='AIRR repertoire metadata file name')
    parser.add_argument('clone_file', type=str, help='clone file')
    args = parser.parse_args()

    if args:
        data = airr.load_repertoire(args.airr_metadata)
        repertoires = { obj['repertoire_id'] : obj for obj in data['Repertoire'] }

        # load clones
        reader = csv.DictReader(open(args.clone_file, 'r'))
        for row in reader:
            if clones.get(row['repertoire_id']) is None:
                clones[row['repertoire_id']] = {}
            if clones[row['repertoire_id']].get(row['clone_id']) is None:
                clones[row['repertoire_id']][row['clone_id']] = row

        for rep_id in repertoires:
            rep = repertoires[rep_id]
            filename = rep['data_processing'][0]['data_processing_files'][0]
            infile = dir + filename.replace('.airr.tsv.gz', '.gene.clone.airr.tsv')
            outname = rep_id + '.gene.clone.airr.tsv'
            print('Processing:', infile,'to:',outname)

            reader = airr.read_rearrangement(infile)
            writer = airr.derive_rearrangement(outname, infile)
            for row in reader:
                if clones.get(row['repertoire_id']):
                    if clones[row['repertoire_id']].get(row['clone_id']):
                        writer.write(row)

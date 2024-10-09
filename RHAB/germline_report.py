#
# Generate summary table for germline reconstruction
#
# Author: Scott Christley
# Date: Dec 15, 2020
#

from __future__ import print_function
import json
import argparse
import os
import sys
import csv
import airr

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Generate germline reconstruction summary.')
    parser.add_argument('repertoire_list', type=str, help='AIRR repertoire metadata file')
    args = parser.parse_args()

    if args:
        data = airr.read_repertoire(args.repertoire_list)
        repertoires = data['Repertoire']

        first = True
        names = ['sample_id']
        usage = {}
        for rep in repertoires:
            sample_id = repertoires[rep]['sample_id']
            usage[sample_id] = {'sample_id': sample_id}
            fields = repertoires[rep]['clones'].split('.clones.')
            print('Processing repertoire:', repertoires[rep]['repertoire_id'], 'sample:', sample_id, 'prefix:', fields[0])

            # non/productive counts
            filename = fields[0] + '.igblast.productive.airr.json'
            try:
                doc = json.load(open(filename, 'r'))
            except:
                print('WARNING: Could not read file', filename, 'skipping...')
                continue
            usage[sample_id]['rearrangement_count'] = doc['RECORDS']
            if first: names.append('rearrangement_count')
            usage[sample_id]['rearrangement_count_productive'] = doc['SELECTED']
            if first: names.append('rearrangement_count_productive')

            # allele counts
            filename = fields[0] + '.igblast.summary.allele.clone.airr.json'
            try:
                doc = json.load(open(filename, 'r'))
            except:
                print('WARNING: Could not read file', filename, 'skipping...')
                continue
            usage[sample_id]['allele_rearrangement_pass'] = doc['PASS']
            if first: names.append('allele_rearrangement_pass')
            usage[sample_id]['allele_rearrangement_fail'] = doc['FAIL']
            if first: names.append('allele_rearrangement_fail')
            usage[sample_id]['allele_clone_count'] = doc['CLONES']
            if first: names.append('allele_clone_count')

            # allele germline counts
            filename = fields[0] + '.igblast.allele.germ.clone.airr.json'
            try:
                doc = json.load(open(filename, 'r'))
            except:
                print('WARNING: Could not read file', filename, 'skipping...')
                continue
            usage[sample_id]['allele_germline_fail'] = doc['FAIL']
            if first: names.append('allele_germline_fail')
            usage[sample_id]['allele_germline_pass'] = doc['PASS']
            if first: names.append('allele_germline_pass')

            # not all clones get germline assignments
            # change-o does not tell us, so we have to count
            clones = {}
            filename = fields[0] + '.igblast.allele.germ.clone.airr.tsv'
            reader = airr.read_rearrangement(filename)
            for row in reader:
                if clones.get(row['clone_id']) is None:
                    clones[row['clone_id']] = row['clone_id']
            usage[sample_id]['allele_germline_count'] = len(clones)
            if first: names.append('allele_germline_count')

            # gene counts
            filename = fields[0] + '.igblast.summary.gene.clone.airr.json'
            try:
                doc = json.load(open(filename, 'r'))
            except:
                print('WARNING: Could not read file', filename, 'skipping...')
                continue
            usage[sample_id]['gene_rearrangement_pass'] = doc['PASS']
            if first: names.append('gene_rearrangement_pass')
            usage[sample_id]['gene_rearrangement_fail'] = doc['FAIL']
            if first: names.append('gene_rearrangement_fail')
            usage[sample_id]['gene_clone_count'] = doc['CLONES']
            if first: names.append('gene_clone_count')

            # gene germline counts
            filename = fields[0] + '.igblast.gene.germ.clone.airr.json'
            try:
                doc = json.load(open(filename, 'r'))
            except:
                print('WARNING: Could not read file', filename, 'skipping...')
                continue
            usage[sample_id]['gene_germline_fail'] = doc['FAIL']
            if first: names.append('gene_germline_fail')
            usage[sample_id]['gene_germline_pass'] = doc['PASS']
            if first: names.append('gene_germline_pass')

            # not all clones get germline assignments
            # change-o does not tell us, so we have to count
            clones = {}
            filename = fields[0] + '.igblast.gene.germ.clone.airr.tsv'
            reader = airr.read_rearrangement(filename)
            for row in reader:
                if clones.get(row['clone_id']) is None:
                    clones[row['clone_id']] = row['clone_id']
            usage[sample_id]['gene_germline_count'] = len(clones)
            if first: names.append('gene_germline_count')

            first = False

        writer = csv.DictWriter(open('clone_report.csv', 'w'), fieldnames=names)
        writer.writeheader()
        for rep in usage:
            writer.writerow(usage[rep])


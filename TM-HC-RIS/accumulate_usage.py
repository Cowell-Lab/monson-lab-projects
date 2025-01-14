from __future__ import print_function
import json
import argparse
import os
import sys
import csv

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Accumulate usage.')
    parser.add_argument('sample_map', type=str, help='sample map file')
    parser.add_argument('usage', type=str, help='usage name to accumulate')
    args = parser.parse_args()

    if args:
        repertoires = {}
        reader = csv.DictReader(open(args.sample_map, 'r'))
        for row in reader:
            repertoires[row['sample_id']] = row

        fieldnames = {}
        usage = {}
        for rep in repertoires:
            sample_id = repertoires[rep]['sample_id']
            usage[sample_id] = {'sample_id': sample_id}
            fields = repertoires[rep]['clones'].split('.clones.')
            filename = fields[0] + '.igblast.' + args.usage + '.tsv'
            reader = csv.DictReader(open(filename, 'r'), dialect='excel-tab')
            for row in reader:
                usage[sample_id][row['gene']] = row['copy_count']
                fieldnames[row['gene']] = 1

        for rep in usage:
            for f in fieldnames:
                if usage[rep].get(f) is None:
                    usage[rep][f] = 0

        names = ['sample_id']
        for f in fieldnames:
            names.append(f)

        writer = csv.DictWriter(open(args.usage + '.csv', 'w'), fieldnames=names)
        writer.writeheader()
        for rep in usage:
            writer.writerow(usage[rep])


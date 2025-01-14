#
# Search for similar CDR3s
#
from __future__ import print_function
import json
import yaml
import argparse
import os
import sys
import csv
import airr
import Levenshtein

seqs = {}
rhabs = {}
totals = {}

def get_duplicate_count(row):
    if row['duplicate_count'] is None:
        return 1
    else:
        x = int(row['duplicate_count'])
        if x == 0:
            return 1
        else:
            return x

def get_gene(gene_call):
    m = gene_call.split(',')
    if len(m) > 1:
        gene_call = m[0]
    m = gene_call.split('*')
    return m[0]

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Search CDR3s.')
    parser.add_argument('airr_metadata', type=str, help='AIRR repertoire metadata file name')
    parser.add_argument('airr_group', type=str, help='AIRR repertoire group file name')
    parser.add_argument('processing_stage', type=str, help='Processing stage')
    parser.add_argument('input_cdr3s', type=str, help='Input AIRR TSV with CDR3s')
    parser.add_argument('group_name', type=str, help='AIRR repertoire group name')
    parser.add_argument('outfile', type=str, help='Output prefix')
    args = parser.parse_args()

    if args:
        data = airr.read_airr(args.airr_metadata)
        repertoires = { obj['repertoire_id'] : obj for obj in data['Repertoire'] }

        data = airr.read_airr(args.airr_group)
        groups = { obj['repertoire_group_id'] : obj for obj in data['RepertoireGroup'] }

        cdr3s = []
        reader = airr.read_rearrangement(args.input_cdr3s)
        for row in reader:
            cdr3s.append(row)
        print(len(cdr3s), 'CDR3s in input file')

        fieldnames = ['repertoire_id', 'distance', 'rhab_sequence_id', 'rhab_junction', 'match_junction', 'sequence_id', 'clone_id', 'v_call', 'j_call', 'duplicate_count']
        writer = csv.DictWriter(open(args.group_name + '.' + args.outfile + '_detail.csv','w'), fieldnames = fieldnames)
        writer.writeheader()

        # search
        score_cutoff = 2
        for rep in groups[args.group_name]['repertoires']:
            rep_id = rep['repertoire_id']
            filename = rep_id + '.' + args.processing_stage + '.airr.tsv'
            reader = airr.read_rearrangement(filename)
            dname = rep_id + '.' + args.processing_stage + '.' + args.outfile + '.airr.tsv'
            dwriter = airr.derive_rearrangement(dname, filename)
            print('Processing:', filename, 'for group', args.group_name)
            for row in reader:
                if totals.get(rep_id) is None:
                    totals[rep_id] = { 'total_seq': 0, 'total_count': 0 }
                totals[rep_id]['total_seq'] += 1
                totals[rep_id]['total_count'] += get_duplicate_count(row)
                found = False
                for c in cdr3s:
                    if get_gene(c['v_call']) != get_gene(row['v_call']):
                        continue
                    if get_gene(c['j_call']) != get_gene(row['j_call']):
                        continue
                    if len(c['junction_aa']) == len(row['junction_aa']):
                        d = Levenshtein.distance(c['junction_aa'], row['junction_aa'], score_cutoff=score_cutoff)
                        if d <= score_cutoff:
                            found = True
                            #print(d, c['junction_aa'], c['v_call'], c['j_call'], row['junction_aa'], row['v_call'], row['j_call'])
                            match = { 'repertoire_id': rep_id }
                            match['distance'] = d
                            match['rhab_sequence_id'] = c['sequence_id']
                            match['rhab_junction'] = c['junction_aa']
                            match['match_junction'] = row['junction_aa']
                            match['sequence_id'] = row['sequence_id']
                            match['clone_id'] = row['clone_id']
                            match['v_call'] = row['v_call']
                            match['j_call'] = row['j_call']
                            match['duplicate_count'] = get_duplicate_count(row)
                            if seqs.get(rep_id) is None:
                                seqs[rep_id] = {}
                            if seqs[rep_id].get(row['sequence_id']) is None:
                                seqs[rep_id][row['sequence_id']] = get_duplicate_count(row)
                            if rhabs.get(rep_id) is None:
                                rhabs[rep_id] = { }
                            if rhabs[rep_id].get(match['rhab_sequence_id']) is None:
                                rhabs[rep_id][match['rhab_sequence_id']] = { 'num_seq':0, 'num_count':0 }
                            rhabs[rep_id][match['rhab_sequence_id']]['num_seq'] += 1
                            rhabs[rep_id][match['rhab_sequence_id']]['num_count'] += get_duplicate_count(row)
                            writer.writerow(match)
                if found:
                    dwriter.write(row)

        fieldnames = ['repertoire_id', 'rhab_sequence_id', 'num_seq', 'num_count', 'total_seq', 'total_count']
        writer = csv.DictWriter(open(args.group_name + '.' + args.outfile + '_summary.csv','w'), fieldnames = fieldnames)
        writer.writeheader()
        for rep_id in rhabs:
            for s in rhabs[rep_id]:
                entry = { 'repertoire_id': rep_id, 'rhab_sequence_id': s, 'num_seq': rhabs[rep_id][s]['num_seq'], 'num_count': rhabs[rep_id][s]['num_count'], 'total_seq': totals[rep_id]['total_seq'], 'total_count': totals[rep_id]['total_count'] }
                writer.writerow(entry)


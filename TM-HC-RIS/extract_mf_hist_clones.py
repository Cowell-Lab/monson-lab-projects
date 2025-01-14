#
# mutation frequency histogram for clones
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

r_freqs = {}
r_counts = {}
s_freqs = {}
s_counts = {}
bins = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
bin_names = []
for i in range(0,len(bins)):
    bin_names.append('bin_' + str(bins[i]))

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='MF histogram for clones.')
    parser.add_argument('airr_metadata', type=str, help='AIRR repertoire metadata file name')
    parser.add_argument('airr_group', type=str, help='AIRR repertoire group file name')
    args = parser.parse_args()

    if args:
        data = airr.load_repertoire(args.airr_metadata)
        repertoires = { obj['repertoire_id'] : obj for obj in data['Repertoire'] }

        data = load_data_file(args.airr_group)
        groups = { obj['repertoire_group_id'] : obj for obj in data['RepertoireGroup'] }

        for group in groups:
            r_freqs[group] = {}
            r_counts[group] = {}
            s_freqs[group] = {}
            s_counts[group] = {}
            for rep in groups[group]['repertoires']:
                rep_id = rep['repertoire_id']
                r_freqs[group][rep_id] = []
                r_counts[group][rep_id] = []
                s_freqs[group][rep_id] = []
                s_counts[group][rep_id] = []
                for i in range(0,len(bins)):
                    r_freqs[group][rep_id].append(0.0)
                    r_counts[group][rep_id].append(0)
                    s_freqs[group][rep_id].append(0.0)
                    s_counts[group][rep_id].append(0)

        # bin by mutation frequency
        reader = csv.DictReader(open('combine_report.clone.csv', 'r'))
        bin1_writer = csv.DictWriter(open('bin_0.05.clone.csv', 'w'), fieldnames=reader.fieldnames)
        bin1_writer.writeheader();
        bin2_writer = csv.DictWriter(open('bin_0.1.clone.csv', 'w'), fieldnames=reader.fieldnames)
        bin2_writer.writeheader();
        bin3_writer = csv.DictWriter(open('bin_0.15.clone.csv', 'w'), fieldnames=reader.fieldnames)
        bin3_writer.writeheader();
        bin4_writer = csv.DictWriter(open('bin_0.2.clone.csv', 'w'), fieldnames=reader.fieldnames)
        bin4_writer.writeheader();
        bin5_writer = csv.DictWriter(open('bin_0.25.clone.csv', 'w'), fieldnames=reader.fieldnames)
        bin5_writer.writeheader();
        bin6_writer = csv.DictWriter(open('bin_0.3.clone.csv', 'w'), fieldnames=reader.fieldnames)
        bin6_writer.writeheader();
        bin7_writer = csv.DictWriter(open('bin_0.35.clone.csv', 'w'), fieldnames=reader.fieldnames)
        bin7_writer.writeheader();
        bin8_writer = csv.DictWriter(open('bin_0.4.clone.csv', 'w'), fieldnames=reader.fieldnames)
        bin8_writer.writeheader();
        bin9_writer = csv.DictWriter(open('bin_0.45.clone.csv', 'w'), fieldnames=reader.fieldnames)
        bin9_writer.writeheader();
        for row in reader:
            for i in range(0,len(bins)):
                if float(row['mu_freq_r']) <= bins[i]:
                    break
                if i == 0:
                    bin1_writer.writerow(row)
                if i == 1:
                    bin2_writer.writerow(row)
                if i == 2:
                    bin3_writer.writerow(row)
                if i == 3:
                    bin4_writer.writerow(row)
                if i == 4:
                    bin5_writer.writerow(row)
                if i == 5:
                    bin6_writer.writerow(row)
                if i == 6:
                    bin7_writer.writerow(row)
                if i == 7:
                    bin8_writer.writerow(row)
                if i == 8:
                    bin9_writer.writerow(row)
            for group in r_freqs:
                for rep in groups[group]['repertoires']:
                    if rep['repertoire_id'] == row['repertoire_id']:
                        rep_id = rep['repertoire_id']
                        for i in range(0,len(bins)):
                            if float(row['mu_freq_r']) <= bins[i]:
                                r_freqs[group][rep_id][i] += float(row['copy_freq'])
                                r_counts[group][rep_id][i] += 1
                                break
            for group in s_freqs:
                for rep in groups[group]['repertoires']:
                    if rep['repertoire_id'] == row['repertoire_id']:
                        rep_id = rep['repertoire_id']
                        for i in range(0,len(bins)):
                            if float(row['mu_freq_s']) <= bins[i]:
                                s_freqs[group][rep_id][i] += float(row['copy_freq'])
                                s_counts[group][rep_id][i] += 1
                                break

        # write output
        for group in r_freqs:
            filename = group + '.rf_hist_freq_clone.csv'
            fieldnames = ['repertoire_id']
            fieldnames.extend(bin_names)
            writer = csv.DictWriter(open(filename,'w'), fieldnames = fieldnames)
            writer.writeheader()
            for rep_id in r_freqs[group]:
                row = r_freqs[group][rep_id]
                entry = { 'repertoire_id': rep_id }
                for i in range(0,len(bins)):
                    entry[bin_names[i]] = row[i]
                writer.writerow(entry)

            filename = group + '.sf_hist_freq_clone.csv'
            fieldnames = ['repertoire_id']
            fieldnames.extend(bin_names)
            writer = csv.DictWriter(open(filename,'w'), fieldnames = fieldnames)
            writer.writeheader()
            for rep_id in s_freqs[group]:
                row = s_freqs[group][rep_id]
                entry = { 'repertoire_id': rep_id }
                for i in range(0,len(bins)):
                    entry[bin_names[i]] = row[i]
                writer.writerow(entry)

            filename = group + '.rf_hist_count_clone.csv'
            fieldnames = ['repertoire_id']
            fieldnames.extend(bin_names)
            writer = csv.DictWriter(open(filename,'w'), fieldnames = fieldnames)
            writer.writeheader()
            for rep_id in r_counts[group]:
                row = r_counts[group][rep_id]
                entry = { 'repertoire_id': rep_id }
                for i in range(0,len(bins)):
                    entry[bin_names[i]] = row[i]
                writer.writerow(entry)

            filename = group + '.sf_hist_count_clone.csv'
            fieldnames = ['repertoire_id']
            fieldnames.extend(bin_names)
            writer = csv.DictWriter(open(filename,'w'), fieldnames = fieldnames)
            writer.writeheader()
            for rep_id in s_counts[group]:
                row = s_counts[group][rep_id]
                entry = { 'repertoire_id': rep_id }
                for i in range(0,len(bins)):
                    entry[bin_names[i]] = row[i]
                writer.writerow(entry)

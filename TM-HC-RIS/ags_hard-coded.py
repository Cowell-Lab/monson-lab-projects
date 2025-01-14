#!/usr/bin/env python3
##########################################################################################
# Author: Scott Christley and Jared L. Ostmeyer
# Date Started: 2022-05-09
# Purpose: Attempt at calculating AGS scores
##########################################################################################

##########################################################################################
# Libraries
##########################################################################################

import csv
import json
import numpy as np

##########################################################################################
# Settings
##########################################################################################

#group_file = '../../../Data/B-cell_Human_Multiple-Sclerosis_Nancy-Monson/Metadata/repertoire_groups.airr.json'
#report_file = '../../../Data/B-cell_Human_Multiple-Sclerosis_Nancy-Monson/Metadata/mutational_report.frequency.repertoire.csv'
group_file = '../../repertoire_groups.airr.json'
report_file = 'mutational_report.frequency.repertoire.csv'

# Verify AGS positions 81 and 89 are 90 and 101, respectively
ags_pos = [ 36, 45, 64, 65, 90, 101 ]

mean = 1.6
std = 0.9

##########################################################################################
# Load RF
##########################################################################################

samples = {}
with open(report_file, 'r') as stream:
  reader = csv.DictReader(stream)
  for row in reader:
    sample = row['repertoire_id']
    values = []
    for pos in ags_pos:
      values.append(float(row['mu_freq_'+str(pos)+'_r']))
    samples[sample] = values

##########################################################################################
# Load groups
##########################################################################################

with open(group_file, 'r', encoding='utf-8') as handle:
  data = json.load(handle)
groups = { obj['repertoire_group_id'] : obj for obj in data['RepertoireGroup'] }

##########################################################################################
# Calculate AGS
##########################################################################################

print('Group', 'Sample', 'AGS', 'Above 6.8', sep=',')
for group in [ 'HC_PB_DNA', 'TM_PB_DNA', 'RIS_PB_DNA' ]:
  for reps in groups[group]['repertoires']:
    sample = reps['repertoire_id']
    f = np.array(samples[sample])
    ags = np.sum((f-mean)/std)
    print(group, sample, ags, 1 if ags > 6.8 else 0, sep=',')


#import airr
import json

rep_file = './repertoires.v2.airr.json'
group_file = './repertoire_groups.airr.json'

with open(group_file, 'r', encoding='utf-8') as handle:
    data = json.load(handle)
groups = { obj['repertoire_group_id'] : obj for obj in data['RepertoireGroup'] }

print(len(groups['HC_PB_DNA']['repertoires']), 'repertoires in HC_PB_DNA')
for reps in groups['HC_PB_DNA']['repertoires']:
    print(reps['repertoire_id'])

print(len(groups['TM_PB_DNA']['repertoires']), 'repertoires in TM_PB_DNA')
for reps in groups['TM_PB_DNA']['repertoires']:
    print(reps['repertoire_id'])

print(len(groups['RIS_PB_DNA']['repertoires']), 'repertoires in RIS_PB_DNA')
for reps in groups['RIS_PB_DNA']['repertoires']:
    print(reps['repertoire_id'])


import airr
import json
import csv

rep_file = '../../repertoires.v2.airr.json'
group_file = '../../repertoire_groups.airr.json'

with open(group_file, 'r', encoding='utf-8') as handle:
    data = json.load(handle)
groups = { obj['repertoire_group_id'] : obj for obj in data['RepertoireGroup'] }

for reps in groups['TM_PB_DNA']['repertoires']:
    fname = reps['repertoire_id'] + '.repertoire.cdr3_aa_sharing.tsv'
    reader = csv.DictReader(open(fname, 'r'), dialect='excel-tab')
    cnt = 0
    for row in reader:
        cnt = cnt + 1
    print(reps['repertoire_id'] + ',' + str(cnt))

#print(len(groups['TM_PB_DNA']['repertoires']), 'repertoires in TM_PB_DNA')
#for reps in groups['TM_PB_DNA']['repertoires']:
#    print(reps['repertoire_id'])

#print(len(groups['RIS_PB_DNA']['repertoires']), 'repertoires in RIS_PB_DNA')
#for reps in groups['RIS_PB_DNA']['repertoires']:
#    print(reps['repertoire_id'])

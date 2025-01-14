
import airr
import json

rep_file = '../../repertoires.v2.airr.json'
group_file = '../../repertoire_groups.airr.json'

with open(group_file, 'r', encoding='utf-8') as handle:
    data = json.load(handle)
groups = { obj['repertoire_group_id'] : obj for obj in data['RepertoireGroup'] }

junctions_HC = {}
junctions_RIS = {}
junctions_TM = {}

#print(len(groups['HC_PB_DNA']['repertoires']), 'repertoires in HC_PB_DNA')
#for reps in groups['HC_PB_DNA']['repertoires']:
#    print(reps['repertoire_id'])

for reps in groups['HC_PB_DNA']['repertoires']:
    fname = reps['repertoire_id'] + '.ighv4.ge3ags.mutations.airr.tsv'
    reader = airr.read_rearrangement(fname)
    for row in reader:
        junctions_HC[row['junction_aa']] = row['junction_aa']

for reps in groups['RIS_PB_DNA']['repertoires']:
    fname = reps['repertoire_id'] + '.ighv4.ge3ags.mutations.airr.tsv'
    reader = airr.read_rearrangement(fname)
    for row in reader:
        junctions_RIS[row['junction_aa']] = row['junction_aa']

for reps in groups['TM_PB_DNA']['repertoires']:
    fname = reps['repertoire_id'] + '.ighv4.ge3ags.mutations.airr.tsv'
    reader = airr.read_rearrangement(fname)
    for row in reader:
        junctions_TM[row['junction_aa']] = row['junction_aa']

for junction in junctions_RIS:
    if junctions_HC.get(junction) is not None:
        continue
    if junctions_TM.get(junction) is not None:
        print(junction)

#print(len(groups['TM_PB_DNA']['repertoires']), 'repertoires in TM_PB_DNA')
#for reps in groups['TM_PB_DNA']['repertoires']:
#    print(reps['repertoire_id'])

#print(len(groups['RIS_PB_DNA']['repertoires']), 'repertoires in RIS_PB_DNA')
#for reps in groups['RIS_PB_DNA']['repertoires']:
#    print(reps['repertoire_id'])

#
# Monson Lab TM/HC/RIS
#
# Create RepertoireGroups for comparative analysis
#

import airr
import json

rep_file = './repertoires.v2.airr.json'
group_file = './repertoire_groups.airr.json'

data = airr.load_repertoire(rep_file)
repertoires = { obj['repertoire_id'] : obj for obj in data['Repertoire'] }

data = json.load(open(group_file,'r'))
groups = { obj['repertoire_group_id'] : obj for obj in data['RepertoireGroup'] }

cohorts = [ "HC_PB_DNA", "TM_PB_DNA", "RIS_PB_DNA"]

for group in cohorts:
    print(group)
    for rep in groups[group]['repertoires']:
        rep_id = rep['repertoire_id']
        r = repertoires[rep_id]
        filename = r['data_processing'][0]['data_processing_files'][0]
        filename = filename.replace('.airr.tsv.gz', '.gene.clone.count.tsv')
        print(filename)

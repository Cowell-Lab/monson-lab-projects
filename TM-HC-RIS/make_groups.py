#
# Monson Lab TM/HC/RIS
#
# Create RepertoireGroups for comparative analysis
#

import airr
import json

rep_file = './repertoires.v2.airr.json'
out_file = './repertoire_groups.airr.json'

# low counts <100
# library8, UTSW23, 6528862801179644396-242ac114-0001-012
exclude_repertoires = ['6528862801179644396-242ac114-0001-012']
# library4, 2983 BC6, 7769631023962395116-242ac114-0001-012
exclude_repertoires.append('7769631023962395116-242ac114-0001-012')
# library4, 1743-7_S7 BC2, 7875974414211355116-242ac114-0001-012
exclude_repertoires.append('7875974414211355116-242ac114-0001-012')
# duplicate subject
# library8, 2929-2_S44, 6531740429267964396-242ac114-0001-012
exclude_repertoires.append('6531740429267964396-242ac114-0001-012')
# excluded by nancy
# library6, 278, 7881214274312475116-242ac114-0001-012
exclude_repertoires.append('7881214274312475116-242ac114-0001-012')
# library6, 279, 7881729674682962412-242ac114-0001-012
exclude_repertoires.append('7881729674682962412-242ac114-0001-012')
# library9, HD281, 7798622053210395116-242ac114-0001-012
exclude_repertoires.append('7798622053210395116-242ac114-0001-012')
# library_10, 1063, 153727805303614995-242ac118-0001-012
exclude_repertoires.append('153727805303614995-242ac118-0001-012')

# right before paper submission, found out that some of the TM subjects were mislabeled
# so we need to exclude them:
# subject 2147
exclude_repertoires.append('7778006210189595116-242ac114-0001-012')
# subject 2567
exclude_repertoires.append('7829459918395675116-242ac114-0001-012')
# subject 2478
exclude_repertoires.append('7796345720543515116-242ac114-0001-012')
# subject 2070
exclude_repertoires.append('7781957580101915116-242ac114-0001-012')
# subject 2929
exclude_repertoires.append('7795572626430235116-242ac114-0001-012')
# subject 2777
exclude_repertoires.append('7829030421666075116-242ac114-0001-012')
# subject 2934
exclude_repertoires.append('7876575709632795116-242ac114-0001-012')
# subject 2619
exclude_repertoires.append('6531439781557244396-242ac114-0001-012')
# subject 2620
exclude_repertoires.append('6528218556085244396-242ac114-0001-012')
# subject 2278
exclude_repertoires.append('152224562455047699-242ac118-0001-012')

data = airr.load_repertoire(rep_file)
reps = data['Repertoire']

cohorts = { "TM": [], "RIS": [], "HC": [],
    "TM_PB": [], "TM_PB_DNA": [], "TM_MB": [], "TM_MB_DNA": [], "TM_NB": [],
    "RIS_PB": [], "RIS_PB_DNA": [], "RIS_MB": [], "RIS_MB_DNA": [], "RIS_NB": [],
    "HC_PB": [], "HC_PB_DNA": [], "HC_MB": [], "HC_MB_DNA": [], "HC_NB": [] }

for r in reps:
    if r['repertoire_id'] in exclude_repertoires:
        continue
    cohorts[r['sample'][0]['disease_state_sample']].append({'repertoire_id':r['repertoire_id']})
    if r['sample'][0]['cell_subset']['id'] == 'CL:0000787':
        cohorts[r['sample'][0]['disease_state_sample'] + '_MB'].append({'repertoire_id':r['repertoire_id']})
        if r['sample'][0]['template_class'] == 'DNA':
            cohorts[r['sample'][0]['disease_state_sample'] + '_MB_DNA'].append({'repertoire_id':r['repertoire_id']})
    if r['sample'][0]['cell_subset']['id'] == 'CL:0000980':
        cohorts[r['sample'][0]['disease_state_sample'] + '_PB'].append({'repertoire_id':r['repertoire_id']})
        if r['sample'][0]['template_class'] == 'DNA':
            cohorts[r['sample'][0]['disease_state_sample'] + '_PB_DNA'].append({'repertoire_id':r['repertoire_id']})
    if r['sample'][0]['cell_subset']['id'] == 'CL:0000788':
        cohorts[r['sample'][0]['disease_state_sample'] + '_NB'].append({'repertoire_id':r['repertoire_id']})


groups = { 'RepertoireGroup': [] }
for c in cohorts:
    groups['RepertoireGroup'].append({'repertoire_group_id': c, 'repertoires': cohorts[c]})

with open(out_file, 'w') as f:
    json.dump(groups, f, indent=2)

for c in cohorts:
    print('group', c, 'has', len(cohorts[c]), 'repertoires.')
    
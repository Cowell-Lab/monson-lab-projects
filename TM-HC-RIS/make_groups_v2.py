#
# Monson Lab TM/HC/RIS
#
# After reviews from Annals of Neurology, we re-organized the comparison groups
#
# Create RepertoireGroups for comparative analysis
#

import airr
import json

rep_file = './repertoires.v2.airr.json'
out_file = './repertoire_groups.airr.json'
out_rep_file = './repertoires.group_only.v2.airr.json'

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

data = airr.read_airr(rep_file)
reps = data['Repertoire']
rep_dict = { obj['repertoire_id'] : obj for obj in data['Repertoire'] }

cohorts = { "RIS_PB_DNA": [], "HC_PB_DNA": [], "CIS_PB_DNA": [],
    "ADVANCING": [], "STABLE": [], "CORD_POS": [], "CORD_NEG": [],
    "RIS_ADVANCING": [], "RIS_STABLE": [], "RIS_POS": [], "RIS_NEG": [],
    "RIS_CSRHI": [], "RIS_CSRLO": [],
    "CIS_ADVANCING": [], "CIS_POS": [], "CIS_NEG": [] }

for r in reps:
    if r['repertoire_id'] in exclude_repertoires:
        continue
    if r['sample'][0]['disease_state_sample'] not in [ 'RIS', 'HC', 'CIS' ]:
        continue
    if r['sample'][0]['cell_subset']['id'] == 'CL:0000980':
        if r['sample'][0]['template_class'] == 'DNA':
            cohorts[r['sample'][0]['disease_state_sample'] + '_PB_DNA'].append({'repertoire_id':r['repertoire_id']})

            if 'CORD-POS' in r['subject']['diagnosis'][0]['study_group_description']:
                cohorts['CORD_POS'].append({'repertoire_id':r['repertoire_id']})
                if r['sample'][0]['disease_state_sample'] == 'RIS':
                    cohorts['RIS_POS'].append({'repertoire_id':r['repertoire_id']})
                if r['sample'][0]['disease_state_sample'] == 'CIS':
                    cohorts['CIS_POS'].append({'repertoire_id':r['repertoire_id']})
            if 'CORD-NEG' in r['subject']['diagnosis'][0]['study_group_description']:
                cohorts['CORD_NEG'].append({'repertoire_id':r['repertoire_id']})
                if r['sample'][0]['disease_state_sample'] == 'RIS':
                    cohorts['RIS_NEG'].append({'repertoire_id':r['repertoire_id']})
                if r['sample'][0]['disease_state_sample'] == 'CIS':
                    cohorts['CIS_NEG'].append({'repertoire_id':r['repertoire_id']})

            if 'CSRHI' in r['subject']['diagnosis'][0]['study_group_description']:
                cohorts['RIS_CSRHI'].append({'repertoire_id':r['repertoire_id']})
            if 'CSRLO' in r['subject']['diagnosis'][0]['study_group_description']:
                cohorts['RIS_CSRLO'].append({'repertoire_id':r['repertoire_id']})

            if r['subject']['diagnosis'][0]['disease_stage'] == 'advancing':
                cohorts['ADVANCING'].append({'repertoire_id':r['repertoire_id']})
                if r['sample'][0]['disease_state_sample'] == 'RIS':
                    cohorts['RIS_ADVANCING'].append({'repertoire_id':r['repertoire_id']})
                if r['sample'][0]['disease_state_sample'] == 'CIS':
                    cohorts['CIS_ADVANCING'].append({'repertoire_id':r['repertoire_id']})
            if r['subject']['diagnosis'][0]['disease_stage'] == 'stable':
                cohorts['STABLE'].append({'repertoire_id':r['repertoire_id']}) 
                if r['sample'][0]['disease_state_sample'] == 'RIS':
                    cohorts['RIS_STABLE'].append({'repertoire_id':r['repertoire_id']}) 
                if r['sample'][0]['disease_state_sample'] == 'CIS':
                    cohorts['CIS_STABLE'].append({'repertoire_id':r['repertoire_id']}) 

# MS2025
exclude_subjects = ["UTSW03", "UTSW06", "UTSW21", "UTSW27", "UTSW28", "1842", "1985", "2229", "2260", "2364", "2405", "3166"]
cohorts['MS2025_PB_DNA'] = []
for rep in cohorts['RIS_PB_DNA']:
    r = rep_dict[rep['repertoire_id']]
    if r['subject']['subject_id'] not in exclude_subjects:
        cohorts['MS2025_PB_DNA'].append({'repertoire_id':r['repertoire_id']})
for rep in cohorts['CIS_PB_DNA']:
    r = rep_dict[rep['repertoire_id']]
    if r['subject']['subject_id'] not in exclude_subjects:
        cohorts['MS2025_PB_DNA'].append({'repertoire_id':r['repertoire_id']})

unique_reps = {}
groups = { 'RepertoireGroup': [] }
for c in cohorts:
    groups['RepertoireGroup'].append({'repertoire_group_id': c, 'repertoires': cohorts[c]})
    for rep in cohorts[c]:
        unique_reps[rep['repertoire_id']] = rep

with open(out_file, 'w') as f:
    json.dump(groups, f, indent=2)

for c in cohorts:
    print('group', c, 'has', len(cohorts[c]), 'repertoires.')

group_reps = []
for rep in unique_reps:
    group_reps.append(rep_dict[rep])

print('total unique repetoires is', len(group_reps))
airr.write_airr(out_rep_file, {'Repertoire': group_reps})

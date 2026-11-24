#
# MOGAD project
# make custom repertoires for cell types
#
# designed to run in docker
# need the source and data directories mounted

import airr
import copy

job_id = "259e13b9-c272-4daf-8e86-4162f1a118a9-007"

# Define data directories
data_dir = "/data/NguyenLab/MOGAD/"
repcalc_dir = data_dir + 'vdjserver/'+job_id+'/'

# take the repertoire AIRR JSON from repcalc
# and generate new repertoire for each cell type

data = airr.read_airr(repcalc_dir + "repertoires.airr.json")

cell_types = ['Naive', 'Memory', 'Plasmablast', 'Atypical']
new_data = { 'Repertoire': [] }
for rep in data['Repertoire']:
    print(rep['subject']['subject_id'])
    for ct in cell_types:
        new_rep = copy.deepcopy(rep)
        new_rep['original_rep_id'] = new_rep['repertoire_id']
        new_rep['repertoire_id'] = "s" + new_rep['subject']['subject_id'] + "_" + ct
        new_data['Repertoire'].append(new_rep)

# make groups

airr.write_airr('repertoires.airr.json', new_data)

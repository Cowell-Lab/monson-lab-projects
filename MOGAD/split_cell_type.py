#
# MOGAD project
# split AIRR TSVs
#
# designed to run in docker
# need the source and data directories mounted

import csv
import airr

job_id = "c885ee98-c9ff-42a0-9079-7a4894754905-007"
processing_stage = "igblast.makedb"

# Define data directories
data_dir = "/data/NguyenLab/MOGAD/"
igblast_dir = data_dir + 'vdjserver/'+job_id+'/'

cell_types = ['Naive', 'Memory', 'Plasmablast', 'Atypical']
cell_map_file = data_dir + "B_cell_barcodes_celltype_v2.csv"
print(f"Reading cell map file {cell_map_file}")
cell_map = {}
reader = csv.DictReader(open(cell_map_file, 'r'))
for row in reader:
    r = row['sample_id'] + '_' + row['subset']
    if cell_map.get(r) is None:
        cell_map[r] = {}
    if row['in_main'] == 'TRUE':
        if cell_map[r].get(row['barcode']) is not None:
            print(f"error: duplicate barcode {row['barcode']} for sample {r}")
        cell_map[r][row['barcode']] = row

total_cnt = 0
for cm in cell_map:
    print(f"{cm} has {len(cell_map[cm])} barcodes")
    total_cnt += len(cell_map[cm])
print(f"Total barcodes is {total_cnt}")

data = airr.read_airr("repertoires.airr.json")
rep_dict = { obj['repertoire_id'] : obj for obj in data['Repertoire'] }

print("Split AIRR TSVs")
for rep in data['Repertoire']:
    print(f"Processing {rep['original_rep_id']} for new {rep['repertoire_id']}")
    input_airr = f"{igblast_dir}{rep['original_rep_id']}.{processing_stage}.airr.tsv"
    output_airr = rep['repertoire_id'] + ".airr.tsv"
    reader = airr.read_rearrangement(input_airr)
    writer = airr.derive_rearrangement(output_airr, input_airr)
    for row in reader:
        if cell_map[rep['repertoire_id']].get(row['adc_annotation_cell_id']):
            row['duplicate_count'] = 1
            writer.write(row)
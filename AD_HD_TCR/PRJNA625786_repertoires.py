#
# Conversion script for Monson Lab ACS TCR study
#
# First run the standard convert_repertoires.py
# Then run this to do specific cleanup
#

import airr
import csv

rep_file = './repertoires.airr.json'
out_file = './PRJNA625786.airr.json'

sra = {}
reader = csv.DictReader(open('SraRunInfo.csv', 'r'))
for row in reader:
    sra[row['LibraryName']] = row

data = airr.load_repertoire(rep_file)
reps = data['Repertoire']

for r in reps:
    # hand coded stuff

    for entry in r['sample']:
        entry['tissue'] = { 'id': 'UBERON:0001359', 'label': 'cerebrospinal fluid' }
        entry['sample_type'] = 'lumbar puncture'
        entry['cell_storage'] = True
        if entry['cell_subset']['label'] == 'CD4+':
            entry['cell_subset']['id'] = 'CL:0000624'
            entry['cell_subset']['id'] = 'CD4-positive, alpha-beta T cell'
        elif entry['cell_subset']['label'] == 'CD8+':
            entry['cell_subset']['id'] = 'CL:0000625'
            entry['cell_subset']['id'] = 'CD8-positive, alpha-beta T cell'
        else:
            print('unknown cell', entry['cell_subset'])

        # samples got renamed after SRA submission
        if entry['sample_id'] == 'AD1':
            entry['sequencing_run_id'] = sra['AD11']['Run']
        elif entry['sample_id'] == 'AD4':
            entry['sequencing_run_id'] = sra['AD14']['Run']
        elif entry['sample_id'] == 'AD5':
            entry['sequencing_run_id'] = sra['AD15']['Run']
        elif entry['sample_id'] == 'AD10':
            entry['sequencing_run_id'] = sra['AD16']['Run']
        elif entry['sample_id'] == 'AD12':
            entry['sequencing_run_id'] = sra['AD17']['Run']
        elif entry['sample_id'] == 'AD15':
            entry['sequencing_run_id'] = sra['AD21']['Run']
        elif entry['sample_id'] == '3093CD8':
            entry['sequencing_run_id'] = 'SRR16206633'
        elif sra.get(entry['sample_id']) is None:
            print('cannot find SRA record', entry['sample_id'])
        else:
            entry['sequencing_run_id'] = sra[entry['sample_id']]['Run']
        entry['sequencing_files']['read_length'] = 300
        entry['sequencing_files']['paired_read_length'] = 300
        entry['sequencing_platform'] = 'Illumina MiSeq'
        entry['sequencing_facility'] = 'UT Southwestern Medical Center'
        del entry['read_length']

    if r['subject']['subject_id'] in ['3091', '3122', '3143', '3312', '3379', '100101', '100112', '120120', '130010', '130076', '150018', '190010', '3093', '3519', '3542', '3675', '80068']:
        r['subject']['age_min'] = int(r['subject']['age'])
        r['subject']['age_max'] = int(r['subject']['age'])
        r['subject']['age_unit']['id'] = 'UO:0000036'
        r['subject']['age_unit']['label'] = 'year'
        r['subject']['diagnosis'][0]['study_group_description'] = 'ACS'
        r['subject']['diagnosis'][0]['disease_diagnosis']['id'] = 'DOID:10652'
        r['subject']['diagnosis'][0]['disease_diagnosis']['label'] = "Alzheimer's disease"
    elif r['subject']['subject_id'] in ['80093', '110039', '120006', '120056', '951309']:
        r['subject']['age_min'] = int(r['subject']['age'])
        r['subject']['age_max'] = int(r['subject']['age'])
        r['subject']['age_unit']['id'] = 'UO:0000036'
        r['subject']['age_unit']['label'] = 'year'
        r['subject']['diagnosis'][0]['study_group_description'] = 'Healthy'
    else:
        print('unknown subject', r['subject'])

    r['subject']['age_event'] = "sampling"
    r['subject']['species']['id'] = "NCBITAXON:9606"
    r['subject']['species']['label'] = "Homo sapiens"
    del r['subject']['age']

    r['data_processing'][0]['data_processing_id'] = '1626eab3-b98f-4140-91eb-b72d92d276c8-007'
    r['data_processing'][0]['primary_annotation'] = True
    r['data_processing'][0]['software_versions'] = 'igblast-stampede2-1.14u5'
    r['data_processing'][0]['germline_database'] = 'VDJServer IMGT 2019.01.23'

    files = []
    fname = ''
    for entry in r['sample']:
        fname = entry['sequencing_files']['filename']
        fname = fname.replace('.gz','.merged.unique.igblast.airr.tsv.gz')
        files.append(fname)
    r['data_processing'][0]['data_processing_files'] = [ ','.join(files) ]

airr.write_repertoire(out_file, reps)

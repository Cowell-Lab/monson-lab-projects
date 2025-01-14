#
# Conversion script for Monson Lab TM-HC-RIS IG study
#
# First run the standard convert_repertoires.py
# Then run this to do specific cleanup
#

import airr
import csv

rep_file = './repertoires.airr.json'
out_file = './IG_repertoires.airr.json'

# Not in SRA yet
#sra = {}
#reader = csv.DictReader(open('SraRunInfo.csv', 'r'))
#for row in reader:
#    sra[row['LibraryName']] = row

data = airr.load_repertoire(rep_file)
reps = data['Repertoire']

for r in reps:
    # hand coded stuff

    for entry in r['sample']:
        entry['tissue'] = { 'id': 'UBERON:0013756', 'label': 'peripheral blood' }
        entry['sample_type'] = 'peripheral venous puncture'
        entry['cell_storage'] = True
        if entry['cell_subset']['label'] == 'PB':
            entry['cell_subset']['id'] = 'CL:0000980'
            entry['cell_subset']['label'] = 'plasmablast'
        elif entry['cell_subset']['label'] == 'MB':
            entry['cell_subset']['id'] = 'CL:0000787'
            entry['cell_subset']['label'] = 'memory B cell'
        elif entry['cell_subset']['label'] == 'NB':
            entry['cell_subset']['id'] = 'CL:0000788'
            entry['cell_subset']['label'] = 'naive B cell'
        else:
            print('unknown cell', entry['cell_subset'])

        if entry['disease_state_sample'] == 'TM':
            r['subject']['diagnosis'][0]['study_group_description'] = 'TM'
            r['subject']['diagnosis'][0]['disease_diagnosis']['id'] = 'DOID:0080743'
            r['subject']['diagnosis'][0]['disease_diagnosis']['label'] = 'transverse myelitis'
        elif entry['disease_state_sample'] == 'HC':
            r['subject']['diagnosis'][0]['study_group_description'] = 'Healthy'
        elif entry['disease_state_sample'] == 'RIS':
            r['subject']['diagnosis'][0]['study_group_description'] = 'RIS'
        else:
            print('unknown disease_state_sample', entry['disease_state_sample'])

        entry['pcr_target'][0]['pcr_target_locus'] = 'IGH'
        entry['complete_sequences'] = 'partial'
        entry['physical_linkage'] = 'none'
        entry['sequencing_files']['read_length'] = 300
        entry['sequencing_files']['paired_read_length'] = 300
        entry['sequencing_platform'] = 'Illumina MiSeq'
        entry['sequencing_facility'] = 'UT Southwestern Medical Center'
        del entry['read_length']

    r['subject']['age_event'] = "sampling"
    r['subject']['species']['id'] = "NCBITAXON:9606"
    r['subject']['species']['label'] = "Homo sapiens"
    del r['subject']['age']

    if r['sample'][0]['sequencing_run_id'] == 'library_1':
        r['data_processing'][0]['data_processing_id'] = '1839ef2c-f8ce-47aa-a5e6-84a88bdffdf3-007'
    elif r['sample'][0]['sequencing_run_id'] == 'library_4':
        r['data_processing'][0]['data_processing_id'] = '8140a58b-9420-4f9e-a591-0f762806d848-007'
    elif r['sample'][0]['sequencing_run_id'] == 'library_6':
        if r['sample'][0]['template_class'] == 'DNA':
            r['data_processing'][0]['data_processing_id'] = 'd4fad922-79fe-4aa7-b0bb-ce6648288a79-007'
        elif r['sample'][0]['template_class'] == 'RNA':
            r['data_processing'][0]['data_processing_id'] = 'abdf135c-6177-4c21-8843-6057fed8b7b5-007'
        else:
            print('unknown template_class', r['sample'][0]['template_class'])
    elif r['sample'][0]['sequencing_run_id'] == 'library_9':
        if r['sample'][0]['template_class'] == 'DNA':
            r['data_processing'][0]['data_processing_id'] = 'b4bd375c-e974-4ebe-8bd6-6bdbe05e05bc-007'
        elif r['sample'][0]['template_class'] == 'RNA':
            r['data_processing'][0]['data_processing_id'] = '38b2524a-e3fc-4d4e-a3fa-25e5bd25e417-007'
        else:
            print('unknown template_class', r['sample'][0]['template_class'])
    elif r['sample'][0]['sequencing_run_id'] == 'library_8':
        r['data_processing'][0]['data_processing_id'] = 'a0908af9-7562-45b3-8fb8-eef18f9dcb74-007'
    else:
        print('unknown sequencing_run_id', r['sample'][0]['sequencing_run_id'])

    r['data_processing'][0]['primary_annotation'] = True
    r['data_processing'][0]['software_versions'] = 'vdj_pipe-stampede2-0.1.7u5 igblast-stampede2-1.17u1'
    r['data_processing'][0]['germline_database'] = 'VDJServer IMGT 2019.01.23'

    if r['sample'][0]['sequencing_run_id'] == 'library_4':
        files = []
        fname = ''
        for entry in r['sample']:
            fname = entry['sequencing_files']['filename']
            fname = fname.replace('.fasta','.igblast.makedb.airr.tsv.gz')
            files.append(fname)
        r['data_processing'][0]['data_processing_files'] = [ ','.join(files) ]
    else:
        files = []
        fname = ''
        for entry in r['sample']:
            fname = entry['sequencing_files']['filename']
            fname = fname.replace('.gz','.merged.unique.igblast.makedb.airr.tsv.gz')
            files.append(fname)
        r['data_processing'][0]['data_processing_files'] = [ ','.join(files) ]

airr.write_repertoire(out_file, reps)

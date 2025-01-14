#
# Monson Lab TM/HC/RIS
#
# Filter sequences based on SHM ratios
#

import argparse
import airr
import csv

# low counts <100
# library8, UTSW22, 6529893593330684396-242ac114-0001-012
# library4, 2983 BC6, 7769631023962395116-242ac114-0001-012
exclude_repertoires = ['6529893593330684396-242ac114-0001-012', '7769631023962395116-242ac114-0001-012']

def get_duplicate_count(fields):
    if fields.get('duplicate_count') is None: return 1
    else: return int(fields['duplicate_count'])

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Split files based on SHM.')
    parser.add_argument('airr_metadata', type=str, help='AIRR repertoire metadata file name')
    args = parser.parse_args()

    if args:
        data = airr.load_repertoire(args.airr_metadata)
        reps = data['Repertoire']
        print('Loaded', len(reps), 'repertoires.')

        for r in reps:
            #filename = r['data_processing'][0]['data_processing_files'][0]
            #infile = filename.replace('.airr.tsv.gz', '.gene.mutations.airr.tsv')
            #outname0 = filename.replace('.airr.tsv.gz', '.0shm.mutations.airr.tsv')
            #outname3 = filename.replace('.airr.tsv.gz', '.gt3shm.mutations.airr.tsv')
            infile = r['repertoire_id'] + '.allele.mutations.airr.tsv'
            outname0 = r['repertoire_id'] + '.0shm.mutations.airr.tsv'
            outname3 = r['repertoire_id'] + '.gt3shm.mutations.airr.tsv'
            print('Processing:', infile,'to:',outname0,'and',outname3)

            try:
                reader = airr.read_rearrangement(infile)
                writer0 = airr.derive_rearrangement(outname0, infile)
                writer3 = airr.derive_rearrangement(outname3, infile)
                for row in reader:
                    if int(row['mu_count_v_r']) + int(row['mu_count_v_s']) == 0:
                        writer0.write(row)
                    if int(row['mu_count_v_r']) + int(row['mu_count_v_s']) > 3:
                        writer3.write(row)

            except:
                print('WARNING: Could not read file:', infile, '... skipping.')


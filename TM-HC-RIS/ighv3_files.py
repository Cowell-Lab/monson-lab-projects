#
# Monson Lab TM/HC/RIS
#
# Filter sequences for just IGHV3
#

import argparse
import airr
import csv

#dir = '../mutations/'
#in_stage = 'gene.mutations'
#out_stage = 'ighv3.mutations'

#dir = '../stats/'

#dir = '../gt3shm/'
#in_stage = 'gt3shm.mutations.aa_properties'
#out_stage = 'gt3shm.ighv3.mutations.aa_properties'

#dir = '../gt3shm/'
#in_stage = 'gt3shm.mutations'
#out_stage = 'gt3shm.ighv3.mutations'

#dir = '../'
#in_stage = 'ighj6.gene.mutations.aa_properties'
#out_stage = 'ighv3.ighj6.gene.mutations.aa_properties'

dir = '../'
in_stage = 'ighj6.gt3shm.mutations.aa_properties'
out_stage = 'ighv3.ighj6.gt3shm.mutations.aa_properties'

def get_duplicate_count(fields):
    if fields.get('duplicate_count') is None: return 1
    else: return int(fields['duplicate_count'])

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Split files for IGHV3.')
    parser.add_argument('airr_metadata', type=str, help='AIRR repertoire metadata file name')
    args = parser.parse_args()

    if args:
        data = airr.load_repertoire(args.airr_metadata)
        reps = data['Repertoire']
        print('Loaded', len(reps), 'repertoires.')

        for r in reps:
            filename = r['repertoire_id'] + '.' + in_stage + '.airr.tsv'
            infile = dir + filename
            outname0 = filename.replace('.' + in_stage + '.airr.tsv', '.' + out_stage + '.airr.tsv')
            print('Processing:', infile,'to:',outname0)

            try:
                reader = airr.read_rearrangement(infile)
                writer0 = airr.derive_rearrangement(outname0, infile)
                for row in reader:
                    if 'IGHV3' in row['v_call']:
                        writer0.write(row)

            except:
                print('WARNING: Could not read file:', infile, '... skipping.')


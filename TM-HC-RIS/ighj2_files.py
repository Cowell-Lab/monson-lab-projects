#
# Monson Lab TM/HC/RIS
#
# Filter sequences for just IGHJ2
#

import argparse
import airr
import csv

dir = '../'
#dir = '../stats/'
#in_stage = 'ighv4.mutations'
#out_stage = 'ighv4.ighj2.mutations'
#in_stage = 'gene.mutations.aa_properties'
#out_stage = 'ighj2.gene.mutations.aa_properties'
in_stage = 'gt3shm.mutations.aa_properties'
out_stage = 'ighj2.gt3shm.mutations.aa_properties'

def get_duplicate_count(fields):
    if fields.get('duplicate_count') is None: return 1
    else: return int(fields['duplicate_count'])

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Split files for IGHJ2.')
    parser.add_argument('airr_metadata', type=str, help='AIRR repertoire metadata file name')
    args = parser.parse_args()

    if args:
        data = airr.read_airr(args.airr_metadata)
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
                    if 'IGHJ2' in row['j_call']:
                        writer0.write(row)

            except:
                print('WARNING: Could not read file:', infile, '... skipping.')


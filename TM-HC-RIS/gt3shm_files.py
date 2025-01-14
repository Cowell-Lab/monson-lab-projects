#
# Monson Lab TM/HC/RIS
#
# Filter sequences for just 3+ SHM
#

import argparse
import airr
import csv

dir = '../new_mutations/'

def get_duplicate_count(fields):
    if fields.get('duplicate_count') is None: return 1
    else: return int(fields['duplicate_count'])

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Split files for IGHV4.')
    parser.add_argument('airr_metadata', type=str, help='AIRR repertoire metadata file name')
    args = parser.parse_args()

    if args:
        data = airr.load_repertoire(args.airr_metadata)
        reps = data['Repertoire']
        print('Loaded', len(reps), 'repertoires.')

        for r in reps:
            filename = r['repertoire_id'] + '.gene.mutations.airr.tsv'
            infile = dir + filename
            outname0 = filename.replace('.gene.mutations.airr.tsv', '.gt3shm.mutations.airr.tsv')
            print('Processing:', infile,'to:',outname0)

            try:
                reader = airr.read_rearrangement(infile)
                writer0 = airr.derive_rearrangement(outname0, infile)
                for row in reader:
                    if int(row['mu_count_v_r']) + int(row['mu_count_v_s']) > 3:
                        writer0.write(row)

            except:
                print('WARNING: Could not read file:', infile, '... skipping.')


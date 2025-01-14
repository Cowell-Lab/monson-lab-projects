#
# Monson Lab TM/HC/RIS
#
# Filter sequences for just IGHV4 and RF at >=3 AGS positions
#

import argparse
import airr
import csv

dir = '../ighv4_mutations/'

# AGS positions in IMGT numbering
ags_pos = [ 36, 45, 64, 65, 90, 101 ]

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Filter files.')
    parser.add_argument('airr_metadata', type=str, help='AIRR repertoire metadata file name')
    args = parser.parse_args()

    if args:
        data = airr.load_repertoire(args.airr_metadata)
        reps = data['Repertoire']
        print('Loaded', len(reps), 'repertoires.')

        for r in reps:
            filename = r['repertoire_id'] + '.ighv4.mutations.airr.tsv'
            infile = dir + filename
            outname0 = filename.replace('.ighv4.mutations.airr.tsv', '.ighv4.2ags.mutations.airr.tsv')
            print('Processing:', infile,'to:',outname0)

            try:
                reader = airr.read_rearrangement(infile)
                writer0 = airr.derive_rearrangement(outname0, infile)
                for row in reader:
                    rf_cnt = 0
                    for pos in ags_pos:
                        if float(row['mu_count_'+str(pos)+'_r_aa']) > 0:
                            rf_cnt += 1
                    if rf_cnt == 2:
                        writer0.write(row)

            except:
                print('WARNING: Could not read file:', infile, '... skipping.')


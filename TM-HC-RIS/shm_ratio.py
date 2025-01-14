#
# Monson Lab TM/HC/RIS
#
# Calculate SHM ratios
#

import argparse
import airr
import csv
import json

def get_duplicate_count(fields):
    if fields.get('duplicate_count') is None: return 1
    else: return int(fields['duplicate_count'])

names = ['repertoire_id', '#seq 0 SHM', '#seq >3 SHM', 'total #seq', 'abundance 0 SHM', 'abundance >3 SHM', 'total abundance']

if (__name__=="__main__"):
    parser = argparse.ArgumentParser(description='Calculate SHM ratios.')
    parser.add_argument('airr_metadata', type=str, help='AIRR repertoire metadata file name')
    parser.add_argument('airr_group', type=str, help='AIRR repertoire group file name')
    args = parser.parse_args()

    if args:
        data = airr.load_repertoire(args.airr_metadata)
        reps_all = { obj['repertoire_id'] : obj for obj in data['Repertoire'] }
        print('Loaded', len(reps_all), 'repertoires.')

        with open(args.airr_group, 'r', encoding='utf-8') as handle:
            data = json.load(handle)
        groups = { obj['repertoire_group_id'] : obj for obj in data['RepertoireGroup'] }
        print('Loaded', len(groups), 'repertoire groups.')

        for group in groups:
            reps = groups[group]['repertoires']
            print(group)

            writer = csv.DictWriter(open(group + '.shm_ratios.csv', 'w'), fieldnames=names)
            writer.writeheader()

            for gr in reps:
                r = reps_all[gr['repertoire_id']]
                #sample_id = r['sample'][0]['sample_id']
                #subject_id = r['subject']['subject_id']
                #filename = r['data_processing'][0]['data_processing_files'][0]
                #filename = filename.replace('.airr.tsv.gz', '.gene.mutations.airr.tsv')
                filename = r['repertoire_id'] + '.gene.mutations.airr.tsv'
                print('Processing:', filename)

                shm0_seq_cnt = 0
                shm3_seq_cnt = 0
                total_seq_cnt = 0
                shm0_ab_cnt = 0
                shm3_ab_cnt = 0
                total_ab_cnt = 0
                try:
                    reader = airr.read_rearrangement(filename)
                    for row in reader:
                        total_ab_cnt += get_duplicate_count(row)
                        total_seq_cnt += 1
                        if int(row['mu_count_v_r']) + int(row['mu_count_v_s']) == 0:
                            shm0_ab_cnt += get_duplicate_count(row)
                            shm0_seq_cnt += 1
                        if int(row['mu_count_v_r']) + int(row['mu_count_v_s']) > 3:
                            shm3_ab_cnt += get_duplicate_count(row)
                            shm3_seq_cnt += 1

                    entry = { 'repertoire_id':r['repertoire_id'] }
                    entry['#seq 0 SHM'] = shm0_seq_cnt
                    entry['#seq >3 SHM'] = shm3_seq_cnt
                    entry['total #seq'] = total_seq_cnt
                    entry['abundance 0 SHM'] = shm0_ab_cnt
                    entry['abundance >3 SHM'] = shm3_ab_cnt
                    entry['total abundance'] = total_ab_cnt
                    writer.writerow(entry)

                except:
                    print('Could not read file:', filename, '... skipping.')


import pandas as pd
import os
import re

#file location and storage
current_dir = os.getcwd()
output_file = os.path.join(current_dir, 'Kreye_2021_vdj_pairmaster_20250721.xlsx')

#emtpy list to hold dataframes
all_data = []

#BR codes
seqID_list = sorted(['113-101', '113-115', '113-175', '113-198', '113-201'], key=len, reverse=True)

#loop through .txt files in the folder
for filename in os.listdir(current_dir):
    if filename.endswith('.tsv'):
      filepath = os.path.join(current_dir, filename)
      df = pd.read_csv(filepath, delimiter='\t')
      df['source_file'] = filename
      all_data.append(df)

#combine df into one
merged_df = pd.concat(all_data, ignore_index=True)

#ensure productive and chain type are strings
merged_df['sequence_id'] = merged_df['sequence_id'].astype(str)
merged_df['productive'] = merged_df['productive'].astype(str).str.lower()
merged_df = merged_df[merged_df['productive'] == 'true']

#extract seq ID
def extract_br_code_from_seqid(seq_id):
    for code in seqID_list:
        if code in seq_id:
            return code
    return None
merged_df['seqID'] = merged_df['sequence_id'].apply(extract_br_code_from_seqid)

#extract chain type from seq id
merged_df['chain_type'] = merged_df['sequence_id'].str.extract(r'\.(heavy|kappa|lambda)\.')

#create a match key, extract well id
merged_df['match_key'] = merged_df['seqID']

#seperate into heavy and light chains
heavy_df = merged_df[merged_df['chain_type'] == 'heavy'].copy()
light_df = merged_df[merged_df['chain_type'].isin(['kappa', 'lambda'])].copy()


#add clear prefixes
heavy_df = heavy_df.add_prefix('HC_')
light_df = light_df.add_prefix('LC_')

#restore original match key for merging
heavy_df['match_key'] = heavy_df['HC_match_key']
light_df['match_key'] = light_df['LC_match_key']

#merge heavy and light chains on match key
paired_df = pd.merge(heavy_df, light_df, on='match_key', how='inner', suffixes=('_HC', '_LC'))

#keep only rows with heavy and light chains
paired_df['seqID'] = paired_df['HC_seqID']

#define desired columns
pairmaster_columns = ['sequence_id', 'productive', 'v_call', 'j_call',
    'fwr1', 'cdr1', 'fwr2', 'cdr2', 'fwr3', 'cdr3', 'fwr4',
    'cdr3_aa_length', 'cdr3_aa_charge', 'mu_aa_count_v_r',
    'v_identity', 'sequence_alignment']

final_columns = ['seqID']
for prefix in ['HC_', 'LC_']:
	final_columns+= [prefix + col for col in pairmaster_columns]

final_df = paired_df[final_columns]

#add order ID columns after sequence_id
final_df.insert(final_df.columns.get_loc('HC_sequence_id') + 1, 'HC_Order ID', '')
final_df.insert(final_df.columns.get_loc('LC_sequence_id') + 1, 'LC_Order ID', '')

#add reference columns after sequence alignment
final_df.insert(final_df.columns.get_loc('HC_sequence_alignment') + 1, 'HC_reference', 'Kreye-JEM-2021')
final_df.insert(final_df.columns.get_loc('LC_sequence_alignment') + 1, 'LC_reference', 'Kreye-JEM-2021')


#save merged DF to excel
final_df.to_excel(output_file, index=False, engine='xlsxwriter')




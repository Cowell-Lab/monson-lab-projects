"""
Create PairMaster Format from nucleotide sequence and summary tsv
"""

import argparse
import pandas as pd

field_names = ['Subject_ID', 'Sample_ID']
heavy_names = ['Sequence_ID_HC', 'VH_Gene', 'JH_Gene',
                'FWR1_HC', 'CDR1_HC', 'FWR2_HC', 'CDR2_HC', 'FWR3_HC', 'CDR3_HC',
                'CDR3_AA_Length_HC', 'CDR3_AA_Charge_HC', 'Replacement_Mutations_HC',
                '%Homology_HC', 'Nucleotide_Sequence_HC', 'Analysis_Reference_HC']
light_names = ['Sequence_ID_LC', 'VL_Gene', 'JL_Gene',
                'FWR1_LC', 'CDR1_LC', 'FWR2_LC', 'CDR2_LC', 'FWR3_LC', 'CDR3_LC',
                'CDR3_AA_Length_LC', 'CDR3_AA_Charge_LC', 'Replacement_Mutations_LC',
                '%Homology_LC', 'Nucleotide_Sequence_LC', 'Analysis_Reference_LC']
field_names += heavy_names
field_names += light_names

def pairmaster(
        summary:str,
        order:str,
        output:str=None,
    ) -> None:
    """
    Creates the PairMaster format

    Parameters
    ----------
    summary : str
        - The file path to the summary.tsv file
    order : str
        - The file path to the to_order.csv file created by ordering.py
    output : str
        - When set, the function will write the pairmaster as a csv to the filepath given.

    Returns
    -------
    df : pd.DataFrame
        - A Pandas DataFrame with the PairMaster format.
    """
    summary_df = pd.read_csv(summary, delimiter='\t')
    order_df = pd.read_csv(order, usecols=['sequence_id', 'to_order'])
    sum_ord_df = summary_df.merge(order_df, how='left', on='sequence_id')

    sum_ord_df['subject_id'] = sum_ord_df['sequence_id'].str.split('.').str[0].str.split('-').str[0]
    sum_ord_df['sample_id'] = sum_ord_df['sequence_id'].str.split('.').str[0].str.split('-').str[1]

    grouped_df = sum_ord_df.groupby(['subject_id', 'sample_id'])

    pairmaster = pd.DataFrame(columns=field_names)
    for sub_samp_id, group_df in grouped_df:
        group_df = group_df.reset_index()

        print(f'Processing subject_id: {sub_samp_id[0]}, sample_id: {sub_samp_id[1]}')

        heavy = []
        light = []
        for i in range(group_df.shape[0]):
            # print(group_df)
            # if group_df.loc[i,'productive']:
            if group_df.loc[i,'sequence_id'].split('.')[1] == 'heavy':
                heavy.append(group_df.loc[i,:])
            else:
                light.append(group_df.loc[i,:])

        for h in heavy:
            for l in light:
                new_row = pd.DataFrame(columns=field_names)

                # both
                new_row.loc[i,'Subject_ID'] = h['subject_id']
                new_row.loc[i,'Sample_ID'] = h['sample_id']

                # HC
                new_row.loc[i,'Sequence_ID_HC'] = h['sequence_id']
                new_row.loc[i,'VH_Gene'] = h['v_call']
                new_row.loc[i,'JH_Gene'] = h['j_call']

                new_row.loc[i,'FWR1_HC'] = h['fwr1']
                new_row.loc[i,'CDR1_HC'] = h['cdr1']
                new_row.loc[i,'FWR2_HC'] = h['fwr2']
                new_row.loc[i,'CDR2_HC'] = h['cdr2']
                new_row.loc[i,'FWR3_HC'] = h['fwr3']
                new_row.loc[i,'CDR3_HC'] = h['cdr3']
                new_row.loc[i,'CDR3_AA_Length_HC'] = h['cdr3_aa_length']
                new_row.loc[i,'CDR3_AA_Charge_HC'] = h['cdr3_aa_charge']
                new_row.loc[i,'Replacement_Mutations_HC'] = h['mu_aa_count_v_r']
                new_row.loc[i,'%Homology_HC'] = h['v_identity']
                new_row.loc[i,'Nucleotide_Sequence_HC'] = h['to_order']
                new_row.loc[i,'Analysis_Reference_HC'] = ''

                # LC
                new_row.loc[i,'Sequence_ID_LC'] = l['sequence_id']
                new_row.loc[i,'VL_Gene'] = l['v_call']
                new_row.loc[i,'JL_Gene'] = l['j_call']

                new_row.loc[i,'FWR1_LC'] = l['fwr1']
                new_row.loc[i,'CDR1_LC'] = l['cdr1']
                new_row.loc[i,'FWR2_LC'] = l['fwr2']
                new_row.loc[i,'CDR2_LC'] = l['cdr2']
                new_row.loc[i,'FWR3_LC'] = l['fwr3']
                new_row.loc[i,'CDR3_LC'] = l['cdr3']
                new_row.loc[i,'CDR3_AA_Length_LC'] = l['cdr3_aa_length']
                new_row.loc[i,'CDR3_AA_Charge_LC'] = l['cdr3_aa_charge']
                new_row.loc[i,'Replacement_Mutations_LC'] = l['mu_aa_count_v_r']
                new_row.loc[i,'%Homology_LC'] = l['v_identity']
                new_row.loc[i,'Nucleotide_Sequence_LC'] = l['to_order']
                new_row.loc[i,'Analysis_Reference_LC'] = ''

                # create dataframe
                pairmaster = pd.concat([pairmaster, new_row], ignore_index=True)

    if output is not None:
        pairmaster.to_csv(output, index=False)

    return pairmaster

def main():
    """
    Runs the pairmaster function with the given arguments
    """
    parser = argparse.ArgumentParser(
        prog='PairMaster Script',
        description="Produce PairMaster Format")
    parser.add_argument('--summary', required=True,
                        help='Path to input summary tsv file; /path/to/<FILE_NAME>.summary.tsv')
    parser.add_argument('--order', required=True,
                        help='Order filepath (from ordering.py); /path/to/<FILE_NAME>.to_order.csv')
    parser.add_argument('--output', required=False, default=None,
                        help='Output filepath and name; /path/to/<OUTPUT_FILE_NAME>.pairmaster.csv')
    args = parser.parse_args()

    pairmaster(args.summary, args.order, args.output)

if __name__ == "__main__":
    main()

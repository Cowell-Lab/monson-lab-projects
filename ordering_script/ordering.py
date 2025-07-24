import pandas as pd
import argparse

consts_df = {'vector' : ['IGH', 'IGK', 'IGL'], # 'HC', 'KLC', 'LLC'
             'start' : ['CTGCAACCGGTGTACATTCC', 'TGTGCTGCAACCGGTGTACAT', 'CTGCTACCGGT'],
             'end' : ['GCGTCGACTTCGCA', 'CGTACGGTGGC', 'GGTCAGCCCAAGGCTGCCCCATCGGTCACTCTGTTCCCGCCCTCGAGTGAGGAG'],
             'restriction sites' : [('ACCGGT','GTCGAC'),('ACCGGT','CGTACG'),('ACCGGT','CTCGAG')],
             'reading frame' : [1, 0, 1]}

codons = {'TTT':'TTC', 'TTC':'TTT', # phe 
          'TTA':'TTG', 'TTG':'TTA', # leu
          'TCT':'TCC', 'TCC':'TCA', 'TCA':'TCG', 'TCG':'TCT', # ser
          'TAT':'TAC', 'TAC':'TAT', # tyr
          'TGT':'TGC', 'TGC':'TGT', # cys
          'CTT':'CTC', 'CTC':'CTA', 'CTA':'CTG', 'CTG':'CTT', # leu
          'CCT':'CCC', 'CCC':'CCA', 'CCA':'CCG', 'CCG':'CCT', # pro
          'CAT':'CAC', 'CAC':'CAT', # his
          'CAA':'CAG', 'CAG':'CAA', # gln
          'CGT':'CGC', 'CGC':'CGA', 'CGA':'CGG', 'CGG':'CGT', # arg
          'ATT':'ATC', 'ATC':'ATA', 'ATA':'ATT', # ile
          'ATG':'ATG', # met
          'ACT': 'ACC', 'ACC':'ACA', 'ACA':'ACG', 'ACG':'ACT', # thr
          'AAT':'AAC', 'AAC':'AAT', # asn
          'AAA':'AAG', 'AAG':'AAA', # lys
          'AGT':'AGC', 'AGC':'AGT', # ser
          'AGA':'AGG', 'AGG':'AGA', # arg
          'GTT':'GTC', 'GTC':'GTA', 'GTA':'GTG', 'GTG':'GTT', # val
          'GCT':'GCC', 'GCC':'GCA', 'GCA':'GCG', 'GCG':'GCT', # ala
          'GAT':'GAC', 'GAC':'GAC', # asp
          'GAA':'GAG', 'GAG':'GAA', # glu
          'GGT':'GGC', 'GGC':'GGA', 'GGA':'GGG', 'GGG':'GGT',  # gly
          'TAA':'TAG', 'TAG':'TGA', 'TGA':'TAA' #stop
         }

stop_codons = ['TAA, TAG, TGA']

def restriction_replacement(seq_id, seq, selector, i):

    '''
    Replaces the first codon in a restriction site that is not within the respective start or end regions.
    The first restriction site should be seen in the start, the second restriction site in the end.
    
    Parameters
    ----------
    seq : string
        - The nucleotide sequence to replace restriction sites
    selector : int
        - Based on the index of consts_df['vector']. e.g. 'IGH'==0, 'IGK'==1, 'IGL'==2
    
    Output
    ------
    new_seq : string
        - The altered nucleotide sequence with replaced codons and no internal restriction sites.
    '''
    new_seq = seq

    rsites = consts_df['restriction sites'][selector]
    ustart = consts_df['start'][selector]
    uend = consts_df['end'][selector]
    
    # find sites
    # first restriction site outside of (or overlapping with) the universal start (rs needed there)
    i1, i2 = 0, 0
    while i1 != -1 and i2 != -1:
        # first restriction site outside of (or overlapping with) the universal start (rs needed there)
        i1 = seq[(len(ustart) - len(rsites[0]) + 1):].find(rsites[0])
        if i1 != -1:
            i1 = i1 + len(ustart) - len(rsites[0]) + 1
            i1 = i1 - (i1 + consts_df['reading frame'][selector]) % 3
            while i1 < 0:
                i1 += 3
            new_seq = seq[0:i1] + codons[seq[i1:i1+3]] + seq[i1+3:]

            if codons[seq[i1:i1+3]].isin(stop_codons):
                print(f'STOP CODON ENCOUNTERED when replacing restriction site 1 @ {i}; {seq_id}')
            
        # second restriction site outside of (or overlapping with) the universal end (rs needed there)
        i2 = seq[:len(new_seq)-len(uend)].find(rsites[1])
        if i2 != -1:
            i2 = i2 - (i2 + consts_df['reading frame'][selector]) % 3
            while i2 < 0:
                i2 += 3
            new_seq = seq[0:i2] + codons[seq[i2:i2+3]] + seq[i2+3:]

            if codons[seq[i2:i2+3]].isin(stop_codons):
                print(f'STOP CODON ENCOUNTERED when replacing restriction site 1 @ {i}; {seq_id}')
        
        

    return new_seq

def main():
    parser = argparse.ArgumentParser(
        prog='Ordering Script',
        description="Produce orderable vectors")
    parser.add_argument('--data', required=True, help='Path to input igblast makedb airr tsv file; /path/to/<FILE_NAME>.igblast.makedb.airr.tsv')
    parser.add_argument('--v_call', required=True, help='Path to v-call reference; /path/to/genes_v_call.csv')
    args = parser.parse_args()

    # read data
    # igblast_makedb = pd.read_csv('./data/04072023_YW.igblast.makedb.airr.tsv', delimiter='\t')
    # var_df = pd.read_csv('./data/genes_v_call.csv').drop(columns=['Unnamed: 0'])
    igblast_makedb = pd.read_csv(args.data, delimiter='\t')
    var_df = pd.read_csv(args.v_call).drop(columns=['Unnamed: 0'])

    # make order df
    order_df = igblast_makedb.loc[:,['sequence_id']]
    order_df['to_order'] = ['']*len(igblast_makedb)
    for i in range(len(igblast_makedb)): # *** TEMP loop through all of igblast_makedb
        # select dictionary index
        locus = igblast_makedb.loc[i,'locus']
        for ii in range(len(consts_df['vector'])):
            if locus == consts_df['vector'][ii]:
                dict_selector = ii
                break

        # select variable region from file
        if locus != 'IGH':
            curr_v_call = igblast_makedb.loc[i, 'v_call']
            # Check if multiple genes provided in v_call
            # If same codon used for each gene then good, otherwise ERROR printed
            if len(curr_v_call.split(',')) > 1:
                v_calls = curr_v_call.split(',')
                codon_list = [var_df[var_df['v_call'] == v_call]['var'].item() for v_call in v_calls]
                codon_set = set(codon_list)

                if len(codon_set) == 1:
                    var = codon_list[0]
                else:
                    # Default behavior to produce an error when two v_call are given with different codons
                    print('ERROR @', i)
                    order_df.loc[i, 'to_order'] = 'ERROR'
                    continue
            else:
                var = var_df[var_df['v_call']==curr_v_call]['var'].values[0]
        else:
            var = ''
        
        # Construct insert
        order_df.loc[i, 'to_order'] = str(consts_df['start'][dict_selector] + \
                                        var + \
                                        igblast_makedb.loc[i,'fwr1'] + igblast_makedb.loc[i,'cdr1'] + \
                                        igblast_makedb.loc[i,'fwr2'] + igblast_makedb.loc[i,'cdr2'] + \
                                        igblast_makedb.loc[i,'fwr3'] + igblast_makedb.loc[i,'cdr3'] + \
                                        igblast_makedb.loc[i,'fwr4'] + \
                                        consts_df['end'][dict_selector]).translate({ord('.'): None})

    # find and replace restriction site
    for i in range(len(order_df)):
        prev_seq = order_df.loc[i, 'to_order']
        order_df.loc[i, 'to_order'] = restriction_replacement(order_df.loc[i, 'sequence_id'], order_df.loc[i, 'to_order'], dict_selector, i)
        if prev_seq != order_df.loc[i, 'to_order']:
            print(f'Restriction site(s) within internal region detected and replaced at row index {i}.')
    
    # save
    order_df.to_csv('./'+args.data.split('/')[-1].split('.')[0]+'.to_order.csv')

if __name__ == '__main__':
    main()

import sys
import pandas as pd


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
          'GGT':'GGC', 'GGC':'GGA', 'GGA':'GGG', 'GGG':'GGT'  # gly
         }

def restriction_replacement(seq, selector):
    sites = consts_df['restriction sites'][selector]
    
    # find sites
    # first restriction site
    i1 = 0
    if seq.find(sites[0], seq.find(sites[0])+1) != -1:
        i1 = seq.find(sites[0], seq.find(sites[0])+1)
        print(sites[0])
    else:
        pass

    # second restriction site
    i2 = 0
    if seq.find(sites[1], seq.find(sites[1])+1) != -1:
        i2 = seq.find(sites[1])
        print(sites[1])
    else:
        pass
    
    if i1 == 0 and i2 == 0:
        return seq

    # replace
    new_seq = seq
    for i in range((3-consts_df['reading frame'][selector]), len(seq),3):
        if i1 >= i and i1 <= i+2:
            #replace codon
            new_seq = new_seq[0:i] + codons[new_seq[i:i+3]] + new_seq[i+3:]
            print(len(seq), len(new_seq))
            # pass
        if i2 >= i and i2 <= i+2:
            #replace codon
            new_seq = new_seq[0:i] + codons[new_seq[i:i+3]] + new_seq[i+3:]
            # print(len(seq), len(new_seq))
            # pass
        if i1 >=i+2 and i2>=i+2:
            return new_seq

def main():
    #read data
    igblast_makedb = pd.read_csv('./data/04072023_YW.igblast.makedb.airr.tsv', delimiter='\t')
    var_df = pd.read_csv('./data/genes_v_call.csv').drop(columns=['Unnamed: 0'])

    # make order df
    order_df = igblast_makedb.loc[:,['sequence_id']]
    order_df['to_order'] = ['']*len(igblast_makedb)
    # print(order_df.head())
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
                    print('ERROR @', i)
                    order_df.loc[i, 'to_order'] = 'ERROR'
                    continue
            else:
                var = var_df[var_df['v_call']==curr_v_call]['var'].values
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
        order_df.loc[i, 'to_order'] = restriction_replacement(order_df.loc[i, 'to_order'], dict_selector)
        
        # if replaced, test for restriction site replacement again
        attempt_counter = 0
        if prev_seq != order_df.loc[i, 'to_order']:
            i -= 1
            attempt_counter += 1
            
            if attempt_counter == 6:
                print('Restriction site error at ' + i + '. Unable to replace duplicate site.')

    # save
    order_df.to_csv('./to_order.csv')


if __name__ == '__main__':
    main()
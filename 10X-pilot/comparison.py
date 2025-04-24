# Author : Sam Wollenburg
# Date : Dec 23, 2024

import pandas as pd
import numpy as np
import copy
import scanpy as sc

# update to add user defined paths as args for read_data()
def read_data():
    '''
    Read in the data necessary to use this. 
    Currently filepaths are hardcoded and will need to be changed.

    Returns
    -------
    vdj : pd.DataFrame()
        - The vdj information.
    gex : pd.DataFrame()
        - The gene expression information.
    '''
    # read in vdj info, make easily searchable
    vdj = pd.read_csv('../../projects/2024_12_23_10X11_14234_0/out/outs/per_sample_outs/2024_12_23_10X11_14234_0_multi/vdj_b/filtered_contig_annotations.csv')

    for i in range(len(vdj['v_gene'])):
        ig = vdj['v_gene'][i][2:]
        chain = ig[0]
        family = ig[2]
        # gene = ig.split('-')[1]
        
        vdj.loc[i, 'chain'] = chain
        vdj.loc[i, 'family'] = family
        # vdj.loc[i, 'gene'] = gene

    # read in mtx data from gex analysis
    mat = sc.read_10x_mtx('../../projects/2024_12_23_10X11_14234_0/out/outs/per_sample_outs/2024_12_23_10X11_14234_0_multi/count/sample_filtered_feature_bc_matrix')
    gex = mat.to_df()
    gex = gex.reset_index(names='cell_id')

    return vdj, gex

def check_pairs(group):
    '''
    rule for groupby's apply filter
    '''
    # has_heavy = 'H' in group[group['productive'] == True]['chain'].values 
    has_heavy = sum(group[group['productive'] == True]['chain'].values == 'H') == 1
    has_light = ('K' in group[group['productive'] == True]['chain'].values) != ('L' in group[group['productive'] == True]['chain'].values)
    return has_heavy and has_light

# create cell selection function
def select_cells(chain, family, vdj, gex):
    '''
    Selects the intersction between airr and gex when using chain and family from vdj

    Parameters
    ----------
    chain : str
        - A string representing which chain to select (H, L, K)
    family : str
        - A string (not an int) representing which family to select.
    vdj : pd.DataFrame()
        - The vdj dataframe read in
    gex : pd.DataFrame()
        - The gene expression dataframe read in

    Returns
    -------
    intersection : pd.DataFrame()
        - The intersection between vdj and gex for cells in vdj that were of the chain and family.
    '''
    # vdj_filtered = vdj[(vdj[vdj['chain']==chain]) & (vdj['family']==family)]
    vdj_chain = vdj[vdj['chain']==chain]
    vdj_fam = vdj_chain[vdj_chain['family'] == family]
    intersection = vdj_fam.merge(gex, how='inner', left_on='barcode', right_on='cell_id')
    return intersection

def select_paired_cells(chain, family, vdj, gex):
    '''
    Selects the intersction between airr and gex when using chain and family from vdj

    Parameters
    ----------
    chain : str
        - A string representing which chain to select (H, L, K)
    family : str
        - A string (not an int) representing which family to select.
    vdj : pd.DataFrame()
        - The vdj dataframe read in
    gex : pd.DataFrame()
        - The gene expression dataframe read in

    Returns
    -------
    intersection : pd.DataFrame()
        - The intersection between vdj and gex for cells in vdj that were of the chain and family.
    '''
    # vdj_filtered = vdj[(vdj[vdj['chain']==chain]) & (vdj['family']==family)]
    # Find pairs and extract barcodes
    paired_cells = vdj.groupby('barcode', group_keys=False).apply(check_pairs, include_groups=False)
    paired_barcodes = paired_cells[paired_cells].index.tolist()

    # Isolate cells with correct chain and family
    vdj_chain = vdj[vdj['chain']==chain]
    vdj_fam = vdj_chain[vdj_chain['family'] == family]

    # Find intersection of gex and vdj_fam, then select for paired barcodes.
    intersection = vdj_fam.merge(gex, how='inner', left_on='barcode', right_on='cell_id')
    intersection = intersection[intersection['barcode'].isin(paired_barcodes)]
    return intersection

def select_gene(df, gene):
    '''
    Provide a gene, receive tuple with negatives or positivies. 
    tuple[0] == false, tuple[1] == true

    Parameters
    ----------
    df : pd.DataFrame()
        - A dataframe of selected cells.
    gene : str
        - The gene to select on
    
    Returns
    -------
    - (neg, pos) : tuple
        - A tuple containing two dataframes. 
        - Index 0 contains cells that were negative for gene expression of 'gene'.
        - Index 1 contains cells that were positive for gene expression of 'gene'.
    '''
    neg = df[df[gene]==0]
    pos = df[df[gene]>0]
    
    return (neg, pos)

def compare(gene:str, vdj, gex, type1='IGHV3', type2='IGHV4'):
    '''
    Compare two selected groups of cells on a gene.  Prints table of result.

    Parameters
    ----------
    gene : str
        - The gene to compare
    vdj : pd.DataFrame()
        - VDJ information to read in
    gex : pd.DataFrame()
        - Gene expression data to read in
    type1 : str
        - The qualifier for cell clonotype
    type2 : str
        - The qualifier for cell clonotype
    '''
    print(f'Selecting {type1} cells...')
    cells1 = select_cells(type1[2], type1[4], vdj, gex)
    print(f'Selecting {type2} cells...')
    cells2 = select_cells(type2[2], type2[4], vdj, gex)

    print(f'Selecting for {gene} in {type1} cells...')
    selected1 = select_gene(cells1, gene)
    print(f'Selecting for {gene} in {type2} cells...')
    selected2 = select_gene(cells2, gene)

    print('='*26)
    print(f'{gene:<8} {"+":<8} {"-":<8}')
    print('-'*26)
    print(f'{type1:<8} {len(selected1[1]):<8} {len(selected1[0]):<8}')
    print(f'{type2:<8} {len(selected2[1]):<8} {len(selected2[0]):<8}')
    print('='*26+'\n')

def compare_selected(gene, cells1, cells2):
    '''
    Compare two selected groups of cells on a gene. Prints table of result.

    Parameters
    ----------
    gene : str
        - The gene to compare
    cells1 : pd.DataFrame()
        - The first group of cells
    cells2 : pd.DataFrame()
        - The second group of cells
    '''
    selected1 = print(f'Selecting for {gene}...')
    selected1 = select_gene(cells1, gene)
    selected2 = select_gene(cells2, gene)

    name1 = 'IG'+cells1['chain'][0]+'V'+cells1['family'][0]
    name2 = 'IG'+cells2['chain'][0]+'V'+cells2['family'][0]

    print('='*26)
    print(f'{gene:<8} {"+":<8} {"-":<8}')
    print('-'*26)
    print(f'{name1:<8} {len(selected1[1]):<8} {len(selected1[0]):<8}')
    print(f'{name2:<8} {len(selected2[1]):<8} {len(selected2[0]):<8}')
    print('='*26+'\n')

def create_clusters(group_name, genes, cells, location='local'):
    '''
    Creates clusters for loupe browser with list of genes and list of cells. 
    Save a csv file containing Loupe Browser clustering data.

    Parameters
    ----------
    group_name : str
        - The name of the group to be used in Loupe Browser
    genes : list
        - A list of the genes to compare
    cells : list
        - A list of groups of cells to compare
    location : list
        - Default is local and generates the file where ran in a ./cluster/ subdirectory. 
        - A file path location for saving the csv must end in a '/'
    '''
    if len(group_name.split('/'))>1:
        group_name_no_slash = '&'.join(group_name.split('/'))
    else:
        group_name_no_slash = group_name

    if location == 'local':
        location = './clusters/' + group_name_no_slash + '.csv'
    else:
        location = location + group_name_no_slash + '.csv'


    # setup df for output
    out = pd.concat(cells, ignore_index=True).loc[:,['cell_id', 'chain', 'family', 'v_gene']]
    out.rename(columns={'cell_id':'Barcode'}, inplace=True)
    out[group_name] = ['negative']*len(out['Barcode'])

    # compare different families and not different genes
    if len(genes) == 0 and len(cells) > 1:
        for i in range(len(out)):
            out.loc[i,group_name] = 'IG'+out.loc[i,'chain']+'V'+out.loc[i,'family']
    # show all genes in family
    elif len(genes) == 0 and len(cells) == 1:
        for i in range(len(out)):
            out.loc[i,group_name] = out.loc[i,'v_gene']
    # otherwise, compare between genes
    else:
        for gene in genes:
            gene_and_cells = [select_gene(cell, gene) for cell in cells]
            bothgenepos = pd.concat([e[1] for e in gene_and_cells], ignore_index=True)
            
            for i in range(len(out)):
                if out.loc[i,'Barcode'] in bothgenepos['cell_id'].values:
                    if out.loc[i,group_name] != 'negative':
                        out.loc[i,group_name] = out.loc[i,group_name] + '/' + gene
                    else:
                        out.loc[i,group_name] = gene
    out.drop(['chain', 'family', 'v_gene'], axis=1, inplace=True)

    # remove all duplicate barcodes and save
    out.drop_duplicates(subset='Barcode', keep=False, inplace=True, ignore_index=True)
    out.to_csv(location, index=False)

def find_families(vdj, chain='H'):
    '''
    Find all of the families in a chain.
    Return a list for create_clusters
    '''
    # find set of familes
    families = vdj[vdj['chain']==chain].loc[:,'family'].unique()
    return ['IG'+chain+'V'+family for family in families]

def main():
    vdj, gex = read_data()
    
    # ighv3 = select_cells('H', '3', vdj, gex)
    # ighv4 = select_cells('H', '4', vdj, gex)

    # create_clusters('IGHV3', [], [ighv3])
    # create_clusters('IGHV4', [], [ighv4])
    # create_clusters('IGHV3/4', [], [ighv3, ighv4])

    # all_IGHV = find_families(vdj, 'H')
    # cells_IGHV = [select_cells(ighv[2], ighv[4], vdj, gex) for ighv in all_IGHV]
    # create_clusters('IGHV', [], cells_IGHV)

    all_IGHV = find_families(vdj, 'H')
    cells_IGHV = [select_paired_cells(ighv[2], ighv[4], vdj, gex) for ighv in all_IGHV]
    create_clusters('IGHV_paired', [], cells_IGHV)

if __name__ == '__main__':
    main()

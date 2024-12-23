# Author : Sam Wollenburg
# Date : Dec 23, 2024

import pandas as pd
import numpy as np
import copy
import scanpy as sc

# update to add user defined paths as args for read_data()
def read_data():
    '''
    Read in the data necessary to use this. Currently paths are hardcoded and will need to be changed.

    Returns
    -------
    vdj : pd.DataFrame()
        The vdj information.
    gex : pd.DataFrame()
        The gene expression information.
    '''
    # read in vdj info, make easily searchable
    vdj = pd.read_csv('../../projects/2024_11_20_10X11_14234_0/out/outs/per_sample_outs/2024_11_20_10X11_14234_0_multi/vdj_b/filtered_contig_annotations.csv')

    for i in range(len(vdj['v_gene'])):
        ig = vdj['v_gene'][i][2:]
        chain = ig[0]
        family = ig[2]
        # gene = ig.split('-')[1]
        
        vdj.loc[i, 'chain'] = chain
        vdj.loc[i, 'family'] = family
        # vdj.loc[i, 'gene'] = gene

    # read in mtx data from gex analysis
    mat = sc.read_10x_mtx('../../projects/2024_11_20_10X11_14234_0/out/outs/per_sample_outs/2024_11_20_10X11_14234_0_multi/count/sample_filtered_feature_bc_matrix')
    gex = mat.to_df()
    gex = gex.reset_index(names='cell_id')

    return vdj, gex

# create cell selection function
def select_cells(chain, family, vdj, gex):
    '''
    Selects the intersction between airr and gex when using chain and family from vdj

    Parameters
    ----------
    chain : char
        A character representing which chain to select (H, L, K)
    family : char
        A character (not an int) representing which family to select.
    vdj : pd.DataFrame()
        The vdj dataframe read in
    gex : pd.DataFrame()
        The gene expression dataframe read in

    Returns
    -------
    intersenction : pd.DataFrame()
        The intersection between vdj and gex for cells in vdj that were of the chain and family.
    '''
    # vdj_filtered = vdj[(vdj[vdj['chain']==chain]) & (vdj['family']==family)]
    vdj_chain = vdj[vdj['chain']==chain]
    vdj_fam = vdj_chain[vdj_chain['family'] == family]
    intersection = vdj_fam.merge(gex, how='inner', left_on='barcode', right_on='cell_id')
    return intersection

def select_gene(df, gene):
    '''
    Provide a gene, receive tuple with negatives or positivies. 
    tuple[0] == false, tuple[1] == true

    Parameters
    ----------
    df : pd.DataFrame()
        A dataframe of selected cells.
    gene : str
        The gene to select on
    
    Returns
    -------
    tuple
        A tuple containing two dataframes. 
        Index 0 contains cells that were negative for gene expression of 'gene'.
        Index 1 contains cells that were positive for gene expression of 'gene'.
    '''
    neg = df[df[gene]==0]
    pos = df[df[gene]>0]

    return (neg, pos)

def compare(gene, vdj, gex, type1='IGHV3', type2='IGHV4'):
    '''
    Compare two selected groups of cells on a gene

    Parameters
    ----------
    gene : str
        The gene to compare
    vdj : pd.DataFrame()
        vdj inforation read in
    gex : pd.DataFrame()
        gene expression data read in
    type1 : str
        The qualifier for cell clonotype
    type2 : str
        The qualifier for cell clonotype

    Returns
    --------
    Printed table of result.
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
    Compare two selected groups of cells on a gene

    Parameters
    ----------
    gene : str
        The gene to compare
    cells1 : pd.DataFrame()
        The first group of cells
    cells2 : pd.DataFrame()
        The second group of cells

    Returns
    --------
    Printed table of result.
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
    Creates clusters for loupe browser with list of genes and list of cells

    Parameters
    ----------
    group_name : str
        The name of the group to be used in Loupe Browser
    genes : list
        A list of the genes to compare
    cells : list
        A list of groups of cells to compare
    location : list
        A file path location where the csv can be saved, must end in a '/'

    Returns
    -------
    A csv file containing Loupe Browser clustering data
    '''
    if len(group_name.split('/'))>1:
        group_name_no_slash = '&'.join(group_name.split('/'))
    if location == 'local':
        location = group_name_no_slash + '.csv'
    else:
        location = location + group_name_no_slash + '.csv'
    
    # genes = ['XBP1', 'LILRB1']

    out = pd.DataFrame()
    out['Barcode'] = pd.concat(cells, ignore_index=True)['cell_id']
    out[group_name] = ['negative']*len(out['Barcode'])
    out = out.drop_duplicates(ignore_index=True)
    # out.head()
    for gene in genes:
        gene_and_cells = [select_gene(cell, gene) for cell in cells]
        bothgenepos = pd.concat([e[1] for e in gene_and_cells], ignore_index=True)
        
        for i in range(len(out)):
            if out.loc[i,'Barcode'] in bothgenepos['cell_id'].values:
                if out.loc[i,group_name] != 'negative':
                    out.loc[i,group_name] = out.loc[i,group_name] + '/' + gene
                else:
                    out.loc[i,group_name] = gene

    out.to_csv(location, index=False)
import pandas as pd
import numpy as np
import copy
import scanpy as sc
from comparison import *

# def check_pairs(group): 
    # has_heavy = 'H' in group[group['productive'] == True]['chain'].values 
    # has_heavy = sum(group[group['productive'] == True]['chain'].values == 'H') == 1
    # has_light = any(chain in group[group['productive'] == True]['chain'].values for chain in ['K','L'])
    # has_light = ('K' in group['chain'].values) != ('L' in group['chain'].values)

    # has_heavy = sum(group[group['productive'] == True]['chain'].values == 'H') == 1
    # has_light = ('K' in group[group['productive'] == True]['chain'].values) != ('L' in group[group['productive'] == True]['chain'].values)

    # return has_heavy and has_light

def main():
    print('Reading data...')
    vdj, gex = read_data()

    # Total pairs
    print('Counting total pairs...')
    print('Total pairs:', vdj.groupby('barcode', group_keys=False).apply(check_pairs, include_groups=False).sum())
    print()

    # TCL1A
    print('Selecting TCL1A...')
    tcl1a = select_gene(gex, 'TCL1A')
    print('Merging TCL1A...')
    tcl1a_neg = vdj.merge(tcl1a[0], how='inner', left_on='barcode', right_on='cell_id')
    tcl1a_pos = vdj.merge(tcl1a[1], how='inner', left_on='barcode', right_on='cell_id')
    print('Counting TCL1A...')
    print('TCL1A- pairs:', tcl1a_neg.groupby('barcode', group_keys=False).apply(check_pairs, include_groups=False).sum())
    print('TCL1A+ pairs:', tcl1a_pos.groupby('barcode', group_keys=False).apply(check_pairs, include_groups=False).sum())
    print()

    # TCL1A-; IRF4
    print('Selecting TCL1A-; IRF4...')
    irf4 = select_gene(tcl1a[0], 'IRF4')
    print('Merging IRF4...')
    irf4_neg = vdj.merge(irf4[0], how='inner', left_on='barcode', right_on='cell_id')
    irf4_pos = vdj.merge(irf4[1], how='inner', left_on='barcode', right_on='cell_id')
    print('Counting IRF4...')
    print('IRF4- pairs:', irf4_neg.groupby('barcode', group_keys=False).apply(check_pairs, include_groups=False).sum())
    print('IRF4+ pairs:', irf4_pos.groupby('barcode', group_keys=False).apply(check_pairs, include_groups=False).sum())
    print()

    # TCL1A-; XBP1
    print('Selecting TCL1A-; XBP1...')
    xbp1 = select_gene(tcl1a[0], 'XBP1')
    print('Merging XBP1...')
    xbp1_neg = vdj.merge(xbp1[0], how='inner', left_on='barcode', right_on='cell_id')
    xbp1_pos = vdj.merge(xbp1[1], how='inner', left_on='barcode', right_on='cell_id')
    print('Counting XBP1...')
    print('XBP1- pairs:', xbp1_neg.groupby('barcode', group_keys=False).apply(check_pairs, include_groups=False).sum())
    print('XBP1+ pairs:', xbp1_pos.groupby('barcode', group_keys=False).apply(check_pairs, include_groups=False).sum())
    print()

    # TCL1A-; XBP1 IRF4
    print('Selecting TCL1A-; IRF4/XBP1...')
    # xbp1 = select_gene(tcl1a[0], 'XBP1')
    # irf4 = select_gene(tcl1a[0], 'IRF4')
    print('Merging IRF4/XBP1...')
    # irf4_neg = vdj.merge(irf4[0], how='inner', left_on='barcode', right_on='cell_id')
    # irf4_pos = vdj.merge(irf4[1], how='inner', left_on='barcode', right_on='cell_id')
    # xbp1_neg = vdj.merge(xbp1[0], how='inner', left_on='barcode', right_on='cell_id')
    # xbp1_pos = vdj.merge(xbp1[1], how='inner', left_on='barcode', right_on='cell_id')

    # Inclusive-or for pos (anything pos works)
    irf4xbp1_pos = pd.concat([irf4_pos, xbp1_pos]).drop_duplicates()
    # Both must be negative. 1) Take inclusive or of negs. Contains some pos. 2) Combine with pos and drop duplicates, keep none.
    irf4xbp1_neg = pd.concat([irf4_neg, xbp1_neg]).drop_duplicates()
    irf4xbp1_neg = pd.concat([irf4xbp1_neg, irf4xbp1_pos, irf4xbp1_pos]).drop_duplicates(keep=False)

    print('Counting IRF4/XBP1...')
    print('IRF4-XBP1- pairs:', irf4xbp1_neg.groupby('barcode', group_keys=False).apply(check_pairs, include_groups=False).sum())
    print('IRF4+/XBP1+ pairs:', irf4xbp1_pos.groupby('barcode', group_keys=False).apply(check_pairs, include_groups=False).sum())
    print()

if __name__ == '__main__':
    main()

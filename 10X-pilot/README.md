# README

Author : Sam Wollenburg

Date : Dec 23, 2024

### Subdirectories

### Files

- `comparison.ipynb` : Contains 10X analysis. Selecting for genes (positive or negative), and select for IG chain and gene family. Can be used, primarily for development of `comparison.py`. Currently there are issues of jupyter kernel crashing. Likely a result of current hardware with limited memory.
  - Must create a chain&family df before `select_gene()`
  - Currently, jupyter kernel is crashing when runing `compare()`. I think this is likely  because of the limited system memory  on my machine but this has not been thouroughly tested. The issue seems to arise from the `select_cells()` function where cells are selected to be from a specific chain and family. To avoid, use `select_cells()` once to get each group, then use `compare_selected()` to see if the groups are positive or negative for the selected gene.
- `comparison.py` : Contains 10X analysis that is recommended to use. Recommendation on use:
  - Read in the data: `vdj, gex = read_data()` Currently, this uses the `filtered_contig_annotations.csv` file for the vdj information and `sample_filtered_feautre_bc_matrix` for the gex information.
  - Select a group of cells:  `ighv3 = select_cells('H', '3', vdj, gex)` As seen, this selects for vdj cells with the H chain and are from family 3. Repeat this to create multiple groups. For example `ighv4 = select_cells('H', '4', vdj, gex)`
  - Create a list of genes to analize: `genes = ['CD19', 'CD27', 'CD38', 'XBP1', 'LILRB1']`
  - Create the clusters to view in loupe browser: `create_clusters('MBC/PB in IGHV3/4', genes, [ighv3, ighv4])` This will generate a group of clusters that can be read and viewed by loupe browser. This will highlight cells that are positive with any of the genes listed. The cells will be seperated by their unique combination of markers present. The group name `MBC/PB in IGHV3/4` will  be changed to `MBC&PB in IGHV3&4` as the forward slashes interfere with the filepath.
- `README.md` : This file.

### Change notes

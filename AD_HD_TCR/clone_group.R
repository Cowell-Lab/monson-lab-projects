# diversity calculations for ACS/HC clone groups

library(alakazam)
library(stringr)

acs_cd4 = read.table('ACS_CD4_cdr3_vj_aa_sharing.tsv', sep='\t')
acs_cd8 = read.table('ACS_CD8_cdr3_vj_aa_sharing.tsv', sep='\t')

hc_cd4 = read.table('HC_cdr3_vj_aa_sharing.tsv', sep='\t')
hc_cd8 = read.table('HC_CD8_cdr3_vj_aa_sharing.tsv', sep='\t')

acs_cd4.tot = sum(acs_cd4$V3)
acs_cd4.len = dim(acs_cd4)[1]
acs_cd8.tot = sum(acs_cd8$V3)
acs_cd8.len = dim(acs_cd8)[1]
hc_cd4.tot = sum(hc_cd4$V3)
hc_cd4.len = dim(hc_cd4)[1]
hc_cd8.tot = sum(hc_cd8$V3)
hc_cd8.len = dim(hc_cd8)[1]

acs_cd4['pct'] = acs_cd4['V3'] / acs_cd4.tot
acs_cd8['pct'] = acs_cd8['V3'] / acs_cd8.tot
hc_cd4['pct'] = hc_cd4['V3'] / hc_cd4.tot
hc_cd8['pct'] = hc_cd8['V3'] / hc_cd8.tot

acs_cd4['clone_id'] = rep(1:acs_cd4.len)
acs_cd8['clone_id'] = rep(1:acs_cd8.len)
hc_cd4['clone_id'] = rep(1:hc_cd4.len)
hc_cd8['clone_id'] = rep(1:hc_cd8.len)

# ACS CD4
acs_cd4.black = subset(acs_cd4, pct >= 0.01)
acs_cd4.grey = subset(acs_cd4, pct < 0.01 & pct >= 0.0001)
acs_cd4.white = subset(acs_cd4, pct < 0.0001)
acs_cd4.black.q = calcDiversity(acs_cd4.black$V3, 1)
acs_cd4.grey.q = calcDiversity(acs_cd4.grey$V3, 1)
acs_cd4.white.q = calcDiversity(acs_cd4.white$V3, 1)

acs_cd4.black.len = dim(acs_cd4.black)[1]
for (i in 1:acs_cd4.black.len) { acs_cd4.black[i, 'v_call'] = unlist(str_split(acs_cd4.black[i, 'V2'], '\\|'))[1] }
acs_cd4.black.gene = countGenes(acs_cd4.black, gene="v_call", mode="gene",copy='V3')
acs_cd4.black.family = countGenes(acs_cd4.black, gene="v_call", mode="family",copy='V3')

acs_cd4.grey.len = dim(acs_cd4.grey)[1]
for (i in 1:acs_cd4.grey.len) { acs_cd4.grey[i, 'v_call'] = unlist(str_split(acs_cd4.grey[i, 'V2'], '\\|'))[1] }
acs_cd4.grey.gene = countGenes(acs_cd4.grey, gene="v_call", mode="gene",copy='V3')
acs_cd4.grey.family = countGenes(acs_cd4.grey, gene="v_call", mode="family",copy='V3')

acs_cd4.white.len = dim(acs_cd4.white)[1]
for (i in 1:acs_cd4.white.len) { acs_cd4.white[i, 'v_call'] = unlist(str_split(acs_cd4.white[i, 'V2'], '\\|'))[1] }
acs_cd4.white.gene = countGenes(acs_cd4.white, gene="v_call", mode="gene",copy='V3')
acs_cd4.white.family = countGenes(acs_cd4.white, gene="v_call", mode="family",copy='V3')

# ACS CD8
acs_cd8.black = subset(acs_cd8, pct >= 0.01)
acs_cd8.grey = subset(acs_cd8, pct < 0.01 & pct >= 0.0001)
acs_cd8.white = subset(acs_cd8, pct < 0.0001)
acs_cd8.black.q = calcDiversity(acs_cd8.black$V3, 1)
acs_cd8.grey.q = calcDiversity(acs_cd8.grey$V3, 1)
acs_cd8.white.q = calcDiversity(acs_cd8.white$V3, 1)

acs_cd8.black.len = dim(acs_cd8.black)[1]
for (i in 1:acs_cd8.black.len) { acs_cd8.black[i, 'v_call'] = unlist(str_split(acs_cd8.black[i, 'V2'], '\\|'))[1] }
#acs_cd8.black.gene = countGenes(acs_cd8.black, gene="v_call", mode="gene",copy='V3')
#acs_cd8.black.family = countGenes(acs_cd8.black, gene="v_call", mode="family",copy='V3')

acs_cd8.grey.len = dim(acs_cd8.grey)[1]
for (i in 1:acs_cd8.grey.len) { acs_cd8.grey[i, 'v_call'] = unlist(str_split(acs_cd8.grey[i, 'V2'], '\\|'))[1] }
acs_cd8.grey.gene = countGenes(acs_cd8.grey, gene="v_call", mode="gene",copy='V3')
acs_cd8.grey.family = countGenes(acs_cd8.grey, gene="v_call", mode="family",copy='V3')

acs_cd8.white.len = dim(acs_cd8.white)[1]
for (i in 1:acs_cd8.white.len) { acs_cd8.white[i, 'v_call'] = unlist(str_split(acs_cd8.white[i, 'V2'], '\\|'))[1] }
acs_cd8.white.gene = countGenes(acs_cd8.white, gene="v_call", mode="gene",copy='V3')
acs_cd8.white.family = countGenes(acs_cd8.white, gene="v_call", mode="family",copy='V3')

# HC CD4
hc_cd4.black = subset(hc_cd4, pct >= 0.01)
hc_cd4.grey = subset(hc_cd4, pct < 0.01 & pct >= 0.0001)
hc_cd4.white = subset(hc_cd4, pct < 0.0001)
hc_cd4.black.q = calcDiversity(hc_cd4.black$V3, 1)
hc_cd4.grey.q = calcDiversity(hc_cd4.grey$V3, 1)
hc_cd4.white.q = calcDiversity(hc_cd4.white$V3, 1)

hc_cd4.black.len = dim(hc_cd4.black)[1]
for (i in 1:hc_cd4.black.len) { hc_cd4.black[i, 'v_call'] = unlist(str_split(hc_cd4.black[i, 'V2'], '\\|'))[1] }
hc_cd4.black.gene = countGenes(hc_cd4.black, gene="v_call", mode="gene",copy='V3')
hc_cd4.black.family = countGenes(hc_cd4.black, gene="v_call", mode="family",copy='V3')

hc_cd4.grey.len = dim(hc_cd4.grey)[1]
for (i in 1:hc_cd4.grey.len) { hc_cd4.grey[i, 'v_call'] = unlist(str_split(hc_cd4.grey[i, 'V2'], '\\|'))[1] }
hc_cd4.grey.gene = countGenes(hc_cd4.grey, gene="v_call", mode="gene",copy='V3')
hc_cd4.grey.family = countGenes(hc_cd4.grey, gene="v_call", mode="family",copy='V3')

hc_cd4.white.len = dim(hc_cd4.white)[1]
for (i in 1:hc_cd4.white.len) { hc_cd4.white[i, 'v_call'] = unlist(str_split(hc_cd4.white[i, 'V2'], '\\|'))[1] }
hc_cd4.white.gene = countGenes(hc_cd4.white, gene="v_call", mode="gene",copy='V3')
hc_cd4.white.family = countGenes(hc_cd4.white, gene="v_call", mode="family",copy='V3')

# HC CD8
hc_cd8.black = subset(hc_cd8, pct >= 0.01)
hc_cd8.grey = subset(hc_cd8, pct < 0.01 & pct >= 0.0001)
hc_cd8.white = subset(hc_cd8, pct < 0.0001)
hc_cd8.black.q = calcDiversity(hc_cd8.black$V3, 1)
hc_cd8.grey.q = calcDiversity(hc_cd8.grey$V3, 1)
hc_cd8.white.q = calcDiversity(hc_cd8.white$V3, 1)

hc_cd8.black.len = dim(hc_cd8.black)[1]
for (i in 1:hc_cd8.black.len) { hc_cd8.black[i, 'v_call'] = unlist(str_split(hc_cd8.black[i, 'V2'], '\\|'))[1] }
hc_cd8.black.gene = countGenes(hc_cd8.black, gene="v_call", mode="gene",copy='V3')
hc_cd8.black.family = countGenes(hc_cd8.black, gene="v_call", mode="family",copy='V3')

hc_cd8.grey.len = dim(hc_cd8.grey)[1]
for (i in 1:hc_cd8.grey.len) { hc_cd8.grey[i, 'v_call'] = unlist(str_split(hc_cd8.grey[i, 'V2'], '\\|'))[1] }
hc_cd8.grey.gene = countGenes(hc_cd8.grey, gene="v_call", mode="gene",copy='V3')
hc_cd8.grey.family = countGenes(hc_cd8.grey, gene="v_call", mode="family",copy='V3')

hc_cd8.white.len = dim(hc_cd8.white)[1]
for (i in 1:hc_cd8.white.len) { hc_cd8.white[i, 'v_call'] = unlist(str_split(hc_cd8.white[i, 'V2'], '\\|'))[1] }
hc_cd8.white.gene = countGenes(hc_cd8.white, gene="v_call", mode="gene",copy='V3')
hc_cd8.white.family = countGenes(hc_cd8.white, gene="v_call", mode="family",copy='V3')

# all trbv families
TRBV = c("TRBV1", "TRBV2", "TRBV3", "TRBV4", "TRBV5", "TRBV6", "TRBV7", "TRBV9", "TRBV10", "TRBV11", "TRBV12", "TRBV13", "TRBV14", "TRBV15", "TRBV16", "TRBV17", "TRBV18", "TRBV19", "TRBV20", "TRBV21", "TRBV23", "TRBV24", "TRBV25", "TRBV26", "TRBV27", "TRBV28", "TRBV29", "TRBV30")
#exclude.TRBV = c("TRBV1", "TRBV3", "TRBV9", "TRBV16", "TRBV17", "TRBV23", "TRBV26")
exclude.TRBV = c()
use.TRBV = TRBV[! TRBV %in% exclude.TRBV]

usage = array(0, c(12,length(use.TRBV)))
colnames(usage) = use.TRBV

genes <- match(acs_cd4.black.family$gene, use.TRBV)
usage[1,genes] = acs_cd4.black.family$copy_freq
genes <- match(acs_cd4.grey.family$gene, use.TRBV)
usage[2,genes] = acs_cd4.grey.family$copy_freq
genes <- match(acs_cd4.white.family$gene, use.TRBV)
usage[3,genes] = acs_cd4.white.family$copy_freq

#genes <- match(acs_cd8.black.family$gene, use.TRBV)
#usage[4,genes] = acs_cd8.black.family$copy_freq
genes <- match(acs_cd8.grey.family$gene, use.TRBV)
usage[5,genes] = acs_cd8.grey.family$copy_freq
genes <- match(acs_cd8.white.family$gene, use.TRBV)
usage[6,genes] = acs_cd8.white.family$copy_freq

genes <- match(hc_cd4.black.family$gene, use.TRBV)
usage[7,genes] = hc_cd4.black.family$copy_freq
genes <- match(hc_cd4.grey.family$gene, use.TRBV)
usage[8,genes] = hc_cd4.grey.family$copy_freq
genes <- match(hc_cd4.white.family$gene, use.TRBV)
usage[9,genes] = hc_cd4.white.family$copy_freq

genes <- match(hc_cd8.black.family$gene, use.TRBV)
usage[10,genes] = hc_cd8.black.family$copy_freq
genes <- match(hc_cd8.grey.family$gene, use.TRBV)
usage[11,genes] = hc_cd8.grey.family$copy_freq
genes <- match(hc_cd8.white.family$gene, use.TRBV)
usage[12,genes] = hc_cd8.white.family$copy_freq

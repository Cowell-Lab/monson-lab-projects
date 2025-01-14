library(gplots)

# B cell library 4/6

# gene usage heatmap

data_dir = '/Users/s166813/Projects/data/immune/MonsonLab/TM-HC-RIS/f4ff01ac-e5ba-4a1e-a1be-5397d995e089-007/'

ighv.seg = read.table(paste(data_dir,'v_gene_usage.csv',sep=''), header=T, sep=',')
rownames(ighv.seg) = t(ighv.seg["sample_id"])

# ighv genes
IGHV = sort(names(ighv.seg))[1:91]
genes <- match(IGHV,names(ighv.seg))
ighv.rel = ighv.seg[,genes] / rowSums(ighv.seg[,genes])

# exclude when below threshold
ighv.avg = colMeans(ighv.rel)
exclude.IGHV = c()
for (n in names(ighv.rel)) {
    if (ighv.avg[n] < 0.01) {
        exclude.IGHV = c(exclude.IGHV, n)
    }
}
use.IGHV = IGHV[! IGHV %in% exclude.IGHV]
genes <- match(use.IGHV,names(ighv.rel))

# gDNA (lib6) vs cDNA (lib4)
# by subject, by cell subset

# TCR library 5

# gene usage heatmap

samples <- match(sort(rownames(ighv.rel)), rownames(ighv.rel))

pdf("gene_heatmap.pdf", width=10)
#heatmap.2(as.matrix(ighv.rel[,genes]), scale="row", dendrogram=c("row"), col = cm.colors(128), trace=c("none"), density.info=c("histogram"), margins=c(7,7), Colv=NA, keysize=2)
heatmap(as.matrix(ighv.rel[samples,genes]), Rowv=NA, Colv=NA, col = cm.colors(128), scale="row", margins=c(5,10))
dev.off()

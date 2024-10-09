library(gplots)

# TCR library 5

# gene usage heatmap

data_dir = '/Users/s166813/Projects/data/immune/MonsonLab/T_cells_pancreatic/7c83e049-bb46-4fe1-8877-93fbdecd1af1-007/segment_counts_data/'

trbv.seg = read.table(paste(data_dir,'trbv_gene.csv',sep=''), header=T, sep=',')

# AD vs HD

AD = c(14, 13, 3, 4, 1, 2, 15, 29, 28, 21, 26, 27)
AD.names = c("AD1","AD4","AD5","AD6","AD7","AD9","AD10","AD11","AD12","AD13","AD14","AD15")
HD = c(7, 8, 5, 6, 9)
HD.names = c("HD2","HD3","HD4","HD6","HD8")

# all trbv genes
TRBV = c("TRBV1", "TRBV2", "TRBV3.1", "TRBV3.2", "TRBV4.1", "TRBV4.2", "TRBV4.3", "TRBV5.1", "TRBV5.3", "TRBV5.4", "TRBV5.5", "TRBV5.6", "TRBV5.7", "TRBV5.8", "TRBV6.1", "TRBV6.2", "TRBV6.3", "TRBV6.4", "TRBV6.5", "TRBV6.6", "TRBV6.7", "TRBV6.8", "TRBV6.9", "TRBV7.1", "TRBV7.2", "TRBV7.3", "TRBV7.4", "TRBV7.6", "TRBV7.7", "TRBV7.8", "TRBV7.9", "TRBV9", "TRBV10.1", "TRBV10.2", "TRBV10.3", "TRBV11.1", "TRBV11.2", "TRBV11.3", "TRBV12.1", "TRBV12.2", "TRBV12.3", "TRBV12.4", "TRBV12.5", "TRBV13", "TRBV14", "TRBV15", "TRBV16", "TRBV17", "TRBV18", "TRBV19", "TRBV20.1", "TRBV20.OR9.2", "TRBV21.1", "TRBV21.OR9.2", "TRBV23.1", "TRBV23.OR9.2", "TRBV24.1", "TRBV24.OR9.2", "TRBV25.1", "TRBV26", "TRBV26.OR9.2", "TRBV27", "TRBV28", "TRBV29.1", "TRBV29.OR9.2", "TRBV30")

exclude.TRBV = c("file", "clones", "TRBV1", "TRBV7.1", "TRBV6.3", "TRBV17", "TRBV26.OR9.2", "TRBV26", "TRBV29.OR9.2", "TRBV23.OR9.2", "TRBV21.OR9.2", "TRBV12.1", "TRBV12.2")
use.TRBV = TRBV[! TRBV %in% exclude.TRBV]

genes = match(use.TRBV,names(trbv.seg))

AD.counts = as.matrix(trbv.seg[AD,genes])
AD.rel = AD.counts / rowSums(AD.counts)
rownames(AD.rel) = AD.names
HD.counts = as.matrix(trbv.seg[HD,genes])
HD.rel = HD.counts / rowSums(HD.counts)
rownames(HD.rel) = HD.names

both.rel = rbind(AD.rel, HD.rel)

#jpeg("AD_HD_gene_heatmap.jpeg", res=300, width=10, height=6,units="in")
pdf("AD_HD_gene_heatmap.pdf", width=10)
heatmap.2(both.rel, scale="column", dendrogram=c("row"), col = cm.colors(128), trace=c("none"), density.info=c("histogram"), margins=c(7,7), Colv=NA, keysize=2)
dev.off()

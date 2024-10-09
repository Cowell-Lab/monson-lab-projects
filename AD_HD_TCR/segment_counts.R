data_dir = '/Users/s166813/Projects/data/immune/MonsonLab/T_cells/7c947479-f938-4489-a53b-1792db0f3d82-007/segment_counts_data/'

trbv.seg = read.table(paste(data_dir,'trbv_gene.csv',sep=''), header=T, sep=',')
trav.seg = read.table(paste(data_dir,'trav_gene.csv',sep=''), header=T, sep=',')

# all trbv families
#genes = c(11, 10, 9, 8, 7, 6, 5, 4, 19, 18, 17, 16, 15, 14, 13, 12, 21, 20, 30, 31, 29, 27, 28, 25, 26, 23, 24, 22)
genes = c(10, 9, 8, 7, 6, 5, 4, 19, 18, 17, 15, 14, 21, 20, 30, 31, 27, 28, 26, 23, 24)

seg.pbmc = trbv.seg[c(1,7),genes]
seg.pbmc.avg = colMeans(seg.pbmc)
seg.pbmc.sd = apply(seg.pbmc,2,sd)
n = dim(seg.pbmc)[1]
seg.pbmc.se = seg.pbmc.sd / sqrt(n)

seg.cd4 = trbv.seg[c(3,5),genes]
seg.cd4.avg = colMeans(seg.cd4)
seg.cd4.sd = apply(seg.cd4,2,sd)
n = dim(seg.cd4)[1]
seg.cd4.se = seg.cd4.sd / sqrt(n)

seg.cd8 = trbv.seg[c(4,6),genes]
seg.cd8.avg = colMeans(seg.cd8)
seg.cd8.sd = apply(seg.cd8,2,sd)
n = dim(seg.cd8)[1]
seg.cd8.se = seg.cd8.sd / sqrt(n)

seg.naive = trbv.seg[c(8,9,10),genes]
seg.naive.avg = colMeans(seg.naive)
seg.naive.sd = apply(seg.naive,2,sd)
n = dim(seg.naive)[1]
seg.naive.se = seg.naive.sd / sqrt(n)

seg.avg = rbind(seg.pbmc.avg, seg.cd4.avg, seg.cd8.avg, seg.naive.avg)
seg.se = rbind(seg.pbmc.se, seg.cd4.se, seg.cd8.se, seg.naive.se)

#pdf(file='trbv_segment_counts.pdf', width=20, height=5)
barcenters = barplot(seg.avg, beside=T, col = c("lightblue", "mistyrose", "lightcyan", "lavender"), ylim=c(0,0.2), legend.text=c('PBMC', 'CD4', 'CD8', 'Naive'))

segments(barcenters, seg.avg - seg.se * 2, barcenters, seg.avg + seg.se * 2, lwd = 1.5)
#dev.off()

# all trav families
#genes = c(13, 14, 15, 16, 17, 18, 19, 20, 21, 8, 9, 6, 7, 12, 10, 11, 4, 5, 25, 24, 23, 22, 29, 28, 27, 26, 30, 36, 33, 34, 35, 37, 38, 32, 31)
genes = c(13, 14, 16, 17, 18, 20, 21, 8, 6, 7, 12, 10, 11, 5, 25, 24, 23, 22, 27, 26, 30, 34, 35, 37, 31)

seg.pbmc = trav.seg[c(11,17),genes]
seg.pbmc.avg = colMeans(seg.pbmc)
seg.pbmc.sd = apply(seg.pbmc,2,sd)
n = dim(seg.pbmc)[1]
seg.pbmc.se = seg.pbmc.sd / sqrt(n)

seg.cd4 = trav.seg[c(13,15),genes]
seg.cd4.avg = colMeans(seg.cd4)
seg.cd4.sd = apply(seg.cd4,2,sd)
n = dim(seg.cd4)[1]
seg.cd4.se = seg.cd4.sd / sqrt(n)

seg.cd8 = trav.seg[c(14,16),genes]
seg.cd8.avg = colMeans(seg.cd8)
seg.cd8.sd = apply(seg.cd8,2,sd)
n = dim(seg.cd8)[1]
seg.cd8.se = seg.cd8.sd / sqrt(n)

seg.naive = trav.seg[c(18,19,20),genes]
seg.naive.avg = colMeans(seg.naive)
seg.naive.sd = apply(seg.naive,2,sd)
n = dim(seg.naive)[1]
seg.naive.se = seg.naive.sd / sqrt(n)

seg.avg = rbind(seg.pbmc.avg, seg.cd4.avg, seg.cd8.avg, seg.naive.avg)
seg.se = rbind(seg.pbmc.se, seg.cd4.se, seg.cd8.se, seg.naive.se)

pdf(file='trav_segment_counts.pdf', width=20, height=5)
barcenters = barplot(seg.avg, beside=T, col = c("lightblue", "mistyrose", "lightcyan", "lavender"), ylim=c(0,0.16), legend.text=c('PBMC', 'CD4', 'CD8', 'Naive'))

segments(barcenters, seg.avg - seg.se * 2, barcenters, seg.avg + seg.se * 2, lwd = 1.5)
dev.off()

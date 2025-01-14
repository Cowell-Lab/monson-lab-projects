# TM/HC/RIS
# libraries 1, 4, 6, 8, 9

# gene family counts

library(tidyr)

data_dir = '/Users/s166813/Projects/data/immune/MonsonLab/TM-HC-RIS/vdjserver/'
lib_dir = paste(data_dir, 'library_all/stats/', sep='')
#lib_dir = paste(data_dir, 'library_all/0shm/', sep='')
#lib_dir = paste(data_dir, 'library_all/gt3shm/', sep='')

# load raw tables
#filename = 'ighv_family_TM-HC-RIS_PB.pdf'
#liball.HC = read.table(paste(lib_dir, 'HC_PB.group.v_call.tsv',sep=''), header=T, sep='\t')
#liball.RIS = read.table(paste(lib_dir, 'RIS_PB.group.v_call.tsv',sep=''), header=T, sep='\t')
#liball.TM = read.table(paste(lib_dir, 'TM_PB.group.v_call.tsv',sep=''), header=T, sep='\t')
#filename = 'ighv_family_TM-HC-RIS_PB_DNA.pdf'
#filename = 'ighv_family_TM-HC-RIS_PB_DNA_0shm.pdf'
#filename = 'ighv_family_TM-HC-RIS_PB_DNA_gt3shm.pdf'
filename = 'seq_Fig2_ighv_family_TM-HC-RIS_PB_DNA.pdf'
#filename = 'seq_Fig2_ighv_family_TM-HC-RIS_PB_DNA_0shm.pdf'
#filename = 'seq_Fig2_ighv_family_TM-HC-RIS_PB_DNA_gt3shm.pdf'
liball.HC = read.table(paste(lib_dir, 'HC_PB_DNA.group.v_call.tsv',sep=''), header=T, sep='\t')
liball.RIS = read.table(paste(lib_dir, 'RIS_PB_DNA.group.v_call.tsv',sep=''), header=T, sep='\t')
liball.TM = read.table(paste(lib_dir, 'TM_PB_DNA.group.v_call.tsv',sep=''), header=T, sep='\t')

# extract the data
#level = 'subgroup'
level = 'subgroup'
mode = 'proportion'
productive = 'TRUE'
ighv.HC = liball.HC[liball.HC$level == level & liball.HC$mode == mode & liball.HC$productive == productive,]
ighv.RIS = liball.RIS[liball.RIS$level == level & liball.RIS$mode == mode & liball.RIS$productive == productive,]
ighv.TM = liball.TM[liball.TM$level == level & liball.TM$mode == mode & liball.TM$productive == productive,]
ighv.all = rbind(ighv.HC, ighv.RIS, ighv.TM)

# spread genes into columns
#ighv.all.seg = ighv.all[,c("repertoire_group_id","gene","duplicate_frequency_avg")] %>% spread(gene, duplicate_frequency_avg, fill=0.0)
#ighv.all.sd = ighv.all[,c("repertoire_group_id","gene","duplicate_frequency_std")] %>% spread(gene, duplicate_frequency_std, fill=0.0)
ighv.all.seg = ighv.all[,c("repertoire_group_id","gene","sequence_frequency_avg")] %>% spread(gene, sequence_frequency_avg, fill=0.0)
ighv.all.sd = ighv.all[,c("repertoire_group_id","gene","sequence_frequency_std")] %>% spread(gene, sequence_frequency_std, fill=0.0)
ighv.all.N = ighv.all[,c("repertoire_group_id","gene","N")] %>% spread(gene, N, fill=0.0)
ighv.all.len = length(ighv.all.seg)

# genes and groups
ighv.all.genes = names(ighv.all.seg[,2:ighv.all.len])
ighv.all.groups = as.vector(ighv.all.seg[,"repertoire_group_id"])

# exclude when below threshold
ighv.all.avg = colMeans(ighv.all.seg[,ighv.all.genes])
exclude.IGHV = c()
for (n in ighv.all.genes) {
    if (ighv.all.avg[n] < 0.01) {
        exclude.IGHV = c(exclude.IGHV, n)
    }
}
ighv.use.genes = ighv.all.genes[! ighv.all.genes %in% exclude.IGHV]
#ighv.use.genes = ighv.all.genes

ighv.all.data = as.matrix(ighv.all.seg[,ighv.use.genes])
ighv.all.se = as.matrix(ighv.all.sd[,ighv.use.genes] / sqrt(ighv.all.N[,ighv.use.genes]))

pdf(file=filename, width=7, height=7)
barcenters = barplot(ighv.all.data, beside=T, ylim=c(0,0.6), col = c("lightblue", "mistyrose", "chocolate"), legend.text=ighv.all.groups)
segments(barcenters, ighv.all.data - ighv.all.se * 2, barcenters, ighv.all.data + ighv.all.se * 2, lwd = 1.5)
dev.off()

# old stuff
if (FALSE) {
genes <- match(use.IGHV,names(ighv.seg))
genes9 <- match(use.IGHV,names(ighv9.seg))

hc.idx = which(ighv.seg[,"diagnosis"] == "HC" & ighv.seg[,"cell.type"] == "PB" )
tm.idx = which(ighv.seg[,"diagnosis"] == "TM" & ighv.seg[,"cell.type"] == "PB" )
ris.idx = which(ighv.seg[,"diagnosis"] == "RIS" & ighv.seg[,"cell.type"] == "PB" )

hc9.idx = which(ighv9.seg[,"diagnosis"] =="HC" & ighv9.seg[,"Cell.type"] == "PB" )
hc9_dna.idx = which((ighv9.seg[,"diagnosis"] == "HC") & (ighv9.seg[,"library.source"] == "DNA" ) & ighv9.seg[,"Cell.type"] == "PB")
hc9_rna.idx = which(ighv9.seg[,"diagnosis"] == "HC" & ighv9.seg[,"library.source"] == "RNA" & ighv9.seg[,"Cell.type"] == "PB")
tm9.idx = which(ighv9.seg[,"diagnosis"] =="TM" & ighv9.seg[,"Cell.type"] == "PB")
ris9.idx = which(ighv9.seg[,"diagnosis"] =="RIS" & ighv9.seg[,"Cell.type"] == "PB")
ris9_dna.idx = which(ighv9.seg[,"diagnosis"] == "RIS" & ighv9.seg[,"library.source"] == "DNA" & ighv9.seg[,"Cell.type"] == "PB")
ris9_rna.idx = which(ighv9.seg[,"diagnosis"] == "RIS" & ighv9.seg[,"library.source"] == "RNA" & ighv9.seg[,"Cell.type"] == "PB")


# TM/RIS/HC comparison
# by cell subset

hc.pb = ighv.seg[hc.idx, genes]
tm.pb = ighv.seg[tm.idx, genes]
ris.pb = ighv.seg[ris.idx, genes]

tm9.pb = ighv9.seg[tm9.idx, genes9]
ris9.pb = ighv9.seg[ris9.idx, genes9]
ris9_dna.pb = ighv9.seg[ris9_dna.idx, genes9]
ris9_rna.pb = ighv9.seg[ris9_rna.idx, genes9]
hc9.pb = ighv9.seg[hc9.idx, genes9]
hc9_dna.pb = ighv9.seg[hc9_dna.idx, genes9]
hc9_rna.pb = ighv9.seg[hc9_rna.idx, genes9]


# plasmablast
tm9.pb = tm9.pb / rowSums(tm9.pb)
tm9.pb.avg = colMeans(tm9.pb)
tm9.pb.sd = apply(tm9.pb,2,sd)
n = dim(tm9.pb)[1]
tm9.pb.se = tm9.pb.sd / sqrt(n)

ris9.pb = ris9.pb / rowSums(ris9.pb)
ris9.pb.avg = colMeans(ris9.pb)
ris9.pb.sd = apply(ris9.pb,2,sd)
n = dim(ris9.pb)[1]
ris9.pb.se = ris9.pb.sd / sqrt(n)

ris9_dna.pb = ris9_dna.pb / rowSums(ris9_dna.pb)
ris9_dna.pb.avg = colMeans(ris9_dna.pb)
ris9_dna.pb.sd = apply(ris9_dna.pb,2,sd)
n = dim(ris9_dna.pb)[1]
ris9_dna.pb.se = ris9_dna.pb.sd / sqrt(n)

ris9_rna.pb = ris9_rna.pb / rowSums(ris9_rna.pb)
ris9_rna.pb.avg = colMeans(ris9_rna.pb)
ris9_rna.pb.sd = apply(ris9_rna.pb,2,sd)
n = dim(ris9_rna.pb)[1]
ris9_rna.pb.se = ris9_rna.pb.sd / sqrt(n)

hc9.pb = hc9.pb / rowSums(hc9.pb)
hc9.pb.avg = colMeans(hc9.pb)
hc9.pb.sd = apply(hc9.pb,2,sd)
n = dim(hc9.pb)[1]
hc9.pb.se = hc9.pb.sd / sqrt(n)

hc9_dna.pb = hc9_dna.pb / rowSums(hc9_dna.pb)
hc9_dna.pb.avg = colMeans(hc9_dna.pb)
hc9_dna.pb.sd = apply(hc9_dna.pb,2,sd)
n = dim(hc9_dna.pb)[1]
hc9_dna.pb.se = hc9_dna.pb.sd / sqrt(n)

hc9_rna.pb = hc9_rna.pb / rowSums(hc9_rna.pb)
hc9_rna.pb.avg = colMeans(hc9_rna.pb)
hc9_rna.pb.sd = apply(hc9_rna.pb,2,sd)
n = dim(hc9_rna.pb)[1]
hc9_rna.pb.se = hc9_rna.pb.sd / sqrt(n)

pb9.avg = rbind(tm9.pb.avg, ris9.pb.avg, hc9.pb.avg)
pb9.se = rbind(tm9.pb.se, ris9.pb.se, hc9.pb.se)

pdf(file='lib9_ighv_family_pb_TM-HC-RIS.pdf', width=10, height=5)
barcenters = barplot(pb9.avg, beside=T, col = c("lightblue", "mistyrose", "forestgreen"), ylim=c(0,1.0), legend.text=c('TM', 'RIS', 'HC'))

segments(barcenters, pb9.avg - pb9.se * 2, barcenters, pb9.avg + pb9.se * 2, lwd = 1.5)
title('Plasmablast Library 9')

if (FALSE) {
points(barcenters[1,], tm.pb[1,], pch=16)
points(barcenters[1,], tm.pb[2,], pch=16, col="red")
points(barcenters[1,], tm.pb[3,], pch=16, col="purple")
points(barcenters[1,], tm.pb[4,], pch=17)
points(barcenters[1,], tm.pb[5,], pch=17, col="red")
points(barcenters[1,], tm.pb[6,], pch=17, col="purple")
points(barcenters[1,], tm.pb[7,], pch=17)
points(barcenters[1,], tm.pb[8,], pch=17, col="red")

points(barcenters[3,], hc.pb[1,], pch=16)
points(barcenters[3,], hc.pb[2,], pch=16, col="red")
points(barcenters[3,], hc.pb[3,], pch=17)
points(barcenters[3,], hc.pb[4,], pch=17, col="red")
points(barcenters[3,], hc.pb[5,], pch=17)
points(barcenters[3,], hc.pb[6,], pch=17, col="red")
}
dev.off()

# DNA vs RNA library 9
dna_rna.pb.avg = rbind(ris9_dna.pb.avg, ris9_rna.pb.avg, hc9_dna.pb.avg, hc9_rna.pb.avg)
dna_rna.pb.se = rbind(ris9_dna.pb.se, ris9_rna.pb.se, hc9_dna.pb.se, hc9_rna.pb.se)

pdf(file='lib9_ighv_family_pb_DNA-RNA.pdf', width=10, height=5)
barcenters = barplot(dna_rna.pb.avg, beside=T, col = c("lightblue", "mistyrose", "forestgreen", "magenta"), ylim=c(0,1.0), legend.text=c('RIS DNA', 'RIS RNA', 'HC DNA', 'HC RNA'))
segments(barcenters, dna_rna.pb.avg - dna_rna.pb.se * 2, barcenters, dna_rna.pb.avg + dna_rna.pb.se * 2, lwd = 1.5)
title('Plasmablast Library 9')
dev.off()

hc.all.pb = rbind(hc.pb, hc9.pb)
tm.all.pb = rbind(tm.pb, tm9.pb)
ris.all.pb = rbind(ris.pb, ris9.pb)

hc.all.pb = hc.all.pb / rowSums(hc.all.pb)
hc.all.pb.avg = colMeans(hc.all.pb)
hc.all.pb.sd = apply(hc.all.pb,2,sd)
n = dim(hc.all.pb)[1]
hc.all.pb.se = hc.all.pb.sd / sqrt(n)

tm.all.pb = tm.all.pb / rowSums(tm.all.pb)
tm.all.pb.avg = colMeans(tm.all.pb)
tm.all.pb.sd = apply(tm.all.pb,2,sd)
n = dim(tm.all.pb)[1]
tm.all.pb.se = tm.all.pb.sd / sqrt(n)

ris.all.pb = ris.all.pb / rowSums(ris.all.pb)
ris.all.pb.avg = colMeans(ris.all.pb)
ris.all.pb.sd = apply(ris.all.pb,2,sd)
n = dim(ris.all.pb)[1]
ris.all.pb.se = ris.all.pb.sd / sqrt(n)

pb.avg = rbind(tm.all.pb.avg, ris.all.pb.avg, hc.all.pb.avg)
pb.se = rbind(tm.all.pb.se, ris.all.pb.se, hc.all.pb.se)

pdf(file='lib_4_6_9_ighv_family_pb_TM-HC-RIS.pdf', width=10, height=5)
barcenters = barplot(pb.avg, beside=T, col = c("lightblue", "mistyrose", "forestgreen"), ylim=c(0,1.0), legend.text=c('TM', 'RIS', 'HC'))

segments(barcenters, pb.avg - pb.se * 2, barcenters, pb.avg + pb.se * 2, lwd = 1.5)
title('Plasmablast Library 4, 6, 9')
dev.off()
}

# TM/HC/RIS
# libraries 1, 4, 6, 8, 9

# gene counts

library(rjson)
library(tidyr)

data_dir = '/Users/s166813/Projects/data/immune/MonsonLab/TM-HC-RIS/vdjserver/'
lib_dir = paste(data_dir, 'library_all/stats/', sep='')
#lib_dir = paste(data_dir, 'library_all/0shm/', sep='')
#lib_dir = paste(data_dir, 'library_all/gt3shm/', sep='')

metadata <- fromJSON(file = paste(data_dir,'repertoires.v2.airr.json',sep=''))

#exclude_repertoires = c()

#nrep = length(metadata$Repertoire)
#for (i in 1:nrep) {
#    rep_id = metadata$Repertoire[[i]]$repertoire_id
#    library = metadata$Repertoire[[i]]$sample[[1]]$sequencing_run_id
#    print(library)
    #ighv.seg = read.table(paste(data_dir,'allele.germ.clone.airr.productive.v_gene_usage.csv',sep=''), header=T, sep=',')
#}

# load raw tables
#filename = 'ighv_gene_TM-HC-RIS.pdf'
#liball.HC = read.table(paste(lib_dir, 'HC.group.v_call.tsv',sep=''), header=T, sep='\t')
#liball.RIS = read.table(paste(lib_dir, 'RIS.group.v_call.tsv',sep=''), header=T, sep='\t')
#liball.TM = read.table(paste(lib_dir, 'TM.group.v_call.tsv',sep=''), header=T, sep='\t')

#filename = 'ighv_gene_TM-HC-RIS_NB.pdf'
#liball.HC = read.table(paste(lib_dir, 'HC_NB.group.v_call.tsv',sep=''), header=T, sep='\t')
#liball.RIS = read.table(paste(lib_dir, 'RIS_NB.group.v_call.tsv',sep=''), header=T, sep='\t')
#liball.TM = read.table(paste(lib_dir, 'TM_NB.group.v_call.tsv',sep=''), header=T, sep='\t')

#filename = 'ighv_gene_TM-HC-RIS_MB.pdf'
#liball.HC = read.table(paste(lib_dir, 'HC_MB.group.v_call.tsv',sep=''), header=T, sep='\t')
#liball.RIS = read.table(paste(lib_dir, 'RIS_MB.group.v_call.tsv',sep=''), header=T, sep='\t')
#liball.TM = read.table(paste(lib_dir, 'TM_MB.group.v_call.tsv',sep=''), header=T, sep='\t')

#filename = 'ighv_gene_TM-HC-RIS_PB_DNA.pdf'
#filename = 'ighv_gene_TM-HC-RIS_PB_DNA_0shm.pdf'
#filename = 'ighv_gene_TM-HC-RIS_PB_DNA_gt3shm.pdf'
filename = 'seq_Fig3_ighv_gene_TM-HC-RIS_PB_DNA.pdf'
#filename = 'seq_Fig3_ighv_gene_TM-HC-RIS_PB_DNA_0shm.pdf'
#filename = 'seq_Fig3_ighv_gene_TM-HC-RIS_PB_DNA_gt3shm.pdf'
liball.HC = read.table(paste(lib_dir, 'HC_PB_DNA.group.v_call.tsv',sep=''), header=T, sep='\t')
liball.RIS = read.table(paste(lib_dir, 'RIS_PB_DNA.group.v_call.tsv',sep=''), header=T, sep='\t')
liball.TM = read.table(paste(lib_dir, 'TM_PB_DNA.group.v_call.tsv',sep=''), header=T, sep='\t')

#filename = 'ighv_gene_TM.pdf'
#liball.HC = read.table(paste(lib_dir, 'TM_MB.group.v_call.tsv',sep=''), header=T, sep='\t')
#liball.RIS = read.table(paste(lib_dir, 'TM_PB.group.v_call.tsv',sep=''), header=T, sep='\t')
#liball.TM = read.table(paste(lib_dir, 'TM_NB.group.v_call.tsv',sep=''), header=T, sep='\t')

#filename = 'ighv_gene_RIS.pdf'
#liball.HC = read.table(paste(lib_dir, 'RIS_MB.group.v_call.tsv',sep=''), header=T, sep='\t')
#liball.RIS = read.table(paste(lib_dir, 'RIS_PB.group.v_call.tsv',sep=''), header=T, sep='\t')
#liball.TM = read.table(paste(lib_dir, 'RIS_NB.group.v_call.tsv',sep=''), header=T, sep='\t')

#filename = 'ighv_gene_HC.pdf'
#liball.HC = read.table(paste(lib_dir, 'HC_MB.group.v_call.tsv',sep=''), header=T, sep='\t')
#liball.RIS = read.table(paste(lib_dir, 'HC_PB.group.v_call.tsv',sep=''), header=T, sep='\t')
#liball.TM = read.table(paste(lib_dir, 'HC_NB.group.v_call.tsv',sep=''), header=T, sep='\t')


# extract the data
#level = 'subgroup'
level = 'gene'
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

pdf(file=filename, width=40, height=5)
#pdf(file=filename, width=7, height=7)
barcenters = barplot(ighv.all.data, beside=T, ylim=c(0,0.25), col = c("lightblue", "mistyrose", "chocolate"), legend.text=ighv.all.groups)
segments(barcenters, ighv.all.data - ighv.all.se * 2, barcenters, ighv.all.data + ighv.all.se * 2, lwd = 1.5)
dev.off()




if (FALSE) {
ighv.seg = read.table(paste(data_dir,'allele.germ.clone.airr.productive.v_gene_usage.csv',sep=''), header=T, sep=',')
rownames(ighv.seg) = t(ighv.seg["sample_id"])

# ighv genes
IGHV = sort(names(ighv.seg))[1:91]
genes <- match(IGHV,names(ighv.seg))
ighv.rel = ighv.seg[,genes] / rowSums(ighv.seg[,genes])

# gDNA (lib6) vs cDNA (lib4)
# by subject, by cell subset

# subject 278
group = c("278-1", "278-2", "278-3", "278-4", "278-5", "278-BC2", "278-BC4", "278-BC6")

# exclude when below threshold
ighv.avg = colMeans(ighv.rel[group,])
exclude.IGHV = c()
for (n in names(ighv.rel)) {
    if (ighv.avg[n] < 0.01) {
        exclude.IGHV = c(exclude.IGHV, n)
    }
}
use.IGHV = IGHV[! IGHV %in% exclude.IGHV]
genes <- match(use.IGHV,names(ighv.seg))

gdna.naive = ighv.seg["278-1",genes]
cdna.naive = ighv.seg["278-BC2",genes]
gdna.mem = ighv.seg["278-2",genes]
cdna.mem = ighv.seg["278-BC4",genes]
cdna.lib6.mem = ighv.seg["278-4",genes]
gdna.pb = ighv.seg["278-3",genes]
cdna.pb = ighv.seg["278-BC6",genes]
cdna.lib6.pb = ighv.seg["278-5",genes]

gdna.naive = gdna.naive / rowSums(gdna.naive)
gdna.mem = gdna.mem / rowSums(gdna.mem)
gdna.pb = gdna.pb / rowSums(gdna.pb)
cdna.naive = cdna.naive / rowSums(cdna.naive)
cdna.mem = cdna.mem / rowSums(cdna.mem)
cdna.pb = cdna.pb / rowSums(cdna.pb)
cdna.lib6.mem = cdna.lib6.mem / rowSums(cdna.lib6.mem)
cdna.lib6.pb = cdna.lib6.pb / rowSums(cdna.lib6.pb)

seg.avg = as.matrix(rbind(gdna.naive, cdna.naive, gdna.mem, cdna.lib6.mem, cdna.mem, gdna.pb, cdna.lib6.pb, cdna.pb))

pdf(file='ighv_gene_278_gDNA_vs_cDNA.pdf', width=40, height=5)
barcenters = barplot(seg.avg, beside=T, col = c("lightblue", "lightblue", "mistyrose", "mistyrose", "mistyrose", "forestgreen", "forestgreen", "forestgreen"), legend.text=c('gDNA naive (lib6)', 'cDNA naive (lib4)', 'gDNA memory (lib6)', 'cDNA memory (lib6)', 'cDNA memory (lib4)', 'gDNA PB (lib6)', 'cDNA PB (lib6)', 'cDNA PB (lib4)'))
title('Subject 278')
dev.off()

# subject 279
group = c("279-1", "279-2", "279-3", "279-4", "279-5", "279-BC2", "279-BC4", "279-BC6")

# exclude when below threshold
ighv.avg = colMeans(ighv.rel[group,])
exclude.IGHV = c()
for (n in names(ighv.rel)) {
    if (ighv.avg[n] < 0.01) {
        exclude.IGHV = c(exclude.IGHV, n)
    }
}
use.IGHV = IGHV[! IGHV %in% exclude.IGHV]
genes <- match(use.IGHV,names(ighv.seg))

gdna.naive = ighv.seg["279-1",genes]
cdna.naive = ighv.seg["279-BC2",genes]
gdna.mem = ighv.seg["279-2",genes]
cdna.mem = ighv.seg["279-BC4",genes]
cdna.lib6.mem = ighv.seg["279-4",genes]
gdna.pb = ighv.seg["279-3",genes]
cdna.pb = ighv.seg["279-BC6",genes]
cdna.lib6.pb = ighv.seg["279-5",genes]

gdna.naive = gdna.naive / rowSums(gdna.naive)
gdna.mem = gdna.mem / rowSums(gdna.mem)
gdna.pb = gdna.pb / rowSums(gdna.pb)
cdna.naive = cdna.naive / rowSums(cdna.naive)
cdna.mem = cdna.mem / rowSums(cdna.mem)
cdna.pb = cdna.pb / rowSums(cdna.pb)
cdna.lib6.mem = cdna.lib6.mem / rowSums(cdna.lib6.mem)
cdna.lib6.pb = cdna.lib6.pb / rowSums(cdna.lib6.pb)

seg.avg = as.matrix(rbind(gdna.naive, cdna.naive, gdna.mem, cdna.lib6.mem, cdna.mem, gdna.pb, cdna.lib6.pb, cdna.pb))

pdf(file='ighv_gene_279_gDNA_vs_cDNA.pdf', width=40, height=5)
barcenters = barplot(seg.avg, beside=T, ylim=c(0,1), col = c("lightblue", "lightblue", "mistyrose", "mistyrose", "mistyrose", "forestgreen", "forestgreen", "forestgreen"), legend.text=c('gDNA naive (lib6)', 'cDNA naive (lib4)', 'gDNA memory (lib6)', 'cDNA memory (lib6)', 'cDNA memory (lib4)', 'gDNA PB (lib6)', 'cDNA PB (lib6)', 'cDNA PB (lib4)'))
title('Subject 279')
dev.off()

# subject 1395
group = c("1395-1", "1395-2", "1395-3", "1395-4", "1395-5", "1395-BC1", "1395-BC3", "1395-BC5")

# exclude when below threshold
ighv.avg = colMeans(ighv.rel[group,])
exclude.IGHV = c()
for (n in names(ighv.rel)) {
    if (ighv.avg[n] < 0.01) {
        exclude.IGHV = c(exclude.IGHV, n)
    }
}
use.IGHV = IGHV[! IGHV %in% exclude.IGHV]
genes <- match(use.IGHV,names(ighv.seg))

gdna.naive = ighv.seg["1395-1",genes]
cdna.naive = ighv.seg["1395-BC1",genes]
gdna.mem = ighv.seg["1395-2",genes]
cdna.mem = ighv.seg["1395-BC3",genes]
cdna.lib6.mem = ighv.seg["1395-4",genes]
gdna.pb = ighv.seg["1395-3",genes]
cdna.pb = ighv.seg["1395-BC5",genes]
cdna.lib6.pb = ighv.seg["1395-5",genes]

gdna.naive = gdna.naive / rowSums(gdna.naive)
gdna.mem = gdna.mem / rowSums(gdna.mem)
gdna.pb = gdna.pb / rowSums(gdna.pb)
cdna.naive = cdna.naive / rowSums(cdna.naive)
cdna.mem = cdna.mem / rowSums(cdna.mem)
cdna.pb = cdna.pb / rowSums(cdna.pb)
cdna.lib6.mem = cdna.lib6.mem / rowSums(cdna.lib6.mem)
cdna.lib6.pb = cdna.lib6.pb / rowSums(cdna.lib6.pb)

seg.avg = as.matrix(rbind(gdna.naive, cdna.naive, gdna.mem, cdna.lib6.mem, cdna.mem, gdna.pb, cdna.lib6.pb, cdna.pb))

pdf(file='ighv_gene_1395_gDNA_vs_cDNA.pdf', width=40, height=5)
barcenters = barplot(seg.avg, beside=T, ylim=c(0,1), col = c("lightblue", "lightblue", "mistyrose", "mistyrose", "mistyrose", "forestgreen", "forestgreen", "forestgreen"), legend.text=c('gDNA naive (lib6)', 'cDNA naive (lib4)', 'gDNA memory (lib6)', 'cDNA memory (lib6)', 'cDNA memory (lib4)', 'gDNA PB (lib6)', 'cDNA PB (lib6)', 'cDNA PB (lib4)'))
title('Subject 1395')
dev.off()

# subject 1842
group = c("1842-1", "1842-2", "1842-3", "1842-4", "1842-5", "1842-BC1", "1842-BC3", "1842-BC5")

# exclude when below threshold
ighv.avg = colMeans(ighv.rel[group,])
exclude.IGHV = c()
for (n in names(ighv.rel)) {
    if (ighv.avg[n] < 0.01) {
        exclude.IGHV = c(exclude.IGHV, n)
    }
}
use.IGHV = IGHV[! IGHV %in% exclude.IGHV]
genes <- match(use.IGHV,names(ighv.seg))

gdna.naive = ighv.seg["1842-1",genes]
cdna.naive = ighv.seg["1842-BC1",genes]
gdna.mem = ighv.seg["1842-2",genes]
cdna.mem = ighv.seg["1842-BC3",genes]
cdna.lib6.mem = ighv.seg["1842-4",genes]
gdna.pb = ighv.seg["1842-3",genes]
cdna.pb = ighv.seg["1842-BC5",genes]
cdna.lib6.pb = ighv.seg["1842-5",genes]

gdna.naive = gdna.naive / rowSums(gdna.naive)
gdna.mem = gdna.mem / rowSums(gdna.mem)
gdna.pb = gdna.pb / rowSums(gdna.pb)
cdna.naive = cdna.naive / rowSums(cdna.naive)
cdna.mem = cdna.mem / rowSums(cdna.mem)
cdna.pb = cdna.pb / rowSums(cdna.pb)
cdna.lib6.mem = cdna.lib6.mem / rowSums(cdna.lib6.mem)
cdna.lib6.pb = cdna.lib6.pb / rowSums(cdna.lib6.pb)

seg.avg = as.matrix(rbind(gdna.naive, cdna.naive, gdna.mem, cdna.lib6.mem, cdna.mem, gdna.pb, cdna.lib6.pb, cdna.pb))

pdf(file='ighv_gene_1842_gDNA_vs_cDNA.pdf', width=40, height=5)
barcenters = barplot(seg.avg, beside=T, ylim=c(0,1), col = c("lightblue", "lightblue", "mistyrose", "mistyrose", "mistyrose", "forestgreen", "forestgreen", "forestgreen"), legend.text=c('gDNA naive (lib6)', 'cDNA naive (lib4)', 'gDNA memory (lib6)', 'cDNA memory (lib6)', 'cDNA memory (lib4)', 'gDNA PB (lib6)', 'cDNA PB (lib6)', 'cDNA PB (lib4)'))
title('Subject 1842')
dev.off()

# subject 2260
group = c("2260-1", "2260-2", "2260-3", "2260-BC1", "2260-BC3", "2260-BC5")

# exclude when below threshold
ighv.avg = colMeans(ighv.rel[group,])
exclude.IGHV = c()
for (n in names(ighv.rel)) {
    if (ighv.avg[n] < 0.01) {
        exclude.IGHV = c(exclude.IGHV, n)
    }
}
use.IGHV = IGHV[! IGHV %in% exclude.IGHV]
genes <- match(use.IGHV,names(ighv.seg))

gdna.naive = ighv.seg["2260-1",genes]
cdna.naive = ighv.seg["2260-BC1",genes]
gdna.mem = ighv.seg["2260-2",genes]
cdna.mem = ighv.seg["2260-BC3",genes]
gdna.pb = ighv.seg["2260-3",genes]
cdna.pb = ighv.seg["2260-BC5",genes]

gdna.naive = gdna.naive / rowSums(gdna.naive)
gdna.mem = gdna.mem / rowSums(gdna.mem)
gdna.pb = gdna.pb / rowSums(gdna.pb)
cdna.naive = cdna.naive / rowSums(cdna.naive)
cdna.mem = cdna.mem / rowSums(cdna.mem)
cdna.pb = cdna.pb / rowSums(cdna.pb)

seg.avg = as.matrix(rbind(gdna.naive, cdna.naive, gdna.mem, cdna.mem, gdna.pb, cdna.pb))

pdf(file='ighv_gene_2260_gDNA_vs_cDNA.pdf', width=40, height=5)
barcenters = barplot(seg.avg, beside=T, ylim=c(0,1), col = c("lightblue", "lightblue", "mistyrose", "mistyrose", "forestgreen", "forestgreen"), legend.text=c('gDNA naive (lib6)', 'cDNA naive (lib4)', 'gDNA memory (lib6)', 'cDNA memory (lib4)', 'gDNA PB (lib6)', 'cDNA PB (lib4)'))
title('Subject 2260')
dev.off()

# TM/RIS/HC comparison
# by cell subset

tm.naive = ighv.seg[c("1395-1", "1842-1", "2260-1", "1395-BC1", "1842-BC1", "2260-BC1"), genes]
ris.naive = ighv.seg[c("UTSW004-1", "UTSW005-1", "UTSW006-1", "UTSW007-1"), genes]
hc.naive = ighv.seg[c("278-1", "279-1", "278-BC2", "279-BC2"), genes]
tm.mem = ighv.seg[c("1395-2", "1842-2", "2260-2", "1395-BC3", "1842-BC3", "2260-BC3"), genes]
ris.mem = ighv.seg[c("UTSW004-2", "UTSW005-2", "UTSW006-2", "UTSW007-2"), genes]
hc.mem = ighv.seg[c("278-2", "279-2", "278-BC4", "279-BC4"), genes]
tm.pb = ighv.seg[c("1395-3", "1842-3", "2260-3", "1395-BC5", "1842-BC5", "2260-BC5"), genes]
ris.pb = ighv.seg[c("UTSW004-3", "UTSW005-3", "UTSW006-3", "UTSW007-3"), genes]
hc.pb = ighv.seg[c("278-3", "279-3", "278-BC6", "279-BC6"), genes]

# naive
tm.naive = tm.naive / rowSums(tm.naive)
tm.naive.avg = colMeans(tm.naive)
tm.naive.sd = apply(tm.naive,2,sd)
n = dim(tm.naive)[1]
tm.naive.se = tm.naive.sd / sqrt(n)

ris.naive = ris.naive / rowSums(ris.naive)
ris.naive.avg = colMeans(ris.naive)
ris.naive.sd = apply(ris.naive,2,sd)
n = dim(ris.naive)[1]
ris.naive.se = ris.naive.sd / sqrt(n)

hc.naive = hc.naive / rowSums(hc.naive)
hc.naive.avg = colMeans(hc.naive)
hc.naive.sd = apply(hc.naive,2,sd)
n = dim(hc.naive)[1]
hc.naive.se = hc.naive.sd / sqrt(n)

naive.avg = rbind(tm.naive.avg, ris.naive.avg, hc.naive.avg)
naive.se = rbind(tm.naive.se, ris.naive.se, hc.naive.se)

pdf(file='ighv_gene_naive_TM-HC-RIS.pdf', width=40, height=5)
barcenters = barplot(naive.avg, beside=T, col = c("lightblue", "mistyrose", "forestgreen"), ylim=c(0,0.5), legend.text=c('TM', 'RIS', 'HC'))

segments(barcenters, naive.avg - naive.se * 2, barcenters, naive.avg + naive.se * 2, lwd = 1.5)
title('Naive')
dev.off()

# memory
tm.mem = tm.mem / rowSums(tm.mem)
tm.mem.avg = colMeans(tm.mem)
tm.mem.sd = apply(tm.mem,2,sd)
n = dim(tm.mem)[1]
tm.mem.se = tm.mem.sd / sqrt(n)

ris.mem = ris.mem / rowSums(ris.mem)
ris.mem.avg = colMeans(ris.mem)
ris.mem.sd = apply(ris.mem,2,sd)
n = dim(ris.mem)[1]
ris.mem.se = ris.mem.sd / sqrt(n)

hc.mem = hc.mem / rowSums(hc.mem)
hc.mem.avg = colMeans(hc.mem)
hc.mem.sd = apply(hc.mem,2,sd)
n = dim(hc.mem)[1]
hc.mem.se = hc.mem.sd / sqrt(n)

mem.avg = rbind(tm.mem.avg, ris.mem.avg, hc.mem.avg)
mem.se = rbind(tm.mem.se, ris.mem.se, hc.mem.se)

pdf(file='ighv_gene_mem_TM-HC-RIS.pdf', width=40, height=5)


barcenters = barplot(mem.avg, beside=T, col = c("lightblue", "mistyrose", "forestgreen"), ylim=c(0,0.8), legend.text=c('TM', 'RIS', 'HC'))

segments(barcenters, mem.avg - mem.se * 2, barcenters, mem.avg + mem.se * 2, lwd = 1.5)
title('Memory')
dev.off()

# plasmablast
tm.pb = tm.pb / rowSums(tm.pb)
tm.pb.avg = colMeans(tm.pb)
tm.pb.sd = apply(tm.pb,2,sd)
n = dim(tm.pb)[1]
tm.pb.se = tm.pb.sd / sqrt(n)

ris.pb = ris.pb / rowSums(ris.pb)
ris.pb.avg = colMeans(ris.pb)
ris.pb.sd = apply(ris.pb,2,sd)
n = dim(ris.pb)[1]
ris.pb.se = ris.pb.sd / sqrt(n)

hc.pb = hc.pb / rowSums(hc.pb)
hc.pb.avg = colMeans(hc.pb)
hc.pb.sd = apply(hc.pb,2,sd)
n = dim(hc.pb)[1]
hc.pb.se = hc.pb.sd / sqrt(n)

pb.avg = rbind(tm.pb.avg, ris.pb.avg, hc.pb.avg)
pb.se = rbind(tm.pb.se, ris.pb.se, hc.pb.se)

pdf(file='ighv_gene_pb_TM-HC-RIS.pdf', width=40, height=5)
barcenters = barplot(pb.avg, beside=T, col = c("lightblue", "mistyrose", "forestgreen"), ylim=c(0,1.0), legend.text=c('TM', 'RIS', 'HC'))

segments(barcenters, pb.avg - pb.se * 2, barcenters, pb.avg + pb.se * 2, lwd = 1.5)
title('Plasmablast')
dev.off()
}

# TCR library 5/7

# gene family counts

data_dir = '/work/data/immune/MonsonLab/AD_HD_TCR/'

trbv.seg = read.table(paste(data_dir,'ACS_HC_family_usage.csv',sep=''), header=T, sep=',')

# all trbv families
TRBV = c("TRBV1", "TRBV2", "TRBV3", "TRBV4", "TRBV5", "TRBV6", "TRBV7", "TRBV9", "TRBV10", "TRBV11", "TRBV12", "TRBV13", "TRBV14", "TRBV15", "TRBV16", "TRBV17", "TRBV18", "TRBV19", "TRBV20", "TRBV21", "TRBV23", "TRBV24", "TRBV25", "TRBV26", "TRBV27", "TRBV28", "TRBV29", "TRBV30")
exclude.TRBV = c("TRBV1", "TRBV3", "TRBV9", "TRBV16", "TRBV17", "TRBV23", "TRBV26")
use.TRBV = TRBV[! TRBV %in% exclude.TRBV]

genes <- match(use.TRBV,names(trbv.seg))

# ACS, Alzheimers CSF
samples.acs.cd4 = trbv.seg[trbv.seg$patient == 'ACS' & trbv.seg$cell.type == 'CD4',]
samples.acs.cd8 = trbv.seg[trbv.seg$patient == 'ACS' & trbv.seg$cell.type == 'CD8',]
# HC, healthy CSF
samples.hc.cd4 = trbv.seg[trbv.seg$patient == 'HC' & trbv.seg$cell.type == 'CD4',]
samples.hc.cd8 = trbv.seg[trbv.seg$patient == 'HC' & trbv.seg$cell.type == 'CD8',]

seg.acs.cd4 = samples.acs.cd4[,genes]
seg.acs.cd4.avg = colMeans(seg.acs.cd4)
seg.acs.cd4.sd = apply(seg.acs.cd4,2,sd)
n = dim(seg.acs.cd4)[1]
seg.acs.cd4.se = seg.acs.cd4.sd / sqrt(n)

seg.acs.cd8 = samples.acs.cd8[,genes]
seg.acs.cd8.avg = colMeans(seg.acs.cd8)
seg.acs.cd8.sd = apply(seg.acs.cd8,2,sd)
n = dim(seg.acs.cd8)[1]
seg.acs.cd8.se = seg.acs.cd8.sd / sqrt(n)

seg.hc.cd4 = samples.hc.cd4[,genes]
seg.hc.cd4.avg = colMeans(seg.hc.cd4)
seg.hc.cd4.sd = apply(seg.hc.cd4,2,sd)
n = dim(seg.hc.cd4)[1]
seg.hc.cd4.se = seg.hc.cd4.sd / sqrt(n)

seg.hc.cd8 = samples.hc.cd8[,genes]
seg.hc.cd8.avg = colMeans(seg.hc.cd8)
seg.hc.cd8.sd = apply(seg.hc.cd8,2,sd)
n = dim(seg.hc.cd8)[1]
seg.hc.cd8.se = seg.hc.cd8.sd / sqrt(n)

seg.cd4.avg = rbind(seg.acs.cd4.avg, seg.hc.cd4.avg)
seg.cd4.se = rbind(seg.acs.cd4.se, seg.hc.cd4.se)

seg.cd8.avg = rbind(seg.acs.cd8.avg, seg.hc.cd8.avg)
seg.cd8.se = rbind(seg.acs.cd8.se, seg.hc.cd8.se)

pdf(file='trbv_ACS_HC_CD4.pdf', width=20, height=5)
barcenters = barplot(seg.cd4.avg, beside=T, col = c("lightblue", "mistyrose"), ylim=c(0,0.3), legend.text=c('ACS CD4', 'HC CD4'))

segments(barcenters, seg.cd4.avg - seg.cd4.se * 2, barcenters, seg.cd4.avg + seg.cd4.se * 2, lwd = 1.5)
dev.off()

pdf(file='trbv_ACS_HC_CD8.pdf', width=20, height=5)
barcenters = barplot(seg.cd8.avg, beside=T, col = c("forestgreen", "chocolate"), ylim=c(0,0.3), legend.text=c('ACS CD8', 'HC CD8'))

segments(barcenters, seg.cd8.avg - seg.cd8.se * 2, barcenters, seg.cd8.avg + seg.cd8.se * 2, lwd = 1.5)
dev.off()

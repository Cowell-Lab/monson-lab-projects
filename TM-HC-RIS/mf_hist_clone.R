# TM/HC/RIS
# libraries 1, 4, 6, 8, 9

# mutation frequency histogram for clones

# library
library(ggplot2)
library(tidyverse)
library(viridis)

data_dir = '/Users/s166813/Projects/data/immune/MonsonLab/TM-HC-RIS/'
lib_dir = paste(data_dir, 'vdjserver/library_all/stats/', sep='')
#lib_dir = paste(data_dir, 'vdjserver/library_all/ighv4_ge3ags_mutations/', sep='')

# replacement
filename = 'rf_hist_freq.HC_PB_DNA.pdf'
mf_hist.hc = read.table(paste(lib_dir, 'HC_PB_DNA.rf_hist_freq_clone.csv',sep=''), header=T, sep=',')
pdf(file=filename, width=9, height=7)
boxplot(mf_hist.hc[,2:10], outline=F, ylim=c(0,1.0))
points(rep(1, length(mf_hist.hc[,"bin_0.05"])), mf_hist.hc[,"bin_0.05"])
points(rep(2, length(mf_hist.hc[,"bin_0.1"])), mf_hist.hc[,"bin_0.1"])
points(rep(3, length(mf_hist.hc[,"bin_0.15"])), mf_hist.hc[,"bin_0.15"])
points(rep(4, length(mf_hist.hc[,"bin_0.2"])), mf_hist.hc[,"bin_0.2"])
points(rep(5, length(mf_hist.hc[,"bin_0.25"])), mf_hist.hc[,"bin_0.25"])
points(rep(6, length(mf_hist.hc[,"bin_0.3"])), mf_hist.hc[,"bin_0.3"])
dev.off()

filename = 'rf_hist_freq.RIS_PB_DNA.pdf'
mf_hist.ris = read.table(paste(lib_dir, 'RIS_PB_DNA.rf_hist_freq_clone.csv',sep=''), header=T, sep=',')
pdf(file=filename, width=9, height=7)
boxplot(mf_hist.ris[,2:10], outline=F, ylim=c(0,1.0))
points(rep(1, length(mf_hist.ris[,"bin_0.05"])), mf_hist.ris[,"bin_0.05"])
points(rep(2, length(mf_hist.ris[,"bin_0.1"])), mf_hist.ris[,"bin_0.1"])
points(rep(3, length(mf_hist.ris[,"bin_0.15"])), mf_hist.ris[,"bin_0.15"])
points(rep(4, length(mf_hist.ris[,"bin_0.2"])), mf_hist.ris[,"bin_0.2"])
points(rep(5, length(mf_hist.ris[,"bin_0.25"])), mf_hist.ris[,"bin_0.25"])
points(rep(6, length(mf_hist.ris[,"bin_0.3"])), mf_hist.ris[,"bin_0.3"])
dev.off()

filename = 'rf_hist_freq.TM_PB_DNA.pdf'
mf_hist.tm = read.table(paste(lib_dir, 'TM_PB_DNA.rf_hist_freq_clone.csv',sep=''), header=T, sep=',')
pdf(file=filename, width=9, height=7)
boxplot(mf_hist.tm[,2:10], outline=F, ylim=c(0,1.0))
points(rep(1, length(mf_hist.tm[,"bin_0.05"])), mf_hist.tm[,"bin_0.05"])
points(rep(2, length(mf_hist.tm[,"bin_0.1"])), mf_hist.tm[,"bin_0.1"])
points(rep(3, length(mf_hist.tm[,"bin_0.15"])), mf_hist.tm[,"bin_0.15"])
points(rep(4, length(mf_hist.tm[,"bin_0.2"])), mf_hist.tm[,"bin_0.2"])
points(rep(5, length(mf_hist.tm[,"bin_0.25"])), mf_hist.tm[,"bin_0.25"])
points(rep(6, length(mf_hist.tm[,"bin_0.3"])), mf_hist.tm[,"bin_0.3"])
dev.off()

filename = 'rf_hist_count.HC_PB_DNA.pdf'
mf_hist.hc = read.table(paste(lib_dir, 'HC_PB_DNA.rf_hist_count_clone.csv',sep=''), header=T, sep=',')
pdf(file=filename, width=9, height=7)
mf_hist.hc.norm = mf_hist.hc[,2:10] / rowSums(mf_hist.hc[,2:10])
boxplot(mf_hist.hc.norm, outline=F, ylim=c(0,0.6))
points(rep(1, length(mf_hist.hc.norm[,"bin_0.05"])), mf_hist.hc.norm[,"bin_0.05"])
points(rep(2, length(mf_hist.hc.norm[,"bin_0.1"])), mf_hist.hc.norm[,"bin_0.1"])
points(rep(3, length(mf_hist.hc.norm[,"bin_0.15"])), mf_hist.hc.norm[,"bin_0.15"])
points(rep(4, length(mf_hist.hc.norm[,"bin_0.2"])), mf_hist.hc.norm[,"bin_0.2"])
points(rep(5, length(mf_hist.hc.norm[,"bin_0.25"])), mf_hist.hc.norm[,"bin_0.25"])
points(rep(6, length(mf_hist.hc.norm[,"bin_0.3"])), mf_hist.hc.norm[,"bin_0.3"])
dev.off()

filename = 'rf_hist_count.RIS_PB_DNA.pdf'
mf_hist.ris = read.table(paste(lib_dir, 'RIS_PB_DNA.rf_hist_count_clone.csv',sep=''), header=T, sep=',')
pdf(file=filename, width=9, height=7)
mf_hist.ris.norm = mf_hist.ris[,2:10] / rowSums(mf_hist.ris[,2:10])
boxplot(mf_hist.ris.norm, outline=F, ylim=c(0,0.6))
points(rep(1, length(mf_hist.ris.norm[,"bin_0.05"])), mf_hist.ris.norm[,"bin_0.05"])
points(rep(2, length(mf_hist.ris.norm[,"bin_0.1"])), mf_hist.ris.norm[,"bin_0.1"])
points(rep(3, length(mf_hist.ris.norm[,"bin_0.15"])), mf_hist.ris.norm[,"bin_0.15"])
points(rep(4, length(mf_hist.ris.norm[,"bin_0.2"])), mf_hist.ris.norm[,"bin_0.2"])
points(rep(5, length(mf_hist.ris.norm[,"bin_0.25"])), mf_hist.ris.norm[,"bin_0.25"])
points(rep(6, length(mf_hist.ris.norm[,"bin_0.3"])), mf_hist.ris.norm[,"bin_0.3"])
dev.off()

filename = 'rf_hist_count.TM_PB_DNA.pdf'
mf_hist.tm = read.table(paste(lib_dir, 'TM_PB_DNA.rf_hist_count_clone.csv',sep=''), header=T, sep=',')
pdf(file=filename, width=9, height=7)
mf_hist.tm.norm = mf_hist.tm[,2:10] / rowSums(mf_hist.tm[,2:10])
boxplot(mf_hist.tm.norm, outline=F, ylim=c(0,0.6))
points(rep(1, length(mf_hist.tm.norm[,"bin_0.05"])), mf_hist.tm.norm[,"bin_0.05"])
points(rep(2, length(mf_hist.tm.norm[,"bin_0.1"])), mf_hist.tm.norm[,"bin_0.1"])
points(rep(3, length(mf_hist.tm.norm[,"bin_0.15"])), mf_hist.tm.norm[,"bin_0.15"])
points(rep(4, length(mf_hist.tm.norm[,"bin_0.2"])), mf_hist.tm.norm[,"bin_0.2"])
points(rep(5, length(mf_hist.tm.norm[,"bin_0.25"])), mf_hist.tm.norm[,"bin_0.25"])
points(rep(6, length(mf_hist.tm.norm[,"bin_0.3"])), mf_hist.tm.norm[,"bin_0.3"])
dev.off()

# silent
if (FALSE) {
filename = 'sf_hist_freq.HC_PB_DNA.pdf'
mf_hist.hc = read.table(paste(lib_dir, 'HC_PB_DNA.sf_hist_freq_clone.csv',sep=''), header=T, sep=',')
pdf(file=filename, width=18, height=7)
boxplot(mf_hist.hc[,2:21])
dev.off()

filename = 'sf_hist_freq.RIS_PB_DNA.pdf'
mf_hist.ris = read.table(paste(lib_dir, 'RIS_PB_DNA.sf_hist_freq_clone.csv',sep=''), header=T, sep=',')
pdf(file=filename, width=18, height=7)
boxplot(mf_hist.ris[,2:21])
dev.off()

filename = 'sf_hist_freq.TM_PB_DNA.pdf'
mf_hist.tm = read.table(paste(lib_dir, 'TM_PB_DNA.sf_hist_freq_clone.csv',sep=''), header=T, sep=',')
pdf(file=filename, width=18, height=7)
boxplot(mf_hist.tm[,2:21])
dev.off()

filename = 'sf_hist_count.HC_PB_DNA.pdf'
mf_hist.hc = read.table(paste(lib_dir, 'HC_PB_DNA.sf_hist_count_clone.csv',sep=''), header=T, sep=',')
pdf(file=filename, width=18, height=7)
boxplot(mf_hist.hc[,2:21] / rowSums(mf_hist.hc[,2:21]))
dev.off()

filename = 'sf_hist_count.RIS_PB_DNA.pdf'
mf_hist.ris = read.table(paste(lib_dir, 'RIS_PB_DNA.sf_hist_count_clone.csv',sep=''), header=T, sep=',')
pdf(file=filename, width=18, height=7)
boxplot(mf_hist.ris[,2:21] / rowSums(mf_hist.ris[,2:21]))
dev.off()

filename = 'sf_hist_count.TM_PB_DNA.pdf'
mf_hist.tm = read.table(paste(lib_dir, 'TM_PB_DNA.sf_hist_count_clone.csv',sep=''), header=T, sep=',')
pdf(file=filename, width=18, height=7)
boxplot(mf_hist.tm[,2:21] / rowSums(mf_hist.tm[,2:21]))
dev.off()
}
# TM/HC/RIS
# libraries 1, 4, 6, 8, 9

# bubbleplot for clone mutation frequencies

# library
library(ggplot2)
library(tidyverse)
library(viridis)

data_dir = '/Users/s166813/Projects/data/immune/MonsonLab/TM-HC-RIS/'
#lib_dir = paste(data_dir, 'vdjserver/analysis/stats/', sep='')
lib_dir = paste(data_dir, 'vdjserver/analysis/gt3shm/', sep='')

if (TRUE) {
#filename = 'Fig7_top_clones_1_mf_PB_DNA_lm.png'
filename = 'Fig7_gt3shm_top_clones_1_mf_PB_DNA_lm.png'
lib.hc = read.table(paste(lib_dir, 'HC_PB_DNA.top_clones_1.csv',sep=''), header=T, sep=',')
lib.hc['group'] = 'HC'
lib.tm = read.table(paste(lib_dir, 'TM_PB_DNA.top_clones_1.csv',sep=''), header=T, sep=',')
lib.tm['group'] = 'TM'
lib.ris = read.table(paste(lib_dir, 'RIS_PB_DNA.top_clones_1.csv',sep=''), header=T, sep=',')
lib.ris['group'] = 'RIS'
}

if (FALSE) {
#filename = 'Fig7_top_clones_5_mf_PB_DNA_lm.png'
filename = 'Fig7_gt3shm_top_clones_5_mf_PB_DNA_lm.png'
lib.hc = read.table(paste(lib_dir, 'HC_PB_DNA.top_clones_5.csv',sep=''), header=T, sep=',')
lib.hc['group'] = 'HC'
lib.tm = read.table(paste(lib_dir, 'TM_PB_DNA.top_clones_5.csv',sep=''), header=T, sep=',')
lib.tm['group'] = 'TM'
lib.ris = read.table(paste(lib_dir, 'RIS_PB_DNA.top_clones_5.csv',sep=''), header=T, sep=',')
lib.ris['group'] = 'RIS'
}

if (FALSE) {
#filename = 'Fig6_clone1_mf_PB_DNA_minus_cdr3s.png'
filename = 'Fig6_clone1_mf_PB_DNA_lm_minus_cdr3s.png'
lib.hc = read.table(paste(lib_dir, 'HC_PB_DNA.top_clones_1_cdr3.csv',sep=''), header=T, sep=',')
lib.hc['group'] = 'HC'
lib.tm = read.table(paste(lib_dir, 'TM_PB_DNA.top_clones_1_cdr3.csv',sep=''), header=T, sep=',')
lib.tm['group'] = 'TM'
lib.ris = read.table(paste(lib_dir, 'RIS_PB_DNA.top_clones_1_cdr3.csv',sep=''), header=T, sep=',')
lib.ris['group'] = 'RIS'
}

if (FALSE) {
#filename = 'Fig6_clone5_mf_PB_DNA_minus_cdr3s.png'
filename = 'Fig6_clone5_mf_PB_DNA_lm_minus_cdr3s.png'
lib.hc = read.table(paste(lib_dir, 'HC_PB_DNA.top_clones_5_cdr3.csv',sep=''), header=T, sep=',')
lib.hc['group'] = 'HC'
lib.tm = read.table(paste(lib_dir, 'TM_PB_DNA.top_clones_5_cdr3.csv',sep=''), header=T, sep=',')
lib.tm['group'] = 'TM'
lib.ris = read.table(paste(lib_dir, 'RIS_PB_DNA.top_clones_5_cdr3.csv',sep=''), header=T, sep=',')
lib.ris['group'] = 'RIS'
}

if (FALSE) {
filename = 'clone_0.3mf_PB_DNA.png'
lib.hc = read.table(paste(lib_dir, 'HC_PB_DNA.mf_clones_0.3.csv',sep=''), header=T, sep=',')
lib.hc['group'] = 'HC'
lib.tm = read.table(paste(lib_dir, 'TM_PB_DNA.mf_clones_0.3.csv',sep=''), header=T, sep=',')
lib.tm['group'] = 'TM'
lib.ris = read.table(paste(lib_dir, 'RIS_PB_DNA.mf_clones_0.3.csv',sep=''), header=T, sep=',')
lib.ris['group'] = 'RIS'
}

color.codes = c("red", "magenta", "blue")
lib.hc['color'] = "red"
lib.ris['color'] = "magenta"
lib.tm['color'] = "blue"

data <- data.frame(rbind(lib.hc, lib.tm, lib.ris))

p <- data %>%
    arrange(desc(copy_freq)) %>%
    ggplot(aes(x=mu_freq_s_aa, y=mu_freq_r_aa, size=copy_freq)) +
    xlim(0, 0.3) +
    ylim(0, 0.5) +
    xlab("silent mutation frequency") +
    ylab("replacement mutation frequency") +
    geom_point(alpha=0.7, aes(color=factor(group))) +
    scale_size(range = c(0, 25), name="Clone Size") +
    scale_colour_manual(name="Cohort", values=c("salmon", "darkolivegreen", "blue")) +
    #scale_color_discrete(name="Cohort") +
    geom_abline(intercept = 0, slope = 1) +
    geom_smooth(data=data[data$group == 'HC',], method=lm , color="salmon", fill="darkgrey", alpha=0.2, se=TRUE, show.legend=FALSE) +
    geom_smooth(data=data[data$group == 'RIS',], method=lm , color="darkolivegreen", fill="darkgrey", alpha=0.2, se=TRUE, show.legend=FALSE) +
    geom_smooth(data=data[data$group == 'TM',], method=lm , color="blue", fill="darkgrey", alpha=0.2, se=TRUE, show.legend=FALSE)

ggsave(p, file=filename, height=6, width=9)

hc.lm = lm(lib.hc$mu_freq_r_aa ~ lib.hc$mu_freq_s_aa)
ris.lm = lm(lib.ris$mu_freq_r_aa ~ lib.ris$mu_freq_s_aa)
tm.lm = lm(lib.tm$mu_freq_r_aa ~ lib.tm$mu_freq_s_aa)

print(summary(hc.lm))
print(sum((coef(hc.lm)[2] * lib.hc$mu_freq_s_aa + coef(hc.lm)[1]) < lib.hc$mu_freq_r_aa))
print(sum((coef(hc.lm)[2] * lib.hc$mu_freq_s_aa + coef(hc.lm)[1]) > lib.hc$mu_freq_r_aa))

print(summary(ris.lm))
print(sum((coef(ris.lm)[2] * lib.ris$mu_freq_s_aa + coef(ris.lm)[1]) < lib.ris$mu_freq_r_aa))
print(sum((coef(ris.lm)[2] * lib.ris$mu_freq_s_aa + coef(ris.lm)[1]) > lib.ris$mu_freq_r_aa))

print(summary(tm.lm))
print(sum((coef(tm.lm)[2] * lib.tm$mu_freq_s_aa + coef(tm.lm)[1]) < lib.tm$mu_freq_r_aa))
print(sum((coef(tm.lm)[2] * lib.tm$mu_freq_s_aa + coef(tm.lm)[1]) > lib.tm$mu_freq_r_aa))

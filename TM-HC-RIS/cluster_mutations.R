# TM/HC/RIS
# libraries 1, 4, 6, 8, 9

# clustering

# library
library(tidyverse)
library(viridis)

# data
data_dir = '/Users/s166813/Projects/data/immune/MonsonLab/TM-HC-RIS/'
#data_dir = '/Users/scottc/Projects/data/immune/MonsonLab/TM-HC-RIS/'

#lib_dir = 'vdjserver/library_all/mutations/'
lib_dir = 'vdjserver/analysis/mutations/'
file_prefix = 'Fig8.all.'
#lib_dir = 'vdjserver/analysis/ighv1_mutations/'
#file_prefix = 'Fig8.ighv1.'
#lib_dir = 'vdjserver/analysis/ighv2_mutations/'
#file_prefix = 'Fig8.ighv2.'
#lib_dir = 'vdjserver/analysis/ighv3_mutations/'
#file_prefix = 'Fig8.ighv3.'
#lib_dir = 'vdjserver/analysis/ighv4_mutations/'
#file_prefix = 'Fig8.ighv4.'
#lib_dir = 'vdjserver/analysis/ighv5_mutations/'
#file_prefix = 'Fig8.ighv5.'
#lib_dir = 'vdjserver/analysis/ighv6_mutations/'
#file_prefix = 'Fig8.ighv6.'
#lib_dir = 'vdjserver/analysis/ighv7_mutations/'
#file_prefix = 'Fig8.ighv7.'

#lib_dir = 'vdjserver/analysis/ighv4_mutations/ighj6/'
#file_prefix = 'Fig8.ighv4.ighj6.'

#lib_dir = 'vdjserver/analysis/AJL10_mutations/'
#lib_dir = 'vdjserver/analysis/AJL10_gt3shm/'
#file_prefix = 'Fig10.AJL10.'

# beginning codons to trim
#trim_codons = 0
trim_codons = 16

#
#group.muts = read.table(paste(data_dir, lib_dir, 'mutational_report.group.csv',sep=''), header=T, sep=',')
#group.muts = read.table(paste(data_dir, lib_dir, 'mutational_report.frequency.group.csv',sep=''), header=T, sep=',')
group.muts = read.table(paste(data_dir, lib_dir, 'gene.mutations.repertoire_group.frequency.mutational_report.csv',sep=''), header=T, sep=',')
#group.muts = read.table(paste(data_dir, lib_dir, 'ighv1.mutations.repertoire_group.frequency.mutational_report.csv',sep=''), header=T, sep=',')
#group.muts = read.table(paste(data_dir, lib_dir, 'ighv2.mutations.repertoire_group.frequency.mutational_report.csv',sep=''), header=T, sep=',')
#group.muts = read.table(paste(data_dir, lib_dir, 'ighv3.mutations.repertoire_group.frequency.mutational_report.csv',sep=''), header=T, sep=',')
#group.muts = read.table(paste(data_dir, lib_dir, 'ighv4.mutations.repertoire_group.frequency.mutational_report.csv',sep=''), header=T, sep=',')
#group.muts = read.table(paste(data_dir, lib_dir, 'ighv5.mutations.repertoire_group.frequency.mutational_report.csv',sep=''), header=T, sep=',')
#group.muts = read.table(paste(data_dir, lib_dir, 'ighv6.mutations.repertoire_group.frequency.mutational_report.csv',sep=''), header=T, sep=',')
#group.muts = read.table(paste(data_dir, lib_dir, 'ighv7.mutations.repertoire_group.frequency.mutational_report.csv',sep=''), header=T, sep=',')
#group.muts = read.table(paste(data_dir, lib_dir, 'ighv4.ighj6.mutations.repertoire_group.frequency.mutational_report.csv',sep=''), header=T, sep=',')
#group.muts = read.table(paste(data_dir, lib_dir, 'gene.mutations.AJL10.repertoire_group.frequency.mutational_report.csv',sep=''), header=T, sep=',')
#group.muts = read.table(paste(data_dir, lib_dir, 'gt3shm.mutations.AJL10.repertoire_group.frequency.mutational_report.csv',sep=''), header=T, sep=',')
rownames(group.muts) = group.muts$repertoire_group_id

#rep.muts = read.table(paste(data_dir, lib_dir, 'mutational_report.repertoire.csv',sep=''), header=T, sep=',')
#rep.muts = read.table(paste(data_dir, lib_dir, 'gene.mutations.repertoire.frequency.mutational_report.csv',sep=''), header=T, sep=',')
#HC.rep.muts = read.table(paste(data_dir, lib_dir, 'HC_PB_DNA.gene.mutations.repertoire.frequency.mutational_report.csv',sep=''), header=T, sep=',')
#RIS.rep.muts = read.table(paste(data_dir, lib_dir, 'RIS_PB_DNA.gene.mutations.repertoire.frequency.mutational_report.csv',sep=''), header=T, sep=',')
#TM.rep.muts = read.table(paste(data_dir, lib_dir, 'TM_PB_DNA.gene.mutations.repertoire.frequency.mutational_report.csv',sep=''), header=T, sep=',')
#rownames(rep.muts) = rep.muts$repertoire_id

# define regions
fwr1_r_nt = c()
fwr1_s_nt = c()
fwr1_r_aa = c()
fwr1_s_aa = c()
for (i in 1:26) {
    fwr1_r_nt = c(fwr1_r_nt, paste("mu_freq_", i, "_r", sep=''))
    fwr1_s_nt = c(fwr1_s_nt, paste("mu_freq_", i, "_s", sep=''))
    fwr1_r_aa = c(fwr1_r_aa, paste("mu_freq_", i, "_r_aa", sep=''))
    fwr1_s_aa = c(fwr1_s_aa, paste("mu_freq_", i, "_s_aa", sep=''))
}
cdr1_r_nt = c()
cdr1_s_nt = c()
cdr1_r_aa = c()
cdr1_s_aa = c()
for (i in 27:38) {
    cdr1_r_nt = c(cdr1_r_nt, paste("mu_freq_", i, "_r", sep=''))
    cdr1_s_nt = c(cdr1_s_nt, paste("mu_freq_", i, "_s", sep=''))
    cdr1_r_aa = c(cdr1_r_aa, paste("mu_freq_", i, "_r_aa", sep=''))
    cdr1_s_aa = c(cdr1_s_aa, paste("mu_freq_", i, "_s_aa", sep=''))
}
fwr2_r_nt = c()
fwr2_s_nt = c()
fwr2_r_aa = c()
fwr2_s_aa = c()
for (i in 39:55) {
    fwr2_r_nt = c(fwr2_r_nt, paste("mu_freq_", i, "_r", sep=''))
    fwr2_s_nt = c(fwr2_s_nt, paste("mu_freq_", i, "_s", sep=''))
    fwr2_r_aa = c(fwr2_r_aa, paste("mu_freq_", i, "_r_aa", sep=''))
    fwr2_s_aa = c(fwr2_s_aa, paste("mu_freq_", i, "_s_aa", sep=''))
}
cdr2_r_nt = c()
cdr2_s_nt = c()
cdr2_r_aa = c()
cdr2_s_aa = c()
for (i in 56:65) {
    cdr2_r_nt = c(cdr2_r_nt, paste("mu_freq_", i, "_r", sep=''))
    cdr2_s_nt = c(cdr2_s_nt, paste("mu_freq_", i, "_s", sep=''))
    cdr2_r_aa = c(cdr2_r_aa, paste("mu_freq_", i, "_r_aa", sep=''))
    cdr2_s_aa = c(cdr2_s_aa, paste("mu_freq_", i, "_s_aa", sep=''))
}
fwr3_r_nt = c()
fwr3_s_nt = c()
fwr3_r_aa = c()
fwr3_s_aa = c()
for (i in 66:104) {
    fwr3_r_nt = c(fwr3_r_nt, paste("mu_freq_", i, "_r", sep=''))
    fwr3_s_nt = c(fwr3_s_nt, paste("mu_freq_", i, "_s", sep=''))
    fwr3_r_aa = c(fwr3_r_aa, paste("mu_freq_", i, "_r_aa", sep=''))
    fwr3_s_aa = c(fwr3_s_aa, paste("mu_freq_", i, "_s_aa", sep=''))
}
pos_cols_r_nt = c(fwr1_r_nt, cdr1_r_nt, fwr2_r_nt, cdr2_r_nt, fwr3_r_nt)
pos_cols_s_nt = c(fwr1_s_nt, cdr1_s_nt, fwr2_s_nt, cdr2_s_nt, fwr3_s_nt)
pos_cols_r_aa = c(fwr1_r_aa, cdr1_r_aa, fwr2_r_aa, cdr2_r_aa, fwr3_r_aa)
pos_cols_s_aa = c(fwr1_s_aa, cdr1_s_aa, fwr2_s_aa, cdr2_s_aa, fwr3_s_aa)

pos_cols_r_nt_avg = c()
pos_cols_r_nt_std = c()
pos_cols_r_nt_N = c()
pos_cols_r_aa_avg = c()
pos_cols_r_aa_std = c()
pos_cols_r_aa_N = c()
for (r in pos_cols_r_nt) {
    pos_cols_r_nt_avg = c(pos_cols_r_nt_avg, paste(r, "_avg", sep=''))
    pos_cols_r_nt_std = c(pos_cols_r_nt_std, paste(r, "_std", sep=''))
    pos_cols_r_nt_N = c(pos_cols_r_nt_N, paste(r, "_N", sep=''))
}
for (r in pos_cols_r_aa) {
    pos_cols_r_aa_avg = c(pos_cols_r_aa_avg, paste(r, "_avg", sep=''))
    pos_cols_r_aa_std = c(pos_cols_r_aa_std, paste(r, "_std", sep=''))
    pos_cols_r_aa_N = c(pos_cols_r_aa_N, paste(r, "_N", sep=''))
}

#lib.rep.muts = rep.muts['6531482731230204396-242ac114-0001-012',pos_cols_r_nt] / rep.muts['6531482731230204396-242ac114-0001-012','total_count']
#lib.rep.muts = rep.muts[,pos_cols_r_nt] / rep.muts[,'total_count']

#lib.muts.avg = group.muts['TM_PB',pos_cols_r_nt_avg]
#lib.muts.se = group.muts['TM_PB',pos_cols_r_nt_std] / sqrt(group.muts['TM_PB','N'])
#lib.muts.avg = group.muts['HC_PB',pos_cols_r_nt_avg]
#lib.muts.se = group.muts['HC_PB',pos_cols_r_nt_std] / sqrt(group.muts['HC_PB','N'])
#lib.muts.avg = group.muts['RIS_PB',pos_cols_r_nt_avg]
#lib.muts.se = group.muts['RIS_PB',pos_cols_r_nt_std] / sqrt(group.muts['RIS_PB','N'])

#group_name = 'HC_PB_DNA'
#group_name = 'RIS_PB_DNA'
group_name = 'TM_PB_DNA'
filename = paste(file_prefix, "mutations_", group_name, "_N=", group.muts[group_name,'N'], ".png", sep='')
lib.muts.avg = group.muts[group_name,pos_cols_r_aa_avg]
lib.muts.se = group.muts[group_name,pos_cols_r_aa_std] / sqrt(group.muts[group_name,pos_cols_r_aa_N])

#filename = paste("mutations_TM_PB_N=", group.muts['TM_PB_DNA','N'], ".png", sep='')
#lib.muts.avg = group.muts['TM_PB_DNA',pos_cols_r_nt_avg]
#lib.muts.se = group.muts['TM_PB_DNA',pos_cols_r_nt_std] / sqrt(group.muts['TM_PB_DNA',pos_cols_r_nt_N])
#filename = paste("mutations_HC_PB_N=", group.muts['HC_PB_DNA','N'], ".png", sep='')
#lib.muts.avg = group.muts['HC_PB_DNA',pos_cols_r_nt_avg]
#lib.muts.se = group.muts['HC_PB_DNA',pos_cols_r_nt_std] / sqrt(group.muts['HC_PB_DNA',pos_cols_r_nt_N])
#filename = paste(file_prefix, "mutations_RIS_PB_N=", group.muts['RIS_PB_DNA','N'], ".png", sep='')
#lib.muts.avg = group.muts['RIS_PB_DNA',pos_cols_r_nt_avg]
#lib.muts.se = group.muts['RIS_PB_DNA',pos_cols_r_nt_std] / sqrt(group.muts['RIS_PB_DNA',pos_cols_r_nt_N])
#lib.muts.avg = group.muts['RIS_PB_DNA',pos_cols_r_aa_avg]
#lib.muts.se = group.muts['RIS_PB_DNA',pos_cols_r_aa_std] / sqrt(group.muts['RIS_PB_DNA',pos_cols_r_aa_N])

# zero out any NaNs
lib.muts.se = rapply(lib.muts.se, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
# trim codons
if (trim_codons > 0) {
    for (i in 1:trim_codons) {
        lib.muts.avg[i] = 0
        lib.muts.se[i] = 0
    }
}

# construct data frame
#r = "_r"
r = "_r_aa"
extra = '_avg'
extra_std = '_std'
empty_bar <- 4
group = c()
place = c()
value = c()
se = c()
se_start = c()
se_end = c()
for (i in 1:26) {
    col = paste("mu_freq_", i, r, extra, sep='')
    group=c(group,"FWR1")
    place = c(place, i)
    avg = as.numeric(lib.muts.avg[col])
    value = c(value, avg)
    col = paste("mu_freq_", i, r, extra_std, sep='')
    stde = as.numeric(lib.muts.se[col])
    se = c(se, stde)
    if (avg - stde < 0) {
        se_start = c(se_start, 0)
    } else {
        se_start = c(se_start, avg - stde)
    }
    se_end = c(se_end, avg + stde)
}
for (i in 1:empty_bar) {
    group=c(group,"FWR1")
    place = c(place, NA)
    value = c(value, 0)
    se = c(se, 0)
    se_start = c(se_start, 0)
    se_end = c(se_end, 0)
}
for (i in 27:38) {
    col = paste("mu_freq_", i, r, extra, sep='')
    group=c(group,"CDR1")
    place = c(place, i)
    avg = as.numeric(lib.muts.avg[col])
    value = c(value, avg)
    col = paste("mu_freq_", i, r, extra_std, sep='')
    stde = as.numeric(lib.muts.se[col])
    se = c(se, stde)
    if (avg - stde < 0) {
        se_start = c(se_start, 0)
    } else {
        se_start = c(se_start, avg - stde)
    }
    se_end = c(se_end, avg + stde)
}
for (i in 1:empty_bar) {
    group=c(group,"CDR1")
    place = c(place, NA)
    value = c(value, 0)
    se = c(se, 0)
    se_start = c(se_start, 0)
    se_end = c(se_end, 0)
}
for (i in 39:55) {
    col = paste("mu_freq_", i, r, extra, sep='')
    group=c(group,"FWR2")
    place = c(place, i)
    avg = as.numeric(lib.muts.avg[col])
    value = c(value, avg)
    col = paste("mu_freq_", i, r, extra_std, sep='')
    stde = as.numeric(lib.muts.se[col])
    se = c(se, stde)
    if (avg - stde < 0) {
        se_start = c(se_start, 0)
    } else {
        se_start = c(se_start, avg - stde)
    }
    se_end = c(se_end, avg + stde)
}
for (i in 1:empty_bar) {
    group=c(group,"FWR2")
    place = c(place, NA)
    value = c(value, 0)
    se = c(se, 0)
    se_start = c(se_start, 0)
    se_end = c(se_end, 0)
}
for (i in 56:65) {
    col = paste("mu_freq_", i, r, extra, sep='')
    group=c(group,"CDR2")
    place = c(place, i)
    avg = as.numeric(lib.muts.avg[col])
    value = c(value, avg)
    col = paste("mu_freq_", i, r, extra_std, sep='')
    stde = as.numeric(lib.muts.se[col])
    se = c(se, stde)
    if (avg - stde < 0) {
        se_start = c(se_start, 0)
    } else {
        se_start = c(se_start, avg - stde)
    }
    se_end = c(se_end, avg + stde)
}
for (i in 1:empty_bar) {
    group=c(group,"CDR2")
    place = c(place, NA)
    value = c(value, 0)
    se = c(se, 0)
    se_start = c(se_start, 0)
    se_end = c(se_end, 0)
}
for (i in 66:104) {
    col = paste("mu_freq_", i, r, extra, sep='')
    group=c(group,"FWR3")
    place = c(place, i)
    avg = as.numeric(lib.muts.avg[col])
    value = c(value, avg)
    col = paste("mu_freq_", i, r, extra_std, sep='')
    stde = as.numeric(lib.muts.se[col])
    se = c(se, stde)
    if (avg - stde < 0) {
        se_start = c(se_start, 0)
    } else {
        se_start = c(se_start, avg - stde)
    }
    se_end = c(se_end, avg + stde)
}
for (i in 1:empty_bar) {
    group=c(group,"FWR3")
    place = c(place, NA)
    value = c(value, 0)
    se = c(se, 0)
    se_start = c(se_start, 0)
    se_end = c(se_end, 0)
}

data <- data.frame(
    place=place,
    group=group,
    value=value,
    se=se,
    se_start=se_start,
    se_end=se_end
)
data$id <- seq(1, nrow(data))

# Get the name and the y place of each label
label_data <- data[data$place %in% c(1,26,27,38,39,55,56,65,66,104),]
number_of_bar <- nrow(data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
#label_data$angle <- ifelse(angle < -90, angle+180, angle)
label_data$angle <- ifelse(angle < -90, angle+270, angle-90)
label_id = label_data$id
label_data$id <- ifelse(label_id < 50, label_id, label_id+0.5)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
 
# prepare a data frame for grid (scales)
ix = sort(base_data$start, index.return=T)
grid_data <- base_data[ix$ix,]
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data[1,]$start = grid_data[1,]$end + 1
grid_data[1,]$end = grid_data[1,]$start + 3
#grid_data <- grid_data[-1,]

# Make the plot
p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  #geom_bar(stat="identity", alpha=0.5) +
  geom_col(alpha=0.5) +
  geom_segment(data=grid_data, aes(x = end, y = 0.05, xend = start, yend = 0.05), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.1, xend = start, yend = 0.1), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.15, xend = start, yend = 0.15), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  annotate("text", x = rep(max(data$id),3) + c(1.1,0.5,0.5), y = c(0.04, 0.09, 0.14), label = c("0.05", "0.1", "0.15") , color="grey", size=3 , angle=5, fontface="bold", hjust=1) +
  geom_segment(aes(x = id, y = se_start, xend = id, yend = se_end), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  ylim(-0.2,0.3) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm"),
    plot.background = element_rect(fill="white")
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=-0.01, label=place, hjust=hjust), color="black",alpha=0.6, size=2.0, angle= label_data$angle, inherit.aes = FALSE ) +
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -0.03, xend = end+0.5, yend = -0.03), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE ) +
  geom_text(data=base_data, aes(x = title, y = -0.08, label=group),  colour = "black", alpha=0.8, size=2, fontface="bold", inherit.aes = FALSE)

ggsave(p, file=filename)

#geom_text(data=label_data, aes(x=id, y=-0.15, label=place, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=1.5, angle= label_data$angle, inherit.aes = FALSE ) +




if (FALSE) {
# Make the plot
p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  #geom_bar(stat="identity", alpha=0.5) +
  geom_col(alpha=0.5) +
  geom_segment(data=grid_data, aes(x = end, y = 0.2, xend = start, yend = 0.2), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.4, xend = start, yend = 0.4), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 0.6, xend = start, yend = 0.6), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  annotate("text", x = rep(max(data$id),3) + c(0.8,0.5,0.3), y = c(0.15, 0.35, 0.55), label = c("0.2", "0.4", "0.6") , color="grey", size=3 , angle=5, fontface="bold", hjust=1) +
  geom_segment(aes(x = id, y = se_start, xend = id, yend = se_end), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  ylim(-0.5,0.65) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm"),
    plot.background = element_rect(fill="white")
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=-0.025, label=place, hjust=hjust), color="black",alpha=0.6, size=1.5, angle= label_data$angle, inherit.aes = FALSE ) +
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -0.05, xend = end+0.5, yend = -0.05), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE ) +
  geom_text(data=base_data, aes(x = title, y = -0.15, label=group),  colour = "black", alpha=0.8, size=2, fontface="bold", inherit.aes = FALSE)

ggsave(p, file=filename)

#geom_text(data=label_data, aes(x=id, y=-0.15, label=place, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=1.5, angle= label_data$angle, inherit.aes = FALSE ) +
}

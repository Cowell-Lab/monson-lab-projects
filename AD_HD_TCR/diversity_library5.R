# TCR library 5

# diversity calculations
library(alakazam)
library(shazam)

data_dir = '/Users/s166813/Projects/data/immune/MonsonLab/AD_HD_TCR/library5/7c83e049-bb46-4fe1-8877-93fbdecd1af1-007/cdr3_sharing_data/'
data_dir2 = '/Users/s166813/Projects/data/immune/MonsonLab/AD_HD_TCR/library7/c17f5818-ae03-4efb-b6de-e31b2345ba54-007/cdr3_sharing_data/'

# AD, Alzheimers CSF
# HD, healthy CSF

clones.ad1 <- read.table(paste(data_dir,'file23_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.ad2 <- read.table(paste(data_dir,'file22_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.ad3 <- read.table(paste(data_dir,'file1_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.ad4 <- read.table(paste(data_dir,'file20_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.ad5 <- read.table(paste(data_dir,'file14_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.ad6 <- read.table(paste(data_dir,'file15_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')

clones.mci1 <- read.table(paste(data_dir,'file0_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.mci2 <- read.table(paste(data_dir,'file3_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.mci3 <- read.table(paste(data_dir,'file2_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.mci4 <- read.table(paste(data_dir,'file19_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.mci5 <- read.table(paste(data_dir,'file13_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.mci6 <- read.table(paste(data_dir,'file16_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')

clones.hd2 <- read.table(paste(data_dir,'file5_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.hd3 <- read.table(paste(data_dir,'file4_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.hd5 <- read.table(paste(data_dir,'file7_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.hd6 <- read.table(paste(data_dir,'file6_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.hd8 <- read.table(paste(data_dir,'file9_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')

clones.cd4.ad1 <- read.table(paste(data_dir2,'sample1_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd4.ad2 <- read.table(paste(data_dir2,'sample44_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd4.ad3 <- read.table(paste(data_dir2,'sample32_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd4.ad4 <- read.table(paste(data_dir2,'sample4_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd4.ad5 <- read.table(paste(data_dir2,'sample15_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd4.ad6 <- read.table(paste(data_dir2,'sample34_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd4.ad7 <- read.table(paste(data_dir2,'sample35_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd4.ad8 <- read.table(paste(data_dir2,'sample47_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd4.ad9 <- read.table(paste(data_dir2,'sample5_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd4.ad10 <- read.table(paste(data_dir2,'sample37_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd4.ad11 <- read.table(paste(data_dir2,'sample51_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')

clones.cd4.mci1 <- read.table(paste(data_dir2,'sample22_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd4.mci2 <- read.table(paste(data_dir2,'sample31_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd4.mci3 <- read.table(paste(data_dir2,'sample39_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd4.mci4 <- read.table(paste(data_dir2,'sample46_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')

clones.cd8.ad1 <- read.table(paste(data_dir2,'sample17_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.ad2 <- read.table(paste(data_dir2,'sample7_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.ad3 <- read.table(paste(data_dir2,'sample14_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.ad4 <- read.table(paste(data_dir2,'sample49_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.ad5 <- read.table(paste(data_dir2,'sample33_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.ad6 <- read.table(paste(data_dir2,'sample9_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.ad7 <- read.table(paste(data_dir2,'sample10_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.ad8 <- read.table(paste(data_dir2,'sample48_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.ad9 <- read.table(paste(data_dir2,'sample20_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.ad10 <- read.table(paste(data_dir2,'sample52_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.ad11 <- read.table(paste(data_dir2,'sample27_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.ad12 <- read.table(paste(data_dir2,'sample24_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.ad13 <- read.table(paste(data_dir2,'sample2_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.ad14 <- read.table(paste(data_dir2,'sample53_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.ad15 <- read.table(paste(data_dir2,'sample12_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.ad16 <- read.table(paste(data_dir2,'sample13_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.ad17 <- read.table(paste(data_dir2,'sample0_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.ad18 <- read.table(paste(data_dir2,'sample41_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')

clones.cd8.mci1 <- read.table(paste(data_dir2,'sample6_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.mci2 <- read.table(paste(data_dir2,'sample23_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.mci3 <- read.table(paste(data_dir2,'sample8_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.mci4 <- read.table(paste(data_dir2,'sample45_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.mci5 <- read.table(paste(data_dir2,'sample3_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.mci6 <- read.table(paste(data_dir2,'sample16_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.mci7 <- read.table(paste(data_dir2,'sample30_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.mci8 <- read.table(paste(data_dir2,'sample19_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.mci9 <- read.table(paste(data_dir2,'sample54_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')

clones.cd8.hd1 <- read.table(paste(data_dir2,'sample38_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.hd2 <- read.table(paste(data_dir2,'sample29_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.hd3 <- read.table(paste(data_dir2,'sample18_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.hd4 <- read.table(paste(data_dir2,'sample43_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.cd8.hd5 <- read.table(paste(data_dir2,'sample28_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')

# Pancreatic tumor
clones.pt2 <- read.table(paste(data_dir,'file24_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pt4 <- read.table(paste(data_dir,'file39_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pt6 <- read.table(paste(data_dir,'file37_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pt8 <- read.table(paste(data_dir,'file12_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pt10 <- read.table(paste(data_dir,'file8_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pt12 <- read.table(paste(data_dir,'file32_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pt14 <- read.table(paste(data_dir,'file34_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pt16 <- read.table(paste(data_dir,'file28_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pt17 <- read.table(paste(data_dir,'file29_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pt18 <- read.table(paste(data_dir,'file30_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pt19 <- read.table(paste(data_dir,'file35_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pt20 <- read.table(paste(data_dir,'file26_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')

# Pancreatic adjacent healthy
clones.pn1 <- read.table(paste(data_dir,'file36_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pn3 <- read.table(paste(data_dir,'file21_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pn5 <- read.table(paste(data_dir,'file10_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pn7 <- read.table(paste(data_dir,'file38_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pn9 <- read.table(paste(data_dir,'file11_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pn11 <- read.table(paste(data_dir,'file31_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pn13 <- read.table(paste(data_dir,'file33_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pn15 <- read.table(paste(data_dir,'file27_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')
clones.pn21 <- read.table(paste(data_dir,'file25_cdr3_vj_aa_sharing.tsv',sep=''), header=T, sep='\t')


# main diversity quotients
q <- c(0, 1, 2)

diversity.ad1 <- calcDiversity(clones.ad1$TOTAL_COUNT, q)
diversity.ad2 <- calcDiversity(clones.ad2$TOTAL_COUNT, q)
diversity.ad3 <- calcDiversity(clones.ad3$TOTAL_COUNT, q)
diversity.ad4 <- calcDiversity(clones.ad4$TOTAL_COUNT, q)
diversity.ad5 <- calcDiversity(clones.ad5$TOTAL_COUNT, q)
diversity.ad6 <- calcDiversity(clones.ad6$TOTAL_COUNT, q)

diversity.mci1 <- calcDiversity(clones.mci1$TOTAL_COUNT, q)
diversity.mci2 <- calcDiversity(clones.mci2$TOTAL_COUNT, q)
diversity.mci3 <- calcDiversity(clones.mci3$TOTAL_COUNT, q)
diversity.mci4 <- calcDiversity(clones.mci4$TOTAL_COUNT, q)
diversity.mci5 <- calcDiversity(clones.mci5$TOTAL_COUNT, q)
diversity.mci6 <- calcDiversity(clones.mci6$TOTAL_COUNT, q)

diversity.hd2 <- calcDiversity(clones.hd2$TOTAL_COUNT, q)
diversity.hd3 <- calcDiversity(clones.hd3$TOTAL_COUNT, q)
diversity.hd5 <- calcDiversity(clones.hd5$TOTAL_COUNT, q)
diversity.hd6 <- calcDiversity(clones.hd6$TOTAL_COUNT, q)
diversity.hd8 <- calcDiversity(clones.hd8$TOTAL_COUNT, q)

diversity.cd4.ad1 <- calcDiversity(clones.cd4.ad1$TOTAL_COUNT, q)
diversity.cd4.ad2 <- calcDiversity(clones.cd4.ad2$TOTAL_COUNT, q)
diversity.cd4.ad3 <- calcDiversity(clones.cd4.ad3$TOTAL_COUNT, q)
diversity.cd4.ad4 <- calcDiversity(clones.cd4.ad4$TOTAL_COUNT, q)
diversity.cd4.ad5 <- calcDiversity(clones.cd4.ad5$TOTAL_COUNT, q)
diversity.cd4.ad6 <- calcDiversity(clones.cd4.ad6$TOTAL_COUNT, q)
diversity.cd4.ad7 <- calcDiversity(clones.cd4.ad7$TOTAL_COUNT, q)
diversity.cd4.ad8 <- calcDiversity(clones.cd4.ad8$TOTAL_COUNT, q)
diversity.cd4.ad9 <- calcDiversity(clones.cd4.ad9$TOTAL_COUNT, q)
diversity.cd4.ad10 <- calcDiversity(clones.cd4.ad10$TOTAL_COUNT, q)
diversity.cd4.ad11 <- calcDiversity(clones.cd4.ad11$TOTAL_COUNT, q)

diversity.cd4.mci1 <- calcDiversity(clones.cd4.mci1$TOTAL_COUNT, q)
diversity.cd4.mci2 <- calcDiversity(clones.cd4.mci2$TOTAL_COUNT, q)
diversity.cd4.mci3 <- calcDiversity(clones.cd4.mci3$TOTAL_COUNT, q)
diversity.cd4.mci4 <- calcDiversity(clones.cd4.mci4$TOTAL_COUNT, q)

diversity.cd8.ad1 <- calcDiversity(clones.cd8.ad1$TOTAL_COUNT, q)
diversity.cd8.ad2 <- calcDiversity(clones.cd8.ad2$TOTAL_COUNT, q)
diversity.cd8.ad3 <- calcDiversity(clones.cd8.ad3$TOTAL_COUNT, q)
diversity.cd8.ad4 <- calcDiversity(clones.cd8.ad4$TOTAL_COUNT, q)
diversity.cd8.ad5 <- calcDiversity(clones.cd8.ad5$TOTAL_COUNT, q)
diversity.cd8.ad6 <- calcDiversity(clones.cd8.ad6$TOTAL_COUNT, q)
diversity.cd8.ad7 <- calcDiversity(clones.cd8.ad7$TOTAL_COUNT, q)
diversity.cd8.ad8 <- calcDiversity(clones.cd8.ad8$TOTAL_COUNT, q)
diversity.cd8.ad9 <- calcDiversity(clones.cd8.ad9$TOTAL_COUNT, q)
diversity.cd8.ad10 <- calcDiversity(clones.cd8.ad10$TOTAL_COUNT, q)
diversity.cd8.ad11 <- calcDiversity(clones.cd8.ad11$TOTAL_COUNT, q)
diversity.cd8.ad12 <- calcDiversity(clones.cd8.ad12$TOTAL_COUNT, q)
diversity.cd8.ad13 <- calcDiversity(clones.cd8.ad13$TOTAL_COUNT, q)
diversity.cd8.ad14 <- calcDiversity(clones.cd8.ad14$TOTAL_COUNT, q)
diversity.cd8.ad15 <- calcDiversity(clones.cd8.ad15$TOTAL_COUNT, q)
diversity.cd8.ad16 <- calcDiversity(clones.cd8.ad16$TOTAL_COUNT, q)
diversity.cd8.ad17 <- calcDiversity(clones.cd8.ad17$TOTAL_COUNT, q)
diversity.cd8.ad18 <- calcDiversity(clones.cd8.ad18$TOTAL_COUNT, q)

diversity.cd8.mci1 <- calcDiversity(clones.cd8.mci1$TOTAL_COUNT, q)
diversity.cd8.mci2 <- calcDiversity(clones.cd8.mci2$TOTAL_COUNT, q)
diversity.cd8.mci3 <- calcDiversity(clones.cd8.mci3$TOTAL_COUNT, q)
diversity.cd8.mci4 <- calcDiversity(clones.cd8.mci4$TOTAL_COUNT, q)
diversity.cd8.mci5 <- calcDiversity(clones.cd8.mci5$TOTAL_COUNT, q)
diversity.cd8.mci6 <- calcDiversity(clones.cd8.mci6$TOTAL_COUNT, q)
diversity.cd8.mci7 <- calcDiversity(clones.cd8.mci7$TOTAL_COUNT, q)
diversity.cd8.mci8 <- calcDiversity(clones.cd8.mci8$TOTAL_COUNT, q)
diversity.cd8.mci9 <- calcDiversity(clones.cd8.mci9$TOTAL_COUNT, q)

diversity.cd8.hd1 <- calcDiversity(clones.cd8.hd1$TOTAL_COUNT, q)
diversity.cd8.hd2 <- calcDiversity(clones.cd8.hd2$TOTAL_COUNT, q)
diversity.cd8.hd3 <- calcDiversity(clones.cd8.hd3$TOTAL_COUNT, q)
diversity.cd8.hd4 <- calcDiversity(clones.cd8.hd4$TOTAL_COUNT, q)
diversity.cd8.hd5 <- calcDiversity(clones.cd8.hd5$TOTAL_COUNT, q)


diversity.pt2 <- calcDiversity(clones.pt2$TOTAL_COUNT, q)
diversity.pt4 <- calcDiversity(clones.pt4$TOTAL_COUNT, q)
diversity.pt6 <- calcDiversity(clones.pt6$TOTAL_COUNT, q)
diversity.pt8 <- calcDiversity(clones.pt8$TOTAL_COUNT, q)
diversity.pt10 <- calcDiversity(clones.pt10$TOTAL_COUNT, q)
diversity.pt12 <- calcDiversity(clones.pt12$TOTAL_COUNT, q)
diversity.pt14 <- calcDiversity(clones.pt14$TOTAL_COUNT, q)
diversity.pt16 <- calcDiversity(clones.pt16$TOTAL_COUNT, q)
diversity.pt17 <- calcDiversity(clones.pt17$TOTAL_COUNT, q)
diversity.pt18 <- calcDiversity(clones.pt18$TOTAL_COUNT, q)
diversity.pt19 <- calcDiversity(clones.pt19$TOTAL_COUNT, q)
diversity.pt20 <- calcDiversity(clones.pt20$TOTAL_COUNT, q)

diversity.pn1 <- calcDiversity(clones.pn1$TOTAL_COUNT, q)
diversity.pn3 <- calcDiversity(clones.pn3$TOTAL_COUNT, q)
diversity.pn5 <- calcDiversity(clones.pn5$TOTAL_COUNT, q)
diversity.pn7 <- calcDiversity(clones.pn7$TOTAL_COUNT, q)
diversity.pn9 <- calcDiversity(clones.pn9$TOTAL_COUNT, q)
diversity.pn11 <- calcDiversity(clones.pn11$TOTAL_COUNT, q)
diversity.pn13 <- calcDiversity(clones.pn13$TOTAL_COUNT, q)
diversity.pn15 <- calcDiversity(clones.pn15$TOTAL_COUNT, q)
diversity.pn21 <- calcDiversity(clones.pn21$TOTAL_COUNT, q)

use_q = 1
diversity.ad <- c(diversity.ad1[use_q], diversity.ad2[use_q], diversity.ad3[use_q], diversity.ad4[use_q], diversity.ad5[use_q], diversity.ad6[use_q])
diversity.mci <- c(diversity.mci1[use_q], diversity.mci2[use_q], diversity.mci3[use_q], diversity.mci4[use_q], diversity.mci5[use_q], diversity.mci6[use_q])
diversity.hd <- c(diversity.hd2[use_q], diversity.hd3[use_q], diversity.hd5[use_q], diversity.hd6[use_q], diversity.hd8[use_q])

diversity.cd4.ad <- c(diversity.cd4.ad1[use_q], diversity.cd4.ad2[use_q], diversity.cd4.ad3[use_q], diversity.cd4.ad4[use_q], diversity.cd4.ad5[use_q], diversity.cd4.ad6[use_q], diversity.cd4.ad7[use_q], diversity.cd4.ad8[use_q], diversity.cd4.ad9[use_q], diversity.cd4.ad10[use_q], diversity.cd4.ad11[use_q])
diversity.cd4.mci <- c(diversity.cd4.mci1[use_q], diversity.cd4.mci2[use_q], diversity.cd4.mci3[use_q], diversity.cd4.mci4[use_q])

diversity.cd8.ad <- c(diversity.cd8.ad1[use_q], diversity.cd8.ad2[use_q], diversity.cd8.ad3[use_q], diversity.cd8.ad4[use_q], diversity.cd8.ad5[use_q], diversity.cd8.ad6[use_q], diversity.cd8.ad7[use_q], diversity.cd8.ad8[use_q], diversity.cd8.ad9[use_q], diversity.cd8.ad10[use_q], diversity.cd8.ad11[use_q], diversity.cd8.ad12[use_q], diversity.cd8.ad13[use_q], diversity.cd8.ad14[use_q], diversity.cd8.ad15[use_q], diversity.cd8.ad16[use_q], diversity.cd8.ad17[use_q], diversity.cd8.ad18[use_q])
diversity.cd8.mci <- c(diversity.cd8.mci1[use_q], diversity.cd8.mci2[use_q], diversity.cd8.mci3[use_q], diversity.cd8.mci4[use_q], diversity.cd8.mci5[use_q], diversity.cd8.mci6[use_q], diversity.cd8.mci7[use_q], diversity.cd8.mci8[use_q], diversity.cd8.mci9[use_q])
diversity.cd8.hd <- c(diversity.cd8.hd1[use_q], diversity.cd8.hd2[use_q], diversity.cd8.hd3[use_q], diversity.cd8.hd4[use_q], diversity.cd8.hd5[use_q])

diversity.acs <- c(diversity.ad, diversity.cd4.ad)

diversity.pt <- c(diversity.pt2[use_q], diversity.pt4[use_q], diversity.pt6[use_q], diversity.pt8[use_q], diversity.pt10[use_q], diversity.pt12[use_q], diversity.pt14[use_q], diversity.pt16[use_q], diversity.pt17[use_q], diversity.pt18[use_q], diversity.pt19[use_q], diversity.pt20[use_q])

diversity.pn <- c(diversity.pn1[use_q], diversity.pn3[use_q], diversity.pn5[use_q], diversity.pn7[use_q], diversity.pn9[use_q], diversity.pn11[use_q], diversity.pn13[use_q], diversity.pn15[use_q], diversity.pn21[use_q])

#C <- list(diversity.ad, diversity.mci, diversity.hd, diversity.cd4.ad, diversity.cd4.mci, diversity.cd8.ad, diversity.cd8.mci, diversity.cd8.hd)
C <- list(diversity.hd, diversity.acs, diversity.cd8.hd, diversity.cd8.ad)

names(C) <- c(paste("HC CD4\n n=" , length(diversity.hd) , sep=""),
    paste("ACS CD4\n n=" , length(diversity.acs) , sep=""),
    paste("HC CD8\n n=" , length(diversity.cd8.hd) , sep=""),
    paste("ACS CD8\n n=" , length(diversity.cd8.ad) , sep=""))

pdf("diversity_AD_CD4_CD8.pdf", width=6)
# Change the mgp argument: avoid text overlaps axis
par(mgp=c(3,2,0))

boxplot(C, ylab="Diversity Q=1", outline=F)

for (i in 1:length(C)) {
    myjitter <- jitter(rep(i, length(C[[i]])), amount=0.1)
    points(myjitter, C[[i]], pch=20, col=rgb(0,0,0,.9))
}
dev.off()

if (FALSE) {
C <- list(diversity.pt, diversity.pn)
names(C) <- c(paste("Pancreatic tumor\n n=" , length(diversity.pt) , sep=""), paste("Pancreatic normal\n n=" , length(diversity.pn) , sep=""))

#pdf("diversity_PT_vs_PN.pdf")
# Change the mgp argument: avoid text overlaps axis
par(mgp=c(3,2,0))

boxplot(C, ylab="Diversity Q=1" )
#dev.off()
}
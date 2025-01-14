# TM/HC/RIS
# libraries 1, 4, 6, 8, 9

# bubbleplot for clone mutation frequencies

# library
library(ggplot2)
library(tidyverse)
library(viridis)

data_dir = '/Users/s166813/Projects/data/immune/MonsonLab/TM-HC-RIS/'
lib_dir = paste(data_dir, 'vdjserver/library_all/stats/', sep='')

hc_files = c('3645-3_S12_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'HD279-3_S30_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'1742_S29_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'HD278-3_S27_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'1055_S39_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'HD281_S26_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'1235_S27_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'1053_S41_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'1764_S40_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'3851_S32_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'1308_S31_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'1768_S30_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'1408_S28_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'3852_S41_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'1540_S38_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'1068_S35_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'3049_S40_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'1065_S42_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'2265_S39_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'1202_S36_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'1420_S37_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv')

tm_files = c('2364_S17_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'2260-3_S9_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'2567_S15_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'2989_S24_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'2405_S14_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'1842-3_S6_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'2229_S45_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'3166_S12_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'1985_S25_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'2777_S22_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'2934_S23_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'2620_S28_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'2147_S13_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'1602_S16_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'2619_S27_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'2478_S18_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'2070_S19_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'1395-3_S3_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'2929-3_S21_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv')

ris_files = c('UTSW17_S9_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW18_S10_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW14_S6_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW15_S7_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW20_S44_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW10_S2_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW22_S30_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW13_S5_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW006-3_S21_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW19_S43_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW09_S1_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW16_S8_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW007-3_S24_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW29_S25_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW27_S23_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW12_S4_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW24_S32_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW11_S3_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW004-3_S15_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW31_S26_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW28_S24_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW32_S42_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW2_S22_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW26_S34_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW005-3_S18_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW21_S29_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW1_S21_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv',
'UTSW25_S33_L001_R1_001.fastq.merged.unique.igblast.makedb.gene.clone.count.tsv')

#filename = 'cumulative_clone.HC_PB_DNA.pdf'
#cohort_files = hc_files
filename = 'cumulative_clone.RIS_PB_DNA.pdf'
cohort_files = ris_files
#filename = 'cumulative_clone.TM_PB_DNA.pdf'
#cohort_files = tm_files

max_len = 0
for (i in 1:length(cohort_files)) {
    clone.hc = read.table(paste(lib_dir, cohort_files[i],sep=''), header=T, sep='\t')
    if (length(clone.hc$copy_freq) > max_len) {
        max_len = length(clone.hc$copy_freq)
    }
}

pdf(file=filename, width=14, height=7)
for (i in 1:length(cohort_files)) {
    clone.hc = read.table(paste(lib_dir, cohort_files[i],sep=''), header=T, sep='\t')
    clone.hc[1,'cumulative_copy_freq'] = clone.hc[1,'copy_freq']
    for (j in 2:length(clone.hc$copy_freq)) {
        clone.hc[j,'cumulative_copy_freq'] = clone.hc[j,'copy_freq'] + clone.hc[j-1,'cumulative_copy_freq']
    }
    if (i == 1) {
        plot(clone.hc$cumulative_copy_freq,pch=20,xlim=c(0,max_len), ylim=c(0,1.1))
    } else {
        points(clone.hc$cumulative_copy_freq,pch=20)
    }
}
dev.off()

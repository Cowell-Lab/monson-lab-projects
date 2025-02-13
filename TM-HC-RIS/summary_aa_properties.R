# TM/HC/RIS
# libraries 1, 4, 6, 8, 9, 10, 11

# Summarize CDR3 AA properties

library(airr)

#data_dir = '/Users/s166813/Projects/data/immune/MonsonLab/TM-HC-RIS/vdjserver/'
#data_dir = '/Volumes/Antigen/Projects/data/immune/MonsonLab/TM-HC-RIS/vdjserver/'
data_dir = '/work/data/immune/MonsonLab/TM-HC-RIS/vdjserver/'

#lib_dir = paste(data_dir, 'analysis_v2/stats/', sep='')
#stage = '.gene.mutations.aa_properties'
#lib_dir = paste(data_dir, 'analysis_v2/stats/ighj6/', sep='')
#stage = '.ighj6.gene.mutations.aa_properties'
#lib_dir = paste(data_dir, 'analysis_v2/stats/ighj6/ighv7/', sep='')
#stage = '.ighv7.ighj6.gene.mutations.aa_properties'

#lib_dir = paste(data_dir, 'analysis_v2/gt3shm/', sep='')
#stage = '.gt3shm.mutations.aa_properties'
lib_dir = paste(data_dir, 'analysis_v2/gt3shm_ighv1/', sep='')
stage = '.gt3shm.ighv1.mutations.aa_properties'
#lib_dir = paste(data_dir, 'analysis_v2/gt3shm/ighj6/', sep='')
#stage = '.ighj6.gt3shm.mutations.aa_properties'
#lib_dir = paste(data_dir, 'analysis_v2/gt3shm/ighj6/ighv7/', sep='')
#stage = '.ighv7.ighj6.gt3shm.mutations.aa_properties'

#g = 'HC_PB_DNA'
g = 'CIS_PB_DNA'
#g = 'RIS_PB_DNA'
#g = 'RIS_ADVANCING'
#g = 'RIS_STABLE'
#g = 'RIS_POS'
#g = 'RIS_NEG'
#v1#g = 'TM_PB_DNA'

#seq_outfile = paste(g, '.IGHV7.IGHJ6.seq.summary.aa_properties.csv', sep='')
#outfile = paste(g, '.IGHV7.IGHJ6.summary.aa_properties.csv', sep='')
#seq_outfile = paste(g, '.IGHJ6.seq.summary.aa_properties.csv', sep='')
#outfile = paste(g, '.IGHJ6.summary.aa_properties.csv', sep='')
seq_outfile = paste(g, '.IGHV1.seq.summary.aa_properties.csv', sep='')
outfile = paste(g, '.IGHV1.summary.aa_properties.csv', sep='')
#seq_outfile = paste(g, '.GT3SHM.seq.summary.aa_properties.csv', sep='')
#outfile = paste(g, '.GT3SHM.summary.aa_properties.csv', sep='')
#seq_outfile = paste(g, '.ALL.seq.summary.aa_properties.csv', sep='')
#outfile = paste(g, '.ALL.summary.aa_properties.csv', sep='')

groups <- read_airr(paste(data_dir, 'repertoire_groups.airr.json',sep=''))
reps <- read_airr(paste(data_dir, 'repertoires.v2.airr.json',sep=''))

ighv.names = c('IGHV1', 'IGHV2', 'IGHV3', 'IGHV4', 'IGHV5', 'IGHV6', 'IGHV7')
ighj.names = c('IGHJ1', 'IGHJ2', 'IGHJ3', 'IGHJ4', 'IGHJ5', 'IGHJ6')
fields = c('cdr3_aa_length', 'cdr3_aa_gravy', 'cdr3_aa_bulk', 'cdr3_aa_aliphatic', 'cdr3_aa_polarity', 'cdr3_aa_charge', 'cdr3_aa_basic', 'cdr3_aa_acidic', 'cdr3_aa_aromatic')

#hc_subject_id = c(1053,1055,1065,1068,1202,1235,1308,1408,1420,1540,1742,1763,1764,1768,2265,3049,3645,3647,3851,3852)

df <- data.frame()
seq_df <- data.frame()
for (i in 1:length(groups$RepertoireGroup)) {
    if (groups$RepertoireGroup[[i]]$repertoire_group_id == g) {
        for (j in 1:length(groups$RepertoireGroup[[i]]$repertoires)) {
            print(groups$RepertoireGroup[[i]]$repertoires[[j]]$repertoire_id)
            seq_df[j, 'repertoire_id'] = groups$RepertoireGroup[[i]]$repertoires[[j]]$repertoire_id
            df[j, 'repertoire_id'] = groups$RepertoireGroup[[i]]$repertoires[[j]]$repertoire_id
            x = read.table(paste(lib_dir, df[j, 'repertoire_id'], stage, '.airr.tsv',sep=''),sep='\t',header=T)
            seq_df[j, 'N'] = dim(x)[1]
            df[j, 'N'] = sum(x[,"duplicate_count"])
            x[,"weight"] = x[,"duplicate_count"] / sum(x[,"duplicate_count"])
            for (k in 1:length(fields)) {
                seq_df[j, paste(fields[k], '_avg', sep='')] = mean(x[,fields[k]])
                seq_df[j, paste(fields[k], '_std', sep='')] = sd(x[,fields[k]])
                xm = weighted.mean(x[,fields[k]], x[,"weight"])
                df[j, paste(fields[k], '_avg', sep='')] = xm
                df[j, paste(fields[k], '_std', sep='')] = sqrt(sum(x[,"weight"] * (x[,fields[k]] - xm)^2) * df[j, 'N'] / (df[j, 'N']-1))
            }
        }
    }
}

write.csv(seq_df, file=paste(lib_dir, seq_outfile, sep=''), row.names=FALSE)
write.csv(df, file=paste(lib_dir, outfile, sep=''), row.names=FALSE)

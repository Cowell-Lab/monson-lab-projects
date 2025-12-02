# TM/HC/RIS
# libraries 1, 4, 6, 8, 9, 10, 11

# Summarize CDR3 AA properties

library(airr)

#data_dir = '/Users/s166813/Projects/data/immune/MonsonLab/TM-HC-RIS/vdjserver/'
#data_dir = '/Volumes/Antigen/Projects/data/immune/MonsonLab/TM-HC-RIS/vdjserver/'
data_dir = '/work/data/immune/MonsonLab/TM-HC-RIS/vdjserver/'

ighv.names = c('IGHV1', 'IGHV2', 'IGHV3', 'IGHV4', 'IGHV5', 'IGHV6', 'IGHV7')
ighj.names = c('IGHJ1', 'IGHJ2', 'IGHJ3', 'IGHJ4', 'IGHJ5', 'IGHJ6')

summaryAAproperties <- function(sub_dir, stage, group_name, out_stage)
{
    lib_dir = paste(data_dir, sub_dir, sep='')
    #stage = '.gt3shm.ighv1.mutations.aa_properties'

    #lib_dir = paste(data_dir, 'analysis_v2/stats/', sep='')
    #stage = '.gene.mutations.aa_properties'
    #lib_dir = paste(data_dir, 'analysis_v2/stats/ighj6/', sep='')
    #stage = '.ighj6.gene.mutations.aa_properties'
    #lib_dir = paste(data_dir, 'analysis_v2/stats/ighj6/ighv7/', sep='')
    #stage = '.ighv7.ighj6.gene.mutations.aa_properties'

    #lib_dir = paste(data_dir, 'analysis_v2/gt3shm/', sep='')
    #stage = '.gt3shm.mutations.aa_properties'
    #lib_dir = paste(data_dir, 'analysis_v2/gt3shm_ighv1/', sep='')
    #stage = '.gt3shm.ighv1.mutations.aa_properties'
    #lib_dir = paste(data_dir, 'analysis_v2/gt3shm/ighj6/', sep='')
    #stage = '.ighj6.gt3shm.mutations.aa_properties'
    #lib_dir = paste(data_dir, 'analysis_v2/gt3shm/ighj6/ighv7/', sep='')
    #stage = '.ighv7.ighj6.gt3shm.mutations.aa_properties'

    #g = 'HC_PB_DNA'
    #g = 'CIS_PB_DNA'
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
    seq_outfile = paste(group_name, out_stage, '.summary.seq.aa_properties.csv', sep='')
    outfile = paste(group_name, out_stage, '.summary.aa_properties.csv', sep='')

    groups <- read_airr(paste(data_dir, 'repertoire_groups.airr.json',sep=''))
    reps <- read_airr(paste(data_dir, 'repertoires.v2.airr.json',sep=''))

    fields = c('cdr3_aa_length', 'cdr3_aa_gravy', 'cdr3_aa_bulk', 'cdr3_aa_aliphatic', 'cdr3_aa_polarity', 'cdr3_aa_charge', 'cdr3_aa_basic', 'cdr3_aa_acidic', 'cdr3_aa_aromatic')

    df <- data.frame()
    seq_df <- data.frame()
    for (i in 1:length(groups$RepertoireGroup)) {
        if (groups$RepertoireGroup[[i]]$repertoire_group_id == group_name) {
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
}

# local_dir = 'analysis_v3/stats/'
# stage = '.gene.mutations.aa_properties'
# summaryAAproperties(local_dir, '.gene.mutations.aa_properties', 'HC_PB_DNA', '.ALL')
# summaryAAproperties(local_dir, '.gene.mutations.aa_properties', 'CIS_PB_DNA', '.ALL')
# summaryAAproperties(local_dir, '.gene.mutations.aa_properties', 'RIS_PB_DNA', '.ALL')
# summaryAAproperties(local_dir, '.gene.mutations.aa_properties', 'ADVANCING', '.ALL')
# summaryAAproperties(local_dir, '.gene.mutations.aa_properties', 'STABLE', '.ALL')
# summaryAAproperties(local_dir, '.gene.mutations.aa_properties', 'CORD_POS', '.ALL')
# summaryAAproperties(local_dir, '.gene.mutations.aa_properties', 'CORD_NEG', '.ALL')
# summaryAAproperties(local_dir, '.gene.mutations.aa_properties', 'RIS_ADVANCING', '.ALL')
# summaryAAproperties(local_dir, '.gene.mutations.aa_properties', 'RIS_STABLE', '.ALL')
# summaryAAproperties(local_dir, '.gene.mutations.aa_properties', 'RIS_POS', '.ALL')
# summaryAAproperties(local_dir, '.gene.mutations.aa_properties', 'RIS_NEG', '.ALL')
# summaryAAproperties(local_dir, '.gene.mutations.aa_properties', 'RIS_CSRHI', '.ALL')
# summaryAAproperties(local_dir, '.gene.mutations.aa_properties', 'RIS_CSRLO', '.ALL')
# summaryAAproperties(local_dir, '.gene.mutations.aa_properties', 'CIS_ADVANCING', '.ALL')
# summaryAAproperties(local_dir, '.gene.mutations.aa_properties', 'CIS_POS', '.ALL')
# summaryAAproperties(local_dir, '.gene.mutations.aa_properties', 'CIS_NEG', '.ALL')

# for (i in 1:length(ighj.names)) {
#     print(ighj.names[i])
#     local_dir = paste('analysis_v3/stats/', tolower(ighj.names[i]), '/', sep='')
#     stage = paste('.', tolower(ighj.names[i]), '.gene.mutations.aa_properties', sep='')
#     outn = paste('.', ighj.names[i], '.ALL', sep='')
# 
#     summaryAAproperties(local_dir, stage, 'HC_PB_DNA', outn)
#     summaryAAproperties(local_dir, stage, 'CIS_PB_DNA', outn)
#     summaryAAproperties(local_dir, stage, 'RIS_PB_DNA', outn)
#     summaryAAproperties(local_dir, stage, 'ADVANCING', outn)
#     summaryAAproperties(local_dir, stage, 'STABLE', outn)
#     summaryAAproperties(local_dir, stage, 'CORD_POS', outn)
#     summaryAAproperties(local_dir, stage, 'CORD_NEG', outn)
#     summaryAAproperties(local_dir, stage, 'RIS_ADVANCING', outn)
#     summaryAAproperties(local_dir, stage, 'RIS_STABLE', outn)
#     summaryAAproperties(local_dir, stage, 'RIS_POS', outn)
#     summaryAAproperties(local_dir, stage, 'RIS_NEG', outn)
#     summaryAAproperties(local_dir, stage, 'RIS_CSRHI', outn)
#     summaryAAproperties(local_dir, stage, 'RIS_CSRLO', outn)
#     summaryAAproperties(local_dir, stage, 'CIS_ADVANCING', outn)
#     summaryAAproperties(local_dir, stage, 'CIS_POS', outn)
#     summaryAAproperties(local_dir, stage, 'CIS_NEG', outn)
# }

local_dir = 'analysis_v3/stats/ighj6/ighv4/'
stage = '.ighv4.ighj6.gene.mutations.aa_properties'
outn = '.IGHV4.IGHJ6.ALL'
summaryAAproperties(local_dir, stage, 'HC_PB_DNA', outn)
summaryAAproperties(local_dir, stage, 'CIS_PB_DNA', outn)
summaryAAproperties(local_dir, stage, 'RIS_PB_DNA', outn)
summaryAAproperties(local_dir, stage, 'ADVANCING', outn)
summaryAAproperties(local_dir, stage, 'STABLE', outn)
summaryAAproperties(local_dir, stage, 'CORD_POS', outn)
summaryAAproperties(local_dir, stage, 'CORD_NEG', outn)
summaryAAproperties(local_dir, stage, 'RIS_ADVANCING', outn)
summaryAAproperties(local_dir, stage, 'RIS_STABLE', outn)
summaryAAproperties(local_dir, stage, 'RIS_POS', outn)
summaryAAproperties(local_dir, stage, 'RIS_NEG', outn)
summaryAAproperties(local_dir, stage, 'RIS_CSRHI', outn)
summaryAAproperties(local_dir, stage, 'RIS_CSRLO', outn)
summaryAAproperties(local_dir, stage, 'CIS_ADVANCING', outn)
summaryAAproperties(local_dir, stage, 'CIS_POS', outn)
summaryAAproperties(local_dir, stage, 'CIS_NEG', outn)


# local_dir = 'analysis_v3/gt3shm/'
# stage = '.gt3shm.mutations.aa_properties'
# outn = '.GT3SHM'
# summaryAAproperties(local_dir, stage, 'HC_PB_DNA', outn)
# summaryAAproperties(local_dir, stage, 'CIS_PB_DNA', outn)
# summaryAAproperties(local_dir, stage, 'RIS_PB_DNA', outn)
# summaryAAproperties(local_dir, stage, 'ADVANCING', outn)
# summaryAAproperties(local_dir, stage, 'STABLE', outn)
# summaryAAproperties(local_dir, stage, 'CORD_POS', outn)
# summaryAAproperties(local_dir, stage, 'CORD_NEG', outn)
# summaryAAproperties(local_dir, stage, 'RIS_ADVANCING', outn)
# summaryAAproperties(local_dir, stage, 'RIS_STABLE', outn)
# summaryAAproperties(local_dir, stage, 'RIS_POS', outn)
# summaryAAproperties(local_dir, stage, 'RIS_NEG', outn)
# summaryAAproperties(local_dir, stage, 'RIS_CSRHI', outn)
# summaryAAproperties(local_dir, stage, 'RIS_CSRLO', outn)
# summaryAAproperties(local_dir, stage, 'CIS_ADVANCING', outn)
# summaryAAproperties(local_dir, stage, 'CIS_POS', outn)
# summaryAAproperties(local_dir, stage, 'CIS_NEG', outn)
# 
# for (i in 1:length(ighj.names)) {
#     print(ighj.names[i])
#     local_dir = paste('analysis_v3/gt3shm/', tolower(ighj.names[i]), '/', sep='')
#     stage = paste('.', tolower(ighj.names[i]), '.gt3shm.mutations.aa_properties', sep='')
#     outn = paste('.', ighj.names[i], '.GT3SHM', sep='')
# 
#     summaryAAproperties(local_dir, stage, 'HC_PB_DNA', outn)
#     summaryAAproperties(local_dir, stage, 'CIS_PB_DNA', outn)
#     summaryAAproperties(local_dir, stage, 'RIS_PB_DNA', outn)
#     summaryAAproperties(local_dir, stage, 'ADVANCING', outn)
#     summaryAAproperties(local_dir, stage, 'STABLE', outn)
#     summaryAAproperties(local_dir, stage, 'CORD_POS', outn)
#     summaryAAproperties(local_dir, stage, 'CORD_NEG', outn)
#     summaryAAproperties(local_dir, stage, 'RIS_ADVANCING', outn)
#     summaryAAproperties(local_dir, stage, 'RIS_STABLE', outn)
#     summaryAAproperties(local_dir, stage, 'RIS_POS', outn)
#     summaryAAproperties(local_dir, stage, 'RIS_NEG', outn)
#     summaryAAproperties(local_dir, stage, 'RIS_CSRHI', outn)
#     summaryAAproperties(local_dir, stage, 'RIS_CSRLO', outn)
#     summaryAAproperties(local_dir, stage, 'CIS_ADVANCING', outn)
#     summaryAAproperties(local_dir, stage, 'CIS_POS', outn)
#     summaryAAproperties(local_dir, stage, 'CIS_NEG', outn)
# }
# 

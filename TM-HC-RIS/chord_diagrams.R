# TM/HC/RIS
# libraries 1, 4, 6, 8, 9

# VJ chord diagrams

library(circlize)
library(rjson)
library(tidyr)

data_dir = '/Users/s166813/Projects/data/immune/MonsonLab/TM-HC-RIS/vdjserver/'
#data_dir = '/Users/scottc/Projects/data/immune/MonsonLab/TM-HC-RIS/vdjserver/'

#lib_dir = paste(data_dir, 'analysis_v2/stats/', sep='')
#file_prefix = 'Fig6.all.'
lib_dir = paste(data_dir, 'analysis_v2/gt3shm/', sep='')
file_prefix = 'Fig6.gt3shm.'

#lib_dir = paste(data_dir, 'analysis/0shm/', sep='')
#lib_dir = paste(data_dir, 'analysis/gt3shm/', sep='')
# old
#lib_dir = paste(data_dir, 'library_all/stats/', sep='')
#lib_dir = paste(data_dir, 'library_all/0shm/', sep='')
#lib_dir = paste(data_dir, 'library_all/gt3shm/', sep='')

metadata <- fromJSON(file = paste(data_dir,'repertoires.v2.airr.json',sep=''))

ighv.names = c('IGHV1', 'IGHV2', 'IGHV3', 'IGHV4', 'IGHV5', 'IGHV6', 'IGHV7')
ighj.names = c('IGHJ1', 'IGHJ2', 'IGHJ3', 'IGHJ4', 'IGHJ5', 'IGHJ6')

#exclude_repertoires = c()

#nrep = length(metadata$Repertoire)
#for (i in 1:nrep) {
#    rep_id = metadata$Repertoire[[i]]$repertoire_id
#    library = metadata$Repertoire[[i]]$sample[[1]]$sequencing_run_id
#    print(library)
    #ighv.seg = read.table(paste(data_dir,'allele.germ.clone.airr.productive.v_gene_usage.csv',sep=''), header=T, sep=',')
#}

# load raw tables
#filename = 'vj_combo_HC.pdf'
#liball = read.table(paste(lib_dir, 'HC.group.vj_combo.tsv',sep=''), header=T, sep='\t')
#filename = 'vj_combo_HC_NB.pdf'
#liball = read.table(paste(lib_dir, 'HC_NB.group.vj_combo.tsv',sep=''), header=T, sep='\t')
#filename = 'vj_combo_HC_MB.pdf'
#liball = read.table(paste(lib_dir, 'HC_MB.group.vj_combo.tsv',sep=''), header=T, sep='\t')
#filename = 'vj_combo_HC_PB.pdf'
#filename = 'vj_combo_HC_PB_0shm.pdf'
#filename = 'vj_combo_HC_PB_gt3shm.pdf'
#liball = read.table(paste(lib_dir, 'HC_PB.group.vj_combo.tsv',sep=''), header=T, sep='\t')
#filename = 'Fig6_HC_PB_DNA_vj_combo_.pdf'
#filename = 'seq_Fig4_vj_combo_HC_PB_DNA.pdf'
#filename = 'seq_Fig4_vj_combo_HC_PB_DNA_0shm.pdf'
#filename = 'seq_Fig4_vj_combo_HC_PB_DNA_gt3shm.pdf'
#filename = 'Fig4_IGHV4_combo_HC_PB_DNA_gt3shm.pdf'
#liball = read.table(paste(lib_dir, 'HC_PB_DNA.group.vj_combo.tsv',sep=''), header=T, sep='\t')

#filename = 'vj_combo_RIS.pdf'
#liball = read.table(paste(lib_dir, 'RIS.group.vj_combo.tsv',sep=''), header=T, sep='\t')
#filename = 'vj_combo_RIS_NB.pdf'
#liball = read.table(paste(lib_dir, 'RIS_NB.group.vj_combo.tsv',sep=''), header=T, sep='\t')
#filename = 'vj_combo_RIS_MB.pdf'
#liball = read.table(paste(lib_dir, 'RIS_MB.group.vj_combo.tsv',sep=''), header=T, sep='\t')
#filename = 'vj_combo_RIS_PB.pdf'
#filename = 'vj_combo_RIS_PB_0shm.pdf'
#filename = 'vj_combo_RIS_PB_gt3shm.pdf'
#liball = read.table(paste(lib_dir, 'RIS_PB.group.vj_combo.tsv',sep=''), header=T, sep='\t')
#filename = 'vj_combo_RIS_PB_DNA.pdf'
#filename = 'seq_Fig4_vj_combo_RIS_PB_DNA.pdf'
#filename = 'seq_Fig4_vj_combo_RIS_PB_DNA_0shm.pdf'
#filename = 'seq_Fig4_vj_combo_RIS_PB_DNA_gt3shm.pdf'
#liball = read.table(paste(lib_dir, 'RIS_PB_DNA.group.vj_combo.tsv',sep=''), header=T, sep='\t')

#filename = 'vj_combo_TM.pdf'
#liball = read.table(paste(lib_dir, 'TM.group.vj_combo.tsv',sep=''), header=T, sep='\t')
#filename = 'vj_combo_TM_NB.pdf'
#liball = read.table(paste(lib_dir, 'TM_NB.group.vj_combo.tsv',sep=''), header=T, sep='\t')
#filename = 'vj_combo_TM_MB.pdf'
#liball = read.table(paste(lib_dir, 'TM_MB.group.vj_combo.tsv',sep=''), header=T, sep='\t')
#filename = 'vj_combo_TM_PB.pdf'
#filename = 'vj_combo_TM_PB_0shm.pdf'
#filename = 'vj_combo_TM_PB_gt3shm.pdf'
#liball = read.table(paste(lib_dir, 'TM_PB.group.vj_combo.tsv',sep=''), header=T, sep='\t')
#filename = 'vj_combo_TM_PB_DNA.pdf'
#filename = 'seq_Fig4_vj_combo_TM_PB_DNA.pdf'
#filename = 'seq_Fig4_vj_combo_TM_PB_DNA_0shm.pdf'
#filename = 'seq_Fig4_vj_combo_TM_PB_DNA_gt3shm.pdf'
#liball = read.table(paste(lib_dir, 'TM_PB_DNA.group.vj_combo.tsv',sep=''), header=T, sep='\t')

#group_name = 'HC_PB_DNA'
group_name = 'CIS_PB_DNA'
#group_name = 'RIS_PB_DNA'
#group_name = 'RIS_ADVANCING'
#group_name = 'RIS_STABLE'
#group_name = 'RIS_POS'
#group_name = 'RIS_NEG'
#v1#group_name = 'TM_PB_DNA'

processing_stage = '.gt3shm.mutations'
#processing_stage = '.gene.mutations'
liball = read.table(paste(lib_dir, group_name, processing_stage, '.group.vj_combo.tsv',sep=''), header=T, sep='\t')

# extract the data
level = 'subgroup|subgroup'
mode = 'proportion'
productive = 'TRUE'
combo = liball[liball$level == level & liball$mode == mode & liball$productive == productive & (liball$v_level %in% ighv.names) & (liball$j_level %in% ighj.names),]
#combo = liball.RIS[liball.RIS$level == level & liball.RIS$mode == mode & liball.RIS$productive == productive,]
#combo = liball.TM[liball.TM$level == level & liball.TM$mode == mode & liball.TM$productive == productive,]

filename = paste(file_prefix, group_name, "_vj_combo", ".pdf", sep='')

# pull out VH4 subsets
if (FALSE) {
    combo = combo[combo$v_level == 'IGHV4',]
    # normalize subset
    combo$sequence_frequency_avg = combo$sequence_frequency_avg / sum(combo$sequence_frequency_avg)
    combo$duplicate_frequency_avg = combo$duplicate_frequency_avg / sum(combo$duplicate_frequency_avg)
}


# abundance-based
df = data.frame(from=combo$v_level, to=combo$j_level, value=combo$duplicate_frequency_avg)
# sequence-based
#df = data.frame(from=combo$v_level, to=combo$j_level, value=combo$sequence_frequency_avg)
#df = data.frame(from=combo$j_level, to=combo$v_level, value=combo$sequence_frequency_avg)

grid.col = c(IGHV1="chocolate", IGHV2="steelblue", IGHV3="forestgreen", IGHV4="cyan", IGHV5="magenta", IGHV6="wheat", IGHV7="salmon")
grid.col = c(grid.col, IGHJ1="chocolate", IGHJ2="steelblue", IGHJ3="forestgreen", IGHJ4="cyan", IGHJ5="magenta", IGHJ6="wheat")

pdf(file=filename)
chordDiagram(df[order(df[,"from"], df[,"to"]),], grid.col=grid.col, annotationTrack = c("grid","axis"), 
    preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(df))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
        facing = "clockwise", niceFacing = TRUE, adj = c(-0.2, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
dev.off()

# VH gene usage statistical tests

data_dir = '/Users/s166813/Projects/data/immune/MonsonLab/TM-HC-RIS/vdjserver/'

sub_dir = 'analysis_v3/stats/'
lib_dir = lib_dir = paste(data_dir, sub_dir, sep='')

g1 = read.table(paste(lib_dir, 'Fig4_HC_PB_DNA', '.VHgene.all.sequence_frequency', '.usage_table.csv',sep=''), header=T, sep=',')
g2 = read.table(paste(lib_dir, 'Fig4_RIS_PB_DNA', '.VHgene.all.sequence_frequency', '.usage_table.csv',sep=''), header=T, sep=',')
g3 = read.table(paste(lib_dir, 'Fig4_CIS_PB_DNA', '.VHgene.all.sequence_frequency', '.usage_table.csv',sep=''), header=T, sep=',')
g4 = read.table(paste(lib_dir, 'Fig4_ADVANCING', '.VHgene.all.sequence_frequency', '.usage_table.csv',sep=''), header=T, sep=',')
g5 = read.table(paste(lib_dir, 'Fig4_STABLE', '.VHgene.all.sequence_frequency', '.usage_table.csv',sep=''), header=T, sep=',')

ighv4family = c("IGHV4.28","IGHV4.30.2","IGHV4.30.4","IGHV4.31","IGHV4.34","IGHV4.38.2","IGHV4.39","IGHV4.4","IGHV4.55","IGHV4.59","IGHV4.61")

# HC vs RIS vs CIS

df1 = data.frame(log(g1[,ighv4family]), cohort='HC')
df2 = data.frame(log(g2[,ighv4family]), cohort='RIS')
df3 = data.frame(log(g3[,ighv4family]), cohort='CIS')
df = rbind(df1, df2, df3)

gr1 = data.frame()
for (i in ighv4family) {
    x = oneway.test(as.formula(paste(i, " ~ cohort")), df)
    gr1 = rbind(gr1, data.frame(gene=i, p=x$p.value, adjusted_p=x$p.value * length(ighv4family), avg=mean(df[,i], na.rm=T), n=sum(!is.na(df[,i]))))
}
print("HC vs RIS vs CIS")
print(gr1)

df = rbind(df2, df3)
res = oneway.test(IGHV4.34 ~ cohort, df)
print(res)
res = t.test(IGHV4.34 ~ cohort, df)
print(res)
res = t.test(IGHV4.39 ~ cohort, df)
print(res)


# HC vs MS2025

e2 = g2
exclude_subjects = c("UTSW03", "UTSW06", "UTSW21", "UTSW27", "UTSW28")
for (i in exclude_subjects) {
    e2 = e2[e2$subject_id != i,]
}
e3 = g3
exclude_subjects = c("1842", "1985", "2229", "2260", "2364", "2405", "3166")
for (i in exclude_subjects) {
    e3 = e3[e3$subject_id != i,]
}

df1 = data.frame(log(g1[,ighv4family]), cohort='HC')
df2 = data.frame(log(e2[,ighv4family]), cohort='MS20205')
df3 = data.frame(log(e3[,ighv4family]), cohort='MS20205')
df = rbind(df1, df2, df3)

gr2 = data.frame()
for (i in ighv4family) {
    x = oneway.test(as.formula(paste(i, " ~ cohort")), df)
    gr2 = rbind(gr2, data.frame(gene=i, p=x$p.value, adjusted_p=x$p.value * length(ighv4family), avg=mean(df[,i], na.rm=T), n=sum(!is.na(df[,i]))))
}
print("")
print("")
print("HC vs MS2025")
print(gr2)


# ADVANCING vs STABLE

df2 = data.frame(log(g4[,ighv4family]), cohort='ADVANCING')
df3 = data.frame(log(g5[,ighv4family]), cohort='STABLE')
df = rbind(df1, df2, df3)

gr3 = data.frame()
for (i in ighv4family) {
    x = oneway.test(as.formula(paste(i, " ~ cohort")), df)
    gr3 = rbind(gr3, data.frame(gene=i, p=x$p.value, adjusted_p=x$p.value * length(ighv4family), avg=mean(df[,i], na.rm=T), n=sum(!is.na(df[,i]))))
}
print("")
print("")
print("ADVANCING vs STABLE")
print(gr3)

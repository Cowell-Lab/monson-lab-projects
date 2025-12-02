# VH|JH gene combo statistical tests

data_dir = '/Users/s166813/Projects/data/immune/MonsonLab/TM-HC-RIS/vdjserver/'

sub_dir = 'analysis_v3/stats/'
lib_dir = lib_dir = paste(data_dir, sub_dir, sep='')

g1 = read.table(paste(lib_dir, 'Fig6_HC_PB_DNA', '.VJcombo.all.sequence_frequency', '.combo_table.csv',sep=''), header=T, sep=',')
g2 = read.table(paste(lib_dir, 'Fig6_RIS_PB_DNA', '.VJcombo.all.sequence_frequency', '.combo_table.csv',sep=''), header=T, sep=',')
g3 = read.table(paste(lib_dir, 'Fig6_CIS_PB_DNA', '.VJcombo.all.sequence_frequency', '.combo_table.csv',sep=''), header=T, sep=',')
g4 = read.table(paste(lib_dir, 'Fig6_ADVANCING', '.VJcombo.all.sequence_frequency', '.combo_table.csv',sep=''), header=T, sep=',')
g5 = read.table(paste(lib_dir, 'Fig6_STABLE', '.VJcombo.all.sequence_frequency', '.combo_table.csv',sep=''), header=T, sep=',')

ighv4family = c("IGHV4.IGHJ1","IGHV4.IGHJ2","IGHV4.IGHJ3","IGHV4.IGHJ4","IGHV4.IGHJ5","IGHV4.IGHJ6")

# RIS within group
comp1 = c("IGHV4.IGHJ1","IGHV4.IGHJ2","IGHV4.IGHJ3","IGHV4.IGHJ4","IGHV4.IGHJ5")
comp2 = c("IGHV4.IGHJ6")
df1 = data.frame(value=log(rowMeans(g2[,comp1], na.rm=T)), cohort='RIS')
df2 = data.frame(value=log(g2[,comp2]), cohort='V6')
df = rbind(df1, df2)
res = oneway.test(value ~ cohort, df)
#print(df1)
#print(df2)
#print(res)
res = t.test(value ~ cohort, df)
print(res)
#data.frame(wg1, avg=rowMeans(wg1[,comp1], na.rm=T))
#m1 = mean(


# HC vs RIS vs CIS

df1 = data.frame(log(g1[,ighv4family]), cohort=as.factor('HC'))
df2 = data.frame(log(g2[,ighv4family]), cohort=as.factor('RIS'))
df3 = data.frame(log(g3[,ighv4family]), cohort=as.factor('CIS'))
df = rbind(df1, df2, df3)

gr1 = data.frame()
for (i in ighv4family) {
    x = oneway.test(as.formula(paste(i, " ~ cohort")), df)
    gr1 = rbind(gr1, data.frame(combo=i, p=x$p.value, adjusted_p=x$p.value * length(ighv4family), avg=mean(df[,i], na.rm=T), n=sum(!is.na(df[,i]))))
}
print("HC vs RIS vs CIS")
print(gr1)

# post-hoc
# IGHV4.IGHJ6
x = aov(IGHV4.IGHJ6 ~ cohort, df)
post_test <- glht(x, linfct = mcp(cohort = "Tukey"))
print(summary(post_test))

df = rbind(df1, df2)
res = t.test(IGHV4.IGHJ6 ~ cohort, df)
print(res)
df = rbind(df1, df3)
res = t.test(IGHV4.IGHJ6 ~ cohort, df)
print(res)
df = rbind(df2, df3)
res = t.test(IGHV4.IGHJ6 ~ cohort, df)
print(res)

# IGHV4.IGHJ4
x = aov(IGHV4.IGHJ4 ~ cohort, df)
post_test <- glht(x, linfct = mcp(cohort = "Tukey"))
print(summary(post_test))


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
    gr2 = rbind(gr2, data.frame(combo=i, p=x$p.value, adjusted_p=x$p.value * length(ighv4family), avg=mean(df[,i], na.rm=T), n=sum(!is.na(df[,i]))))
}
print("HC vs MS2025")
print(gr2)


# ADVANCING vs STABLE

df2 = data.frame(log(g4[,ighv4family]), cohort='ADVANCING')
df3 = data.frame(log(g5[,ighv4family]), cohort='STABLE')
df = rbind(df1, df2, df3)

gr3 = data.frame()
for (i in ighv4family) {
    x = oneway.test(as.formula(paste(i, " ~ cohort")), df)
    gr3 = rbind(gr3, data.frame(combo=i, p=x$p.value, adjusted_p=x$p.value * length(ighv4family), avg=mean(df[,i], na.rm=T), n=sum(!is.na(df[,i]))))
}
print("ADVANCING vs STABLE")
print(gr3)

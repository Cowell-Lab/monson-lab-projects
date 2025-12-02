# Mutation statistical tests

data_dir = '/Users/s166813/Projects/data/immune/MonsonLab/TM-HC-RIS/vdjserver/'

#sub_dir = 'analysis_v3/stats/'
#stage = '.gene.mutations'
sub_dir = 'analysis_v3/stats_ighv4/'
stage = '.ighv4.mutations.aa_properties'
lib_dir = lib_dir = paste(data_dir, sub_dir, sep='')

g1 = read.table(paste(lib_dir, 'Fig8_HC_PB_DNA', stage, '.repertoire.frequency', '.mutational_report.csv',sep=''), header=T, sep=',')
g2 = read.table(paste(lib_dir, 'Fig8_RIS_PB_DNA', stage, '.repertoire.frequency', '.mutational_report.csv',sep=''), header=T, sep=',')
g3 = read.table(paste(lib_dir, 'Fig8_CIS_PB_DNA', stage, '.repertoire.frequency', '.mutational_report.csv',sep=''), header=T, sep=',')
g4 = read.table(paste(lib_dir, 'Fig8_ADVANCING', stage, '.repertoire.frequency', '.mutational_report.csv',sep=''), header=T, sep=',')
g5 = read.table(paste(lib_dir, 'Fig8_STABLE', stage, '.repertoire.frequency', '.mutational_report.csv',sep=''), header=T, sep=',')

origfields = c("mu_freq","mu_freq_r","mu_freq_fwr1_s","mu_freq_fwr1_r","mu_freq_fwr2_s","mu_freq_fwr2_r","mu_freq_fwr3_s","mu_freq_fwr3_r","mu_freq_cdr1_s","mu_freq_cdr1_r","mu_freq_cdr2_s","mu_freq_cdr2_r","mu_cdrrsratio","mu_fwrrsratio")

# HC vs RIS vs CIS

df1 = data.frame(g1[,origfields], cohort='HC')
df2 = data.frame(g2[,origfields], cohort='RIS')
df3 = data.frame(g3[,origfields], cohort='CIS')
df = rbind(df1, df2, df3)

mutfields = c("mu_freq", "mu_freq_r", "mu_freq_fwr", "mu_freq_fwr_r", "mu_freq_cdr", "mu_freq_cdr_r")

df['mu_freq'] = log(df[,'mu_freq'])
df['mu_freq_r'] = log(df[,'mu_freq_r'])
df['mu_freq_fwr'] = log((df[,'mu_freq_fwr1_s'] + df[,'mu_freq_fwr1_r'] + df[,'mu_freq_fwr2_s'] + df[,'mu_freq_fwr2_r'] + df[,'mu_freq_fwr3_s'] + df[,'mu_freq_fwr3_r']) / 6)
df['mu_freq_fwr_r'] = log((df[,'mu_freq_fwr1_r'] + df[,'mu_freq_fwr2_r'] + df[,'mu_freq_fwr3_r']) / 3)
df['mu_freq_cdr'] = log((df[,'mu_freq_cdr1_s'] + df[,'mu_freq_cdr1_r'] + df[,'mu_freq_cdr2_s'] + df[,'mu_freq_cdr2_r']) / 4)
df['mu_freq_cdr_r'] = log((df[,'mu_freq_cdr1_r'] + df[,'mu_freq_cdr2_r']) / 2)

gr1 = data.frame()
for (i in mutfields) {
    x = oneway.test(as.formula(paste(i, " ~ cohort")), df)
    gr1 = rbind(gr1, data.frame(mut=i, p=x$p.value, adjusted_p=x$p.value * length(mutfields), log_avg=mean(df[,i], na.rm=T), avg=mean(exp(df[,i]), na.rm=T), n=sum(!is.na(df[,i]))))
}
print("HC vs RIS vs CIS")
print(gr1)

ratiofields = c("mu_fwrrsratio", "mu_cdrrsratio")

gr2 = data.frame()
for (i in ratiofields) {
    x = oneway.test(as.formula(paste(i, " ~ cohort")), df)
    gr2 = rbind(gr2, data.frame(mut=i, p=x$p.value, adjusted_p=x$p.value * length(ratiofields), avg=mean(df[,i], na.rm=T), n=sum(!is.na(df[,i]))))
}
print("HC vs RIS vs CIS")
print(gr2)

for (i in mutfields) {
    gr = df[df$cohort == 'RIS' | df$cohort == 'HC',]
    res = oneway.test(as.formula(paste(i, " ~ cohort")), gr)
    print("HC vs RIS")
    print(res)
    print(res$p.value * 3)

    gr = df[df$cohort == 'CIS' | df$cohort == 'HC',]
    res = oneway.test(as.formula(paste(i, " ~ cohort")), gr)
    print("HC vs CIS")
    print(res)
    print(res$p.value * 3)

    gr = df[df$cohort == 'RIS' | df$cohort == 'CIS',]
    res = oneway.test(as.formula(paste(i, " ~ cohort")), gr)
    print("RIS vs CIS")
    print(res)
    print(res$p.value * 3)
}

for (i in ratiofields) {
    gr = df[df$cohort == 'RIS' | df$cohort == 'HC',]
    res = oneway.test(as.formula(paste(i, " ~ cohort")), gr)
    print("HC vs RIS")
    print(res)
    print(res$p.value * 3)

    gr = df[df$cohort == 'CIS' | df$cohort == 'HC',]
    res = oneway.test(as.formula(paste(i, " ~ cohort")), gr)
    print("HC vs CIS")
    print(res)
    print(res$p.value * 3)

    gr = df[df$cohort == 'RIS' | df$cohort == 'CIS',]
    res = oneway.test(as.formula(paste(i, " ~ cohort")), gr)
    print("RIS vs CIS")
    print(res)
    print(res$p.value * 3)
}

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

df1 = data.frame(g1[,origfields], cohort='HC')
df2 = data.frame(e2[,origfields], cohort='MS2025')
df3 = data.frame(e3[,origfields], cohort='MS2025')
df = rbind(df1, df2, df3)

df['mu_freq'] = log(df[,'mu_freq'])
df['mu_freq_r'] = log(df[,'mu_freq_r'])
df['mu_freq_fwr'] = log((df[,'mu_freq_fwr1_s'] + df[,'mu_freq_fwr1_r'] + df[,'mu_freq_fwr2_s'] + df[,'mu_freq_fwr2_r'] + df[,'mu_freq_fwr3_s'] + df[,'mu_freq_fwr3_r']) / 6)
df['mu_freq_fwr_r'] = log((df[,'mu_freq_fwr1_r'] + df[,'mu_freq_fwr2_r'] + df[,'mu_freq_fwr3_r']) / 3)
df['mu_freq_cdr'] = log((df[,'mu_freq_cdr1_s'] + df[,'mu_freq_cdr1_r'] + df[,'mu_freq_cdr2_s'] + df[,'mu_freq_cdr2_r']) / 4)
df['mu_freq_cdr_r'] = log((df[,'mu_freq_cdr1_r'] + df[,'mu_freq_cdr2_r']) / 2)

gr1 = data.frame()
for (i in mutfields) {
    x = oneway.test(as.formula(paste(i, " ~ cohort")), df)
    gr1 = rbind(gr1, data.frame(mut=i, p=x$p.value, adjusted_p=x$p.value * length(mutfields), log_avg=mean(df[,i], na.rm=T), avg=mean(exp(df[,i]), na.rm=T), n=sum(!is.na(df[,i]))))
}
print("")
print("")
print("HC vs MS2025")
print(gr1)

ratiofields = c("mu_fwrrsratio", "mu_cdrrsratio")

gr2 = data.frame()
for (i in ratiofields) {
    x = oneway.test(as.formula(paste(i, " ~ cohort")), df)
    gr2 = rbind(gr2, data.frame(mut=i, p=x$p.value, adjusted_p=x$p.value * length(ratiofields), avg=mean(df[,i], na.rm=T), n=sum(!is.na(df[,i]))))
}
print("HC vs MS2025")
print(gr2)

gr = df[df$cohort == 'MS2025',]
print(mean(exp(gr$mu_freq)))
print(mean(exp(gr$mu_freq_cdr)))
print(mean(exp(gr$mu_freq_fwr)))
print(mean(exp(gr$mu_freq_r)))
print(mean(exp(gr$mu_freq_cdr_r)))
print(mean(exp(gr$mu_freq_fwr_r)))
print(mean(gr$mu_cdrrsratio))
print(mean(gr$mu_fwrrsratio))

# ADVANCING vs STABLE

df1 = data.frame(g4[,origfields], cohort='ADVANCING')
df2 = data.frame(g5[,origfields], cohort='STABLE')
df = rbind(df1, df2)

mutfields = c("mu_freq", "mu_freq_r", "mu_freq_fwr", "mu_freq_fwr_r", "mu_freq_cdr", "mu_freq_cdr_r")

df['mu_freq'] = log(df[,'mu_freq'])
df['mu_freq_r'] = log(df[,'mu_freq_r'])
df['mu_freq_fwr'] = log((df[,'mu_freq_fwr1_s'] + df[,'mu_freq_fwr1_r'] + df[,'mu_freq_fwr2_s'] + df[,'mu_freq_fwr2_r'] + df[,'mu_freq_fwr3_s'] + df[,'mu_freq_fwr3_r']) / 6)
df['mu_freq_fwr_r'] = log((df[,'mu_freq_fwr1_r'] + df[,'mu_freq_fwr2_r'] + df[,'mu_freq_fwr3_r']) / 3)
df['mu_freq_cdr'] = log((df[,'mu_freq_cdr1_s'] + df[,'mu_freq_cdr1_r'] + df[,'mu_freq_cdr2_s'] + df[,'mu_freq_cdr2_r']) / 4)
df['mu_freq_cdr_r'] = log((df[,'mu_freq_cdr1_r'] + df[,'mu_freq_cdr2_r']) / 2)

gr1 = data.frame()
for (i in mutfields) {
    x = oneway.test(as.formula(paste(i, " ~ cohort")), df)
    gr1 = rbind(gr1, data.frame(mut=i, p=x$p.value, adjusted_p=x$p.value * length(mutfields), log_avg=mean(df[,i], na.rm=T), avg=mean(exp(df[,i]), na.rm=T), n=sum(!is.na(df[,i]))))
}
print("ADVANCING vs STABLE")
print(gr1)

ratiofields = c("mu_fwrrsratio", "mu_cdrrsratio")

gr2 = data.frame()
for (i in ratiofields) {
    x = oneway.test(as.formula(paste(i, " ~ cohort")), df)
    gr2 = rbind(gr2, data.frame(mut=i, p=x$p.value, adjusted_p=x$p.value * length(ratiofields), avg=mean(df[,i], na.rm=T), n=sum(!is.na(df[,i]))))
}
print("ADVANCING vs STABLE")
print(gr2)

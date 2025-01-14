#!/usr/bin/env Rscript

# Alakazam gene usage
#
# Author: Scott Christley
# Date: Sep 3, 2020
# 

# based upon this script for parsing args with optparse
# https://bitbucket.org/kleinstein/immcantation/src/master/pipelines/shazam-threshold.R

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("alakazam"))
suppressPackageStartupMessages(library("airr"))

# Define commmandline arguments
opt_list <- list(make_option(c("-d", "--db"), dest="DB",
                             help="Tabulated data file, in AIRR format (TSV)."),
                 make_option(c("-o", "--output"), dest="OUTFILE",
                             help="output filename prefix"))

# Parse arguments
opt <- parse_args(OptionParser(option_list=opt_list))

# Check input file
if (!("DB" %in% names(opt))) {
    stop("You must provide a database file with the -d option.")
}

# Check output file
if (!("OUTFILE" %in% names(opt))) {
    stop("You must provide an output filename prefix with the -o option.")
}

# Read rearrangement data
db <- airr::read_rearrangement(opt$DB)

# allele
# v call
genes <- countGenes(db, gene='v_call', mode='allele', copy='duplicate_count')
write.table(genes, row.names=F, sep='\t', file=paste(opt$OUTFILE, '.all.v_allele_usage.tsv', sep=''))
genes <- countGenes(db, gene='v_call', mode='allele', copy='duplicate_count', group='productive')
write.table(genes[genes$productive,], row.names=F, sep='\t', file=paste(opt$OUTFILE, '.productive.v_allele_usage.tsv', sep=''))
write.table(genes[!genes$productive,], row.names=F, sep='\t', file=paste(opt$OUTFILE, '.unproductive.v_allele_usage.tsv', sep=''))

# d call
genes <- countGenes(db, gene='d_call', mode='allele', copy='duplicate_count')
write.table(genes, row.names=F, sep='\t', file=paste(opt$OUTFILE, '.all.d_allele_usage.tsv', sep=''))
genes <- countGenes(db, gene='d_call', mode='allele', copy='duplicate_count', group='productive')
write.table(genes[genes$productive,], row.names=F, sep='\t', file=paste(opt$OUTFILE, '.productive.d_allele_usage.tsv', sep=''))
write.table(genes[!genes$productive,], row.names=F, sep='\t', file=paste(opt$OUTFILE, '.unproductive.d_allele_usage.tsv', sep=''))

# j call
genes <- countGenes(db, gene='j_call', mode='allele', copy='duplicate_count')
write.table(genes, row.names=F, sep='\t', file=paste(opt$OUTFILE, '.all.j_allele_usage.tsv', sep=''))
genes <- countGenes(db, gene='j_call', mode='allele', copy='duplicate_count', group='productive')
write.table(genes[genes$productive,], row.names=F, sep='\t', file=paste(opt$OUTFILE, '.productive.j_allele_usage.tsv', sep=''))
write.table(genes[!genes$productive,], row.names=F, sep='\t', file=paste(opt$OUTFILE, '.unproductive.j_allele_usage.tsv', sep=''))

# TODO: Alakazam throws an error and stops execution if no data in c_call field
# We need to figure out to check for this
#genes <- countGenes(db, gene='c_call', group='repertoire_id', mode='allele', copy='duplicate_count', fill=T)
#write.table(genes, row.names=F, sep='\t', file='c_allele_usage.tsv')

# gene
genes <- countGenes(db, gene='v_call', mode='gene', copy='duplicate_count')
write.table(genes, row.names=F, sep='\t', file=paste(opt$OUTFILE, '.all.v_gene_usage.tsv', sep=''))
genes <- countGenes(db, gene='v_call', mode='gene', copy='duplicate_count', group='productive')
write.table(genes[genes$productive,], row.names=F, sep='\t', file=paste(opt$OUTFILE, '.productive.v_gene_usage.tsv', sep=''))
write.table(genes[!genes$productive,], row.names=F, sep='\t', file=paste(opt$OUTFILE, '.unproductive.v_gene_usage.tsv', sep=''))

genes <- countGenes(db, gene='d_call', mode='gene', copy='duplicate_count')
write.table(genes, row.names=F, sep='\t', file=paste(opt$OUTFILE, '.all.d_gene_usage.tsv', sep=''))
genes <- countGenes(db, gene='d_call', mode='gene', copy='duplicate_count', group='productive')
write.table(genes[genes$productive,], row.names=F, sep='\t', file=paste(opt$OUTFILE, '.productive.d_gene_usage.tsv', sep=''))
write.table(genes[!genes$productive,], row.names=F, sep='\t', file=paste(opt$OUTFILE, '.unproductive.d_gene_usage.tsv', sep=''))

genes <- countGenes(db, gene='j_call', mode='gene', copy='duplicate_count')
write.table(genes, row.names=F, sep='\t', file=paste(opt$OUTFILE, '.all.j_gene_usage.tsv', sep=''))
genes <- countGenes(db, gene='j_call', mode='gene', copy='duplicate_count', group='productive')
write.table(genes[genes$productive,], row.names=F, sep='\t', file=paste(opt$OUTFILE, '.productive.j_gene_usage.tsv', sep=''))
write.table(genes[!genes$productive,], row.names=F, sep='\t', file=paste(opt$OUTFILE, '.unproductive.j_gene_usage.tsv', sep=''))

#genes <- countGenes(db, gene='c_call', mode='gene', copy='duplicate_count')
#write.table(genes, row.names=F, sep='\t', file='c_gene_usage.tsv')

# family/subgroup
genes <- countGenes(db, gene='v_call', mode='family', copy='duplicate_count')
write.table(genes, row.names=F, sep='\t', file=paste(opt$OUTFILE, '.all.v_subgroup_usage.tsv', sep=''))
genes <- countGenes(db, gene='v_call', mode='family', copy='duplicate_count', group='productive')
write.table(genes[genes$productive,], row.names=F, sep='\t', file=paste(opt$OUTFILE, '.productive.v_subgroup_usage.tsv', sep=''))
write.table(genes[!genes$productive,], row.names=F, sep='\t', file=paste(opt$OUTFILE, '.unproductive.v_subgroup_usage.tsv', sep=''))

genes <- countGenes(db, gene='d_call', mode='family', copy='duplicate_count')
write.table(genes, row.names=F, sep='\t', file=paste(opt$OUTFILE, '.all.d_subgroup_usage.tsv', sep=''))
genes <- countGenes(db, gene='d_call', mode='family', copy='duplicate_count', group='productive')
write.table(genes[genes$productive,], row.names=F, sep='\t', file=paste(opt$OUTFILE, '.productive.d_subgroup_usage.tsv', sep=''))
write.table(genes[!genes$productive,], row.names=F, sep='\t', file=paste(opt$OUTFILE, '.unproductive.d_subgroup_usage.tsv', sep=''))

genes <- countGenes(db, gene='j_call', mode='family', copy='duplicate_count')
write.table(genes, row.names=F, sep='\t', file=paste(opt$OUTFILE, '.all.j_subgroup_usage.tsv', sep=''))
genes <- countGenes(db, gene='j_call', mode='family', copy='duplicate_count', group='productive')
write.table(genes[genes$productive,], row.names=F, sep='\t', file=paste(opt$OUTFILE, '.productive.j_subgroup_usage.tsv', sep=''))
write.table(genes[!genes$productive,], row.names=F, sep='\t', file=paste(opt$OUTFILE, '.unproductive.j_subgroup_usage.tsv', sep=''))


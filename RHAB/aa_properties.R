#!/usr/bin/env Rscript

# Alakazam AA properties
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

aa_db <- aminoAcidProperties(db, seq="cdr3_aa", nt=FALSE, label="cdr3")
#aa_db <- aminoAcidProperties(aa_db, seq="junction_aa", trim=TRUE, nt=FALSE, label="junction")

airr::write_rearrangement(aa_db, file=paste(opt$OUTFILE, '.aa_properties.airr.tsv', sep=''))

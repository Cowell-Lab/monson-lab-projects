#
# This is MonsonLab specific for "single cell" Sanger
# sequencing data analysis workflow. Takes the output from
# IgBlast and performs some additional calculations and
# generates summary TSV.
#
# Author: Scott Christley
# Date: Oct 13, 2021
#
# This script assumes the docker environment and has hard-coded paths

# input file
NAME=$1
if [[ "x$NAME" == "x" ]] ; then
    echo "Need to specify NAME as parameter."
    exit 1
fi

# calculate CDR3 properties
/work/vdjserver/irplus/tapis/vdjserver/common/aa_properties.R -d ${NAME}.igblast.makedb.airr.tsv -o ${NAME}.igblast.makedb.aa_properties.airr.tsv

# germline alignments with IMGT gaps and D mask
CreateGermlines.py -d ${NAME}.igblast.makedb.aa_properties.airr.tsv -g dmask --failed -r /work/vdjserver/db/db.2019.01.23/germline/human/ReferenceDirectorySet/IG_VDJ.fna

# mutations
/work/vdjserver/irplus/tapis/vdjserver/common/mutational_analysis.R -d ${NAME}.igblast.makedb.aa_properties.airr_germ-pass.tsv -o ${NAME}

# summary TSV
python3 /work/research/immune/code/MonsonLab/RHAB/rhab_summary_v3.py ${NAME}.mutations.airr.tsv ${NAME}.igblast.airr.tsv ${NAME}.summary.tsv

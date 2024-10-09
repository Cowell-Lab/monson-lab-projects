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

# assign repertoire_id
repertoire_id=MonsonLab_RHAB_${NAME}
data_processing_id=MonsonLab_RHAB_${NAME}
metadata_file=${repertoire_id}.airr.json
python3 /work/vdjserver/irplus/tapis/vdjserver/common/assign_repertoire_id.py ${repertoire_id} ${data_processing_id} ${NAME}.igblast.makedb.airr.tsv ${NAME}.igblast.makedb.new.airr.tsv
# generate fake metadata file
python3 /work/vdjserver/irplus/tapis/vdjserver/common/create_single_metadata.py ${repertoire_id} ${metadata_file}

# calculate CDR3 properties
/work/vdjserver/irplus/tapis/vdjserver/common/aa_properties.R -d ${NAME}.igblast.makedb.new.airr.tsv -o ${NAME}.igblast.makedb.aa_properties.airr.tsv

# germline alignments with IMGT gaps and D mask
CreateGermlines.py -d ${NAME}.igblast.makedb.aa_properties.airr.tsv -g dmask --failed -r /work/vdjserver/db/db.2019.01.23/germline/human/ReferenceDirectorySet/IG_VDJ.fna

# mutations
/work/vdjserver/irplus/tapis/vdjserver/common/mutational_analysis.R -d ${NAME}.igblast.makedb.aa_properties.airr_germ-pass.tsv -o ${NAME}

# repcalc annotations
# generate config
python3 /work/vdjserver/irplus/tapis/vdjserver/common/repcalc_create_config.py --init /work/vdjserver/irplus/tapis/vdjserver/common/mutational_template.json ${metadata_file} ${repertoire_id}.mutational_config.json
python3 /work/vdjserver/irplus/tapis/vdjserver/common/repcalc_create_config.py --rearrangementFile ${NAME}.mutations.airr.tsv --repertoireID $repertoire_id ${repertoire_id}.mutational_config.json

# run it
repcalc ${repertoire_id}.mutational_config.json

# summary TSV
python3 /work/research/immune/code/MonsonLab/RHAB/rhab_summary_v3.py ${repertoire_id}.mutations.airr.tsv ${NAME}.igblast.airr.tsv ${NAME}.summary.tsv
python3 /work/research/immune/code/MonsonLab/RHAB/rhab_summary_v4.py ${repertoire_id}.mutations.airr.tsv ${NAME}.igblast.airr.tsv ${repertoire_id}.summary.tsv

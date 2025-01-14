
# simple script to walk through the clone/rearrangements files and calculate statistics
# Assumes running in the immcantation docker
# Need to install RepCalc V2 in the docker

metadata_file=../../repertoires.v2.airr.json
#group_file=../../repertoire_groups.airr.json
germline_db=/work/vdjserver/db/db.2019.01.23/germline/vdjserver_human_germline.airr.json
common_dir=/work/vdjserver/irplus/tapis/vdjserver/common

# The Tapis script assumes some programs are in the current directory so copy them
cp ${common_dir}/common_functions.sh .
cp ${common_dir}/clonal_abundance.R .
cp ${common_dir}/repcalc_create_config.py .
cp ${common_dir}/gene_usage_template.json .
cp ${common_dir}/gene_combo_template.json .

#fileList=(278-5_S5_R1_001.fastq.merged.total-BC6.unique.igblast.airr.tsv)
fileList=($(ls *.gene.clone.airr.tsv))
#fileList=($(ls *.ighv4.ge3ags.mutations.airr.tsv))

count=0
while [ "x${fileList[count]}" != "x" ]
do
    file=${fileList[count]}
    echo $file
    fileOutname="${file##*/}"
    #fileBasename="${fileOutname%%.*}"
    # assume and strip airr.tsv
    fileBasename="${fileOutname%.*}"
    fileBasename="${fileBasename%.*}"
    echo $fileBasename

    # Clonal abundance
    Rscript ./clonal_abundance.R -d $file -o $fileBasename

    count=$(( $count + 1 ))
done

fileList=($(ls *.allele.clone.airr.tsv))
if [ "x${group_file}" != "x" ]; then
    python3 repcalc_create_config.py --init gene_usage_template.json ${metadata_file} --groups ${group_file} --germline ${germline_db} gene_usage_config.json
    python3 repcalc_create_config.py --init gene_combo_template.json ${metadata_file} --groups ${group_file} --germline ${germline_db} gene_combo_config.json
else
    python3 repcalc_create_config.py --init gene_usage_template.json ${metadata_file} --germline ${germline_db} gene_usage_config.json
    python3 repcalc_create_config.py --init gene_combo_template.json ${metadata_file} --germline ${germline_db} gene_combo_config.json
fi

count=0
while [ "x${fileList[count]}" != "x" ]
do
    file=${fileList[count]}
    echo $file
    fileOutname="${file##*/}"
    #fileBasename="${fileOutname%%.*}"
    # assume and strip airr.tsv
    fileBasename="${fileOutname%.*}"
    fileBasename="${fileBasename%.*}"
    echo $fileBasename

    # Clonal abundance
    Rscript ./clonal_abundance.R -d $file -o $fileBasename

    # Repcalc config
    fileOrig="${fileBasename%.*}"
    fileOrig="${fileOrig%.*}"
    fileOrig=${fileOrig}.airr.tsv.gz
    echo $fileOrig
    repertoire_id=$(python3 ${common_dir}/get_repertoire_id_for_file.py ${metadata_file} $fileOrig)
    echo $repertoire_id
    python3 repcalc_create_config.py --rearrangementFile $file --repertoireID $repertoire_id gene_usage_config.json
    python3 repcalc_create_config.py --rearrangementFile $file --repertoireID $repertoire_id gene_combo_config.json

    count=$(( $count + 1 ))
done

repcalc gene_usage_config.json
repcalc gene_combo_config.json

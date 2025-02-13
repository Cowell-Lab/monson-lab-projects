
# simple script to walk through the clone files and calculation mutations
# Assumes running in the immcantation docker

fileList=($(ls *.clone.airr.tsv))
#fileList=($(ls 2070_S19_L001_R1_001*.clone.airr.tsv))
#fileList=(1068_S35_L001_R1_001.fastq.merged.unique.igblast.makedb.allele.clone.airr.tsv)

metadata_file=../../repertoires.v2.airr.json
#metadata_file=../repertoires.airr.json
germline_fasta=/work/vdjserver/db/db.2019.01.23/germline/human/ReferenceDirectorySet/IG_VDJ.fna
common_dir=/work/vdjserver/vdjserver-agave/common
#common_dir=/work/vdjserver/vdjserver-tapis/scripts


# The Tapis script assumes some programs are in the current directory so copy them
cp ${common_dir}/common_functions.sh .
cp ${common_dir}/create_germlines.sh .
cp ${common_dir}/mutational_analysis.R .
cp ${common_dir}/mutational_analysis.sh .
cp ${common_dir}/mutational_template.json .
cp ${common_dir}/get_repertoire_id_for_file.py .
cp ${common_dir}/repcalc_create_config.py .
cp ${common_dir}/parse_changeo.py .

    count=0
    while [ "x${fileList[count]}" != "x" ]
    do
        file=${fileList[count]}
        echo $file
        fileOutname="${file##*/}"
        echo $fileOutname

        # save filenames for later processing
        # assume and strip airr.tsv
        fileExtname="${file%.*}" # file.clone.airr.tsv -> file.clone.airr
        airrFilename="${fileExtname%.*}" # file.clone.airr -> file.clone
        baseFilename="${airrFilename%.*}" # file.clone -> file
        echo $baseFilename

        fileOrig="${baseFilename%.*}" # strip the allele/gene
        fileOrig=${fileOrig}.airr.tsv.gz
        echo $fileOrig
        repertoire_id=$(python3 ./get_repertoire_id_for_file.py ${metadata_file} $fileOrig)
        echo $repertoire_id
        if [ "x${repertoire_id}" == "x" ]; then
            echo "Cannot determine repertoire_id for $file"
            exit 1
        fi

        bash ./create_germlines.sh ${fileOutname} ${airrFilename} ${germline_fasta}
        germFilename="${airrFilename}.germ.airr.tsv"

        bash mutational_analysis.sh ${metadata_file} ${repertoire_id} ${germFilename} ${baseFilename}

        count=$(( $count + 1 ))
    done

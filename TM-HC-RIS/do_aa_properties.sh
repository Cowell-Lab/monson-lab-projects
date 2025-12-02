
# simple script to walk through the clone/rearrangements files and calculate CDR3 AA properties
# Assumes running in the immcantation docker
# Need to install RepCalc V2 in the docker

metadata_file=../../repertoires.v2.airr.json
#group_file=../../repertoire_groups.airr.json
germline_db=/work/vdjserver/db/db.2019.01.23/germline/vdjserver_human_germline.airr.json
common_dir=/work/vdjserver/vdjserver-tapis/scripts

# Create environment similar to Tapis app
cp ${common_dir}/common_functions.sh .
cp ${common_dir}/aa_properties.R .

processing_stage=gt3shm.mutations
#processing_stage=gene.mutations
fileList=($(ls *.${processing_stage}.airr.tsv))
#fileList=($(ls *.allele.mutations.airr.tsv))

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

    # calculate CDR3 properties
    Rscript ./aa_properties.R -d ${file} -o ${fileBasename}.aa_properties.airr.tsv

    count=$(( $count + 1 ))
done

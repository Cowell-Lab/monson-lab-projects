
# simple script to walk through the clone/rearrangements files and calculate statistics
# Assumes running in the immcantation docker
# Need to install RepCalc V2 in the docker

metadata_file=../../repertoires.v2.airr.json
group_file=../../repertoire_groups.airr.json
germline_db=/work/vdjserver/db/db.2019.01.23/germline/vdjserver_human_germline.airr.json
common_dir=/work/vdjserver/irplus/tapis/vdjserver/common

# The Tapis script assumes some programs are in the current directory so copy them
cp ${common_dir}/common_functions.sh .
cp ${common_dir}/clonal_abundance.R .

processing_stage=gt3shm.mutations
fileList=($(ls *.${processing_stage}.airr.tsv))

count=0
while [ "x${fileList[count]}" != "x" ]
do
    file=${fileList[count]}
    echo $file
    fileOutname="${file##*/}"
    # assume and strip airr.tsv
    fileBasename="${fileOutname%.*}"
    fileBasename="${fileBasename%.*}"
    echo $fileBasename

    # Clonal abundance
    Rscript ./clonal_abundance.R -d $file -o $fileBasename

    count=$(( $count + 1 ))
done


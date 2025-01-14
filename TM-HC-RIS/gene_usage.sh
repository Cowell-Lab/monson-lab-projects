
# simple script to walk through the clone/rearrangements files and calculate statistics
# Assumes running in the immcantation docker

common_dir=/work/vdjserver/irplus/tapis/vdjserver/common

# The Tapis script assumes some programs are in the current directory so copy them
cp ${common_dir}/common_functions.sh .
cp ${common_dir}/gene_usage.R .
cp ${common_dir}/clonal_abundance.R .
cp ${common_dir}/clone_report.py .

fileList=($(ls *.clone.airr.tsv))
#fileList=(278-5_S5_R1_001.fastq.merged.total-BC6.unique.igblast.airr.tsv)

count=0
while [ "x${fileList[count]}" != "x" ]
do
    file=${fileList[count]}
    echo $file
    fileOutname="${file##*/}"
    fileBasename="${fileOutname%%.*}"
    echo $fileBasename

    # Clonal abundance
    Rscript ./clonal_abundance.R -d $file -o $fileBasename

    count=$(( $count + 1 ))
done

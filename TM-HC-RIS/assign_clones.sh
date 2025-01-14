
# simple script to walk through the IgBlast MakeDB files from VDJServer and assign clone
# Assumes running in the immcantation docker

fileList=($(ls *.igblast.makedb.airr.tsv))
#fileList=(2070_S19_L001_R1_001.fastq.merged.unique.igblast.makedb.airr.tsv)

common_dir=/work/vdjserver/vdjserver-agave/common

# The Tapis script assumes some programs are in the current directory so copy them
cp ${common_dir}/changeo_clones.sh .
cp ${common_dir}/common_functions.sh .
cp ${common_dir}/find_threshold.R .
cp ${common_dir}/parse_changeo.py .
cp ${common_dir}/clone_report.py .

    cloneFileList=()
    count=0
    while [ "x${fileList[count]}" != "x" ]
    do
        file=${fileList[count]}
        echo $file
        fileOutname="${file##*/}"
        echo $fileOutname
        #noArchive $fileOutname

        # save filenames for later processing
        # assume and strip airr.tsv
        fileBasename="${fileOutname%.*}"
        fileBasename="${fileBasename%.*}"
        echo $fileBasename
        cloneFileList[${#cloneFileList[@]}]=$fileBasename

        bash ./changeo_clones.sh $file $fileBasename 1

        count=$(( $count + 1 ))
    done

    #generate clone report
    python3 ./clone_report.py *.igblast.makedb.airr.tsv

    echo ${cloneFileList}

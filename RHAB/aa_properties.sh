    
    fileList=($(ls *.igblast.airr.tsv))

    count=0
    while [ "x${fileList[count]}" != "x" ]
    do
        file=${fileList[count]}
        echo $file
        fileOutname="${file##*/}"
        #echo $fileOutname
        #noArchive $fileOutname

        # save filenames for later processing
        fileBasename="${fileOutname%.*}"
        fileBasename="${fileBasename%.*}"
        echo $fileBasename

        /work/research/immune/code/MonsonLab/RHAB/aa_properties.R -d $file -o $fileBasename

        count=$(( $count + 1 ))
    done

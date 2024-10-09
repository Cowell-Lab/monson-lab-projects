    
    fileList=($(ls *.germ.airr.tsv))

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

        # mutations
        /work/research/immune/code/MonsonLab/RHAB/mutational_analysis.R -d $file -o $fileBasename

        count=$(( $count + 1 ))
    done

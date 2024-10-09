    
    fileList=($(ls *.seq))
    outfile=$1
    if [ -z "$1" ]; then
        echo Provide output file name
        exit 1
    fi

    rm -f $outfile
    count=0
    while [ "x${fileList[count]}" != "x" ]
    do
        file=${fileList[count]}
        #echo $file
        fileOutname="${file#*_}"
        #echo $fileOutname
        seqid="${fileOutname%%_*}"
        echo $seqid

        echo ">${seqid}" >> $outfile
        cat $file >> $outfile

        count=$(( $count + 1 ))
    done

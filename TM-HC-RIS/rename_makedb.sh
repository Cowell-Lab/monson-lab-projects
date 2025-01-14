
# manually assigning repertoire_id creates a new.airr.tsv, this script swaps new and original files

fileList=($(ls *.igblast.makedb.new.airr.tsv))
#fileList=(1395-2_S2_R1_001.fastq.merged.unique.igblast.makedb.new.airr.tsv)

    count=0
    while [ "x${fileList[count]}" != "x" ]
    do
        file=${fileList[count]}
        echo $file
        fileOutname="${file##*/}"
        echo $fileOutname
        #noArchive $fileOutname

        # save filenames for later processing
        # strip new.airr.tsv
        fileBasename="${fileOutname%.*}"
        fileBasename="${fileBasename%.*}"
        fileBasename="${fileBasename%.*}"
        echo $fileBasename
        mv ${fileBasename}.airr.tsv ${fileBasename}.orig.airr.tsv
        mv $file ${fileBasename}.airr.tsv

        count=$(( $count + 1 ))
    done

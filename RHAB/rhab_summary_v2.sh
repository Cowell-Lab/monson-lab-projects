    
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
        fileBasename="${fileBasename%.*}"
        echo $fileBasename
        aa_file=${fileBasename}.igblast.aa_properties.airr.tsv
        mut_file=${fileBasename}.igblast.germ.mutations.airr.tsv
        out_file=${fileBasename}.summary.airr.tsv

        python3 /work/research/immune/code/MonsonLab/RHAB/rhab_summary_v2.py $mut_file $aa_file $out_file

        count=$(( $count + 1 ))
    done

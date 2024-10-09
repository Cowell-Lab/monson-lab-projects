    
    fileList=($(ls *.igblast.makedb.airr.tsv))

    organism=human
    seq_type=IG
    vdj_db='/work/vdjserver/db/db.2019.01.23/germline/'${organism}'/ReferenceDirectorySet/'${seq_type}'_VDJ.fna'

    count=0
    while [ "x${fileList[count]}" != "x" ]
    do
        file=${fileList[count]}
        echo $file
        fileOutname="${file##*/}"
        echo $fileOutname
        #noArchive $fileOutname

        # save filenames for later processing
        fileBasename="${fileOutname%.*}"
        fileBasename="${fileBasename%.*}"
        fileName="${fileBasename}.airr_germ-pass.tsv"
        failName="${fileBasename}.airr_germ-fail.tsv"

        fileBasename="${fileBasename%.*}"
        echo $fileBasename
        echo $fileName

        # create germlines
        CreateGermlines.py -d ${file} -r $vdj_db -g dmask --failed | tee germ.log
        newPass="${fileBasename}.germ.airr.tsv"
        newFail="${fileBasename}.germ-fail.airr.tsv"
        mv $fileName $newPass
        mv $failName $newFail
        summaryName="${fileBasename}.germ.airr.json"
        python3 /work/vdjserver/vdjserver-agave/common/parse_changeo.py germ.log $summaryName

        count=$(( $count + 1 ))
    done

    echo ${cloneFileList}

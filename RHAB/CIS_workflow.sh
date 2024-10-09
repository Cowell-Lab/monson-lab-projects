    
    fileList=($(ls *.igblast.airr.tsv))
    #fileList=(278-5_S5_R1_001.fastq.merged.total-BC6.unique.igblast.airr.tsv)
    
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
        fileBasename="${fileOutname%.*}"
        fileBasename="${fileBasename%.*}"
        fileName="${fileBasename}.airr_clone-pass.tsv"
        failName="${fileBasename}.airr_clone-fail.tsv"
        echo $fileBasename
        echo $fileName
        echo $failName
        cloneFileList[${#cloneFileList[@]}]=$fileName

        # find the threshold that DefineClones needs
        rm -f threshold.dat
        (time -p /work/vdjserver/irplus/tapis/apps/changeo/0.4/common/find_threshold.R -d $file -o threshold.dat) 2> timing.dat
        cat threshold.dat
        threshold="0.16"
        if [ -f threshold.dat ]; then
            threshold=$(cat threshold.dat)
            rm -f threshold.dat
        fi

        # allele mode
        tmpLogFile=clone.allele.log
        echo DefineClones.py -d $file --mode allele --act set --model ham --norm len --dist $threshold --failed
        DefineClones.py -d $file --mode allele --act set --model ham --norm len --dist $threshold --failed | tee $tmpLogFile
        newPass="${fileBasename}.allele.clone.airr.tsv"
        newFail="${fileBasename}.allele.clone-fail.airr.tsv"
        mv $fileName $newPass
        if [ -f $failName ]; then
            mv $failName $newFail
        fi
        summaryName="${fileBasename}.allele.summary.clone.airr.json"
        python3 /work/vdjserver/irplus/tapis/apps/changeo/0.4/common/parse_changeo.py $tmpLogFile $summaryName
        rm -f $tmpLogFile
        
        # gene mode
        tmpLogFile=clone.gene.log
        echo DefineClones.py -d $file --mode gene --act set --model ham --norm len --dist $threshold --failed
        DefineClones.py -d $file --mode gene --act set --model ham --norm len --dist $threshold --failed | tee $tmpLogFile
        newPass="${fileBasename}.gene.clone.airr.tsv"
        newFail="${fileBasename}.gene.clone-fail.airr.tsv"
        mv $fileName $newPass
        if [ -f $failName ]; then
            mv $failName $newFail
        fi
        summaryName="${fileBasename}.gene.summary.clone.airr.json"
        python3 /work/vdjserver/irplus/tapis/apps/changeo/0.4/common/parse_changeo.py $tmpLogFile $summaryName
        rm -f $tmpLogFile

        count=$(( $count + 1 ))
    done

    echo ${cloneFileList}

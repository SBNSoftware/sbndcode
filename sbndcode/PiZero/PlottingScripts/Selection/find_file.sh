#!/bin/bash

FILE_LIST=/pnfs/sbnd/scratch/users/hlay/ncpizero/NCPiZeroAv2/rockbox/reco1/files.list

RUN_NUM=$1
SUBRUN_NUM=$2

echo "Looking for run $RUN_NUM subrun $SUBRUN_NUM"

for FILE in `cat $FILE_LIST`
do
    if [[ "$FILE" == *"4709159_0"* ]];
    then
        continue
    fi

    echo $FILE

    META=`sam_metadata_dumper $FILE`

    while read -r METALINE;
    do
        if [[ "$METALINE" == *"runs"* ]];
        then
            echo $METALINE
            continue
            SETS=`echo $METALINE | cut -d : -f 2`
            SETS=${SETS//[[:blank:]]/}
            SETS=${SETS#"["}
            SETS=${SETS%"],"}
            for SET in ${SETS//],/] }
            do
                SET=${SET#"["}
                SET=${SET%",\"physics\"]"}
                FILE_RUN_NUM=`echo $SET | cut -d ',' -f 1`
                FILE_SUBRUN_NUM=`echo $SET | cut -d ',' -f 2`

                if [[ $FILE_RUN_NUM -eq $RUN_NUM && $FILE_SUBRUN_NUM -eq $SUBRUN_NUM ]];
                then
                    echo $FILE
                    return
                fi
            done
        fi
    done < <(echo "$META")
done

h#!/bin/bash

mkdir tmp
cd tmp

export DISPLAY=localhost:1

FILE_LIST=/pnfs/sbnd/scratch/users/hlay/ncpizero/NCPiZeroAv6/rockbox/reco2/files.list

while read LINE;
do
    RUN_NUM=`echo $LINE | cut -d ' ' -f 1`
    SUBRUN_NUM=`echo $LINE | cut -d ' ' -f 2`
    EV_NUM=`echo $LINE | cut -d ' ' -f 3`

    echo "Looking for run $RUN_NUM subrun $SUBRUN_NUM event $EV_NUM"

    for FILE in `cat $FILE_LIST`
    do
        META=`sam_metadata_dumper $FILE`

        while read -r METALINE;
        do
            if [[ "$METALINE" == *"runs"* ]];
            then
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
                        lar -c event_display_and_dump.fcl -s $FILE -n 1 -e $RUN_NUM:$SUBRUN_NUM:$EV_NUM
                        break 3
                    fi
                done
            fi
        done < <(echo "$META")
    done
done < <(cat $1)

cd ../

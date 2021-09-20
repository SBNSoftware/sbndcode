#!/usr/bin/env bash

echo "%%%------------------------------%%%"
echo "%  Executing fcl_checks.sh script  %"
echo "%%%------------------------------%%%"

#########################################################################
### Establish environment and read in the arguments from the cfg file ###
#########################################################################

CWD="$(pwd)"
WORK_DIR="${CWD}/fcl_tests"

echo -e "\nWorking Directory: ${WORK_DIR}"

declare x=$@

while :
do
    case "x$1" in
	x--refdir)               REF_DIR="${2}";                         shift; shift;;
	x--update-ref-files)     UPDATE_REF_FILE_ON=1;                   shift;;
	x)                                                               break;;
    esac
done


################################################
### Create working directory and output file ###
################################################

rm -Rf "$WORK_DIR"
mkdir "$WORK_DIR"
cd "$WORK_DIR"


######################################
### Loop over fcl files for checks ###
######################################

ACCESS_REF_DIR=${REF_DIR///pnfs//cvmfs/sbnd.osgstorage.org/pnfs/fnal.gov/usr}
exit_code_parsing=0
exit_code_dump=0

echo -e "\nReference Directory: ${ACCESS_REF_DIR}"
echo "ls ${ACCESS_REF_DIR}"
eval ifdh ls $ACCESS_REF_DIR


echo -e "\nList of files from: ${SBNDCODE_DIR}/test/fcl_file_checks.list"
echo "$(cat ${SBNDCODE_DIR}/test/fcl_file_checks.list)"
echo

export datestamp=$(date +"%Y%m%d%H%M")

if [[ ${UPDATE_REF_FILE_ON} -gt 0 ]]; then
    echo -e "Updating ref files for fcl checks\n"
fi

for fcl in `cat ${SBNDCODE_DIR}/test/fcl_file_checks.list`

do
    echo -e "\nTesting fcl file ${fcl}"

    fcl_dir="${WORK_DIR}/${fcl%.fcl}"
    mkdir "$fcl_dir"
    cd "$fcl_dir"


    #######################################
    ### Perform debug check on fcl file ###
    #######################################

    fclout=${fcl%.fcl}_fhicl_dump.out
    larout=lar.out
    larerr=lar.err
    lar -c ${SBNDCODE_DIR}/fcl/${fcl} --debug-config $fclout > $larout 2> $larerr


    ##############################################
    ### Check fcl could be parsed successfully ###
    ##############################################

    stat=$?

    if [[ $stat -ne 0 && $stat -ne 1 ]]; then
        echo "Error parsing ${fcl} --- Exited with status ${stat}"
        echo "==================================================================================================================================="
        cat ${larerr}
        echo "==================================================================================================================================="
        let exit_code_parsing=201
        continue
    else
        echo "Successful parsing of ${fcl} --- Exited with status ${stat}"
    fi

    ############################################
    ### Update reference file (if requested) ###
    ############################################

    if [[ ${UPDATE_REF_FILE_ON} -gt 0 ]]; then
        echo "Removing ${REF_DIR}/${fcl%.fcl}_fhicl_dump.out"
        ifdh rm ${REF_DIR}/${fcl%.fcl}_fhicl_dump.out
        echo "Copying: ${fclout} to ${REF_DIR}/old/${fcl%.fcl}_fhicl_dump_${datestamp}.out"
        ifdh cp ${fclout} ${REF_DIR}/old/${fcl%.fcl}_fhicl_dump_${datestamp}.out
        echo "Copying: ${fclout} to ${REF_DIR}/${fcl%.fcl}_fhicl_dump.out"
        ifdh cp ${fclout} ${REF_DIR}/${fcl%.fcl}_fhicl_dump.out
        continue
    fi

    ###############################################################################
    ### Check the diff between the fcl configs for the ref and current versions ###
    ###############################################################################

    if [[ ! -f ${ACCESS_REF_DIR}/${fcl%.fcl}_fhicl_dump.out ]]
    then
	echo "Cannot open reference file: ${ACCESS_REF_DIR}/${fcl%.fcl}_fhicl_dump.out it does not exist"
	let exit_code_dump=201
	continue
    fi

    fcl_diff=$(diff ${ACCESS_REF_DIR}/${fcl%.fcl}_fhicl_dump.out ${fclout})

    if [[ -n $fcl_diff ]]; then
        echo "Non-zero diff from file: ${fcl}"
        echo "==================================================================================================================================="
        echo "$fcl_diff"
        echo "==================================================================================================================================="
        let exit_code_dump=201
    else
        echo "Zero diff from file: ${fcl}"
    fi
    echo
done

##################################
### Prepare exit code and exit ###
##################################

exit_code=0

if [[ $exit_code_parsing -ne 0 && $exit_code_dump -ne 0 ]]; then
    ERRORSTRING="Error parsing fcl file and differences in fhicl dump outputs"
    let exit_code=201
elif [[ $exit_code_parsing -ne 0 ]]; then
    ERRORSTRING="Error parsing fcl file"
    let exit_code=201
elif [[ $exit_code_dump -ne 0 ]]; then
    ERRORSTRING="Differences in fhicl dump outputs"
    let exit_code=201
else
    ERRORSTRING="All checks successful"
fi

echo "Exiting with exit code: $exit_code"
echo "$ERRORSTRING"

exit $exit_code

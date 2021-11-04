#!/usr/bin/env bash

echo "%%%------------------------------%%%"
echo "%  Executing fcl_checks.sh script  %"
echo "%%%------------------------------%%%"

#########################################################################
### Establish environment and read in the arguments from the cfg file ###
#########################################################################

CWD="$(pwd)"
WORK_DIR="${CWD}/fcl_tests"
LOCAL_REF_DIR="${WORK_DIR}/references"

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


######################################################
### Create working and reference store directories ###
######################################################

rm -rf "${WORK_DIR}"
mkdir "${WORK_DIR}"
mkdir "${LOCAL_REF_DIR}"
cd "${WORK_DIR}"


######################################
### Loop over fcl files for checks ###
######################################

ACCESS_REF_DIR=${REF_DIR///pnfs//cvmfs/sbnd.osgstorage.org/pnfs/fnal.gov/usr}
exit_code_parsing=0
exit_code_dump=0

echo -e "\nReference Directory: ${ACCESS_REF_DIR}"
echo "ls ${ACCESS_REF_DIR}"
eval ifdh ls $ACCESS_REF_DIR

if [[ ! -f ${ACCESS_REF_DIR}/fhicl_dump_references.tar.gz ]]
then
    echo -e "\nCannot find reference tar: ${ACCESS_REF_DIR}/fhicl_dump_references.tar.gz"
    echo "Exiting with exit code: $exit_code"
    exit 1
else
    echo -e "\nFound reference tar: ${ACCESS_REF_DIR}/fhicl_dump_references.tar.gz"
    echo "cp ${ACCESS_REF_DIR}/fhicl_dump_references.tar.gz ${LOCAL_REF_DIR}/fhicl_dump_references.tar.gz"
    eval ifdh cp ${ACCESS_REF_DIR}/fhicl_dump_references.tar.gz ${LOCAL_REF_DIR}/fhicl_dump_references.tar.gz
    echo -e "\nExtract tar to local references directory"
    echo "tar -xzvf references/fhicl_dump_references.tar.gz -C ${LOCAL_REF_DIR}"
    eval tar -xzvf references/fhicl_dump_references.tar.gz -C ${LOCAL_REF_DIR}
fi

echo -e "\nList of reference files in: ${LOCAL_REF_DIR}"
echo "$(ls ${LOCAL_REF_DIR}/*)"

echo -e "\nList of files from: ${SBNDCODE_DIR}/test/fcl_file_checks.list"
echo "$(cat ${SBNDCODE_DIR}/test/fcl_file_checks.list)"

export datestamp=$(date +"%Y%m%d%H%M")

if [[ ${UPDATE_REF_FILE_ON} -gt 0 ]]; then
    echo -e "\nUpdating ref files for fcl checks"
fi

for fcl in `cat ${SBNDCODE_DIR}/test/fcl_file_checks.list`

do
    echo -e "\nTesting fcl file ${fcl}"


    #######################################
    ### Perform debug check on fcl file ###
    #######################################

    fclout=${fcl%.fcl}_fhicl_dump.out
    larout=${fcl%.fcl}_lar.out
    larerr=${fcl%.fcl}_lar.err
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


    ###############################################################################
    ### Check the diff between the fcl configs for the ref and current versions ###
    ###############################################################################

    if [[ ${UPDATE_REF_FILE_ON} -gt 0 ]]; then
	continue
    fi

    if [[ ! -f ${LOCAL_REF_DIR}/${fcl%.fcl}_fhicl_dump.out ]]
    then
	echo "Cannot open reference file: ${LOCAL_REF_DIR}/${fcl%.fcl}_fhicl_dump.out it does not exist"
	let exit_code_dump=201
	continue
    fi

    fcl_diff=$(diff ${LOCAL_REF_DIR}/${fcl%.fcl}_fhicl_dump.out ${fclout})

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


############################################
### Update reference file (if requested) ###
############################################

if [[ ${UPDATE_REF_FILE_ON} -gt 0 ]]; then

    cd ${WORK_DIR}
    echo -e "\nMaking tar of new output files"
    echo "tar -czvf fhicl_dump_references.tar.gz *_fhicl_dump.out"
    eval tar -czvf fhicl_dump_references.tar.gz *_fhicl_dump.out
    echo -e "Copy reference tar to pnfs"
    echo "ifdh cp fhicl_dump_references.tar.gz ${REF_DIR}/fhicl_dump_references_${datestamp}.tar.gz"
    eval ifdh cp fhicl_dump_references.tar.gz ${REF_DIR}/fhicl_dump_references_${datestamp}.tar.gz
    echo -e "Moving old reference tar to backup"
    echo "ifdh rename ${REF_DIR}/fhicl_dump_references.tar.gz ${REF_DIR}/fhicl_dump_references.tar.gz.bak"
    eval ifdh rename ${REF_DIR}/fhicl_dump_references.tar.gz ${REF_DIR}/fhicl_dump_references.tar.gz.bak
    echo -e "Copy new tar to reference name"
    echo "ifdh cp ${REF_DIR}/fhicl_dump_references_${datestamp}.tar.gz ${REF_DIR}/fhicl_dump_references.tar.gz"
    eval ifdh cp ${REF_DIR}/fhicl_dump_references_${datestamp}.tar.gz ${REF_DIR}/fhicl_dump_references.tar.gz
    
    continue
fi


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

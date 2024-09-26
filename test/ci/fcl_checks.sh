#!/usr/bin/env bash

echo "%%%------------------------------%%%"
echo "%  Executing $(basename ${0}) script  %"
echo "%%%------------------------------%%%"

echo "args: ${@}"

#########################################################################
### Establish environment and read in the arguments from the cfg file ###
#########################################################################

CWD="$(pwd)"

WORKSPACE=${WORKSPACE:-$PWD}
ci_cur_exp_name=${ci_cur_exp_name:-${EXPERIMENT}}

WORK_DIR="${CWD}/fcl_tests"
LOCAL_REF_DIR="${WORK_DIR}/references"

echo -e "\nWorking Directory: ${WORK_DIR}"

while :
do
    case "x$1" in
        x--refdir)               REF_DIR="${2}";                         shift; shift;;
        x--output-file)          CUROUTPUT="${2}";                       shift; shift;;
        x--input-file)           INPUT_FILE="${2}";                      shift; shift;;
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

ACCESS_REF_DIR=${REF_DIR///pnfs/http://fndcadoor.fnal.gov:8000/pnfs/fnal.gov/usr}
REF_FILE=${INPUT_FILE}

exit_code_parsing=0
exit_code_dump=0

echo -e "\nRemote Reference Directory: ${ACCESS_REF_DIR}"

if [[ ${UPDATE_REF_FILE_ON} -gt 0 ]]; then
    echo -e "\nUpdating ref files for fcl checks"
    export datestamp=$(date +"%Y%m%d%H%M")

else
    #if IFDH_DEBUG=0 ifdh ll ${ACCESS_REF_DIR}/${REF_FILE}
    if IFDH_PROXY_ENABLE=0 IFDH_TOKEN_ENABLE=0 BEARER_TOKEN_FILE="" ifdh ll ${ACCESS_REF_DIR}/${REF_FILE}
    then
	echo -e "\nFound reference tar: ${ACCESS_REF_DIR}/${REF_FILE}"
	echo "ifdh cp ${ACCESS_REF_DIR}/${REF_FILE} ${LOCAL_REF_DIR}/${REF_FILE}"
	IFDH_PROXY_ENABLE=0 IFDH_TOKEN_ENABLE=0 BEARER_TOKEN_FILE="" ifdh cp ${ACCESS_REF_DIR}/${REF_FILE} ${LOCAL_REF_DIR}/${REF_FILE}
	echo -e "\nExtract tar to local references directory"
	echo "tar -xzvf references/${REF_FILE} -C ${LOCAL_REF_DIR}"
	tar -xzvf references/${REF_FILE} -C ${LOCAL_REF_DIR}
    else
	exit_code=${?}
	echo -e "\nCannot find reference tar: ${ACCESS_REF_DIR}/${REF_FILE}"
	echo "Exiting with exit code: ${exit_code}"
	EXITSTATUS=F
	ERRORSTRING="Failure accessing reference tar file~Check the log"
	echo "`basename $CWD`~${exit_code}~${EXITSTATUS}~$ERRORSTRING" >> $WORKSPACE/data_production_stats${ci_cur_exp_name}.log
	exit $exit_code
    fi

    echo -e "\nList of reference files in: ${LOCAL_REF_DIR}"
    echo "$(ls ${LOCAL_REF_DIR}/*)"
fi

echo -e "\nList of files from: ${SBNDCODE_DIR}/test/fcl_file_checks.list"
echo "$(cat ${SBNDCODE_DIR}/test/fcl_file_checks.list)"

for fcl in `cat ${SBNDCODE_DIR}/test/fcl_file_checks.list`

do
    echo -e "\nTesting fcl file ${fcl}"


    #######################################
    ### Perform debug check on fcl file ###
    #######################################

    fclout=${fcl%.fcl}_fhicl_dump.out
    larout=${fcl%.fcl}_lar.out
    larerr=${fcl%.fcl}_lar.err
    lar -c ${fcl} --debug-config $fclout > $larout 2> $larerr


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
    echo "tar -czvf ${CWD}/${CUROUTPUT} *_fhicl_dump.out"
    tar -czvf ${CWD}/${CUROUTPUT} *_fhicl_dump.out
fi


##################################
### Prepare exit code and exit ###
##################################

exit_code=0

if [[ $exit_code_parsing -ne 0 && $exit_code_dump -ne 0 ]]; then
    ERRORSTRING="Error parsing fcl file and differences in fhicl dump outputs~Check error in log"
    let exit_code=201
elif [[ $exit_code_parsing -ne 0 ]]; then
    ERRORSTRING="Error parsing fcl file~Check error in log"
    let exit_code=201
elif [[ $exit_code_dump -ne 0 ]]; then
    ERRORSTRING="Differences in fhicl dump outputs~Update reference files"
    let exit_code=201
else
    ERRORSTRING="All checks successful"
fi

echo -e "\nExiting with exit code: $exit_code"
echo $ERRORSTRING | cut -d~ -f 1


if [[ $exit_code -ne 0 ]]; then
    EXITSTATUS=W
    echo "`basename $CWD`~${exit_code}~${EXITSTATUS}~$ERRORSTRING" >> $WORKSPACE/data_production_stats${ci_cur_exp_name}.log
fi

exit $exit_code

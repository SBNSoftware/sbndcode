#!/usr/bin/env bash


function usage {
      cat <<EOF
   usage: $0 [options]
      running CI tests for ${proj_PREFIX}_ci.
   options:
      --executable               Define the executable to run
      --nevents                  Define the number of events to process
      --stage                    Define the stage number used to parse the right testmask column number
      --fhicl                    Set the FHiCl file to use to run the test
      --input-file               Set the file on which you want to run the test
      --reference-files          Set the reference file used by the test, this can be a comma separated list
      --outputs                  Define a list of pairs <output_stream>:<output_filename> using "," as separator
      --stage-name               Define the name of the test
      --testmask                 Define the bit-mask to enable the different test phases
                                 (currently there are 3 test phases: data_production; compare_data_products; compare_data_product_size.)
      --update-ref-files         Flag to activate the "Update Reference Files" mode
      --input-files-to-fetch     List of input files to be downloaded before to execute the data production
      --reference-files-to-fetch List of reference files to be downloaded before the product comparison
      --extra-function           Define and extra function to run with list of required arguments; the elements need to be comma separated
      --extra-options            Define and extra options/arguments for the executable; the elements need to be comma separated
EOF
}

function initialize
{
    TASKSTRING="initialize"
    ERRORSTRING="F~Error initializing the test~Check the log"
    trap 'LASTERR=$?; FUNCTION_NAME=${FUNCNAME[0]:-main};  exitstatus ${LASTERR} trap ${LINENO}; exit ${LASTERR}' ERR

    echo "running CI tests for ${proj_PREFIX}_ci."
    echo "ci_cur_exp_name: ${ci_cur_exp_name}"
    echo
    echo "initialize $@"

    #~~~~~~~~~~~~~~~ DEFAULT VALUES ~~~~~~~~~~~~~~~~
    EXECUTABLE_NAME=no_executable_defined
    NEVENTS=1
    UPDATE_REF_FILE_ON=0
    INPUT_FILE=""
    INPUT_FILES_TO_FETCH=""
    REFERENCE_FILES=""
    REFERENCE_FILES_TO_FETCH=""
    WORKSPACE=${WORKSPACE:-$PWD}
    EXTRA_OPTIONS=""
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    #~~~~~~~~~~~~~~~~~~~~~~GET VALUE FROM THE CI_TESTS.CFG ARGS SECTION~~~~~~~~~~~~~~~
    while :
    do
      case "x$1" in
      x-h|x--help)                 usage;                                                       exit;;
      x--executable)               EXECUTABLE_NAME="${2}";                                      shift; shift;;
      x--nevents)                  NEVENTS="${2}";                                              shift; shift;;
      x--stage)                    STAGE="${2}";                                                shift; shift;;
      x--fhicl)                    FHiCL_FILE="${2}";                                           shift; shift;;
      x--input-file)               INPUT_FILE="${2}";                                           shift; shift;;
      x--reference-files)          REFERENCE_FILES="${2}";                                       shift; shift;;
      x--outputs)                  OUTPUT_LIST="${2}"; OUTPUT_STREAM="${OUTPUT_LIST//,/ -o }";  shift; shift;;
      x--stage-name)               STAGE_NAME="${2}";                                           shift; shift;;
      x--testmask)                 TESTMASK="${2}";                                             shift; shift;;
      x--update-ref-files)         UPDATE_REF_FILE_ON=1;                                        shift;;
      x--input-files-to-fetch)     INPUT_FILES_TO_FETCH="${2}";                                 shift; shift;;
      x--reference-files-to-fetch) REFERENCE_FILES_TO_FETCH="${2}";                             shift; shift;;
      x--extra-function)           EXTRA_FUNCTION="${2}";                                       shift; shift;;
      x--extra-options)            EXTRA_OPTIONS="${2//,/ }";                                   shift; shift;;
      x)                                                                                break;;
      x*)            echo "Unknown argument $1"; usage; exit 1;;
      esac
    done

    if [ ${UPDATE_REF_FILE_ON} -gt 0 ]; then
        echo -e "\n***************************************************"
        echo "This CI build is running to update reference files:"
        echo "- data product comparison is disabled"
        echo "- number of events is set to 1"
        echo "- existing reference files will not be used"
        echo -e "***************************************************\n"
        TESTMASK=""
        if [[ "$(basename ${0})" != *"lariatsoft"* ]]; then
            NEVENTS=1
        fi
        REFERENCE_FILES=""
        REFERENCE_FILES_TO_FETCH=""
    fi

    if [ -n "${INPUT_FILES_TO_FETCH}" ]; then
        fetch_files input ${INPUT_FILES_TO_FETCH}
    fi
    if [ -n "${REFERENCE_FILES_TO_FETCH}" ]; then
        fetch_files reference ${REFERENCE_FILES_TO_FETCH}
    fi

    #~~~~~~~~~~~~~~~~~~~~~PARSE THE TESTMASK FILE TO UNDERSTAND WHICH FUNCTION TO RUN ~~~~~~~~~~~~
    if [ -n "${TESTMASK}" ];then
        check_data_production=${TESTMASK:0:1}
        check_compare_names=$((${TESTMASK:1:1}&&${check_data_production}))
        check_compare_size=$((${TESTMASK:2:1}&&${check_compare_names}))
    else
        check_data_production=1
        check_compare_names=0
        check_compare_size=0
    fi

    echo "Input file:  ${INPUT_FILE}"
    echo "Output files: ${OUTPUT_LIST}"
    echo "FHiCL file:  ${FHiCL_FILE}"
    echo "Testmask: ${TESTMASK}"
    echo
    echo -e "\nRunning\n `basename $0` $@"

    exitstatus $?
}


function fetch_files
{
    old_taskstring="$TASKSTRING"
    old_errorstring="$ERRORSTRING"
    TASKSTRING="fetching $1 files"

    ERRORSTRING="F~Error in fetching $1 files~Check if the $1 files are available"

    echo "fetching $1 files for ${proj_PREFIX}_ci."
    echo
    echo "fetch_files $@"
    echo

    maxretries_backup=$IFDH_CP_MAXRETRIES
    debug_backup=$IFDH_DEBUG

    export IFDH_DEBUG=1
    export IFDH_CP_MAXRETRIES=0

    for file in ${2//,/ }
    do
        echo "Command: ifdh cp -D $file ./"
        ifdh cp -D $file ./ > fetch_inputs.log  2>&1
        local copy_exit_code=$?

        if [[ $copy_exit_code -ne 0 ]]; then
            echo "Failed to fetch $file"
            exitstatus 211
        fi
    done

    export IFDH_DEBUG=$debug_backup
    export IFDH_CP_MAXRETRIES=$maxretries_backup
    exitstatus $copy_exit_code
    TASKSTRING="$old_taskstring"
    ERRORSTRING="$old_errorstring"
}


function data_production
{
    TASKSTRING="data_production"
    ERRORSTRING="F~Error in data production~Check the log"
    trap 'LASTERR=$?; FUNCTION_NAME=${FUNCNAME[0]:-main};  exitstatus ${LASTERR} trap ${LINENO}; exit ${LASTERR}' ERR

    export TMPDIR=${PWD} #Temporary directory used by IFDHC

    #~~~~~~~~~~~~~IF THE TESTMASK VALUE IS SET TO 1 THEN RUN THE PRODUCTION OF THR DATA~~~~~~~~~~~~~~~~~~
    if [[ "${1}" -eq 1 ]]
    then
        # This is LArIAT specific.
        # The slicer stage has output filename in the form file_%#.root
        if [[ "${STAGE_NAME}" == "slicer" && "$(basename ${0})" == *"lariatsoft"* ]]; then  # This is LArIAT specific.
            OUTPUT_STREAM=${OUTPUT_STREAM//.root/_%#.root}
        fi

        echo -e "\nNumber of events for ${STAGE_NAME} stage: $NEVENTS\n"
        echo ${EXECUTABLE_NAME} --rethrow-all -n ${NEVENTS} ${EXTRA_OPTIONS} ${OUTPUT_STREAM:+-o ${OUTPUT_STREAM}} --config ${FHiCL_FILE} ${INPUT_FILE}
        echo

        (
            local counter=0
            local expcode_exitcode=20
            until [[ ${expcode_exitcode} -ne 20 || ${counter} -gt 5 ]]; do
                ${EXECUTABLE_NAME} --rethrow-all -n ${NEVENTS} ${EXTRA_OPTIONS} ${OUTPUT_STREAM:+-o "$OUTPUT_STREAM"} --config ${FHiCL_FILE} ${INPUT_FILE}
                expcode_exitcode=$?
                if [[ ${expcode_exitcode} -ne 0 ]]; then
                    let $((counter++))
                    echo -e "\n\n*** ${EXECUTABLE_NAME} can not access the input file, wait 30 s, then retry #${counter}\n\n"
                    sleep 30
                fi
            done
 	    exit $expcode_exitcode
        )

    else
        echo -e "\nCI MSG BEGIN\n Stage: ${STAGE_NAME}\n Task: ${TASKSTRING}\n skipped\nCI MSG END\n"
    fi
    exitstatus $?

    # This is uBooNE specific.
    # The mergeana stage doesn't produce an artRoot file,
    # but the CI expect to have it.
    if [[ "${STAGE_NAME}" == "mergeana" && "$(basename ${0})" == *"uboonecode"* ]]; then
        for CUR_OUT in ${OUTPUT_STREAM//-o/}; do
            touch ${CUR_OUT}
        done
    fi

    # This is LArIAT specific.
    # The slicer stage has output filename in the form file_%#.root
    if [[ "${OUTPUT_STREAM}" == *"_%#"* &&  "${STAGE_NAME}" == "slicer" && "$(basename ${0})" == *"lariatsoft"* ]]; then

        for CUR_OUT in ${OUTPUT_STREAM//-o/} # this is LArIAT specific
        do
            CUR_OUT=${CUR_OUT//*:/}
            CUR_OUT2=$(echo $CUR_OUT | sed -e "s/_%#// ; s/_1.root/.root/")
            ln -fv ${CUR_OUT//_%#.root/_1.root} ${CUR_OUT2} && rm ${CUR_OUT//_%#.root/_1.root}
        done
    fi

}

function generate_data_dump
{
    TASKSTRING="generate_data_dump for ${file_stream} output stream"
    ERRORSTRING="W~Error during dump Generation~Check the log"

    trap 'LASTERR=$?; FUNCTION_NAME=${FUNCNAME[0]:-main};  exitstatus ${LASTERR} trap ${LINENO}; exit ${LASTERR}' ERR

#     local NEVENTS=1

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PRINT THE COMMAND TO LOG AND THEN GENERATE THE DUMP FOR THE REFERENCE FILE ~~~~~~~~~~~~~~~~~~~
    echo -e "\nGenerating Dump for ${reference_file}"
    REF_DUMP_FILE=$(basename ${reference_file} | sed -e 's/.root/.dump/')
    echo "${EXECUTABLE_NAME} --rethrow-all -n ${NEVENTS} --config eventdump.fcl ${reference_file} 2>&1 | tee ${REF_DUMP_FILE}"

    (
        set -o pipefail

        local counter=0
        local expcode_exitcode=20
        until [[ ${expcode_exitcode} -ne 20 || ${counter} -gt 5 ]]; do
            ${EXECUTABLE_NAME} --rethrow-all -n ${NEVENTS} --config eventdump.fcl "${reference_file}" 2>&1 | tee ${REF_DUMP_FILE}
            expcode_exitcode=$?
            if [[ ${expcode_exitcode} -ne 0 ]]; then
                let $((counter++))
                echo -e "\n\n*** ${EXECUTABLE_NAME} can not access the input file, wait 30 s, then retry #${counter}\n\n"
                sleep 30
            fi
        done
        exit ${expcode_exitcode}
    )

    #~~~~~~~~~~~~~~~~~~~~~~~~~SAVE IN A VARIABLE THE PARSED REFERENCE DUMP FILE ~~~~~~~~~~~~~~~~~~~~~~
    OUTPUT_REFERENCE=$(cat "${REF_DUMP_FILE}" | sed -e  '/PRINCIPAL TYPE:/,/^\s*$/!d ; s/PRINCIPAL TYPE:.*$// ; /^\s*$/d' )

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PRINT THE COMMAND TO LOG AND THEN GENERATE THE DUMP FOR THE CURRENT FILE ~~~~~~~~~~~~~~~~~~~
    echo -e "\nGenerating Dump for ${current_file}"
    echo "${EXECUTABLE_NAME} --rethrow-all -n ${NEVENTS} --config eventdump.fcl ${current_file} 2>&1 | tee ${current_file//.root}.dump"

    (
        set -o pipefail

        local counter=0
        local expcode_exitcode=20
        until [[ ${expcode_exitcode} -ne 20 || ${counter} -gt 5 ]]; do
            ${EXECUTABLE_NAME} --rethrow-all -n ${NEVENTS} --config eventdump.fcl "${current_file}" 2>&1 | tee "${current_file//.root}".dump
            expcode_exitcode=$?
            if [[ ${expcode_exitcode} -ne 0 ]]; then
                let $((counter++))
                echo -e "\n\n*** ${EXECUTABLE_NAME} can not access the input file, wait 30 s, then retry #${counter}\n\n"
                sleep 30
            fi
        done
        exit ${expcode_exitcode}
    )
    #~~~~~~~~~~~~~~~~~~~~~~~~~SAVE IN A VARIABLE THE PARSED CURRENT DUMP FILE ~~~~~~~~~~~~~~~~~~~~~~
    OUTPUT_CURRENT=$(cat "${current_file//.root}".dump | sed -e  '/PRINCIPAL TYPE:/,/^\s*$/!d ; s/PRINCIPAL TYPE:.*$// ; /^\s*$/d' )

    echo -e "\nReference files for ${file_stream} output stream:"
    echo -e "\n${REF_DUMP_FILE}\n"
    echo "$OUTPUT_REFERENCE"
    echo -e "\nCurrent files for ${file_stream} output stream:"
    echo -e "\n${current_file//.root}.dump\n"
    echo "$OUTPUT_CURRENT"

    exitstatus $?
}

function compare_products_names
{
    TASKSTRING="compare_products_names for ${file_stream} output stream"
    ERRORSTRING="W~Error comparing products names~check the log"

    if [[ "$1" -eq 1 ]]
    then
        REF_DUMP_FILE=$(basename ${reference_file} | sed -e 's/.root/.dump/')
        echo -e "\nCompare products names for ${file_stream} output stream."
        #~~~~~~~~~~~~~~~~CHECK IF THERE'S A DIFFERENCE BEETWEEN THE TWO DUMP FILES IN THE FIRST FOUR COLUMNS~~~~~~~~~~~~~~
        DIFF=$(diff  <(sed 's/\.//g ; /PRINCIPAL TYPE:/,/^\s*$/!d ; s/PRINCIPAL TYPE:.*$// ; /^\s*$/d' ${REF_DUMP_FILE} | cut -d "|" -f -4 ) <(sed 's/\.//g ; /PRINCIPAL TYPE:/,/^\s*$/!d ; s/PRINCIPAL TYPE:.*$// ; /^\s*$/d' ${current_file//.root/.dump} | cut -d "|" -f -4 ) )
        STATUS=$?

        echo -e "\nCheck for added/removed data products"
        echo -e "difference(s)\n"
        #~~~~~~~~~~~~~~~IF THERE'S A DIFFERENCE EXIT WITH ERROR CODE 201~~~~~~~~~~~~~~~
        if [[ "${STATUS}" -ne 0  ]]; then
            echo "${DIFF}"
            ERRORSTRING="W~Differences in products names~Request new reference files"
            exitstatus 201
        else
            echo -e "none\n\n"
        fi
    else
        echo -e "\nCI MSG BEGIN\n Stage: ${STAGE_NAME}\n Task: ${TASKSTRING}\n skipped\nCI MSG END\n"
        exitstatus $?
    fi
}

function compare_products_sizes
{
    TASKSTRING="compare_products_sizes for ${file_stream} output stream"
    ERRORSTRING="W~Error comparing product sizes~Check the log"


    if [[ "${1}" -eq 1 ]]
    then

        REF_DUMP_FILE=$(basename ${reference_file} | sed -e 's/.root/.dump/')
        echo -e "\nCompare products sizes for ${file_stream} output stream.\n"
        #~~~~~~~~~~~~~~~~CHECK IF THERE'S A DIFFERENCE BEETWEEN THE TWO DUMP FILES,IN ALL THE COLUMNS~~~~~~~~~~~~~~
        DIFF=$(diff  <(sed 's/\.//g ; /PRINCIPAL TYPE:/,/^\s*$/!d ; s/PRINCIPAL TYPE:.*$// ; /^\s*$/d' ${REF_DUMP_FILE}) <(sed 's/\.//g ; /PRINCIPAL TYPE:/,/^\s*$/!d ; s/PRINCIPAL TYPE:.*$// ; /^\s*$/d' ${current_file//.root/.dump}) )
        STATUS=$?
        echo -e "\nCheck for differences in the size of data products"
        echo -e "difference(s)\n"

        #~~~~~~~~~~~~~~~IF THERE'S A DIFFERENCE EXIT WITH ERROR CODE 202 ~~~~~~~~~~~~~~~~~~~~~~~
        if [[ "${STATUS}" -ne 0 ]]; then
            echo "${DIFF}"
            ERRORSTRING="W~Differences in products sizes~Request new reference files"
            exitstatus 202
        else
            echo -e "none\n\n"
        fi
    else
        echo -e "\nCI MSG BEGIN\n Stage: ${STAGE_NAME}\n Task: ${TASKSTRING}\n skipped\nCI MSG END\n"
        exitstatus $?
    fi
}


#~~~~~~~~~~~~~~~~~~~~~~~PRINT AN ERROR MESSAGE IN THE PROGRAM EXIT WITH AN ERROR CODE~~~~~~~~~~~~~~~~
function exitstatus
{
    EXITSTATUS="$1"
    if [ "$2" == "trap" ];then
        echo -e "\nCI MSG BEGIN\n Script: `basename $0`\n Function: ${FUNCTION_NAME} - error at line ${3}\n Stage: ${STAGE_NAME}\n Task: ${TASKSTRING}\n exit status: $EXITSTATUS\nCI MSG END\n"
    else
        echo -e "\nCI MSG BEGIN\n Stage: ${STAGE_NAME}\n Task: ${TASKSTRING}\n exit status: ${EXITSTATUS}\nCI MSG END\n"
    fi
    #don't exit if the fetch of the reference failed,because we need to produce one and then upload it
    if [[ "${EXITSTATUS}" -ne 0 ]]; then
        if [[ -n "$ERRORSTRING" ]];then
            echo "`basename $PWD`~${EXITSTATUS}~$ERRORSTRING" >> $WORKSPACE/data_production_stats${ci_cur_exp_name}.log
        fi
        exit "${EXITSTATUS}"
    fi
}


#~~~~~~~~~~~~~~~~~~~~~~~EXTRA FUNCTIONS TO RUN SPECIAL TESTS~~~~~~~~~~~~~~~~
function compare_anatree
{

    source ${GENERIC_CI_DIR}/bin/reporter_functions.sh

    THISCIDIR=$(eval "echo \${${PROJ_PREFIX}_CI_DIR}")

    root -l -b -q ${THISCIDIR}/test/compare_anatree.C\(\"${1}\",\"${2}\"\)



    if [ -n "${BUILD_ID}" ] # not in a jenkins build, don't send plots to CI web app
    then
        BASE_JOB=`echo $JOB_NAME | sed -e 's;/.*;;' -e 's/_jenkins.*//'`
        export report_fullname="${BASE_JOB}/${BUILD_NUMBER}"
        export report_serverurl="${HUDSON_URL}"

        for f in *.gif
        do
            bf=`basename $f`
            hist_desc="${bf//.gif/} plot"
            hist_name="${bf//.gif/}"
            report_img "ci_tests" "${ci_cur_exp_name}" "" "$(basename $PWD)" "$hist_name" "$f" "$hist_desc"
            # report_img "ci_tests" "" "end" "$hist_name" "$f" "$hist_desc"
        done

    fi

}


#~~~~~~~~~~~~~~~~~~~~~~~~MAIN OF THE SCRIPT~~~~~~~~~~~~~~~~~~
trap 'LASTERR=$?; FUNCTION_NAME=${FUNCNAME[0]:-main};  exitstatus ${LASTERR} trap ${LINENO}; exit ${LASTERR}' ERR

initialize $@

data_production "${check_data_production}"

#~~~~~~~~~~~~~~~~PROCESS ALL THE FILES DECLARED INTO THE OUTPUT LIST~~~~~~~~~~~~~~~~~
for filename in ${OUTPUT_LIST//,/ }
do
    file_stream=$(echo "${filename}" | cut -d ':' -f 1)
    current_file=$(echo "${filename}" | cut -d ':' -f 2)

    echo "filename: ${filename}"
    echo "file_stream: ${file_stream}"
    echo "current_file: ${current_file}"

    reference_file=""
    if [[ "${REFERENCE_FILES}" = *"xroot"* ]]; then
        for temp_file in ${REFERENCE_FILES//,/ }
        do
            if [[ "${temp_file}" = *"${current_file//Current/Reference}"* ]]; then
                reference_file=${temp_file}
            fi
        done

        ### XROOTD_PATH=$(echo "${reference_file}" | sed -e 's#\(^.*\)/'${STAGE_NAME}'/.*#\1/'${STAGE_NAME}/'#')

    else
        ### if [ -n "${build_platform}" ]; then
        ###     reference_file=$(echo "${current_file%`echo ${build_platform}`*}${build_platform}${current_file#*`echo ${build_platform}`}")
        ### else
        ###     reference_file=$(echo "${current_file}")
        ### fi
        ### reference_file="${reference_file//Current/Reference}"
        reference_file=$(echo ${current_file//Current/Reference} | sed -e 's#'${build_identifier}'##' )
    fi

    if [[ "${check_compare_names}" -eq 1  || "${check_compare_size}" -eq 1 ]]; then
        generate_data_dump
    else
        break
    fi

    compare_products_names "${check_compare_names}"

    compare_products_sizes "${check_compare_size}"

done

#~~~~~~~~~~~~~~~~RUN EXTRA FUNCTION~~~~~~~~~~~~~~~~~
if [ -n "${EXTRA_FUNCTION}" ]; then
    ${EXTRA_FUNCTION//,/ }
fi

#!/usr/bin/env bash
#
# Runs a fast simulation + reconstruction chain for test purposes
#

#############################################################################
###  Tests to be executed in chain:
###
declare -ar TestNames=(
    'prodsingle_mu_bnblike_newflux'
    'g4_sce'
    'detsim_sce_lite'
    'standard_reco1_sbnd'
    'reco2_sce'
)

#############################################################################


# Set WireCell env vars
export WIRECELL_PATH=${SBNDCODE_DIR}/sbndcode/WireCell/cfg/:${SBND_DATA_DIR}/WireCell


function ExecBase() {
	local LogFile="$1"
	shift
	local -a Cmd=( "$@" )
	
	local -i res=0
	if [[ -z "${FAKE//0}" ]]; then
		if [[ -n "$LogFile" ]]; then
			"${Cmd[@]}" >& "$LogFile"
			res=$?
			if [[ $res != 0 ]]; then
				cat "$LogFile"
			fi
		else
			"${Cmd[@]}"
			res=$?
		fi
	else
		echo "  (skipped in dry run mode)"
	fi
	return $res
} # ExecBase()


function Exec() {
	local LogFile="$1"
	shift
	local -a Cmd=( "$@" )
	
	echo "\$ ${Cmd[@]}${LogFile:+ >& ${LogFile}}"
	ExecBase "$LogFile" "${Cmd[@]}"
} # Exec()


function ExecCommand() {
	local LogFile="$1"
	local CommandFile="$2"
	
	cat <<-EOH
	\$ $(< "$CommandFile")${LogFile:+" >& ${LogFile}"}
	---------------------------------------------------------------------
	EOH
	ExecBase "$LogFile" 'sh' "$CommandFile"
} # ExecCommand()


function GetWorkDir() {
	local TestName="$1"
	echo "$TestName"
} # GetWorkDir()


function chdir() {
	local DirName="$1"
	mkdir "$DirName" && cd "$DirName"
} # chdir()

function Link() {
	local Source="$1"
	local Target="$2"
	
	if [[ -d "$Target" ]]; then
		Target="${Target%/}/$(basename "$Source")"
	fi
	
	[[ -h "$Target" ]] && rm "$Target"
	ln -s "$Source" "$Target"
} # Link()


function MakeCommandFile() {
	local CommandFile="$1"
	shift
	rm -f "$CommandFile"
	local Command
	local Token
	for Token in "$@" ; do
		[[ -n "$Command" ]] && Command+=' '
		Command+="'${Token}'"
	done
	echo "$Command" > "$CommandFile"
	chmod a+x "$CommandFile"
} # MakeCommandFile()


# number of events to be generated
: ${NEvents:=1}

#
# clean up from previous tests
#
for TestName in "${TestNames[@]}" ; do
	WorkDir="$(GetWorkDir "$TestName")"
	[[ -d "$WorkDir" ]] || continue
	Exec '' rm -r "${WorkDir%/}/"
done

InputFile=''
OutputTreeFile=''
BaseDir="$(pwd)"
declare -i iTest=0
declare -ir NTests="${#TestNames[@]}"
for TestName in "${TestNames[@]}" ; do
	#
	# For each test:
	# - test name:         <TestName>
	# - working directory: <TestName>
	# - configuration:     <TestName>.fcl
	# - input file:        from the previous job
	# - output files:      <TestName>-art.root  (art events)
	#                      <TestName>-hist.root (TFileService)
	# - log file:          <TestName>.out
	# - command run:       run_<TestName>.cmd
	#
	# A link to the output file is placed into the base directory
	#

	let ++iTest
	
	#
	# move into a directory specific to the job
	#
	WorkDir="$(GetWorkDir "$TestName")"
	Exec '' chdir "$WorkDir" || exit $?
	
	cat <<-EOH
	=====================================================================
	   [${iTest}/${NTests}]  Starting test:  ${TestName}
	=====================================================================
	EOH
	
	###
	### prepare the command
	###
	declare -a Cmd=( 'lar' '--rethrow-all' )
	
	ConfigFile="${TestName}.fcl"
	ConfigDumpFile="${ConfigFile%.fcl}.cfg"
	Cmd=( "${Cmd[@]}" '--config' "$ConfigFile" '--config-out' "$ConfigDumpFile" )

	# prepare the input	
	if [[ -n "$OutputTreeFile" ]]; then
		InputFile="$OutputTreeFile"
		Cmd=( "${Cmd[@]}" '--source' "../${InputFile}" )
	else
		Cmd=( "${Cmd[@]}" '--nevts' "$NEvents" )
	fi
	
	OutputTreeFile="${TestName}-art.root"
	OutputHistFile="${TestName}-hist.root"
	Cmd=( "${Cmd[@]}" '--output' "$OutputTreeFile" '--TFileName' "$OutputHistFile" )
	
	
	###
	### run the test
	###
	LogFile="${TestName}.out"
	CommandFile="run_${TestName}.cmd"
	MakeCommandFile "$CommandFile" "${Cmd[@]}"
	ExecCommand "$LogFile" "$CommandFile"
	res=$?

	if [[ $res != 0 ]]; then
		cat <<-EOM >&2
		*********************************************************************
		***  [${iTest}/${NTests}]  Test '${TestName}' FAILED!!! (exit code: ${res})
		*********************************************************************
		EOM
		exit $res
	fi	
	
	#
	# place a link to the output in the base directory
	#
	Exec '' Link "${WorkDir}/${OutputTreeFile}" "${BaseDir%/}/"
	
	# go back to base dir
	cd "$BaseDir"
	
done
cat <<-EOM
=====================================================================
***  All ${NTests} tests succeeded.
=====================================================================
EOM
exit 0


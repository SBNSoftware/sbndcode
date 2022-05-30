#!/usr/bin/env bash
#
# Checks all FHiCL files under the specified directory.
#
# Usage:
#     
#     testFHiCLfiles.sh  [file or directory ...]
#     
# By default all FHiCL files under the current directory are tested.
#
#

SCRIPTNAME="$(basename "$0")"
SCRIPTDIR="$(dirname "$0")"

: ${DefaultTest:='dump'}
: ${DEBUG:=0}

################################################################################
declare -A ColorCodes
function SetupColors() {
  local UseColors="${1:-1}"
  local -r Escape=$'\e['
  if [[ $UseColors != 0 ]]; then
    ColorCodes=(
      ['reset']="${Escape}0m"
      ['red']="${Escape}31m"
      ['yellow']="${Escape}1;33m"
      ['bright red']="${Escape}1;31m"
      ['green']="${Escape}32m"
      ['bright green']="${Escape}1;32m"
      ['cyan']="${Escape}36m"
      )
  fi
} # SetupColors()

function GetColor() {
  local Key="$1"
  echo -e "${ColorCodes[$Key]}"
}

function ApplyColor() {
  local Color="$1"
  shift
  echo -e "$(GetColor "$Color")$*$(GetColor 'reset')"
} # ApplyColor()

function isDebugging() {
  local -i Level="${1:-1}"
  [[ "$DEBUG" -ge "$Level" ]]
}
function DBGN() {
  local -i Level="$1"
  isDebugging "$Level" || return
  shift
  STDERR 'green' "DBG[${Level}]| $*" ;
}
function DBG() { DBGN 1 "$@" ; }
function INFO() { [[ "$Quiet" == 0 ]] && ApplyColor 'cyan' "$*" ; }
function SUCCESS() { ApplyColor 'bright green' "$*" ; }
function STDERR() {
  local Color="$1"
  shift
  ApplyColor "$Color" "$*" >&2
} # STDERR()
function ERROR() { STDERR 'red' "ERROR: $*" ; }
function FATAL() {
  local -i Code="$1"
  shift
  STDERR 'bright red' "FATAL (${Code}): $*"
  exit $Code
} # FATAL()
function LASTFATAL() {
  local -i res="$?"
  [[ $res == 0 ]] || FATAL "$res" "$@"
} # LASTFATAL()


function DebugTest() {
  local -ir MaxLevel="${1:-10}"
  local -i Level
  for (( Level = 1 ; Level < 10 ; ++Level )); do
    DBGN "$Level" "Debug level ${Level} shown"
  done
} # DebugTest()


################################################################################
function printHelp() {
  cat <<EOH
Performs checks on the specified FHiCL files.

Usage:  ${SCRIPTNAME}  [options] [-|--] [Input ...]

A test is performed on every file specified as "Input" on the command line.
If an input specification is a file, that file is tested; if the input is a
directory, all files with suffix '.fcl' under that directory and its
subdirectories are tested.
If no input is specified, the current directory is used as starting point.

Options:
--listtests , -L
    lists all the supported tests
--test=TEST ["${DefaultTest}"]
    choose which test to perform; use \`--listtests\` to see the supported tests
--dump
    shortcut for \`--test=dump\`
--validate
    shortcut for \`--test=validate\`
--exclude=PATTERN
    exclude all files and directories that match this pattern (Bash regex)
--exclude-from=FILE
    equivalent to \`--exclude=PATTERN\` for all patterns in FILE; each line
    in FILE is interpreted as a pattern, except for the lines starting with
    a '#' (no spaces allowed before it)
--quiet
    will print only errors or a single success message
--debug[=LEVEL]
    sets the verbosity to the specified LEVEL for debugging purposes
--color | --no-color  [default: enabled]
    enables or disables ANSI color codes in the output
--help , -h , -?
    shows this usage information

EOH
  
} # printHelp()


################################################################################
declare -Ar TestDescriptions=(
  ['dump']='FHiCL syntax check (`fhicl-dump`)'
  ['validate']='ask `art` to perform validation (`lar --validate-config`)'
)

function printTests() {
  local TestKey
  for TestKey in "${!TestDescriptions[@]}" ; do
    local TestDescription="${TestDescriptions[$TestKey]}"
    cat <<EOM
${TestKey}
    ${TestDescription}
EOM
  done
} # printTests()


function isTestType() {
  local Key="$1"
  for TestKey in "${!TestDescriptions[@]}" ; do
    [[ "$TestKey" == "$Key" ]] && return 0
  done
  return 1
} # isTestType()


function SelectTest() {
  local Type="${1:-${TestType}}"
  
  case "$Type" in
    # add non-standard types here
    ( * )
      if isTestType "$Type" ; then
        TestProc="Test_${Type}"
        TestPrep="${TestProc}_prep"
        DBG "Selected test '${Type}'."
        return 0
      fi
      return 1
      ;;
  esac
  
} # SelectTest()


function isExcluded() {
  
  local Path="$1"
  
  local ExcludePattern
  for ExcludePattern in "${ExcludedPatterns[@]}" ; do
    [[ "$Path" =~ $ExcludePattern ]] && return 0
  done
  return 1

} # isExcluded()


function ExcludeFromFile() {
  
  local File="$1"
  
  [[ -r "$File" ]] || FATAL 2 "Can't read exclusion file '${File}'!"
  
  local -i nExcluded=0
  local Pattern
  while read Pattern ; do
    [[ -z "$Pattern" ]] && continue
    [[ "${Pattern:0:1}" == '#' ]] && continue
    DBGN 2 "Exclude pattern: '${Pattern}'"
    ExcludedPatterns+=( "$Pattern" )
    let ++nExcluded
  done < "$File"
  
  DBGN 1 "${nExcluded} exclusion patterns loaded from '${File}'."
  
} # ExcludeFromFile()


################################################################################
function Test_dump_prep() {
  
  if [[ -z "$fhicldump" ]]; then
    fhicldump="$(which 'fhicl-dump' 2> /dev/null)"
    [[ -x "$fhicldump" ]] || FATAL 1 "Couldn't find 'fhicl-dump' executable!"
  fi
  
} # Test_dump_prep()


function Test_dump() {
  
  local FHiCLpath="$1"
  
  INFO "Testing: '${FHiCLpath}'"
  
  local WithPath=0
  [[ "$FHiCLpath" =~ / ]] && WithPath=1
  
  local -a Options
  [[ "$WithPath" == 1 ]] && Options+=( '-l' 'after1' )
  
  $fhicldump --quiet "${Options[@]}" --config "$FHiCLpath"
  local -i res=$?
  
  [[ $res == 0 ]] || ERROR "File '${FHiCLpath}' failed verification (code: ${res})."
  return $res
} # Test_dump()


################################################################################
function Test_validate_prep() {
  
  if [[ -z "$lar" ]]; then
    lar="$(which 'lar' 2> /dev/null)"
    [[ -x "$lar" ]] || FATAL 1 "Couldn't find 'lar' executable!"
  fi
  
} # Test_validate_prep()


function Test_validate() {
  
  local FHiCLpath="$1"
  
  INFO "Validating: '${FHiCLpath}'"
  
  local WithPath=0
  [[ "$FHiCLpath" =~ / ]] && WithPath=1
  
  local -a Options
  
  local Cwd="$(pwd)"
  local -a Cmd=( $lar --validate-config '/dev/null' "${Options[@]}" --config "$FHiCLpath" )
  "${Cmd[@]}" > /dev/null
  local -i res=$?
  
  if [[ $res != 0 ]]; then
    ERROR "File '${FHiCLpath}' failed validation (code: ${res})."
    INFO "Command was:\n${Cwd}\$ ${Cmd[@]}"
  fi
  return $res
} # Test_validate()


################################################################################
function TestFHiCLfile() {
  
  local FilePath="$1"
  
  if isExcluded "$FilePath" ; then
    echo "Skipping '${FilePath}' (excluded)"
    return
  fi
  
  "$TestProc" "$FilePath"
  local -i res=$?
  if [[ $res != 0 ]]; then
    Errors["$FilePath"]="$res"
    return 1
  fi
  return 0
} # TestFHiCLfile()


function TestFHiCLdirectory() {
  
  local DirPath="$1"
  
  if [[ ! -d "$DirPath" ]]; then
    ERROR "Path '${DirPath}' is not a directory."
    return 1
  fi
  
  if isExcluded "$DirPath" ; then
    echo "Skipping '${DirPath}' and subdirectories (excluded)"
    return
  fi
  
  local -i nErrors=0
  local -i res
  local FilePath
  echo "Testing FHiCL files in '${DirPath}' and subdirectories"
  while read FilePath ; do
    TestFHiCLfile "$FilePath" || let ++nErrors
  done < <( find "$DirPath" -name "*.fcl" )
  
  [[ $nErrors -gt 0 ]] # return value
} # TestFHiCLfile()


function TestFHiCLpath() {
  local Path="$1"
  
  if [[ -d "$Path" ]]; then
    TestFHiCLdirectory "$Path"
  else
    TestFHiCLfile "$Path"
  fi
} # TestFHiCLpath()

################################################################################
###
### argument parsing
###
declare -a InputPaths
declare -a ExcludedPatterns
declare TestType="$DefaultTest"
declare -i DoHelp=0 DoListTests=0 UseColors=1 Quiet=0
declare -i NoMoreOptions=0
for (( iParam = 1 ; iParam <= $# ; ++iParam )); do
  Param="${!iParam}"
  if [[ $NoMoreOptions == 0 ]] && [[ "${Param:0:1}" == '-' ]]; then
    case "$Param" in
      ( '--test='* )             TestType="${Param#--*=}" ;;
      ( '--exclude='* )          ExcludedPatterns+=( "${Param#--*=}" ) ;;
      ( '--exclude-from='* )     ExcludedPatternFiles+=( "${Param#--*=}" ) ;;
      ( '--listtests' | '-L' )   DoListTests=1 ;;
      ( '--color' )              UseColors=1 ;;
      ( '--no-color' )           UseColors=0 ;;
      ( '--quiet' | '-q' )       Quiet=1 ;;
      ( '--debug' | '-d' )       DEBUG=1 ;;
      ( '--debug='* )            DEBUG="${Param#--*=}" ;;
      ( '--help' | '-h' | '-?' ) DoHelp=1 ;;
      ( '--' | '-' )             NoMoreOptions=1 ;;
      ( '--'* )
        if isTestType "${Param#--}" ; then
          TestType="${Param#--}"
        else
          FATAL 1 "Unsupported option: '${Param}'."
        fi
    esac
    [[ $DoHelp == 1 ]] && break
  else
    InputPaths+=( "$Param" )
  fi
done

for File in "${ExcludedPatternFiles[@]}" ; do
  ExcludeFromFile "$File"
done

################################################################################

SetupColors "$UseColors"
DebugTest

if [[ $DoHelp == 1 ]]; then
  printHelp
  exit
fi

if [[ $DoListTests == 1 ]]; then
  echo "Supported tests are:"
  printTests
  exit
fi

#
# test selection and preparation
#
if ! SelectTest "$TestType" ; then
  ERROR "Unrecognised test: '${TestType}'. Supported tests are:"
  printTests
  FATAL 1 "Unrecognised test: '${TestType}'."
fi


[[ "${#InputPaths[@]}" == 0 ]] && InputPaths=( '.' )

if [[ -n "$TestPrep" ]]; then
  "$TestPrep"
  LASTFATAL 1 "Preparation to test '${TestType}' failed."
fi

#
# the test
#
declare -A Errors
for InputPath in "${InputPaths[@]}" ; do
  
  TestFHiCLpath "$InputPath"
  
done

if [[ "${#Errors[@]}" -gt 0 ]]; then
  STDERR 'bright red' "${#Errors[@]} errors found in test '${TestType}':"
  for Failed in "${!Errors[@]}" ; do
    STDERR 'bright red' "- '${Failed}' (code: ${Errors[$Failed]})"
  done
  exit 1
else
  SUCCESS "All tested files passed the verification."
  exit 0
fi

################################################################################

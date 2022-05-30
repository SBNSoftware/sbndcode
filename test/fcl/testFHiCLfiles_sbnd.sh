#!/usr/bin/env bash

SCRIPTDIR="$(dirname "$0")"

#
# platform compatibility test
#
if [[ "${BASH_VERSINFO[0]}" -lt 4 ]]; then
  echo "This test script requires BASH 4.0 or newer." >&2
  exit # the test is declared a success
fi

#
# find the actual test script
# currently located alongside this file
#
declare -r TestScriptName='testFHiCLfiles.sh'
declare TestScript
TestScript="$(which "$TestScriptName" 2> /dev/null)"
[[ $? != 0 ]] && TestScript="${SCRIPTDIR}/${TestScriptName}"
[[ -r "$TestScript" ]] || TestScript="./${TestScriptName}"

if [[ ! -r "$TestScript" ]]; then
  echo "FATAL: test script '${TestScriptName}' not found!" >&2
  exit 2
fi

#
# decide what to test
#
declare -a TestDirs

if [[ -n "$SBNDCODE_DIR" ]]; then
  DistributedFHiCLdir="${SBNDCODE_DIR}/fcl"
  if [[ -d "$DistributedFHiCLdir" ]]; then
    echo "Will test the distributed job configuration directory '${DistributedFHiCLdir}'"
    TestDirs+=( "$DistributedFHiCLdir" )
  fi
fi

if [[ -n "$MRB_BUILDDIR" ]]; then
  InstalledFHiCLdir="${MRB_BUILDDIR}/sbndcode/fcl"
  if [[ -d "$InstalledFHiCLdir" ]]; then
    echo "Will test the installed job configuration directory in MRB build area ('${InstalledFHiCLdir}')"
    TestDirs+=( "$InstalledFHiCLdir" )
  fi
fi

if [[ -n "$MRB_SOURCE" ]]; then
  SourceFHiCLdir="${MRB_SOURCE}/sbndcode/sbndcode/JobConfigurations"
  if [[ -d "$SourceFHiCLdir" ]]; then
    echo "Will test the job configuration directory in MRB source area ('${SourceFHiCLdir}')"
    TestDirs+=( "$SourceFHiCLdir" )
  fi
fi

if [[ "${#TestDirs[@]}]" == 0 ]]; then
  echo "FATAL: no suitable FHiCL directory found to be tested!" >&2
  exit 1
fi

#
# run the test
#
"$SHELL" "$TestScript" "$@" "${TestDirs[@]}"


#!/usr/bin/env bash
#
# Purpose: exercise compilation of gallery in a compatible environment.
# This script will return a non-zero exit code on compilation error.
#

# the first simple command exiting with an error will cause the shell to exit;
# sometimes this interferes with the commands (e.g. `setup`); in that case,
# some trick is needed.
set -e

declare -r CWD="$(pwd)"
declare -r WorkDir="${CWD}/galleryTest"

#
# check that `setup` function is available
#
if ! declare -f setup >& /dev/null ; then
  if [[ -z "$SETUP_UPS" ]]; then
    echo "UPS not correctly set up!!"
  else
    # provide setup and unsetup functions; compared to the ones from UPS,
    # the `ups` executable is searched in PATH rather than hard-coded
    function setup() { source "$(ups setup "$@")" ; }
    function unsetup() { source "$(ups unsetup "$@")" ; }
  fi
fi

#
# create a working area
#
declare -r SourceDir="${WorkDir}/source"
declare -r BuildDir="${WorkDir}/build"

declare -r SBNDgalleryBaseDir="${SBNDCODE_DIR}/sbndcode/gallery"

rm -Rf "$WorkDir"
mkdir "$WorkDir"
cd "$WorkDir"

cp -a "${SBNDgalleryBaseDir}/galleryAnalysis" "$SourceDir"

rm -Rf "$BuildDir"
mkdir "$BuildDir"

#
# set up gallery
#
ups active
source "${SBNDgalleryBaseDir}/helpers/sbnd_gallery_setup"

#
# proceed with compilation
#
cmake "$SourceDir"
make


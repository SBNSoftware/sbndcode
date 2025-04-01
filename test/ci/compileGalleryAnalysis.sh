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
# create a working area
#
declare -r SourceDir="${WorkDir}/source"
declare -r BuildDir="${WorkDir}/build"

declare -r SBNDgalleryBaseDir="${SBNDCODE_DIR}/examples/gallery"

rm -Rf "$WorkDir"
mkdir "$WorkDir"
cd "$WorkDir"

cp -a "${SBNDgalleryBaseDir}/galleryAnalysis" "$SourceDir"

rm -Rf "$BuildDir"
mkdir "$BuildDir"

#
# set up
#
# set up some UPS version of CMake only if none is set up already:
[[ -n "$CMAKE_DIR" ]] || source "$(ups setup cmake 'v3_27_4')"
[[ -n "$CETMODULES_DIR" ]] || source "$(ups setup cetmodules 'v3_24_00')"

#
# proceed with compilation
#
cd "$BuildDir"
cmake "$SourceDir"
make


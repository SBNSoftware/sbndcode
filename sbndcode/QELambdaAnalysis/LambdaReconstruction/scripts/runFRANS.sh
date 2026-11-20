#!/bin/bash

# Default values
default_file=""
default_debug=0
default_mode=-1
default_n=-1
default_ev=-1
default_sr=-1
default_nskip=-1


# Parse the optional flags using getopts
while getopts ":s:d:m:n:e:r:k:" opt; do
  case $opt in
    s) argFile="$OPTARG";;
    d) argDebug="$OPTARG";;
    m) argMode="$OPTARG";;
    n) argN="$OPTARG";;
    e) argE="$OPTARG";;
    r) argSr="$OPTARG";;
    k) argNskip="$OPTARG";;
    \?) echo "Invalid option: -$OPTARG" >&2
        exit 1;;
  esac
done
shift $((OPTIND-1))

# Set the required arguments
argFile="${argFile:-$default_file}"
argDebug="${argDebug:-$default_debug}"
argMode="${argMode:-$default_mode}"
argN="${argN:-$default_n}"
argE="${argE:-$default_ev}"
argSr="${argSr:-$default_sr}"
argNskip="${argNskip:-$default_nskip}"


# Launch the ROOT interactive session
root -l <<EOF
.L $TPCLINES_DIR/RunAlgoFRANS.cpp


RunAlgoFRANS($argDebug, $argMode, $argN, $argNskip, $argE, $argSr, "$argFile", "", "")

.q
EOF

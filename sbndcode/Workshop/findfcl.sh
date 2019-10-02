#!/bin/bash
 
if [ $# -ne 1 ]; then
    echo "Error: please pass a fcl file name (can include wildcards)"
  exit 1
fi

if [ -z ${FHICL_FILE_PATH+x} ]; then
  echo "Error: FHICL_FILE_PATH has not been set!"
  exit 2
fi

SEARCH_PATHS=`echo $FHICL_FILE_PATH | sed 's/:/\n/g'`
for THIS_PATH in $SEARCH_PATHS; do
  if [ -d $THIS_PATH ]; then
    find $THIS_PATH -name $1
  fi
done

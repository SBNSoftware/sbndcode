#!/bin/bash

FileNameList=$1
FindFilePath=$2


if [ $# -ne 2 ]; then
  echo "Usage: $0 <FileNameList> <FindFilePath>"
  exit 1
fi
# Check if both files exist
if [ ! -f "$FileNameList" ]; then
  echo "File $FileNameList not found."
  exit 1
fi

cp $FileNameList $FileNameList.copy

# Use a while loop to process both files simultaneously
while IFS= read -r line1 <&3; do
    # Process lines from both files here
    find $FindFilePath*/$line1 | pnfsToXRootD >> $FileNameList.output | grep -v "Reading input URLs from standard input."
done 3< "$FileNameList"

mv $FileNameList $FileNameList.old
mv $FileNameList.output $FileNameList
rm $FileNameList.old


#find /pnfs/sbnd/persistent/users/jiaoyang/dirt_neutrino/v09_59_00/dirt_neutrino_crt_beam_tel_1004/00/reco/outdir/*/$line1 | pnfsToXRootD >> $FileNameList.output | grep -v "Reading input URLs from standard input."

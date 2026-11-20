#!/bin/bash

# Check if three filenames are provided as command-line arguments
if [ $# -ne 3 ]; then
  echo "Usage: $0 <EventFileList> <FindFilePath> <OutputDirectory>"
  return 1
fi

EventFileList=$1
FindFilePath=$2
OutDir=$3

# Check if the input file exists
if [ ! -f "$EventFileList" ]; then
  echo "File $EventFileList not found."
  return 1
fi

# Check if the output directory exists, create it if not
mkdir -p $OutDir
rm -f $OutDir/crtEventDisplay_*.root

echo "Input file $EventFileList"

# Use a while loop to process the file
while IFS=' ' read -r line1 line2; do
    # Process lines from the file
    fullXRootDPath=$(find $FindFilePath*/$line2 | pnfsToXRootD | grep -v "Reading input URLs from standard input.")
    echo "Full xRootD Path: $fullXRootDPath"

    # Process lines from the second file
    echo "lar -c run_tpcana_hypana_MC.fcl $fullXRootDPath -e $line1 -n 1"
    lar -c run_tpcana_hypana_MC.fcl $fullXRootDPath -e $line1 -n 1
    mv analyze*.root $OutDir/analyzeTPC_$line1.root
done < "$EventFileList"

# Calculate and display elapsed time
start_time=$(date +%s)
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
minutes=$((elapsed_time / 60))
seconds=$((elapsed_time % 60))
echo "Elapsed time: $minutes minutes and $seconds seconds"

#!/bin/bash

# Check if three filenames are provided as command-line arguments
if [ $# -ne 3 ]; then
  echo "Usage: $0 <EventFileList> <FindFilePath> <OutputDirectory>"
  return 1
fi

EventFileList=$1
FindFilePath=$2
OutDir=$3
FhiclFile=$4

# Check if the input file exists
if [ ! -f "$EventFileList" ]; then
  echo "File $EventFileList not found."
  return 1
fi

# Initialize counter
counter=1

echo "Input file $EventFileList"

# remove output path
rm -r "$OutDir"
mkdir "$OutDir"

# Use a while loop to process the file
while IFS=' ' read -r line1 line2; do
    # Check if line2 exists in the external file
    if grep -q "$line2" "$FindFilePath"; then
        # Save the full xRootD path in the variable
        fullXRootDPath=$(grep "$line2" "$FindFilePath")
        echo "Found a match: $fullXRootDPath"
        # Add your processing logic here
    else
        echo "No match found for line 2: $line2"
        # Add your logic for lines with no match
    fi

    # Process lines from the second file
    #echo "lar -c run_tpcana_hypana_MC.fcl $fullXRootDPath -e $line1 -n 1"
    #lar -c run_tpcana_hypana_MC.fcl $fullXRootDPath -e $line1 -n 1
    echo "lar -c $FhiclFile $fullXRootDPath -e $line1 -n 1"
    lar -c $FhiclFile -s $fullXRootDPath -e $line1 -n 1 -T debugArtFile.root
    mv debugArtFile.root "$OutDir/debug_${counter}_${line1}.root"

    # Increment counter
    ((counter++))
done < "$EventFileList"

# Calculate and display elapsed time
start_time=$(date +%s)
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
minutes=$((elapsed_time / 60))
seconds=$((elapsed_time % 60))
echo "Elapsed time: $minutes minutes and $seconds seconds"

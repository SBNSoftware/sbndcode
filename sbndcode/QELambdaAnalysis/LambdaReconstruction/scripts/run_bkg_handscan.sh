#!/bin/bash

start_time=$(date +%s)  # Get the start time in seconds since epoch
# Check if two filenames are provided as command-line arguments
if [ $# -ne 3 ]; then
  echo "Usage: $0 <FileNameList> <EventNumberList> <OutDir>"
  exit 1
fi

FileNameList=$1
EventNumberList=$2
OutDir=$3

mkdir -p $OutDir
rm -f $OutPutDir/crtEventDisplay_*.root

# Check if both files exist
if [ ! -f "$FileNameList" ]; then
  echo "File $FileNameList not found."
  exit 1
fi

if [ ! -f "$EventNumberList" ]; then
  echo "File $EventNumberList not found."
  exit 1
fi
# Use a while loop to process both files simultaneously
while IFS= read -r line1 <&3 && IFS= read -r line2 <&4; do
    # Process lines from both files here
    echo "lar -c run_tpcana_hypana_MC.fcl $line1 -e $line2 -n 1"
    lar -c run_tpcana_hypana_MC.fcl $line1 -e $line2 -n 1
    mv analyze*.root $OutDir/analyzeTPC_$line2.root
done 3< "$FileNameList" 4< "$EventNumberList"

end_time=$(date +%s)  # Get the end time in seconds since epoch
elapsed_time=$((end_time - start_time))  # Calculate the elapsed time in seconds

# Calculate minutes and seconds
minutes=$((elapsed_time / 60))
seconds=$((elapsed_time % 60))

echo "Elapsed time: $minutes minutes and $seconds seconds"
#cp eventdisplay_crt_bt_sbnd.fcl.og eventdisplay_crt_bt_sbnd.fcl
#echo "physics.analyzers.crtEventDisplay.Run_SubRun_Event: [ $1, $2, $3 ]"
#echo "physics.analyzers.crtEventDisplay.Run_SubRun_Event: [ $1, $2, $3 ]" >> eventdisplay_crt_bt_sbnd.fcl
#
#lar -c eventdisplay_crt_bt_sbnd.fcl -e $1:$2:$3 -n 1 -T crt_event_display_$1_$2_$3.root

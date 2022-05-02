#!/bin/bash
#
# Author: bjpjones@fnal.gov from echurch@fnal.gov from dbox@fnal.gov
#
# Small subset of a script to run the optical library building job on the grid within larbatch infrastructure, 
#
# To run this job:
#
# 
# You will get outputs in the area specified by the "outstage" variable 
# which is specified below.
#
# The form of the output is one file for each few voxels. These then need 
# stitching together, which is done after all jobs are done, with a
# dedicated stitching script.
#

umask 0002
verbose=T

echo $FHICL_FILE_PATH "fhicl file path"



offset=`echo "($NEVT/($NJOBS))" | bc`   
fOFFSET=`echo "($offset* ($PROCESS))" | bc`
fFIRSTEVENT=`echo "($fOFFSET+1)" | bc`     #add +1 to start on the right event.




# And then tell the user about it:
echo "physics.producers.generator.Offset: $fOFFSET " >> $FCL
echo "source.firstEvent: $fFIRSTEVENT " >> $FCL


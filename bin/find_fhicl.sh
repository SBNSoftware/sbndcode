#!/bin/bash

if [ x$FHICL_FILE_PATH = x ]; then
  echo "***************************************"
  echo "Variable FHICL_FILE_PATH not found."
  echo "You porbably haven't set up larsoft yet,"
  echo "Try 'setup uboonecode vXX_XX_XX -q e10:prof"
  echo "OR 'mrbsetenv'"
  echo "***************************************"
  exit 1

fi

FHICL_SEARCH_FILE=$1

if [ x$FHICL_SEARCH_FILE = x ]; then
  echo "***************************************"
  echo "USAGE: find_fhicl <fhicl file name>."
  echo "Note that \$FHICL_FILE_PATH must be defined."
  echo "Try 'setup uboonecode vXX_XX_XX -q e10:prof"
  echo "OR 'mrbsetenv' if it isn't"
  echo "***************************************"
  exit 1


fi

SEARCH_PATHS=(`awk '{split($0,array,":"); for (a in array)  printf "%s ", array[a]; printf "\n";}' <<< $FHICL_FILE_PATH`)

if [ ! -d "srcs" ]; then
  echo "***************************************"
  echo "NOTE: I do not see a 'srcs' directory in $PWD"
  echo "I will continue to search \$FHICL_FILE_PATH,"
  echo "but you will unlikely have write access to the file."
  echo "Check out a version of uboonecode using the instructions"
  echo "Provided here: " 
  echo "https://cdcvs.fnal.gov/redmine/projects/uboonecode/wiki/Uboone_guide"
  echo "***************************************"



#if srcs directory exists, add it to the search path
else
  SEARCH_PATHS=("${SEARCH_PATHS[@]}" "srcs") 
fi


for elt in ${SEARCH_PATHS[*]};
do
  
  #skip local dirs autmoatically added to the path but do not exist
 #echo $CHECKED_WORKING_DIR
 if [ ! -d "$elt" ]; then
   continue
 fi
 
 # also, skip the current working dir (".")
 if [ "$elt" == "." ]; then
   continue
 fi
 
 
 #echo $elt 
 FOUND_FHICL=`find $elt -name $FHICL_SEARCH_FILE`
 
 if [ -n "$FOUND_FHICL" ]; then
   echo "=========================="
   echo "Found fhicl file(s):"
    
   awk -F:" " '{printf "%s \n", $1}' <<< $FOUND_FHICL
 fi
    
done

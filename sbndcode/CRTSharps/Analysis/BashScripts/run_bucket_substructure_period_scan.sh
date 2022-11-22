#!/bin/bash

periods=(18.50 18.60 18.70 18.80 18.90 18.92 18.93 18.94 18.95 18.96 19.00 19.10 19.20 19.30 19.40 19.50)
for ((i=0; i<${#periods[@]}; ++i))
do
    export period=${periods[i]}
    root -l -b -q "../PlottingScripts/BucketStructureT1.C(${period})"
    root -l -b -q "../PlottingScripts/BucketStructureT0TDC.C(${period})"
done

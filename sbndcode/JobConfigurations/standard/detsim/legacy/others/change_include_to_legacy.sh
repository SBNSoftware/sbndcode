#! /bin/bash

for f in legacy_detsim*.fcl; do
    echo sed -i s/"standard_detsim_sbnd.fcl"/"legacy_detsim_sbnd.fcl"/ $f
    sed -i s/"standard_detsim_sbnd.fcl"/"legacy_detsim_sbnd.fcl"/ $f

    echo sed -i s/"detsim_/"legacy_detsim_/ $f
    sed -i s/"detsim_/"legacy_detsim_/ $f  
done


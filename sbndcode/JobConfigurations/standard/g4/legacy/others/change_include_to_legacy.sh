#! /bin/bash

for f in legacy_g4*.fcl; do
    echo sed -i s/"standard_g4_sbnd.fcl"/"legacy_g4_sbnd.fcl"/ $f
    sed -i s/"standard_g4_sbnd.fcl"/"legacy_g4_sbnd.fcl"/ $f

    echo sed -i s/"g4_/"legacy_g4_/ $f
    sed -i s/"g4_/"legacy_g4_/ $f  
done


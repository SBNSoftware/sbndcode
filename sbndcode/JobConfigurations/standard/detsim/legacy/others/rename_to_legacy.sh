#! /bin/bash

for f in detsim_*.fcl; do
    echo mv $f legacy_${f}
    mv $f legacy_${f}
done


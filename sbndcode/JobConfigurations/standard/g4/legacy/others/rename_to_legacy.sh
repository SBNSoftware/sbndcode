#! /bin/bash

for f in g4*.fcl; do
    echo mv $f legacy_${f}
    mv $f legacy_${f}
done


#!/bin/bash

#Generate geometry without wires -- t is tpb percent coverage; b corresponds to scintillation bars
./generate_gdml.pl -w 0 -t 1 -b 0 -i lar1-gdml-parameters.xml -o lar1-gdml-fragments.xml
./make_gdml.pl -i lar1-gdml-fragments.xml -o lar1ndvX_nowires.gdml

#Generate geometry with wires
./generate_gdml.pl -w 1 -t 1 -b 0 -i lar1-gdml-parameters.xml -o lar1-gdml-fragments.xml
./make_gdml.pl -i lar1-gdml-fragments.xml -o lar1ndvX.gdml

#Copy both back up one level to be seen by LArSoft jobs
cp lar1ndvX_nowires.gdml ..
cp lar1ndvX.gdml ..

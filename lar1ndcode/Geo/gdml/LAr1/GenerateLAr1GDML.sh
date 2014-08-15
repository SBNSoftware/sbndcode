#!/bin/bash

#Generate geometry without wires
./generate_gdml.pl -w 0 -i lar1-gdml-parameters.xml -o lar1-gdml-fragments.xml
./make_gdml.pl -i lar1-gdml-fragments.xml -o lar1_nowires.gdml

#Generate geometry with wires
./generate_gdml.pl -w 1 -i lar1-gdml-parameters.xml -o lar1-gdml-fragments.xml
./make_gdml.pl -i lar1-gdml-fragments.xml -o lar1.gdml

#Copy both back up one level to be seen by LArSoft jobs
cp lar1_nowires.gdml ..
cp lar1.gdml ..

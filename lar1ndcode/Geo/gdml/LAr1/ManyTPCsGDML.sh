#!/bin/bash

#Generate geometry without wires
./manyTPCs_gdml.pl -w 0 -i lar1-gdml-parameters.xml -o lar1-gdml-fragments.xml
./make_gdml.pl -i lar1-gdml-fragments.xml -o manyTPClar1nd_nowires.gdml

#Generate geometry with wires
./manyTPCs_gdml.pl -w 1 -i lar1-gdml-parameters.xml -o lar1-gdml-fragments.xml
./make_gdml.pl -i lar1-gdml-fragments.xml -o manyTPClar1nd.gdml

#Copy both back up one level to be seen by LArSoft jobs
cp manyTPClar1nd_nowires.gdml ..
cp manyTPClar1nd.gdml ..

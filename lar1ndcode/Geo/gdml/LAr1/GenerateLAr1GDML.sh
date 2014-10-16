#!/bin/bash

#Generate geometry without wires
./generate_gdml.pl -w 0 -i lar1-gdml-parameters.xml -o lar1-gdml-fragments.xml
./make_gdml.pl -i lar1-gdml-fragments.xml -o lar1nd_nowires.gdml

#Generate geometry with wires
./generate_gdml.pl -w 1 -i lar1-gdml-parameters.xml -o lar1-gdml-fragments.xml
./make_gdml.pl -i lar1-gdml-fragments.xml -o lar1nd.gdml

#Generate geometry with 16 pmts, behind wires
#./generate_gdml.pl -w 1 -p16 -i lar1-gdml-parameters.xml -o lar1-gdml-fragments.xml
#./make_gdml.pl -i lar1-gdml-fragments.xml -o lar1nd_16pmts.gdml

#Generate geometry with 32 pmts, behind wires and on sides
#./generate_gdml.pl -w 1 -p32 -i lar1-gdml-parameters.xml -o lar1-gdml-fragments.xml
#./make_gdml.pl -i lar1-gdml-fragments.xml -o lar1nd_32pmts.gdml


#Copy both back up one level to be seen by LArSoft jobs
cp lar1nd_nowires.gdml ..
cp lar1nd.gdml ..

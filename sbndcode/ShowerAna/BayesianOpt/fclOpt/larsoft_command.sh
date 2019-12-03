#!/bin/bash

lar -c /sbnd/data/users/dbarker/showers/Shower/SBNShower/BayesTune/SBNRecoValidation.fcl `for file in $(cat /sbnd/data/users/etyley/showers/reconstruction/MC/tracs/filter/files.txt); do samweb -e sbnd get-file-access-url --schema=root $file; done` -n 100

root -b -q '/sbnd/data/users/dbarker/showers/Shower/SBNShower/ShowerValidation.C({"showervalidationGraphs.root"},"",0,{"SBNShowerNew"})'

return awk '/./{line=$0} END{print line}' histcomp.txt


source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
setup root v6_26_06 -q e26:p3913:prof
g++ -o preparseGDML preparseGDML.cpp `root-config --cflags --glibs`

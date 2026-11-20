#!/bin/bash

# Get the directory of the script
LAMBDAANA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

LAMBDAANA_DIR="$LAMBDAANA_DIR"
LAMBDAANA_SRC="$LAMBDAANA_DIR""/src"
LAMBDAANA_BUILD="$LAMBDAANA_DIR""/build"

# Set an environment variable with the script's directory
export LAMBDAANA_DIR="$LAMBDAANA_DIR"
export LAMBDAANA_SRC="$LAMBDAANA_SRC"
export LAMBDAANA_BUILD="$LAMBDAANA_BUILD"

echo "Top directory: $LAMBDAANA_DIR"
echo "Source directory: $LAMBDAANA_SRC"
echo "Build directory: $LAMBDAANA_BUILD"

alias runtpclines='$LAMBDAANA_BUILD/src/RunAlgoTPCLines'
alias runfrans='$LAMBDAANA_BUILD/src/RunAlgoFRANS'
alias runHandScan='sh $LAMBDAANA_DIR/scripts/runHandScanEvents.sh'
alias runHandScan2='sh $LAMBDAANA_DIR/scripts/runHandScanEvents2.sh'
alias loadHypana='root $LAMBDAANA_DIR/src/Analysis/LoadAnalysisMacros.C'
alias loadLambdaAna='root $LAMBDAANA_DIR/src/Analysis/LoadAnalysisMacros.C'
alias loadLambdaAna2D='root $LAMBDAANA_DIR/src/Analysis/LambdaPlot2DDistributions.C'
alias loadLambdaCalo='root $LAMBDAANA_DIR/src/Analysis/PlotCalorimetryDisplay.C'
alias loadLambdaAnaKinematics='root $LAMBDAANA_DIR/src/Analysis/PlotRecoKinematics.C'

alias loadFRANSAna='root $LAMBDAANA_DIR/src/FRANSAna/LoadFRANSAna.C'


#Search Paths
export LAMBDAANA_SEARCHPATH=$LAMBDAANA_DIR/scripts:.:$LAMBDAANA_DIR/src/config/

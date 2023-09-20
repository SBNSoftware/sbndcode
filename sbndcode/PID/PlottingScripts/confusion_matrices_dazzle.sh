#!/bin/bash

root -l -b -q 'ConfusionMatrixDazzle.C++(true, false, true, -1, -1)'
root -l -b -q 'ConfusionMatrixDazzle.C(true, false, false, -1, -1)'
root -l -b -q 'ConfusionMatrixDazzle.C(false, true, true, -1, -1)'
root -l -b -q 'ConfusionMatrixDazzle.C(false, true, false, -1, -1)'
root -l -b -q 'ConfusionMatrixDazzle.C(true, false, true, .5, .5)'
root -l -b -q 'ConfusionMatrixDazzle.C(true, false, false, .5, .5)'
root -l -b -q 'ConfusionMatrixDazzle.C(false, true, true, .5, .5)'
root -l -b -q 'ConfusionMatrixDazzle.C(false, true, false, .5, .5)'
root -l -b -q 'ConfusionMatrixDazzle.C(true, false, true, .8, .8)'
root -l -b -q 'ConfusionMatrixDazzle.C(true, false, false, .8, .8)'
root -l -b -q 'ConfusionMatrixDazzle.C(false, true, true, .8, .8)'
root -l -b -q 'ConfusionMatrixDazzle.C(false, true, false, .8, .8)'

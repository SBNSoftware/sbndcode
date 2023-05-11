#!/bin/bash

root -l -b -q 'ConfusionMatrixRazzle.C++(true, false, true, -1, -1)'
root -l -b -q 'ConfusionMatrixRazzle.C(true, false, false, -1, -1)'
root -l -b -q 'ConfusionMatrixRazzle.C(false, true, true, -1, -1)'
root -l -b -q 'ConfusionMatrixRazzle.C(false, true, false, -1, -1)'
root -l -b -q 'ConfusionMatrixRazzle.C(true, false, true, .8, .8)'
root -l -b -q 'ConfusionMatrixRazzle.C(true, false, false, .8, .8)'
root -l -b -q 'ConfusionMatrixRazzle.C(false, true, true, .8, .8)'
root -l -b -q 'ConfusionMatrixRazzle.C(false, true, false, .8, .8)'

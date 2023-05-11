#!/bin/bash

root -l -b -q 'ScoreDistributionsDazzle.C++(false, -1, -1)'
root -l -b -q 'ScoreDistributionsDazzle.C(true, -1, -1)'
root -l -b -q 'ScoreDistributionsDazzle.C(false, .8, .8)'
root -l -b -q 'ScoreDistributionsDazzle.C(true, .8, .8)'

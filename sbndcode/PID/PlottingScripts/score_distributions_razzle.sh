#!/bin/bash

root -l -b -q 'ScoreDistributionsRazzle.C++(false, -1, -1)'
root -l -b -q 'ScoreDistributionsRazzle.C(true, -1, -1)'
root -l -b -q 'ScoreDistributionsRazzle.C(false, .8, .8)'
root -l -b -q 'ScoreDistributionsRazzle.C(true, .8, .8)'

root -l -b -q 'ScoreDistributionsRazzle.C++(false, -1, -1, true)'

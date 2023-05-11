#!/bin/bash

root -l -b -q 'ConfusionMatrix.C++("Razzled_Standard", "BDTG", false, true, false, false, -1, -1, false, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_Standard", "BDTG", false, true, false, true, -1, -1, false, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_Standard", "BDTG", false, false, true, true, -1, -1, false, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_Standard", "BDTG", false, true, false, true, -1, -1, true, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_Standard", "BDTG", false, true, false, true, -1, -1, false, true)'
root -l -b -q 'ConfusionMatrix.C("Razzled_Standard", "BDTG", false, false, true, true, -1, -1, true, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_Standard", "BDTG", false, false, true, true, -1, -1, false, true)'

root -l -b -q 'ConfusionMatrix.C("Razzled_Standard", "BDTG", false, true, false, false, .8, .8, false, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_Standard", "BDTG", false, true, false, true, .8, .8, false, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_Standard", "BDTG", false, false, true, true, .8, .8, false, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_Standard", "BDTG", false, true, false, true, .8, .8, true, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_Standard", "BDTG", false, true, false, true, .8, .8, false, true)'
root -l -b -q 'ConfusionMatrix.C("Razzled_Standard", "BDTG", false, false, true, true, .8, .8, true, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_Standard", "BDTG", false, false, true, true, .8, .8, false, true)'

root -l -b -q 'ConfusionMatrix.C("Razzled_with_other_category", "BDTG", true, true, false, false, -1, -1, false, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_with_other_category", "BDTG", true, true, false, true, -1, -1, false, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_with_other_category", "BDTG", true, false, true, true, -1, -1, false, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_with_other_category", "BDTG", true, true, false, true, -1, -1, true, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_with_other_category", "BDTG", true, true, false, true, -1, -1, false, true)'
root -l -b -q 'ConfusionMatrix.C("Razzled_with_other_category", "BDTG", true, false, true, true, -1, -1, true, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_with_other_category", "BDTG", true, false, true, true, -1, -1, false, true)'

root -l -b -q 'ConfusionMatrix.C("Razzled_with_other_category", "BDTG", true, true, false, false, .8, .8, false, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_with_other_category", "BDTG", true, true, false, true, .8, .8, false, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_with_other_category", "BDTG", true, false, true, true, .8, .8, false, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_with_other_category", "BDTG", true, true, false, true, .8, .8, true, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_with_other_category", "BDTG", true, true, false, true, .8, .8, false, true)'
root -l -b -q 'ConfusionMatrix.C("Razzled_with_other_category", "BDTG", true, false, true, true, .8, .8, true, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_with_other_category", "BDTG", true, false, true, true, .8, .8, false, true)'

root -l -b -q 'ConfusionMatrix.C("Razzled_no_track_score_inputs", "BDTG", false, true, false, false, -1, -1, false, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_no_track_score_inputs", "BDTG", false, true, false, true, -1, -1, false, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_no_track_score_inputs", "BDTG", false, false, true, true, -1, -1, false, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_no_track_score_inputs", "BDTG", false, true, false, true, -1, -1, true, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_no_track_score_inputs", "BDTG", false, true, false, true, -1, -1, false, true)'
root -l -b -q 'ConfusionMatrix.C("Razzled_no_track_score_inputs", "BDTG", false, false, true, true, -1, -1, true, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_no_track_score_inputs", "BDTG", false, false, true, true, -1, -1, false, true)'

root -l -b -q 'ConfusionMatrix.C("Razzled_no_track_score_inputs", "BDTG", false, true, false, false, .8, .8, false, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_no_track_score_inputs", "BDTG", false, true, false, true, .8, .8, false, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_no_track_score_inputs", "BDTG", false, false, true, true, .8, .8, false, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_no_track_score_inputs", "BDTG", false, true, false, true, .8, .8, true, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_no_track_score_inputs", "BDTG", false, true, false, true, .8, .8, false, true)'
root -l -b -q 'ConfusionMatrix.C("Razzled_no_track_score_inputs", "BDTG", false, false, true, true, .8, .8, true, false)'
root -l -b -q 'ConfusionMatrix.C("Razzled_no_track_score_inputs", "BDTG", false, false, true, true, .8, .8, false, true)'


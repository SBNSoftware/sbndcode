#!/bin/bash

root -l -b -q 'ScoreDistributions.C++("Razzled_no_track_score", "BDT::BDTG", false, false, -1, -1)'
root -l -b -q 'ScoreDistributions.C("Razzled_no_track_score", "BDT::BDTG", false, true, -1, -1)'
root -l -b -q 'ScoreDistributions.C("Razzled_no_track_score", "BDT::BDTG", false, false, .8, .8)'
root -l -b -q 'ScoreDistributions.C("Razzled_no_track_score", "BDT::BDTG", false, true, .8, .8)'
return
root -l -b -q 'ScoreDistributions.C++("Razzled_standard", "BDT::BDTG", false, false, -1, -1)'
root -l -b -q 'ScoreDistributions.C("Razzled_standard", "BDT::BDTG", false, true, -1, -1)'
root -l -b -q 'ScoreDistributions.C("Razzled_standard", "BDT::BDTG", false, false, .8, .8)'
root -l -b -q 'ScoreDistributions.C("Razzled_standard", "BDT::BDTG", false, true, .8, .8)'

root -l -b -q 'ScoreDistributions.C("Razzled_with_other_category", "BDT::BDTG", true, false, -1, -1)'
root -l -b -q 'ScoreDistributions.C("Razzled_with_other_category", "BDT::BDTG", true, true, -1, -1)'
root -l -b -q 'ScoreDistributions.C("Razzled_with_other_category", "BDT::BDTG", true, false, .8, .8)'
root -l -b -q 'ScoreDistributions.C("Razzled_with_other_category", "BDT::BDTG", true, true, .8, .8)'

root -l -b -q 'ScoreDistributions.C("Razzled_no_track_score_inputs", "BDT::BDTG", false, false, -1, -1)'
root -l -b -q 'ScoreDistributions.C("Razzled_no_track_score_inputs", "BDT::BDTG", false, true, -1, -1)'
root -l -b -q 'ScoreDistributions.C("Razzled_no_track_score_inputs", "BDT::BDTG", false, false, .8, .8)'
root -l -b -q 'ScoreDistributions.C("Razzled_no_track_score_inputs", "BDT::BDTG", false, true, .8, .8)'

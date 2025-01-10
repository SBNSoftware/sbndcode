/**
 * @file    TrackAnalysis.cpp
 * @brief   Does something with the tracks (implementation file).
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    October 23, 2017
 * @see     TrackAnalysis.h
 * 
 */

#include "TrackAnalysis.h"

// LArSoft libraries

// canvas libraries

// ROOT libraries

// C/C++ standard libraries
#include <algorithm> // std::count_if()


TrackAnalysis::TrackAnalysis(fhicl::ParameterSet const& config)
  : fMinLength(config.get<double>("MinLength", 0.0))
  {}
  
void TrackAnalysis::prepare() {
  if (fDir) {
    fHNTracks = std::make_unique<TH1F>
      ("HNTracks", "Number of tracks;number of tracks;events", 50, 0., 50.);
    fHNTracks->SetDirectory(fDir);
  }
  else fHNTracks.reset();
} // TrackAnalysis::prepare()


void TrackAnalysis::processTracks(std::vector<recob::Track> const& tracks) {
  
  unsigned int nTracksAboveThreshold = std::count_if(
    tracks.begin(), tracks.end(),
    [this](auto const& track){ return track.Length() > fMinLength; }
    );
  
  if (fHNTracks) fHNTracks->Fill(nTracksAboveThreshold);
  
} // TrackAnalysis::processTracks()


void TrackAnalysis::finish() {
  if (fHNTracks) {
    fDir->cd();
    fHNTracks->Write();
  }
} // TrackAnalysis::finish()


/**
 * @file    TrackAnalysis.h
 * @brief   Does something with the tracks.
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    October 23, 2017
 * @see     galleryAnalysis.cpp
 * 
 */

#ifndef TRACKANALYSIS_H
#define TRACKANALYSIS_H

// LArSoft libraries
#include "lardataobj/RecoBase/Track.h"
#include "larcorealg/Geometry/GeometryCore.h"

// canvas libraries
#include "fhiclcpp/ParameterSet.h"

// ROOT libraries
#include "TDirectory.h"
#include "TH1.h"

// C/C++ standard libraries
#include <vector>
#include <memory> // std::unique_ptr<>


/**
 * @brief Track analysis example.
 * 
 * Configuration
 * --------------
 * 
 * * *MinLength* (real, default: 0): minimum track length, in centimetres
 * 
 */
class TrackAnalysis {
  
  geo::GeometryCore const* fGeom = nullptr;
  TDirectory* fDir = nullptr;
  double fMinLength = 0.0;                  ///< Minimum length, in centimetres.
  std::unique_ptr<TH1> fHNTracks;
  
    public:
  
  TrackAnalysis(fhicl::ParameterSet const& config);
  
  void setup(geo::GeometryCore const& geom, TDirectory* outDir)
    {
      fGeom = &geom; 
      fDir = outDir;
    }
  
  void prepare();
  
  void processTracks(std::vector<recob::Track> const& tracks);
  
  void finish();
  
}; // class TrackAnalysis



#endif // TRACKANALYSIS_H

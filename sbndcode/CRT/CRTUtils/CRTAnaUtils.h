#ifndef CRTANAUTILS_H_SEEN
#define CRTANAUTILS_H_SEEN


///////////////////////////////////////////////
// CRTAnaUtils.h
//
// Reco utilities for doing CRT reconstruction in ana modules
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"

// c++
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <cmath> 
#include <memory>

// ROOT
#include "TVector3.h"
#include "TTree.h"


namespace CRTAnaUtils{

  std::vector<std::vector<sbnd::crt::CRTHit>> CreateCRTTzeros(std::vector<sbnd::crt::CRTHit> crtHits, double fTimeLimit);

  std::vector<sbnd::crt::CRTTrack> CreateCRTTracks(std::vector<std::vector<sbnd::crt::CRTHit>> crtTzeros, 
                                                   double fAverageHitDist, bool fUseTopPlane, double fDistanceLimit);

  std::vector<sbnd::crt::CRTTrack> CreateCRTTracks(std::vector<sbnd::crt::CRTHit> crtHits, double fTimeLimit, 
                                                   double fAverageHitDist, bool fUseTopPlane, double fDistanceLimit);

  double T0FromCRTHits(recob::Track tpcTrack, std::vector<sbnd::crt::CRTHit> crtHits, int tpc, double fMinTrackLength, 
                       double fTrackDirectionFrac, double fDistanceLimit);

  double T0FromCRTTracks(recob::Track tpcTrack, std::vector<sbnd::crt::CRTTrack> crtTracks, int tpc, double fMaxAngleDiff, 
                         double fMaxDistance);

}

#endif

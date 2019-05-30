#ifndef COSMICREMOVALUTILS_H_SEEN
#define COSMICREMOVALUTILS_H_SEEN


///////////////////////////////////////////////
// CosmicRemovalUtils.h
//
// Reco utilities for doing cosmic removal in ana modules
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
#include "art_root_io/TFileService.h" 
#include "art_root_io/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

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
#include "TGraph.h"
#include "TF1.h"


namespace sbnd{
namespace CosmicRemovalUtils{

  bool InFiducial(geo::Point_t point, double fiducial);

  bool InFiducial(geo::Point_t point, double fiducial, double fiducialTop);

  bool InFiducial(geo::Point_t point, double minXCut, double minYCut, double minZCut, double maxXCut, double maxYCut, double maxZCut);

  int DetectedInTPC(std::vector<art::Ptr<recob::Hit>> hits);

  std::pair<std::vector<double>, std::vector<double>> FakeTpcFlashes(std::vector<simb::MCParticle> particles);

  bool BeamFlash(std::vector<double> flashes, double beamTimeLimit);
  
}
}

#endif

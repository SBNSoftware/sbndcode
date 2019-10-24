#ifndef CRT_TPCGEOUTIL_H_
#define CRT_TPCGEOUTIL_H_
// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Hit.h"

// c++
#include <vector>

namespace sbnd {
namespace TPCGeoUtil {
  int DetectedInTPC(std::vector<art::Ptr<recob::Hit>> hits);
  // Work out the drift limits for a collection of hits
  std::pair<double, double> XLimitsFromHits(const geo::GeometryCore *GeometryService, std::vector<art::Ptr<recob::Hit>> hits);
  // Is point inside given TPC
  bool InsideTPC(geo::Point_t point, const geo::TPCGeo& tpc, double buffer);
  int DriftDirectionFromHits(const geo::GeometryCore *GeometryService, std::vector<art::Ptr<recob::Hit>> hits);
} // namespace TPCGeoUtil
} // namespace sbnd
#endif

#ifndef CRT_TPCGEOUTIL_H_
#define CRT_TPCGEOUTIL_H_
// framework
#include "canvas/Persistency/Common/Ptr.h"

// LArSoft
#include "larcorealg/Geometry/fwd.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Hit.h"

// c++
#include <vector>

namespace sbnd::TPCGeoUtil {
  int DetectedInTPC(std::vector<art::Ptr<recob::Hit>> const& hits);
  // Work out the drift limits for a collection of hits
  std::pair<double, double> XLimitsFromHits(const geo::GeometryCore *GeometryService,
                                            std::vector<art::Ptr<recob::Hit>> const& hits);
  // Is point inside given TPC
  bool InsideTPC(geo::Point_t point, const geo::TPCGeo& tpc, double buffer);
  int DriftDirectionFromHits(const geo::GeometryCore *GeometryService,
                             std::vector<art::Ptr<recob::Hit>> const& hits);
} // namespace sbnd::TPCGeoUtil
#endif

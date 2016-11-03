///////////////////////////////////////////////////////////////////////////////
/// \file AuxDetChannelMapSBNDAlg.h
/// \brief Algorithm class for SBND auxiliary detector channel mapping
///
/// Ported from AuxDetChannelMapLArIATAlg.h (Author: brebel@fnal.gov)
///
/// \version $Id:  $
/// \author mastbaum@uchicago.edu
///////////////////////////////////////////////////////////////////////////////

#ifndef SBND_AuxDetChannelMapSBNDAlg_h
#define SBND_AuxDetChannelMapSBNDAlg_h

#include <vector>
#include <set>
#include "larcore/Geometry/AuxDetChannelMapAlg.h"
#include "sbndcode/Geo/AuxDetGeoObjectSorterSBND.h"
#include "fhiclcpp/ParameterSet.h"
#include "TVector3.h"

namespace geo {

  class AuxDetChannelMapSBNDAlg : public AuxDetChannelMapAlg {
  public:
    AuxDetChannelMapSBNDAlg(fhicl::ParameterSet const& p);

    void Initialize(AuxDetGeometryData_t& geodata) override;

    void Uninitialize();

    uint32_t PositionToAuxDetChannel(
        double const worldLoc[3],
        std::vector<geo::AuxDetGeo*> const& auxDets,
        size_t& ad,
        size_t& sv) const;

    const TVector3 AuxDetChannelToPosition(
        uint32_t const& channel,
        std::string const& auxDetName,
        std::vector<geo::AuxDetGeo*> const& auxDets) const;

  private:
    geo::AuxDetGeoObjectSorterSBND fSorter; ///< Class to sort geo objects
  };

}  // namespace geo

#endif  // SBND_AuxDetChannelMapSBNDAlg_h


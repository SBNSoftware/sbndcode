///////////////////////////////////////////////////////////////////////////////
/// \file AuxDetGeoObjectSorterSBND.h
/// \brief Interface to algorithm class for sorting of AuxDetGeo objects
///
/// Ported from AuxDetGeoObjectSorterLArIAT.h (Author: brebel@fnal.gov)
///
/// \version $Id:  $
/// \author mastbaum@uchicago.edu
///////////////////////////////////////////////////////////////////////////////

#ifndef SBND_AuxDetGeoObjectSorterSBND_h
#define SBND_AuxDetGeoObjectSorterSBND_h

#include <vector>

#include "larcore/Geometry/AuxDetGeoObjectSorter.h"

namespace geo {

  class AuxDetGeoObjectSorterSBND : public AuxDetGeoObjectSorter {
  public:

    AuxDetGeoObjectSorterSBND(fhicl::ParameterSet const& p);

    ~AuxDetGeoObjectSorterSBND();

    void SortAuxDets (std::vector<geo::AuxDetGeo*>& adgeo) const;
    void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*>& adsgeo) const;

  };

}

#endif  // SBND_AuxDetGeoObjectSorterSBND_h


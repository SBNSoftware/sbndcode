///////////////////////////////////////////////////////////////////////////////
/// \file CRTGeoObjectSorter.h
/// \brief Interface to algorithm class for sorting of AuxDetGeo objects
///
/// Ported from AuxDetGeoObjectSorterLArIAT.h (Author: brebel@fnal.gov)
///
/// \version $Id:  $
/// \author mastbaum@uchicago.edu
///////////////////////////////////////////////////////////////////////////////

#ifndef SBND_CRTGeoObjectSorter_h
#define SBND_CRTGeoObjectSorter_h

#include <vector>

#include "larcore/Geometry/AuxDetGeoObjectSorter.h"

namespace geo {

  class CRTGeoObjectSorter : public AuxDetGeoObjectSorter {
  public:

    CRTGeoObjectSorter(fhicl::ParameterSet const& p);

    ~CRTGeoObjectSorter();

    void SortAuxDets (std::vector<geo::AuxDetGeo*>& adgeo) const;
    void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*>& adsgeo) const;

  };

}

#endif  // SBND_CRTGeoObjectSorter_h


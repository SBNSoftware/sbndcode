////////////////////////////////////////////////////////////////////////
/// \file CRTGeoObjectSorter.cxx
/// \brief Interface to algorithm class for sorting of AuxDetGeo objects
///
/// Ported from AuxDetGeoObjectSorterLArIAT.cxx (Author: brebel@fnal.gov)
///
/// \version $Id:  $
/// \author mastbaum@uchicago.edu
////////////////////////////////////////////////////////////////////////

#include "sbndcode/CRT/CRTGeoObjectSorter.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"

namespace geo {

  //----------------------------------------------------------------------------
  // Define sort order for AuxDets in standard configuration
  static bool sortAuxDetSBND(const AuxDetGeo* ad1, const AuxDetGeo* ad2) {
    // Sort using the center of the detector - primary ordering by z,
    // then y, and x
    double c1[3] = {0, 0, 0};
    double c2[3] = {0, 0, 0};
    ad1->GetCenter(c1);
    ad2->GetCenter(c2);

    for (int i=2; i>0; i--) {
      if (c1[i] != c2[i]){
	return c1[i] < c2[i];
      }
    }

    return c1[0] < c2[0];

  }

  //----------------------------------------------------------------------------
  // Define sort order for AuxDetSensitives in standard configuration
  static bool sortAuxDetSensitiveSBND(const AuxDetSensitiveGeo* ad1,
                                      const AuxDetSensitiveGeo* ad2)
  {
    // Sort using the center of the detector - primary ordering by z,
    // then y, and x
    double c1[3] = {0, 0, 0};
    double c2[3] = {0, 0, 0};
    ad1->GetCenter(c1);
    ad2->GetCenter(c2);

    for (int i=2; i>0; i--) {
      if (c1[i] != c2[i]) {
	return c1[i] < c2[i];
      }
    }

    return c1[0] < c2[0];
  }

  //----------------------------------------------------------------------------
  CRTGeoObjectSorter::CRTGeoObjectSorter(
      fhicl::ParameterSet const&) {}

  //----------------------------------------------------------------------------
  CRTGeoObjectSorter::~CRTGeoObjectSorter() {}

  //----------------------------------------------------------------------------
  void CRTGeoObjectSorter::SortAuxDets(
      std::vector<geo::AuxDetGeo*> & adgeo) const {
    std::sort(adgeo.begin(), adgeo.end(), sortAuxDetSBND);
  }

  //----------------------------------------------------------------------------
  void CRTGeoObjectSorter::SortAuxDetSensitive(
      std::vector<geo::AuxDetSensitiveGeo*> & adsgeo) const {
    std::sort(adsgeo.begin(), adsgeo.end(), sortAuxDetSensitiveSBND);
  }

}  // namespace geo


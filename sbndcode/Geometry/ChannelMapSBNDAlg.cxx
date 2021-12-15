/**
 * @file ChannelMappingSBNDAlg.cxx
 *
 * Overrides ChannelMapAlg in larcorealg to change exceptions to messages
 */

#include "sbndcode/Geometry/ChannelMapSBNDAlg.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace geo {
  
  size_t ChannelMapSBNDAlg::NearestAuxDet(const double* point, std::vector<geo::AuxDetGeo> const& auxDets, double tolerance) const
  {
    double HalfCenterWidth = 0.;
    double localPoint[3] = {0.};


    for(size_t a = 0; a < auxDets.size(); ++a) {

      auxDets[a].WorldToLocal(point, localPoint);

      HalfCenterWidth = 0.5 * (auxDets[a].HalfWidth1() + auxDets[a].HalfWidth2());
      
      if( localPoint[2] >= - (auxDets[a].Length()/2 + tolerance) &&
          localPoint[2] <=   (auxDets[a].Length()/2 + tolerance) &&
          localPoint[1] >= - auxDets[a].HalfHeight() - tolerance &&
          localPoint[1] <=   auxDets[a].HalfHeight() + tolerance &&
          // if AuxDet a is a box, then HalfSmallWidth = HalfWidth
          localPoint[0] >= - HalfCenterWidth + localPoint[2]*(HalfCenterWidth - auxDets[a].HalfWidth2())/(0.5 * auxDets[a].Length()) - tolerance &&
          localPoint[0] <=   HalfCenterWidth - localPoint[2]*(HalfCenterWidth - auxDets[a].HalfWidth2())/(0.5 * auxDets[a].Length()) + tolerance
          ) return a;

    }// for loop over AudDet a

    // log a message because we couldn't find the aux det volume, exception in base class
    mf::LogDebug("ChannelMapSBND") << "Can't find AuxDet for position ("
                                   << point[0] << ","
                                   << point[1] << ","
                                   << point[2] << ")\n";
    
    return UINT_MAX;

  }

  //----------------------------------------------------------------------------
  size_t ChannelMapSBNDAlg::NearestSensitiveAuxDet(const double* point, std::vector<geo::AuxDetGeo> const& auxDets, double tolerance) const
  {
    double HalfCenterWidth = 0.;
    double localPoint[3] = {0.};

    size_t auxDetIdx = this->NearestAuxDet(point, auxDets, tolerance);

    geo::AuxDetGeo const& adg = auxDets[auxDetIdx];

    for(size_t a = 0; a < adg.NSensitiveVolume(); ++a) {

      geo::AuxDetSensitiveGeo const& adsg = adg.SensitiveVolume(a);
      adsg.WorldToLocal(point, localPoint);

      HalfCenterWidth = 0.5 * (adsg.HalfWidth1() + adsg.HalfWidth2());

      if( localPoint[2] >= - (adsg.Length()/2 + tolerance) &&
          localPoint[2] <=   (adsg.Length()/2 + tolerance) &&
          localPoint[1] >= - adsg.HalfHeight() - tolerance &&
          localPoint[1] <=   adsg.HalfHeight() + tolerance &&
          // if AuxDet a is a box, then HalfSmallWidth = HalfWidth
          localPoint[0] >= - HalfCenterWidth + localPoint[2]*(HalfCenterWidth - adsg.HalfWidth2())/(0.5 * adsg.Length()) - tolerance &&
          localPoint[0] <=   HalfCenterWidth - localPoint[2]*(HalfCenterWidth - adsg.HalfWidth2())/(0.5 * adsg.Length()) + tolerance 
          ) return a;
    }// for loop over AuxDetSensitive a

    // log a message because we couldn't find the sensitive aux det volume, exception in base class
    mf::LogDebug("ChannelMapSBND") << "Can't find AuxDetSensitive for position ("
                                   << point[0] << ","
                                   << point[1] << ","
                                   << point[2] << ")\n";
    
    return UINT_MAX;
  }

} // namespace geo

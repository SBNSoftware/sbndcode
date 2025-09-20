/**
 * @file ChannelMappingSBNDAlg.cxx
 *
 * Overrides ChannelMapAlg in larcorealg to change exceptions to messages
 */

#include "sbndcode/Geometry/ChannelMapSBNDAlg.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace geo {
  
  size_t ChannelMapSBNDAlg::NearestAuxDet(Point_t const& point, std::vector<geo::AuxDetGeo> const& auxDets, double tolerance) const
  {
    double HalfCenterWidth = 0.;


    for(size_t a = 0; a < auxDets.size(); ++a) {

      auto const localPoint = auxDets[a].toLocalCoords(point);

      HalfCenterWidth = 0.5 * (auxDets[a].HalfWidth1() + auxDets[a].HalfWidth2());
      
      if( localPoint.Z() >= - (auxDets[a].Length()/2 + tolerance) &&
          localPoint.Z() <=   (auxDets[a].Length()/2 + tolerance) &&
          localPoint.Y() >= - auxDets[a].HalfHeight() - tolerance &&
          localPoint.Y() <=   auxDets[a].HalfHeight() + tolerance &&
          // if AuxDet a is a box, then HalfSmallWidth = HalfWidth
          localPoint.X() >= - HalfCenterWidth + localPoint.Z()*(HalfCenterWidth - auxDets[a].HalfWidth2())/(0.5 * auxDets[a].Length()) - tolerance &&
          localPoint.X() <=   HalfCenterWidth - localPoint.Z()*(HalfCenterWidth - auxDets[a].HalfWidth2())/(0.5 * auxDets[a].Length()) + tolerance
          ) return a;

    }// for loop over AudDet a

    // log a message because we couldn't find the aux det volume, exception in base class
    mf::LogDebug("ChannelMapSBND") << "Can't find AuxDet for position ("
                                   << point.X() << ","
                                   << point.Y() << ","
                                   << point.Z() << ")\n";
    
    return UINT_MAX;

  }

  //----------------------------------------------------------------------------
  size_t ChannelMapSBNDAlg::NearestSensitiveAuxDet(Point_t const& point, std::vector<geo::AuxDetGeo> const& auxDets, double tolerance) const
  {
    double HalfCenterWidth = 0.;

    size_t auxDetIdx = this->NearestAuxDet(point, auxDets, tolerance);

    if(auxDetIdx == UINT_MAX)
      return UINT_MAX;

    geo::AuxDetGeo const& adg = auxDets[auxDetIdx];

    for(size_t a = 0; a < adg.NSensitiveVolume(); ++a) {

      geo::AuxDetSensitiveGeo const& adsg = adg.SensitiveVolume(a);
      auto const localPoint = adsg.toLocalCoords(point);

      HalfCenterWidth = 0.5 * (adsg.HalfWidth1() + adsg.HalfWidth2());

      if( localPoint.Z() >= - (adsg.Length()/2 + tolerance) &&
          localPoint.Z() <=   (adsg.Length()/2 + tolerance) &&
          localPoint.Y() >= - adsg.HalfHeight() - tolerance &&
          localPoint.Y() <=   adsg.HalfHeight() + tolerance &&
          // if AuxDet a is a box, then HalfSmallWidth = HalfWidth
          localPoint.X() >= - HalfCenterWidth + localPoint.Z()*(HalfCenterWidth - adsg.HalfWidth2())/(0.5 * adsg.Length()) - tolerance &&
          localPoint.X() <=   HalfCenterWidth - localPoint.Z()*(HalfCenterWidth - adsg.HalfWidth2())/(0.5 * adsg.Length()) + tolerance
          ) return a;
    }// for loop over AuxDetSensitive a

    // log a message because we couldn't find the sensitive aux det volume, exception in base class
    mf::LogDebug("ChannelMapSBND") << "Can't find AuxDetSensitive for position ("
                                   << point.X() << ","
                                   << point.Y() << ","
                                   << point.Z() << ")\n";
    
    return UINT_MAX;
  }

} // namespace geo

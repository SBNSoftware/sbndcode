#ifndef __FLASHFINDERFMWKINTERFACE_CXX__
#define __FLASHFINDERFMWKINTERFACE_CXX__

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "larcore/Geometry/Geometry.h"
#include "FlashFinderFMWKInterface.h"

namespace lightana {

  std::vector<size_t> ListOpChannels(int cryostat) {
    std::vector<size_t> res;
    ::art::ServiceHandle<geo::Geometry> geo;
    if(cryostat<0) {
      for(size_t opch=0; opch<geo->MaxOpChannel(); ++opch) {
        if(geo->IsValidOpChannel(opch)) continue;
        res.push_back(opch);
      }
    }else{
      auto const& bbox = geo->Cryostat(cryostat).Boundaries();
      for(size_t opch=0; opch<geo->MaxOpChannel(); ++opch) {
        if(geo->IsValidOpChannel(opch)) continue;
        auto const& pt = geo->OpDetGeoFromOpChannel(opch).GetCenter();
        if(!bbox.ContainsPosition(pt)) continue;
        res.push_back(opch);
      }
    }
    return res;
  }

  std::vector<size_t> ListOpDets(int cryostat) {
    std::vector<size_t> res;
    ::art::ServiceHandle<geo::Geometry> geo;
    if(cryostat<0) {
      for(size_t opdet=0; opdet<geo->NOpDets(); ++opdet) {
        res.push_back(opdet);
      }
    }else{
      auto const& bbox = geo->Cryostat(cryostat).Boundaries();
      for(size_t opdet=0; opdet<geo->NOpDets(); ++opdet) {
        auto const& pt = geo->OpDetGeoFromOpDet(opdet).GetCenter();
        if(!bbox.ContainsPosition(pt)) continue;
        res.push_back(opdet);
      }
    }
    return res;
  }

  size_t NOpDets(int cryostat) {
    ::art::ServiceHandle<geo::Geometry> geo;
    if(cryostat<0)
      return geo->NOpDets();
    else
      return geo->Cryostat(cryostat).NOpDet();
  }

  size_t OpDetFromOpChannel(size_t opch) {
    ::art::ServiceHandle<geo::Geometry> geo;
    return geo->OpDetFromOpChannel(opch);
  }

  void OpDetCenterFromOpChannel(size_t opch, double *xyz) {
    ::art::ServiceHandle<geo::Geometry> geo;
    geo->OpDetGeoFromOpChannel(opch).GetCenter(xyz); 
  }

}
#endif


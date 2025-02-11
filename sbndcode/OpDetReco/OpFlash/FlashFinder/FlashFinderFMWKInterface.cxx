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
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "FlashFinderFMWKInterface.h"

namespace lightana {

  std::vector<size_t> ListOpChannels(int cryostat) {
    std::vector<size_t> res;
    ::art::ServiceHandle<geo::Geometry> geo;
    auto const& channelMapAlg =
      art::ServiceHandle<geo::WireReadout const>()->Get();
    if(cryostat<0) {
      for(size_t opch=0; opch<channelMapAlg.MaxOpChannel(); ++opch) {
        if(channelMapAlg.IsValidOpChannel(opch)) continue;
        res.push_back(opch);
      }
    }else{
      auto const& bbox = geo->Cryostat(geo::CryostatID(cryostat)).Boundaries();
      for(size_t opch=0; opch<channelMapAlg.MaxOpChannel(); ++opch) {
        if(channelMapAlg.IsValidOpChannel(opch)) continue;
        auto const& pt = channelMapAlg.OpDetGeoFromOpChannel(opch).GetCenter();
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
      auto const& bbox = geo->Cryostat(geo::CryostatID(cryostat)).Boundaries();
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
      return geo->Cryostat(geo::CryostatID(cryostat)).NOpDet();
  }

  std::vector<size_t> ListOpChannelsByTPC(int tpc) {
    std::vector<size_t> res;
    ::art::ServiceHandle<geo::Geometry> geo;
    auto const& channelMapAlg =
      art::ServiceHandle<geo::WireReadout const>()->Get();
    if(tpc<0) {
      for(size_t opch=0; opch<channelMapAlg.MaxOpChannel(); ++opch) {
        if(channelMapAlg.IsValidOpChannel(opch)) continue;
        res.push_back(opch);
      }
    }else{
      // auto const& bbox = geo->TPC(tpc).BoundingBox();
      for(size_t opch=0; opch<channelMapAlg.MaxOpChannel(); ++opch) {
        if(!channelMapAlg.IsValidOpChannel(opch)) continue;
        auto const& pt = channelMapAlg.OpDetGeoFromOpChannel(opch).GetCenter();
        // std::cout << "pt: " << pt.X() << ", " << pt.Y() << ", " << pt.Z() << std::endl;
        // if(!bbox.ContainsPosition(pt)) continue;
        if(pt.X() < 0 && tpc == 0) {   
          res.push_back(opch);
        }
        if(pt.X() > 0 && tpc == 1) {   
          res.push_back(opch);
        }
      }
    }
    return res;
  }

  std::vector<int> PDNamesToList(std::vector<std::string> pd_names) {

    std::vector<int> out_ch_v;

    opdet::sbndPDMapAlg pds_map;

    for (auto name : pd_names) {
      auto ch_v = pds_map.getChannelsOfType(name);
      out_ch_v.insert(out_ch_v.end(), ch_v.begin(), ch_v.end());
    }

    return out_ch_v;

  }

  size_t OpDetFromOpChannel(size_t opch) {
    auto const& channelMapAlg =
      art::ServiceHandle<geo::WireReadout const>()->Get();
    return channelMapAlg.OpDetFromOpChannel(opch);
  }

  void OpDetCenterFromOpChannel(size_t opch, double *xyz) {
    auto const& channelMapAlg =
      art::ServiceHandle<geo::WireReadout const>()->Get();
    auto const tmp = channelMapAlg.OpDetGeoFromOpChannel(opch).GetCenter();
    xyz[0] = tmp.X();
    xyz[1] = tmp.Y();
    xyz[2] = tmp.Z();
  }

}
#endif

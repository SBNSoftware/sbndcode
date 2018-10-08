///////////////////////////////////////////////////////////////////////////////
/// \file CRTChannelMapAlg.cxx
/// \brief Algorithm class for SBND auxiliary detector channel mapping
///
/// Ported from AuxDetChannelMapLArIATAlg.cxx (Author: brebel@fnal.gov)
///
/// \version $Id:  $
/// \author mastbaum@uchicago.edu
///////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "sbndcode/CRT/CRTChannelMapAlg.h"
#include "TVector3.h"
#include <iostream>
#include <ostream>

namespace geo {

  //---------------------------------------------------------------------------
  CRTChannelMapAlg::CRTChannelMapAlg(
      fhicl::ParameterSet const& p)
    : fSorter(geo::CRTGeoObjectSorter(p)) {}

  //---------------------------------------------------------------------------
  void CRTChannelMapAlg::Initialize(AuxDetGeometryData_t& geodata) {
    Uninitialize();

    std::vector<geo::AuxDetGeo*>& adgeo = geodata.auxDets;

    // Sort the AuxDetGeo objects and map them to names of the detectors
    //fSorter.SortAuxDets(adgeo);

    // Map the AuxDetGeo names to their position in the sorted vector
    //
    // Each tagger is composed of scintillator modules, composed of 16 strips.
    // In the geometry, CRTStripArrays are AuxDets and CRTStrips are the
    // AuxDetSensitives. Each strip has two SiPM channels, one per optical
    // fiber (not in the geometry).
    //
    // 2x TaggerTop:  (5x2  +5x2)x16 = 320 strips,  640 channels
    // 2x TaggerSide: (4x2  +5x2)x16 = 288 strips,  576 channels
    // 1x TaggerFace: (4x2  +4x2)x16 = 256 strips,  512 channels
    // 1x TaggerFace: (4x2-1+4x2)x16 = 240 strips,  480 channels
    // 1x TaggerBot:           34x16 = 544 strips, 1088 channels

    fADGeoToName.clear();
    fADGeoToChannelAndSV.clear();

    for (size_t a=0; a<adgeo.size(); a++){
      std::string volName(adgeo[a]->TotalVolume()->GetName());

      size_t nsv = adgeo[a]->NSensitiveVolume();
      if (nsv != 16) {
        throw cet::exception("CRTChannelMap")
        << "Wrong number of sensitive volumes for CRT volume "
        << volName << " (got " << nsv << ", expected 16)" << std::endl;
      }

      fADGeoToName[a] = volName;
      fNameToADGeo[volName] = a;

      if (volName.find("CRTStripArray") != std::string::npos) {
        for (size_t svID=0; svID<16; svID++) {
          for (size_t ich=0; ich<2; ich++) {
            size_t chID = 2 * svID + ich;
            fADGeoToChannelAndSV[a].push_back(std::make_pair(chID, svID));
          }
        }
      }
    }
  }

  //----------------------------------------------------------------------------
  void CRTChannelMapAlg::Uninitialize() {}

  //----------------------------------------------------------------------------
  uint32_t CRTChannelMapAlg::PositionToAuxDetChannel(
      double const worldLoc[3],
      std::vector<geo::AuxDetGeo*> const& auxDets,
      size_t& ad,
      size_t& sv) const {

    // Set the default to be that we don't find the position in any AuxDet
    uint32_t channel = UINT_MAX;

    // Figure out which detector we are in
    ad = 0;
    sv = this->NearestSensitiveAuxDet(worldLoc, auxDets, ad);

    // Get the origin of the sensitive volume in the world coordinate system
    double svOrigin[3] = {0, 0, 0};
    double localOrigin[3] = {0, 0, 0};

    auxDets[ad]->SensitiveVolume(sv).LocalToWorld(localOrigin, svOrigin);

    // Check to see which AuxDet this position corresponds to
    auto gnItr = fADGeoToName.find(ad);
    if (gnItr != fADGeoToName.end()){
      // Get the vector of channel and sensitive volume pairs
      auto csvItr = fADGeoToChannelAndSV.find(ad);

      if (csvItr == fADGeoToChannelAndSV.end()) {
        throw cet::exception("CRTChannelMapAlg")
        << "No entry in channel and sensitive volume map for AuxDet index "
        << ad;
      }

      // N.B. This is the ID on the nth channel, and the strip has n and n+1
      channel = 2 * sv + 0;
    }

    if (channel == UINT_MAX) {
      throw cet::exception("CRTChannelMapAlg")
      << "position ("
      << worldLoc[0] << "," << worldLoc[1] << "," << worldLoc[2]
      << ") does not correspond to any AuxDet";
    }

    return channel;
  }

  //----------------------------------------------------------------------------
  const TVector3 CRTChannelMapAlg::AuxDetChannelToPosition(
      uint32_t const& channel,
      std::string const& auxDetName,
      std::vector<geo::AuxDetGeo*> const& auxDets) const {
    double x = 0;
    double y = 0;
    double z = 0;

    // Figure out which detector we are in
    size_t ad = UINT_MAX;
    if (fNameToADGeo.count(auxDetName) > 0) {
      ad = fNameToADGeo.find(auxDetName)->second;
    }
    else {
      throw cet::exception("CRTChannelMapAlg")
      << "No AuxDetGeo with name " << auxDetName;
    }

    // Get the vector of channel and sensitive volume pairs
    auto csvItr = fADGeoToChannelAndSV.find(ad);

    if (csvItr == fADGeoToChannelAndSV.end()) {
      throw cet::exception("CRTChannelMapAlg")
      << "No entry in channel and sensitive volume"
      << " map for AuxDet index " << ad << " bail";
    }

    // Loop over the vector of channel and sensitive volumes to determine the
    // sensitive volume for this channel. Then get the origin of the sensitive
    // volume in the world coordinate system.
    double svOrigin[3] = {0, 0, 0};
    double localOrigin[3] = {0, 0, 0};
    for (auto csv : csvItr->second) {
      if (csv.first == channel) {
        // Get the center of the sensitive volume for this channel
        auxDets[ad]->SensitiveVolume(csv.second).LocalToWorld(localOrigin,
                                                              svOrigin);

        x = svOrigin[0];
        y = svOrigin[1];
        z = svOrigin[2];

        break;
      }
    }

    return TVector3(x, y, z);
  }

}  // namespace geo


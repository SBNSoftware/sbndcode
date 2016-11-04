///////////////////////////////////////////////////////////////////////////////
/// \file AuxDetChannelMapSBNDAlg.cxx
/// \brief Algorithm class for SBND auxiliary detector channel mapping
///
/// Ported from AuxDetChannelMapLArIATAlg.cxx (Author: brebel@fnal.gov)
///
/// \version $Id:  $
/// \author mastbaum@uchicago.edu
///////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetSensitiveGeo.h"
#include "larcore/Geometry/AuxDetGeometryCore.h"
#include "sbndcode/Geo/AuxDetChannelMapSBNDAlg.h"
#include <iostream>
#include <ostream>

namespace geo {

  void throwBadSDCountException(std::string& volName, size_t n, size_t ne) {
    throw cet::exception("AuxDetChannelMapSBND")
    << "Wrong number of sensitive volumes for CRT volume "
    << volName << " (got " << n << ", expected " << ne << ")" << std::endl;
  }

  //---------------------------------------------------------------------------
  AuxDetChannelMapSBNDAlg::AuxDetChannelMapSBNDAlg(
      fhicl::ParameterSet const& p)
    : fSorter(geo::AuxDetGeoObjectSorterSBND(p)) {}

  //---------------------------------------------------------------------------
  void AuxDetChannelMapSBNDAlg::Initialize(AuxDetGeometryData_t& geodata) {
    Uninitialize();

    std::vector<geo::AuxDetGeo*>& adgeo = geodata.auxDets;

    // Sort the AuxDetGeo objects and map them to names of the detectors
    fSorter.SortAuxDets(adgeo);

    // Map the AuxDetGeo names to their position in the sorted vector
    //
    // Each tagger is composed of scintillator modules, composed of 16 strips.
    // Each strip has two SiPM channels, one per optical fiber.
    //
    // 2x TaggerTop:  (5x2+5x2)x16 = 320 strips,  640 channels
    // 2x TaggerSide: (4x2+5x2)x16 = 288 strips,  576 channels
    // 2x TaggerFace: (4x2+4x2)x16 = 256 strips,  512 channels
    // 1x TaggerBot:         34x16 = 544 strips, 1088 channels
    //
    fADGeoToName.clear();
    fADGeoToChannelAndSV.clear();

    for (size_t a=0; a<adgeo.size(); a++){
      std::string volName(adgeo[a]->TotalVolume()->GetName());

      std::cout << ">> " << volName
                << " " << adgeo[a]->NSensitiveVolume() << std::endl;

      fADGeoToName[a] = volName;
      fNameToADGeo[volName] = a;

      double origin[3] = {0,0,0};
      adgeo[a]->GetCenter(origin);
      std::cout << "origin: "
                << origin[0] << " " << origin[1] << " " << origin[2]
                << std::endl;

      // For now, one channel per module
      if (volName.find("TaggerTop") != std::string::npos) {
        for (size_t c=0; c<20; c++) {
          fADGeoToChannelAndSV[a].push_back(std::make_pair(c, c));
        }
      }
      else if (volName.find("TaggerSide") != std::string::npos) {
        for (size_t c=0; c<18; c++) {
          fADGeoToChannelAndSV[a].push_back(std::make_pair(c, c));
        }
      }
      else if (volName.find("TaggerFace") != std::string::npos) {
        for (size_t c=0; c<16; c++) {
          fADGeoToChannelAndSV[a].push_back(std::make_pair(c, c));
        }
      }
      else if (volName.find("TaggerBot") != std::string::npos) {
        for (size_t c=0; c<34; c++) {
          fADGeoToChannelAndSV[a].push_back(std::make_pair(c, c));
        }
      }
    }
  }

  //----------------------------------------------------------------------------
  void AuxDetChannelMapSBNDAlg::Uninitialize() {}

  //----------------------------------------------------------------------------
  uint32_t AuxDetChannelMapSBNDAlg::PositionToAuxDetChannel(
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
      std::cout << "NAME " << gnItr->second << std::endl;
      if (csvItr == fADGeoToChannelAndSV.end()) {
        throw cet::exception("AuxDetChannelMapSBNDAlg")
        << "No entry in channel and sensitive volume"
        << " map for AuxDet index " << ad << " bail";
      }

      if (gnItr->second.find("Tagger") != std::string::npos) {
        channel = sv;
      }
    }

    if (channel == UINT_MAX) {
      throw cet::exception("AuxDetChannelMapSBNDAlg")
      << "position ("
      << worldLoc[0] << "," << worldLoc[1] << "," << worldLoc[2]
      << ") does not correspond to any AuxDet";
    }

    return channel;
  }

  //----------------------------------------------------------------------------
  const TVector3 AuxDetChannelMapSBNDAlg::AuxDetChannelToPosition(
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
      throw cet::exception("AuxDetChannelMapSBNDAlg")
      << "No AuxDetGeo with name " << auxDetName;
    }

    // Get the vector of channel and sensitive volume pairs
    auto csvItr = fADGeoToChannelAndSV.find(ad);
    if (csvItr == fADGeoToChannelAndSV.end()) {
      throw cet::exception("AuxDetChannelMapSBNDAlg")
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

        if (auxDetName.find("Tagger") != std::string::npos) {
          x = svOrigin[0];
          y = svOrigin[1];
          z = svOrigin[2];
        }

        break;
      }
    }

    return TVector3(x, y, z);
  }

}  // namespace geo


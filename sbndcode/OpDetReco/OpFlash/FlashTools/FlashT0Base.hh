///////////////////////////////////////////////////////////////////////
/// File: FlashT0Base.h
///
/// Interface class for a tool to calculate the recob::OpFlash t0
/// from the associated recob::OpHits
///
/// Created by Fran Nicolas, June 2022
////////////////////////////////////////////////////////////////////////

#ifndef SBND_FLASHT0BASE_H
#define SBND_FLASHT0BASE_H

#include "sbndcode/OpDetReco/OpFlash/FlashFinder/FlashFinderTypes.h"
#include "larcore/Geometry/WireReadout.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include <vector>
#include <string>
#include <cstddef>

namespace lightana
{
  class FlashT0Base {

  public:

    // Constructor
    FlashT0Base()
      : fWireReadout(
          art::ServiceHandle<geo::WireReadout>()->Get()
        )
      , fNOpChannels(fWireReadout.NOpChannels())
      , fChannelWeights(fNOpChannels, 0.0)
    {}

    // Default destructor
    virtual ~FlashT0Base() noexcept = default;


    // Method to calculate the OpFlash t0
    virtual double GetFlashT0(
        double flash_peaktime,
        LiteOpHitArray_t ophit_list
    ) = 0;


    // Method to calculate the weighted average magnitude
    double GetAverageMagnitude(
        std::vector<double> const& MagnitudeVector
    ) {

      // Check vector sizes
      if (MagnitudeVector.size() != fChannelWeights.size())
        throw art::Exception(art::errors::LogicError)
          << "Flash T0 Base error. "
          << "Size of magnitude vector and weights are not equal! with size " << MagnitudeVector.size() << " and " << fChannelWeights.size() << "!";

      double sum = 0.0;
      double weight_sum = 0.0;
      int n = 0;

      for (size_t i = 0; i < fChannelWeights.size(); i++) {

        sum += MagnitudeVector[i] * fChannelWeights[i];
        weight_sum += fChannelWeights[i];
        n++;
      }

      double average_magnitude = weight_sum ? sum / weight_sum : 0.0;

      return average_magnitude;
    }


    // Method to calculate the number of PMTs contributing to the flash
    int GetFlashNPMTs(std::string pdType)
    {
      int nChannels = 0;

      for (size_t i = 0; i < fChannelWeights.size(); i++) {

        if (fPDSMap.pdType(i) == pdType &&
            fChannelWeights[i] > 0) {

          nChannels++;
        }
      }

      return nChannels;
    }


  protected:

    geo::WireReadoutGeom const& fWireReadout =
        art::ServiceHandle<geo::WireReadout>()->Get();

    size_t fNOpChannels;

    std::vector<double> fChannelWeights =
        std::vector<double>(fWireReadout.NOpChannels(), 0.0);

    opdet::sbndPDMapAlg fPDSMap;

  };
}

#endif
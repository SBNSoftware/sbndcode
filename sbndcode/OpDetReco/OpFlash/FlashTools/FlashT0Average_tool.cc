///////////////////////////////////////////////////////////////////////
/// File: FlashT0SelectedChannels_tool.cc
///
/// Base class: FlashT0Base
///
/// Algorithm description: it averages the OpHit times
/// for the photon detectors using a PE-dependent weight
///
/// Weight:
///     w = PE / (PE + PE0)
///
/// Only the highest-PE OpHit from each channel is used.
/// Only OpHits above MinHitPE are considered.
///
/// Created by Fran Nicolas, June 2022
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/types/Atom.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Utilities/ToolConfigTable.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

#include "FlashT0Base.hh"

namespace lightana{

  class FlashT0Average : FlashT0Base
  {

  public:

    // Configuration parameters
    struct Config {

      fhicl::Atom<double> PreWindow {
        fhicl::Name("PreWindow"),
        fhicl::Comment(
          "Consider OpHits in the interval "
          "[FlashPeakTime-PreWindow, FlashPeakTime+PostWindow]"
        )
      };

      fhicl::Atom<double> PostWindow {
        fhicl::Name("PostWindow"),
        fhicl::Comment(
          "Consider OpHits in the interval "
          "[FlashPeakTime-PreWindow, FlashPeakTime+PostWindow]"
        )
      };

      fhicl::Atom<double> MinHitPE {
        fhicl::Name("MinHitPE"),
        fhicl::Comment(
          "Minimum number of reconstructed PE to consider the OpHit "
          "for the t0 calculation"
        )
      };

      fhicl::Atom<double> PE0 {
        fhicl::Name("PE0"),
        fhicl::Comment(
          "PE scale controlling the saturation of the PMT weight: "
          "weight = PE / (PE + PE0)"
        )
      };

      fhicl::Atom<std::string> PDType {
        fhicl::Name("PDType"),
        fhicl::Comment(
          "Type of PD to use: pmt_coated or pmt_uncoated"
        )
      };

    };

    // Constructor
    explicit FlashT0Average(art::ToolConfigTable<Config> const& config);

    // Method to calculate the OpFlash t0
    double GetFlashT0(
      double flash_peaktime,
      LiteOpHitArray_t ophit_list
    ) override;

  private:

    double fPreWindow;
    double fPostWindow;
    double fMinHitPE;
    double fPE0;
    std::string fPDType;

    opdet::sbndPDMapAlg fPDSMap;

  };


  FlashT0Average::FlashT0Average(
    art::ToolConfigTable<Config> const& config
  )
    : fPreWindow  { config().PreWindow()  },
      fPostWindow { config().PostWindow() },
      fMinHitPE   { config().MinHitPE()   },
      fPE0        { config().PE0()        },
      fPDType     { config().PDType()     }
  {
  }


  double FlashT0Average::GetFlashT0(
    double flash_time,
    LiteOpHitArray_t ophit_list
  ){
    // Store the highest-PE hit for each channel
    // channel -> (PE, peak_time)
    std::map<int, std::pair<double, double>> best_hit_per_channel;

    for(auto const& hit : ophit_list) {

      int channel = hit.channel;

      // Basic selection
      if( hit.peak_time < flash_time + fPostWindow &&
          hit.peak_time > flash_time - fPreWindow &&
          hit.pe > fMinHitPE &&
          fPDSMap.pdType(channel) == fPDType ) {

        // Keep only the hit with the largest PE for each channel
        auto it = best_hit_per_channel.find(channel);

        if(it == best_hit_per_channel.end() ||
           hit.pe > it->second.first) {

          best_hit_per_channel[channel] =
            std::make_pair(hit.pe, hit.peak_time);
        }
      }
    }


    // Calculate weighted average
    //
    // Weight:
    //     w = PE / (PE + PE0)
    //
    fChannelWeights.assign(fWireReadout.NOpChannels(), 0.0);
    double weighted_time_sum = 0.0;
    double weight_sum = 0.0;

    for(const auto& [channel, hit] : best_hit_per_channel) {


      double pe   = hit.first;
      double time = hit.second;

      double weight = pe / (pe + fPE0);

      weighted_time_sum += weight * time;
      weight_sum += weight;
      fChannelWeights[channel] = weight;
    }


    if(weight_sum > 0.0)
      return weighted_time_sum / weight_sum;
    else
      return flash_time;
  }

}

DEFINE_ART_CLASS_TOOL(lightana::FlashT0Average)
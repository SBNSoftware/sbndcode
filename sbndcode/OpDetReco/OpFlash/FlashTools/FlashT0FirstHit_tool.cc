///////////////////////////////////////////////////////////////////////
/// File: FlashT0FirstHit_tool.cc
///
/// Base class: FlashT0Base.hh
///
/// Algorithm description: The OpFlash t0 is set to the OpHit
/// with the largest signal (#PE). Only OpHits above the
/// MinHitPE cut are considered
///
/// Created by Fran Nicolas, June 2022
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/types/Atom.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Utilities/ToolConfigTable.h"

#include "FlashT0Base.hh"

namespace lightana{

  class FlashT0FirstHit : FlashT0Base
  {

  public:

    //Configuration parameters
    struct Config {

      fhicl::Atom<double> MinHitPE {
        fhicl::Name("MinHitPE"),
        fhicl::Comment("Minimum number of reconstructed PE to consider the OpHit for the t0 calculation")
      };

    };

    // Default constructor
    explicit FlashT0FirstHit(art::ToolConfigTable<Config> const& config);

    // Method to calculate the OpFlash t0
    double GetFlashT0(double flash_peaktime, LiteOpHitArray_t ophit_list) override;

  private:

    double fMinHitPE;

  };

  FlashT0FirstHit::FlashT0FirstHit(art::ToolConfigTable<Config> const& config)
    : fMinHitPE{ config().MinHitPE() }
  {
  }

  double FlashT0FirstHit::GetFlashT0(double flash_time, LiteOpHitArray_t ophit_list){

    double start_time = flash_time;

    for(auto & ophit : ophit_list){
      if(ophit.peak_time<start_time && ophit.pe>fMinHitPE) start_time = ophit.peak_time;
    }

    return start_time;
  }

}

DEFINE_ART_CLASS_TOOL(lightana::FlashT0FirstHit)

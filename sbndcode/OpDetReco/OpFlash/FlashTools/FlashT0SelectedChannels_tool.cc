///////////////////////////////////////////////////////////////////////
/// File: FlashT0SelectedChannels_tool.cc
///
/// Base class: FlashT0Base.hh
///
/// Algorithm description: it averages the OpHit times
/// for the photon detectors (PDs) containing a certain fraction (PDFraction)
/// of the integrated signal. Only OpHit in the interval
/// [FlashPeakTime-PreWindow, FlashPeakTime-PostWindow] and above
/// MinHitPE are taken into account
///
/// Created by Fran Nicolas, June 2022
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/types/Atom.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Utilities/ToolConfigTable.h"

#include "FlashT0Base.hh"

namespace lightana{

  class FlashT0SelectedChannels : FlashT0Base
  {

  public:

    //Configuration parameters
    struct Config {

      fhicl::Atom<double> PDFraction {
        fhicl::Name("PDFraction"),
        fhicl::Comment("Consider PDs containg a PDFraction of the signal.")
      };

      fhicl::Atom<double> PreWindow {
        fhicl::Name("PreWindow"),
        fhicl::Comment("Consider OpHits in the interval [FlashPeakTime-PreWindow, FlashPeakTime-PostWindow]")
      };

      fhicl::Atom<double> PostWindow {
        fhicl::Name("PostWindow"),
        fhicl::Comment("Consider OpHits in the interval [FlashPeakTime-PreWindow, FlashPeakTime-PostWindow]")
      };

      fhicl::Atom<double> MinHitPE {
        fhicl::Name("MinHitPE"),
        fhicl::Comment("Minimum number of reconstructed PE to consider the OpHit for the t0 calculation")
      };

    };

    // Default constructor
    explicit FlashT0SelectedChannels(art::ToolConfigTable<Config> const& config);

    // Method to calculate the OpFlash t0
    double GetFlashT0(double flash_peaktime, LiteOpHitArray_t ophit_list) override;

  private:

    double fPDFraction;
    double fPreWindow;
    double fPostWindow;
    double fMinHitPE;

  };

  FlashT0SelectedChannels::FlashT0SelectedChannels(art::ToolConfigTable<Config> const& config)
    : fPDFraction { config().PDFraction() },
    fPreWindow  { config().PreWindow()  },
    fPostWindow { config().PostWindow() },
    fMinHitPE   { config().MinHitPE()   }
  {
  }

  double FlashT0SelectedChannels::GetFlashT0(double flash_time, LiteOpHitArray_t ophit_list){

    std::vector< std::pair<double, double> > selected_hits;
    double pe_sum = 0;

    // fill vector with selected hits in the specified window
    for(auto const& hit : ophit_list) {
      if( hit.peak_time<flash_time+fPostWindow && hit.peak_time>flash_time-fPreWindow && hit.pe>fMinHitPE){
        selected_hits.push_back( std::make_pair(hit.pe, hit.peak_time));
        pe_sum += hit.pe ;
      }
    }

    if(pe_sum>0){
      // sort vector by number of #PE (ascending order)
      std::sort( selected_hits.begin(), selected_hits.end(), std::greater< std::pair<double, double> >() );

      double flasht0_mean=0, pe_count=0;
      int nophits=0;

      // loop over selected ophits
      for (size_t ix=0; ix<selected_hits.size(); ix++) {
        pe_count += selected_hits[ix].first;
        flasht0_mean += selected_hits[ix].second;
        nophits++;
        if( pe_count/pe_sum>fPDFraction ) break;
      }

      return flasht0_mean/nophits;
    }
    else
      return flash_time;
  }

}

DEFINE_ART_CLASS_TOOL(lightana::FlashT0SelectedChannels)

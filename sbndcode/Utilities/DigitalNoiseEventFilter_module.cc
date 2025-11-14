////////////////////////////////////////////////////////////////////////
// Class:       DigitalNoiseEventFilter
// Plugin Type: filter (Unknown Unknown)
// File:        DigitalNoiseEventFilter_module.cc
//
// Generated at Tue Oct  7 14:28:27 2025 by Thomas Junk using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "sbnobj/Common/Analysis/TPCChannelInfo.hh"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"

#include <memory>

class DigitalNoiseEventFilter;


class DigitalNoiseEventFilter : public art::EDFilter {
public:
  explicit DigitalNoiseEventFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DigitalNoiseEventFilter(DigitalNoiseEventFilter const&) = delete;
  DigitalNoiseEventFilter(DigitalNoiseEventFilter&&) = delete;
  DigitalNoiseEventFilter& operator=(DigitalNoiseEventFilter const&) = delete;
  DigitalNoiseEventFilter& operator=(DigitalNoiseEventFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:

  std::string   fTag;                          // input tag for reading TPCChannelInfo data
  float         fEvenFractionHighCut;          // upper cut on the fraction of ADC samples that are even
  float         fEvenFractionLowCut;           // lower cut on the fraction of ADC samples that are even
  int           fEvenFractionNumChannelCut;    // max # of channels with suspect even fractions to keep event. Otherwise suspect CE noise
  float         fRMSCut;                       // how noisy a channel has to be before we think it might have CE noise
  int           fRMSNumChannelCut;             // keep event if this many channels or fewer have RMS>fRMSCut.  Otherwise suspect CE noise
  bool          fKeepIfNoChannelInfo;          // If we do not find TPCChannelInfo, do we keep or reject the event?
  bool          fUseBadChanStatus;             // true if we skip over bad channels in calculating numbers of suspect CE noise channels

};


DigitalNoiseEventFilter::DigitalNoiseEventFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}  // ,
    // More initializers here.
{
  fTag                       = p.get<std::string>("TPCChannelInfoDataLabel","daq");
  fEvenFractionHighCut       = p.get<float>("EvenFractionHighCut",0.6);
  fEvenFractionLowCut        = p.get<float>("EvenFractionLowCut",0.4);
  fEvenFractionNumChannelCut = p.get<int>("EvenFractionNumChannelCut",10);
  fRMSCut                    = p.get<float>("RMSCut",1100.0);
  fRMSNumChannelCut         = p.get<int>("RMSNumChannelCut",50);
  fKeepIfNoChannelInfo       = p.get<bool>("KeepIfNoChannelInfo",true);
  fUseBadChanStatus          = p.get<bool>("UseBadChanStatus",true);
  
  consumes<std::vector<anab::TPCChannelInfo>>(fTag);
}

bool DigitalNoiseEventFilter::filter(art::Event& e)
{
  lariov::ChannelStatusProvider const& channelStatus(art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider());
  
  auto channelinfos = e.getHandle<std::vector<anab::TPCChannelInfo>>(fTag);
  if (!channelinfos || channelinfos->empty()) return fKeepIfNoChannelInfo;

  int neven_suspect = 0;
  int nrms_suspect = 0;
  double fevdiff = (fEvenFractionHighCut - fEvenFractionLowCut)/2.0;
  for (const auto &ci : *channelinfos) {
    if (fUseBadChanStatus && channelStatus.IsBad(ci.channel)) continue;
    if (ci.even_fraction > fEvenFractionHighCut || ci.even_fraction < fEvenFractionLowCut) ++neven_suspect;
    if (ci.rms > fRMSCut) {
      ++nrms_suspect;
    }
    else { // try a diagonal cut in case a waveform is noisy and then becomes stuck partway through
      if (fRMSCut > 0 && fevdiff > 0) {
	if ( (ci.rms/fRMSCut + std::abs(ci.even_fraction - 0.5)/fevdiff) > 1) ++nrms_suspect;
      }
    }
  }
  return (neven_suspect <= fEvenFractionNumChannelCut && nrms_suspect <= fRMSNumChannelCut);
}

DEFINE_ART_MODULE(DigitalNoiseEventFilter)

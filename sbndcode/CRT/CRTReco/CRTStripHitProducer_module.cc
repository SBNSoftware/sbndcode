////////////////////////////////////////////////////////////////////////
// Class:       CRTStripHitProducer
// Plugin Type: producer
// File:        CRTStripHitProducer_module.cc
//
// Author:      Henry Lay (h.lay@lancaster.ac.uk)
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTStripHit.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

#include <memory>

namespace sbnd::crt {
  class CRTStripHitProducer;
}


class sbnd::crt::CRTStripHitProducer : public art::EDProducer {
public:
  explicit CRTStripHitProducer(fhicl::ParameterSet const& p);

  CRTStripHitProducer(CRTStripHitProducer const&) = delete;
  CRTStripHitProducer(CRTStripHitProducer&&) = delete;
  CRTStripHitProducer& operator=(CRTStripHitProducer const&) = delete;
  CRTStripHitProducer& operator=(CRTStripHitProducer&&) = delete;

  void produce(art::Event& e) override;

  std::vector<CRTStripHit> CreateStripHits(art::Ptr<FEBData> &data);

private:

  CRTGeoAlg           fCRTGeoAlg;
  std::string         fFEBDataModuleLabel;
  uint16_t            fADCThreshold;
  std::vector<double> fErrorCoeff;
};


sbnd::crt::CRTStripHitProducer::CRTStripHitProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fCRTGeoAlg(p.get<fhicl::ParameterSet>("CRTGeoAlg", fhicl::ParameterSet()))
  , fFEBDataModuleLabel(p.get<std::string>("FEBDataModuleLabel"))
  , fADCThreshold(p.get<uint16_t>("ADCThreshold"))
  , fErrorCoeff(p.get<std::vector<double>>("ErrorCoeff"))
  {
    produces<std::vector<CRTStripHit>>();
    produces<art::Assns<FEBData, CRTStripHit>>();
  }

void sbnd::crt::CRTStripHitProducer::produce(art::Event& e)
{
  auto stripHitVec      = std::make_unique<std::vector<CRTStripHit>>();
  auto stripHitDataAssn = std::make_unique<art::Assns<FEBData, CRTStripHit>>();
  
  art::Handle<std::vector<FEBData>> FEBDataHandle;
  e.getByLabel(fFEBDataModuleLabel, FEBDataHandle);
  
  std::vector<art::Ptr<FEBData>> FEBDataVec;
  art::fill_ptr_vector(FEBDataVec, FEBDataHandle);

  for(auto data : FEBDataVec)
    {
      std::vector<CRTStripHit> newStripHits = CreateStripHits(data);
      
      for(auto hit : newStripHits)
	{
	  stripHitVec->push_back(hit);
	  util::CreateAssn(*this, e, *stripHitVec, data, *stripHitDataAssn);
	}
    }

  e.put(std::move(stripHitVec));
  e.put(std::move(stripHitDataAssn));
}

std::vector<sbnd::crt::CRTStripHit> sbnd::crt::CRTStripHitProducer::CreateStripHits(art::Ptr<FEBData> &data)
{
  std::vector<CRTStripHit> stripHits;

  const uint32_t mac5  = data->Mac5();
  const uint32_t unixs = data->UnixS();

  // Only consider "real data" readouts, not clock resets etc
  if(data->Flags() != 3)
    return stripHits;
  
  const CRTModuleGeo module = fCRTGeoAlg.GetModule(mac5 * 32);

  // Correct for FEB readout cable length
  // (time is FEB-by-FEB not channel-by-channel)
  const uint32_t t0 = data->Ts0() + module.t0CableDelayCorrection;
  const uint32_t t1 = data->Ts1() + module.t1CableDelayCorrection;

  // Iterate via strip (2 SiPMs per strip)
  const auto &sipm_adcs = data->ADC();
  for(unsigned adc_i = 0; adc_i < 32; adc_i+=2)
    {
      // Calculate SiPM channel number
      const uint16_t channel = mac5 * 32 + adc_i;

      const CRTStripGeo strip = fCRTGeoAlg.GetStrip(channel);
      const CRTSiPMGeo sipm1  = fCRTGeoAlg.GetSiPM(channel);
      const CRTSiPMGeo sipm2  = fCRTGeoAlg.GetSiPM(channel+1);

      // Subtract channel pedestals
      const uint16_t adc1 = sipm1.pedestal < sipm_adcs[adc_i]   ? sipm_adcs[adc_i] - sipm1.pedestal   : 0;
      const uint16_t adc2 = sipm2.pedestal < sipm_adcs[adc_i+1] ? sipm_adcs[adc_i+1] - sipm2.pedestal : 0;

      // Keep hit if both SiPMs above threshold
      if(adc1 > fADCThreshold && adc2 > fADCThreshold)
	{
	  // Access width of strip from the geometry algorithm
	  const double width = strip.width;

	  // Use light ratio to infer lateral position
	  const double pos = width / 2. * tanh(log(1. * adc2/adc1)) + width / 2.;
	  double err       = fErrorCoeff[0] * width + fErrorCoeff[1] * pos + fErrorCoeff[2] * pos * pos;

	  // Ensure error does not allow positions outside the strip geometry
	  if(pos + err > width)
	    err = width - pos;

	  if(pos - err < 0)
	    err = pos;

	  stripHits.emplace_back(channel, t0, t1, unixs, pos, err, adc1, adc2);
	}
    }

  return stripHits;
}

DEFINE_ART_MODULE(sbnd::crt::CRTStripHitProducer)

#include "CRTHitRecoAlg.h"

namespace sbnd{

  CRTHitRecoAlg::CRTHitRecoAlg(const Config& config)
  {
    this->reconfigure(config);
  }

  CRTHitRecoAlg::CRTHitRecoAlg() {}

  CRTHitRecoAlg::~CRTHitRecoAlg() {}

  void CRTHitRecoAlg::reconfigure(const Config& config) 
  {
    fADCThreshold = config.ADCThreshold();
    return;
  }

  std::vector<CRTStripHit> CRTHitRecoAlg::ProduceStripHits(std::vector<art::Ptr<sbnd::crt::FEBData>> &dataVec) 
  {
    std::vector<CRTStripHit> stripHits;

    // Iterate through each FEBData
    for(auto const &data : dataVec)
      {
	uint32_t mac5 = data->Mac5();

	CRTModuleGeo module = fCrtGeo.GetModule(mac5 * 32);

	// Correct for FEB readout cable length
	uint32_t t0 = data->Ts0() + module.cableDelayCorrection;
	uint32_t t1 = data->Ts1() + module.cableDelayCorrection;

	// Iterate via strip (2 SiPMs per strip)
	const auto &sipm_adcs = data->ADC();
	for(unsigned adc_i = 0; adc_i < 32; adc_i+=2)
	  {
	    // Add an offset for the SiPM channel number
	    uint16_t channel = mac5 * 32 + adc_i;

	    CRTStripGeo strip = fCrtGeo.GetStrip(channel);
	    CRTSiPMGeo sipm1  = fCrtGeo.GetSiPM(channel);
	    CRTSiPMGeo sipm2  = fCrtGeo.GetSiPM(channel+1);

	    // Subtract channel pedestals
	    uint16_t adc1 = sipm_adcs[adc_i]   - sipm1.pedestal;
	    uint16_t adc2 = sipm_adcs[adc_i+1] - sipm2.pedestal;

	    // Keep hit if both SiPMs above threshold
	    if(adc1 > fADCThreshold && adc2 > fADCThreshold)
	      {
		// Access width of strip from the geometry algorithm
		// === TO-DO === //
		// AMEND THE CODE AND IMPROVE CALCULATION
		double width = strip.width;
		double x     = width/2.;
		double ex    = width/2.;

		// Create hit
		stripHits.emplace_back(channel, t0, t1, x, ex, adc1, adc2);
	      }
	  }
      }
    return stripHits;
  }
}

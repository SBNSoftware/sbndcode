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
	// The FEB Mac5s in data don't correspond to the IDs used in the gdml
	uint32_t mac5 = /*fMac5ToGeometryIDMap[*/data->Mac5()/*]*/;
	
	// Iterate via strip (2 SiPMs per strip)
	const auto &sipm_adcs = data->ADC();
	for(unsigned adc_i = 0; adc_i < 32; adc_i+=2)
	  {
	    // Add an offset for the SiPM channel number
	    uint16_t channel = mac5 * 32 + adc_i;

	    // Subtract channel pedestals
	    uint16_t adc1 = sipm_adcs[adc_i]   /*- fStripPedestalMap[channel]*/;
	    uint16_t adc2 = sipm_adcs[adc_i+1] /*- fStripPedestalMap[channel+1]*/;

	    // Keep hit if both SiPMs above threshold
	    if(adc1 > fADCThreshold && adc2 > fADCThreshold)
	      {
		// Correct for FEB readout cable length
		uint32_t t0 = data->Ts0() /*+ fCableLengthCorrectionMap[mac5]*/;
		uint32_t t1 = data->Ts1() /*+ fCableLengthCorrectionMap[mac5]*/;
		
		// Access width of strip from the geometry algorithm
		// === TO-DO === //
		// AMEND THE CODE AND IMPROVE CALCULATION
		std::string stripName = fCrtGeo.ChannelToStripName(channel);
		if (stripName.empty()) {
		  throw cet::exception("CRTHitRecoAlg")
		    << "Cannot find strip name for channel " << channel << std::endl;
		}
		double width = fCrtGeo.GetStrip(stripName).width;
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

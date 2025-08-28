////////////////////////////////////////////////////////////////////////
// Class:       SBNDOpT0FinderAna
// Plugin Type: analyzer (art v3_05_01)
// File:        SBNDOpT0FinderAna_module.cc
//
// Generated at Mon Oct 12 13:52:37 2020 by Marco Del Tutto using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "sbnobj/SBND/OpFlashTiming/CorrectedOpFlashTiming.h"

#include <numeric>

class LightPropagationCorrectionAna;


class LightPropagationCorrectionAna : public art::EDAnalyzer {
public:
  explicit LightPropagationCorrectionAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LightPropagationCorrectionAna(LightPropagationCorrectionAna const&) = delete;
  LightPropagationCorrectionAna(LightPropagationCorrectionAna&&) = delete;
  LightPropagationCorrectionAna& operator=(LightPropagationCorrectionAna const&) = delete;
  LightPropagationCorrectionAna& operator=(LightPropagationCorrectionAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  std::string fLightPropCorrectionLabel; ///< The light propagation correction label producer

};


LightPropagationCorrectionAna::LightPropagationCorrectionAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  fLightPropCorrectionLabel = p.get<std::string>("LightPropCorrectionLabel");
}

void LightPropagationCorrectionAna::analyze(art::Event const& e)
{

  std::cout << "Run: " << e.id().run() << " Sub: " << e.id().subRun() << " Evt: " << e.id().event() << std::endl;

  // Get all the T0 objects
  ::art::Handle<std::vector<sbnd::OpFlashTiming::CorrectedOpFlashTiming>> correctedopflash_h;
  e.getByLabel(fLightPropCorrectionLabel, correctedopflash_h);
  if(!correctedopflash_h.isValid() || correctedopflash_h->empty()) {
    mf::LogWarning("fLightPropCorrectionLabel") << "No LightPropCorrectionLabel objects with label " << fLightPropCorrectionLabel << std::endl;
    return;
  }


  std::vector<art::Ptr<sbnd::OpFlashTiming::CorrectedOpFlashTiming>> correctedopflash_v;
  art::fill_ptr_vector(correctedopflash_v, correctedopflash_h);

  // Get the T0->Slice association
  art::FindManyP<recob::Slice> correctedflash_to_slices(correctedopflash_h, e, fLightPropCorrectionLabel);
  art::FindManyP<recob::OpFlash> correctedflash_to_flashes(correctedopflash_h, e, fLightPropCorrectionLabel);

  for (size_t i = 0; i < correctedopflash_v.size(); i++) {

    // The T0 object
    auto const correctedopflash =  correctedopflash_v[i];

    // The associations from T0 to both slices and flashes
    std::vector<art::Ptr<recob::Slice>> slice_v = correctedflash_to_slices.at(correctedopflash.key());
    std::vector<art::Ptr<recob::OpFlash>> flash_v = correctedflash_to_flashes.at(correctedopflash.key());

    // One T0 object is always associated to one and only one Slice and Flash
    assert(slice_v.size() == 1);
    assert(flash_v.size() == 1);
    
    std::cout << " Corrected flash time is " << correctedopflash->OpFlashT0 << std::endl;
    std::cout << " Associated with slice id " << slice_v[0]->ID() << std::endl;
    std::cout << "Corrected flash time light only  " << correctedopflash->UpstreamTime_lightonly << std::endl;
    std::cout << "Corrected flash time tpc z corr " << correctedopflash->UpstreamTime_tpczcorr << std::endl;
    std::cout << "Corrected flash time prop corr tpc z corr " << correctedopflash->UpstreamTime_propcorr_tpczcorr << std::endl;
  }
}

DEFINE_ART_MODULE(LightPropagationCorrectionAna)
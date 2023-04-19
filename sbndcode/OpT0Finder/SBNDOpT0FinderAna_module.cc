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
#include "lardataobj/AnalysisBase/T0.h"
#include "sbnobj/Common/Reco/OpT0FinderResult.h"

#include <numeric>

class SBNDOpT0FinderAna;


class SBNDOpT0FinderAna : public art::EDAnalyzer {
public:
  explicit SBNDOpT0FinderAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SBNDOpT0FinderAna(SBNDOpT0FinderAna const&) = delete;
  SBNDOpT0FinderAna(SBNDOpT0FinderAna&&) = delete;
  SBNDOpT0FinderAna& operator=(SBNDOpT0FinderAna const&) = delete;
  SBNDOpT0FinderAna& operator=(SBNDOpT0FinderAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  std::string _t0_producer; ///< The T0 producer (to be set)

};


SBNDOpT0FinderAna::SBNDOpT0FinderAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  _t0_producer = p.get<std::string>("T0Producer");
}

void SBNDOpT0FinderAna::analyze(art::Event const& e)
{

  // Get all the T0 objects
  ::art::Handle<std::vector<sbn::OpT0Finder>> opt0_h;
  e.getByLabel(_t0_producer, opt0_h);
  if(!opt0_h.isValid() || opt0_h->empty()) {
    mf::LogWarning("SBNDOpT0FinderAna") << "Don't have good T0s." << std::endl;
    return;
  }
  std::vector<art::Ptr<sbn::OpT0Finder>> opt0_v;
  art::fill_ptr_vector(opt0_v, opt0_h);

  // Get the T0->Slice association
  art::FindManyP<recob::Slice> opt0_to_slices(opt0_h, e, _t0_producer);
  art::FindManyP<recob::OpFlash> opt0_to_flashes(opt0_h, e, _t0_producer);

  for (size_t n_t0 = 0; n_t0 < opt0_v.size(); n_t0++) {

    // The T0 object
    auto const opt0 =  opt0_v[n_t0];

    // The associations from T0 to both slices and flashes
    std::vector<art::Ptr<recob::Slice>> slice_v = opt0_to_slices.at(opt0.key());
    std::vector<art::Ptr<recob::OpFlash>> flash_v = opt0_to_flashes.at(opt0.key());

    // One T0 object is always associated to one and only one Slice and Flash
    assert(slice_v.size() == 1);
    assert(flash_v.size() == 1);

    int nopdet = std::accumulate((opt0->opch).begin(),  (opt0->opch).end(), 0);

    std::cout << "T0 obj: " << n_t0 << ", score: " << opt0->score << std::endl;
    std::cout << "\t is associated with slice ID " << slice_v[0]->ID() << std::endl;
    std::cout << "\t is associated with flash with time " << opt0->time 
                                         << ", total meas PE " << opt0->measPE  
                                         << ", total hypo PE " << opt0->hypoPE
                                         << ", using " <<  nopdet << " optical detectors for scoring."
                                         << std::endl;
  }

}

DEFINE_ART_MODULE(SBNDOpT0FinderAna)

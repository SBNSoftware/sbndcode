////////////////////////////////////////////////////////////////////////
// Class:       FullOscNuFilter
// Module Type: filter
// File:        FullOscNuFilter_module.cc
//
// Copy of GenNuFilter for the full-osc workflow. Reads a single labelled
// MCTruth collection and can restrict the selection to the first
// (oscillated) neutrino only, via OscOnly.
////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <string>
#include <vector>

#include "TVector3.h"

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "sbndcode/RecoUtils/RecoUtils.h"

namespace filt {

  class FullOscNuFilter : public art::EDFilter {
  public:
    explicit FullOscNuFilter(fhicl::ParameterSet const& pset);
    bool filter(art::Event& e) override;

  private:
    bool fVtxInTPC;
    std::vector<int> fLepPDGs;
    bool fCC;
    bool fNC;
    bool fOscOnly;
    std::string fGenieModuleLabel;
  };


  FullOscNuFilter::FullOscNuFilter(fhicl::ParameterSet const& pset)
    : EDFilter(pset)
    , fVtxInTPC(pset.get<bool>("VtxInTPC"))
    , fLepPDGs(pset.get<std::vector<int>>("LepPDGs"))
    , fCC(pset.get<bool>("CC"))
    , fNC(pset.get<bool>("NC"))
    , fOscOnly(pset.get<bool>("OscOnly", false))
    , fGenieModuleLabel(pset.get<std::string>("GenieModuleLabel", "generator"))
  {}


  bool FullOscNuFilter::filter(art::Event& e)
  {
    art::Handle<std::vector<simb::MCTruth>> mclist;
    e.getByLabel(fGenieModuleLabel, mclist);

    // Safety: missing or empty collection -> reject event instead of crashing
    if (!mclist.isValid()) {
      mf::LogWarning("FullOscNuFilter")
        << "No std::vector<simb::MCTruth> with label '" << fGenieModuleLabel
        << "' in event " << e.id() << "; rejecting event.";
      return false;
    }
    if (mclist->empty()) return false;

    const size_t end = fOscOnly ? 1 : mclist->size();

    for (size_t i = 0; i < end; ++i) {
      const art::Ptr<simb::MCTruth> mc_truth(mclist, i);
      if (mc_truth->Origin() != simb::kBeamNeutrino) continue;

      const simb::MCNeutrino& nu = mc_truth->GetNeutrino();

      // Check vertex in TPC
      TVector3 vtx(nu.Nu().Vx(), nu.Nu().Vy(), nu.Nu().Vz());
      if (fVtxInTPC && !RecoUtils::IsInsideTPC(vtx, 0)) continue;

      // Check if CC or NC specified
      if (!(fCC && fNC)) {
        if (fCC && nu.CCNC() == simb::kNC) continue;
        if (fNC && nu.CCNC() == simb::kCC) continue;
      }

      // Check which flavour of interaction (outgoing lepton PDG)
      const int lepPdg = nu.Lepton().PdgCode();
      if (!fLepPDGs.empty() && fLepPDGs[0] != 0) {
        if (std::find(fLepPDGs.begin(), fLepPDGs.end(), lepPdg) == fLepPDGs.end()) continue;
      }

      return true;
    }
    return false;
  }


  DEFINE_ART_MODULE(FullOscNuFilter)

}

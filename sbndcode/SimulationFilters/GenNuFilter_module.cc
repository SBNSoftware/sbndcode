#include <iostream>
#include <algorithm>

#include "TGeoManager.h"

#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "sbndcode/RecoUtils/RecoUtils.h"

namespace filt{

  class GenNuFilter : public art::EDFilter {
    public:
      explicit GenNuFilter(fhicl::ParameterSet const & pset);
      virtual bool filter(art::Event& e) override;
      void reconfigure(fhicl::ParameterSet const& pset);
      virtual void beginJob() override;

    private:

      bool fVtxInTPC;
      std::vector<int> fLepPDGs;
      bool fCC;
      bool fNC;

  };


  GenNuFilter::GenNuFilter(fhicl::ParameterSet const & pset)
  : EDFilter(pset)
  {
    this->reconfigure(pset);
    
  }


  void GenNuFilter::reconfigure(fhicl::ParameterSet const& pset){
    fVtxInTPC = pset.get<bool>("VtxInTPC");
    fLepPDGs = pset.get<std::vector<int>>("LepPDGs");
    fCC = pset.get<bool>("CC");
    fNC = pset.get<bool>("NC");
  }


  bool GenNuFilter::filter(art::Event & e){
    std::vector< art::Handle< std::vector<simb::MCTruth> > > mclists;
    e.getManyByType(mclists);
    for (unsigned int i = 0; i < mclists.size() ; i++){
      for (unsigned int j = 0; j < mclists[i]->size(); j++){
        //Should have the truth record for the event now
        const art::Ptr<simb::MCTruth> mc_truth(mclists[i],j);
        if(mc_truth->Origin() == simb::kBeamNeutrino){
          simb::MCNeutrino nu = mc_truth->GetNeutrino();
          // Check vertex in TPC
          double vtxX = nu.Nu().Vx();
          double vtxY = nu.Nu().Vy();
          double vtxZ = nu.Nu().Vz();
          TVector3 vtx(vtxX, vtxY, vtxZ);
          if(fVtxInTPC && !RecoUtils::IsInsideTPC(vtx, 0)) continue;
          // Check if CC or NC specified
          if(!(fCC && fNC)){
            if(fCC && nu.CCNC() == simb::kNC) continue;
            if(fNC && nu.CCNC() == simb::kCC) continue;
          }
          // Check which flavour of interaction
          int lepPdg = nu.Lepton().PdgCode();
          if(fLepPDGs.size() > 0){
            if(fLepPDGs[0] != 0){
              if(std::find(fLepPDGs.begin(), fLepPDGs.end(), lepPdg) == fLepPDGs.end()) continue;
            }
          }
          return true;
        }
      }
    }

    return false;
  }


  void GenNuFilter::beginJob() {
  }


  DEFINE_ART_MODULE(GenNuFilter)

}

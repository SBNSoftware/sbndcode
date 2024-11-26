////////////////////////////////////////////////////////////////////////
// Class:       SBNDCVNTest
// Plugin Type: analyzer (Unknown Unknown)
// File:        SBNDCVNTest_module.cc
//
// Generated at Fri Nov 22 10:10:07 2024 by Tingjun Yang using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larrecodnn/CVN/func/Result.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "TTree.h"

using namespace std;

namespace lcvn {
  class SBNDCVNTest;
}

class lcvn::SBNDCVNTest : public art::EDAnalyzer {
public:
  explicit SBNDCVNTest(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SBNDCVNTest(SBNDCVNTest const&) = delete;
  SBNDCVNTest(SBNDCVNTest&&) = delete;
  SBNDCVNTest& operator=(SBNDCVNTest const&) = delete;
  SBNDCVNTest& operator=(SBNDCVNTest&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Declare member data here.
  art::InputTag fSliceLabel;
  art::InputTag fCVNLabel;
  art::InputTag fHitLabel;
  art::InputTag fHitMatchLabel;

  TTree *anatree;
  int run;
  int event;
  vector<int> nhits;
  vector<float> numuscore;
  vector<float> nuescore;
  vector<float> cosmicscore;
  vector<float> ncscore;
};


lcvn::SBNDCVNTest::SBNDCVNTest(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fSliceLabel(p.get<art::InputTag>("SliceLabel"))
  , fCVNLabel(p.get<art::InputTag>("CVNLabel"))
  , fHitLabel(p.get<art::InputTag>("HitLabel"))
  , fHitMatchLabel(p.get<art::InputTag>("HitMatchLabel"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void lcvn::SBNDCVNTest::analyze(art::Event const& e)
{
  run = e.run();
  event = e.id().event();
  nhits.clear();
  numuscore.clear();
  nuescore.clear();
  cosmicscore.clear();
  ncscore.clear();
  
  auto slcHandle = e.getHandle< std::vector<recob::Slice> >(fSliceLabel);
  if (!slcHandle){
    cout<<"slcHandle invalid"<<endl;
    return;
  }
  //auto mcHandle = e.getHandle< std::vector<simb::MCTruth>>(fGenieLabel);
  auto hitHandle = e.getHandle< std::vector<recob::Hit>>(fHitLabel);
  
  art::FindManyP<recob::Hit> findManyHits(slcHandle, e, fSliceLabel);
  art::FindOne<lcvn::Result> findOneCVN(slcHandle, e, fCVNLabel);
  art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> fmhitmc(hitHandle, e, fHitMatchLabel);
  // ParticleInventoryService
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
    
  for (size_t i = 0; i<slcHandle->size(); ++i){
    if (findOneCVN.at(i).isValid()){
      auto const & cvn = findOneCVN.at(i).ref().fOutput;
//      for (size_t j = 0; j<cvn.size(); ++j){
//        for (size_t k = 0; k<cvn[j].size(); ++k){
//          std::cout<<cvn[j][k]<<" ";
//          if (k==cvn[j].size()-1) std::cout<<std::endl;
//        }
//      }
      numuscore.push_back(cvn[0][0]);
      nuescore.push_back(cvn[0][1]);
      cosmicscore.push_back(cvn[0][2]);
      ncscore.push_back(cvn[0][3]);
    }
    if(findManyHits.isValid()){
      auto const & slice_hits = findManyHits.at(i);
      nhits.push_back(slice_hits.size());
    }
    /*
    if(findManyHits.isValid()){
      auto const & slice_hits = findManyHits.at(i);
      cout<<"slice "<<i<<" nhits "<<slice_hits.size()<<endl;
      double avgtpc[3] = {0};
      double avgwire[3] = {0};
      double avgtick[3] = {0};
      int nhits[3] ={0};
      for (auto const & hit : slice_hits){
        auto const & wid = hit->WireID();
        avgtpc[wid.Plane] += wid.TPC;
        avgwire[wid.Plane] += wid.Wire;
        avgtick[wid.Plane] += hit->PeakTime();
        ++nhits[wid.Plane];
      }
      for (int ipl = 0; ipl <3; ++ipl){
        std::cout<<"plane="<<ipl<<" tpc="<<avgtpc[ipl]/nhits[ipl]<<" wire="<<avgwire[ipl]/nhits[ipl]<<" tick="<<avgtick[ipl]/nhits[ipl]<<std::endl;
      }
    }
    */     
    /*
    if (!e.isRealData()){
      if(findManyHits.isValid()){
        auto const & slice_hits = findManyHits.at(i);
        double tot_slice_eng = 0;
        double tot_slice_nu_eng = 0;
        double tot_slice_cos_eng = 0;
        
//        int fNhits_tpc_0_pl_0 = Get_Hit_Count(0,0,slice_hits);
//        int fNhits_tpc_0_pl_1 = Get_Hit_Count(0,1,slice_hits);
//        int fNhits_tpc_0_pl_2 = Get_Hit_Count(0,2,slice_hits);
//        int fNhits_tpc_1_pl_0 = Get_Hit_Count(1,0,slice_hits);
//        int fNhits_tpc_1_pl_1 = Get_Hit_Count(1,1,slice_hits);
//        int fNhits_tpc_1_pl_2 = Get_Hit_Count(1,2,slice_hits);
//        int fNhits_total = slice_hits.size();
        
        //std::vector<double> mc_truth_eng(mcHandle->size(),0.);
        std::cout<<"slice "<<i<<" total hits "<<slice_hits.size()<<std::endl;
        for(auto const & hit : slice_hits){
          auto particles = fmhitmc.at(hit.key());
          auto hitmatch = fmhitmc.data(hit.key());
          for(size_t e = 0; e<particles.size(); ++e){
            if (!particles[e]) continue;
            if (!particles[e]->TrackId()) continue;
            //if (!pi_serv->TrackIdToMotherParticle_P(particles[e]->TrackId())) continue;
            //size_t trkid = (pi_serv->TrackIdToMotherParticle_P(particles[e]->TrackId()))->TrackId();
            tot_slice_eng += hitmatch[e]->energy;
            auto & mctruth = pi_serv->TrackIdToMCTruth_P(particles[e]->TrackId());
            if (mctruth){
              std::cout<<"origin= "<<mctruth->Origin()<<" pdg= "<<particles[e]->PdgCode()<<std::endl;
              if (mctruth->Origin() == simb::kBeamNeutrino){
                tot_slice_nu_eng += hitmatch[e]->energy;
              }
              else if (mctruth->Origin() == simb::kCosmicRay){
                tot_slice_cos_eng += hitmatch[e]->energy;
              }
              //mc_truth_eng[mctruth.key()] += hitmatch[e]->energy;
            }
          }
        } // loop over hits in the selected slice
	
        //ftotsliceE = tot_slice_eng;
    	
        if(tot_slice_eng > 0){
          std::cout << "Total energy : " << tot_slice_eng << "\n";
          std::cout << "Total cosmic energy : " << tot_slice_cos_eng << "\n";
          std::cout << "Total nu energy : " << tot_slice_nu_eng << "\n";
          std::cout << "Cosmic fraction : " << double(tot_slice_cos_eng)/tot_slice_eng << "\n";
          std::cout << "Neutrino fraction : " << double(tot_slice_nu_eng)/tot_slice_eng << "\n";
          
//          ftotsliceNuE = tot_slice_nu_eng;
//          ftotsliceCosE = tot_slice_cos_eng;
//          ftotsliceOthE = ftotsliceE - ftotsliceNuE - ftotsliceCosE;
//          fsliceNuEfrac = double(ftotsliceNuE)/ftotsliceE;
//          fsliceCosEfrac = double(ftotsliceCosE)/ftotsliceE;
//          fsliceOthEfrac = double(ftotsliceOthE)/ftotsliceE;
        }
      }
    }
    */
  }
  anatree->Fill();
}

void lcvn::SBNDCVNTest::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  anatree = tfs->make<TTree>("anatree", "anatree");
  anatree->Branch("run", &run);
  anatree->Branch("event", &event);
  anatree->Branch("nhits", &nhits);
  anatree->Branch("numuscore", &numuscore);
  anatree->Branch("nuescore", &nuescore);
  anatree->Branch("cosmicscore", &cosmicscore);
  anatree->Branch("ncscore", &ncscore);
}

DEFINE_ART_MODULE(lcvn::SBNDCVNTest)

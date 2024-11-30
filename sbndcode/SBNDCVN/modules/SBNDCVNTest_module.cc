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
  vector<int> nhits0;
  vector<int> nhits1;
  vector<int> nhits2;
  vector<float> numuscore;
  vector<float> nuescore;
  vector<float> cosmicscore;
  vector<float> ncscore;
  vector<float> edep;
  vector<float> nufrac;
  vector<float> enu;
  vector<float> y;
  vector<int> nupdg;
  vector<int> ccnc;
  vector<int> type;
  vector<double> nuvtxx;
  vector<double> nuvtxy;
  vector<double> nuvtxz;
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
  nhits0.clear();
  nhits1.clear();
  nhits2.clear();
  numuscore.clear();
  nuescore.clear();
  cosmicscore.clear();
  ncscore.clear();
  edep.clear();
  nufrac.clear();
  enu.clear();
  y.clear();
  nupdg.clear();
  ccnc.clear();
  type.clear();
  nuvtxx.clear();
  nuvtxy.clear();
  nuvtxz.clear();

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
      int this_nhits0 = 0;
      int this_nhits1 = 0;
      int this_nhits2 = 0;
      for (auto const & hit: slice_hits){
        int plane = hit->WireID().Plane;
        if (plane == 0){
          ++this_nhits0;
        }
        else if (plane == 1){
          ++this_nhits1;
        }
        else if (plane == 2){
          ++this_nhits2;
        }
      }
      nhits0.push_back(this_nhits0);
      nhits1.push_back(this_nhits1);
      nhits2.push_back(this_nhits2);
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

    if (!e.isRealData()){
      double tot_slice_eng = 0;
      double tot_slice_nu_eng = 0;
      double tot_slice_cos_eng = 0;
      double this_nufrac = 0;
      double this_enu = 0;
      double this_y = 0;
      int    this_nupdg = 0;
      int    this_ccnc = 0;
      int    this_type = 0;
      double this_nuvtxx = -1000;
      double this_nuvtxy = -1000;
      double this_nuvtxz = -1000;
      if(findManyHits.isValid()){
        auto const & slice_hits = findManyHits.at(i);
        
//        int fNhits_tpc_0_pl_0 = Get_Hit_Count(0,0,slice_hits);
//        int fNhits_tpc_0_pl_1 = Get_Hit_Count(0,1,slice_hits);
//        int fNhits_tpc_0_pl_2 = Get_Hit_Count(0,2,slice_hits);
//        int fNhits_tpc_1_pl_0 = Get_Hit_Count(1,0,slice_hits);
//        int fNhits_tpc_1_pl_1 = Get_Hit_Count(1,1,slice_hits);
//        int fNhits_tpc_1_pl_2 = Get_Hit_Count(1,2,slice_hits);
//        int fNhits_total = slice_hits.size();
        
        //std::vector<double> mc_truth_eng(mcHandle->size(),0.);
        //std::cout<<"slice "<<i<<" total hits "<<slice_hits.size()<<std::endl;
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
              //std::cout<<"origin= "<<mctruth->Origin()<<" pdg= "<<particles[e]->PdgCode()<<std::endl;
              if (mctruth->Origin() == simb::kBeamNeutrino){
                tot_slice_nu_eng += hitmatch[e]->energy;
                this_enu = mctruth->GetNeutrino().Nu().E();
                this_y = mctruth->GetNeutrino().Y();
                this_nupdg = mctruth->GetNeutrino().Nu().PdgCode();
                this_ccnc = mctruth->GetNeutrino().CCNC();
                this_type = mctruth->GetNeutrino().InteractionType();
                this_nuvtxx = mctruth->GetNeutrino().Nu().Vx();
                this_nuvtxy = mctruth->GetNeutrino().Nu().Vy();
                this_nuvtxz = mctruth->GetNeutrino().Nu().Vz();
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
//          std::cout << "Total energy : " << tot_slice_eng << "\n";
//          std::cout << "Total cosmic energy : " << tot_slice_cos_eng << "\n";
//          std::cout << "Total nu energy : " << tot_slice_nu_eng << "\n";
//          std::cout << "Cosmic fraction : " << double(tot_slice_cos_eng)/tot_slice_eng << "\n";
//          std::cout << "Neutrino fraction : " << double(tot_slice_nu_eng)/tot_slice_eng << "\n";
          this_nufrac = tot_slice_nu_eng/tot_slice_eng;
//          ftotsliceNuE = tot_slice_nu_eng;
//          ftotsliceCosE = tot_slice_cos_eng;
//          ftotsliceOthE = ftotsliceE - ftotsliceNuE - ftotsliceCosE;
//          fsliceNuEfrac = double(ftotsliceNuE)/ftotsliceE;
//          fsliceCosEfrac = double(ftotsliceCosE)/ftotsliceE;
//          fsliceOthEfrac = double(ftotsliceOthE)/ftotsliceE;
        }
      }
      edep.push_back(tot_slice_eng);
      nufrac.push_back(this_nufrac);
      enu.push_back(this_enu);
      y.push_back(this_y);
      nupdg.push_back(this_nupdg);
      ccnc.push_back(this_ccnc);
      type.push_back(this_type);
      nuvtxx.push_back(this_nuvtxx);
      nuvtxy.push_back(this_nuvtxy);
      nuvtxz.push_back(this_nuvtxz);
    }
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
  anatree->Branch("nhits0", &nhits0);
  anatree->Branch("nhits1", &nhits1);
  anatree->Branch("nhits2", &nhits2);
  anatree->Branch("numuscore", &numuscore);
  anatree->Branch("nuescore", &nuescore);
  anatree->Branch("cosmicscore", &cosmicscore);
  anatree->Branch("ncscore", &ncscore);
  anatree->Branch("edep", &ncscore);
  anatree->Branch("nufrac", &nufrac);
  anatree->Branch("enu", &enu);
  anatree->Branch("y", &y);
  anatree->Branch("nupdg", &nupdg);
  anatree->Branch("ccnc", &ccnc);
  anatree->Branch("type", &type);
  anatree->Branch("nuvtxx", &nuvtxx);
  anatree->Branch("nuvtxy", &nuvtxy);
  anatree->Branch("nuvtxz", &nuvtxz);

}

DEFINE_ART_MODULE(lcvn::SBNDCVNTest)

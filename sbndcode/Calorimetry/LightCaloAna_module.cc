////////////////////////////////////////////////////////////////////////
// Class:       LightCaloAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        LightCaloAna_module.cc
//
// Generated at Fri Oct 14 13:43:49 2022 by Lynn Tung using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Additional framework includes
#include "art_root_io/TFileService.h"

#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

// LArSoft includes 
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

// SBND includes
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"

#include "TFile.h"
#include "TTree.h"


namespace sbnd {
  class LightCaloAna;
}


class sbnd::LightCaloAna : public art::EDAnalyzer {
public:
  explicit LightCaloAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  LightCaloAna(LightCaloAna const&) = delete;
  LightCaloAna(LightCaloAna&&) = delete;
  LightCaloAna& operator=(LightCaloAna const&) = delete;
  LightCaloAna& operator=(LightCaloAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // define functions 
  // art::Ptr<recob::Slice> SelectSlice(std::vector<art::Ptr<recob::Slice>> slice_v, std::vector<float> nuscore_v, std::vector<float> fmscore_v
  //                                    std::vector<art::Ptr<recob::OpFlash>> opflash_v);

  std::vector<art::Ptr<recob::OpFlash>> MatchOpFlash(std::vector<art::Ptr<sbn::SimpleFlashMatch>> fm_v,
                                                     std::vector<art::Ptr<recob::OpFlash>> flash_v);
  std::vector<std::string> _opflash_producer_v; ///< The OpFlash producers (to be set)
  std::string _slice_producer; ///< The Slice producer (to be set)
  std::string _flashmatch_producer;
  float _nuscore_cut; 
  float _fmscore_cut;

  TTree* _tree;
  int _run, _subrun, _event;
  std::vector<int> _slc_pfpid;
  std::vector<float> _slc_nuscore;
  std::vector<float> _slc_fmscore; 
  std::vector<float> _slc_fmtime;
};


sbnd::LightCaloAna::LightCaloAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  _opflash_producer_v = p.get<std::vector<std::string>>("OpFlashProducers");
  _slice_producer = p.get<std::string>("SliceProducer");
  _flashmatch_producer = p.get<std::string>("FlashMatchProducer");
  _nuscore_cut = p.get<float>("nuScoreCut");
  _fmscore_cut = p.get<float>("fmScoreCut");

  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("slice_tree","");
  _tree->Branch("run",             &_run,                             "run/I");
  _tree->Branch("subrun",          &_subrun,                          "subrun/I");
  _tree->Branch("event",           &_event,                           "event/I");
  _tree->Branch("slc_pfpid", "std::vector<int>", &_slc_pfpid);
  _tree->Branch("slc_nuscore", "std::vector<float>", &_slc_nuscore);
  _tree->Branch("slc_fmscore", "std::vector<float>", &_slc_fmscore);
}

void sbnd::LightCaloAna::analyze(art::Event const& e)
{
  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();
  std::cout << "run: " << _run <<  ", subrun: " << _subrun  << ", event: " <<  _event << std::endl;

  // get slices 
  ::art::Handle<std::vector<recob::Slice>> slice_h;
  e.getByLabel(_slice_producer, slice_h);
  if(!slice_h.isValid() || slice_h->empty()){
    std::cout << "dont have good slices!" << std::endl;
    return;
  }

  ::art::Handle<std::vector<recob::PFParticle>> pfp_h;
  e.getByLabel(_slice_producer, pfp_h);
  if(!pfp_h.isValid() || pfp_h->empty()) {
    std::cout << "don't have good PFParticle!" << std::endl;
    return;
  }

  auto const & flash0_h = e.getValidHandle<std::vector<recob::OpFlash>>(_opflash_producer_v[0]);
  auto const & flash1_h = e.getValidHandle<std::vector<recob::OpFlash>>(_opflash_producer_v[1]);
  if( (!flash0_h.isValid() || flash0_h->empty()) && (!flash1_h.isValid() || flash1_h->empty())) {
    std::cout << "don't have good flashes from producer " << _opflash_producer_v[0] << " or "  << _opflash_producer_v[1] << std::endl;
    return;
  }

  // Construct the vector of Slices
  std::vector<art::Ptr<recob::Slice>> slice_v;
  art::fill_ptr_vector(slice_v, slice_h);

  // Get associations 
  art::FindManyP<recob::PFParticle> slice_to_pfps (slice_h, e, _slice_producer);
  art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_to_meta(pfp_h, e, _slice_producer);
  art::FindManyP<sbn::SimpleFlashMatch> pfps_to_sfm (pfp_h, e, _flashmatch_producer);

  _slc_pfpid.clear();
  _slc_nuscore.clear();
  _slc_fmscore.clear(); 

  std::vector<art::Ptr<recob::Slice>> match_slices_v; 
  std::vector<art::Ptr<sbn::SimpleFlashMatch>> match_fm_v;

  for (size_t n_slice=0; n_slice < slice_v.size(); n_slice++){
    float nu_score = -9999;
    float fm_score = -9999;
    float fm_time  = -9999;
    auto slice = slice_v[n_slice];
    bool found_fm = false;
    std::vector<art::Ptr<recob::PFParticle>> pfp_v = slice_to_pfps.at(n_slice);
    for (size_t n_pfp=0; n_pfp < pfp_v.size(); n_pfp++){
      auto pfp = pfp_v[n_pfp];

      // only select the PRIMARY pfp 
      if(!pfp->IsPrimary() && !(abs(pfp->PdgCode()) == 12 || abs(pfp->PdgCode()) == 14|| abs(pfp->PdgCode()) == 16))
        continue;
      _slc_pfpid.push_back(pfp->Self());

      // if primary, get nu-score 
      const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpmeta_v = pfp_to_meta.at(pfp->Self());
      const art::Ptr<larpandoraobj::PFParticleMetadata> pfpmeta = pfpmeta_v.front();
      larpandoraobj::PFParticleMetadata::PropertiesMap propmap = pfpmeta->GetPropertiesMap();
      if (propmap.count("NuScore")) nu_score = propmap.at("NuScore");
      else nu_score = -1;
      _slc_nuscore.push_back(nu_score);
      
      // get fm-score 
      std::vector<art::Ptr<sbn::SimpleFlashMatch>> fm_v = pfps_to_sfm.at(pfp.key());
      if (fm_v.empty()){
        _slc_fmscore.push_back(-999);
        continue;
      }
      if (fm_v.size() > 1)
        std::cout << "more than one match for one pfp?" << std::endl;
      for (size_t n_fm=0; n_fm < fm_v.size(); n_fm++){
        auto fm = fm_v.at(n_fm);
        fm_score = fm->score.total;
        fm_time  = fm->time;
        if (nu_score > _nuscore_cut && fm_score < _fmscore_cut && fm_score > 0){
          found_fm = true;
          match_fm_v.push_back(fm);
        }
      } // end flashmatch loop
      if (found_fm ==true) match_slices_v.push_back(slice);
      _slc_fmscore.push_back(fm_score);
      _slc_fmtime.push_back(fm_time);
    } // end pfp loop
  } // end slice loop
  if (match_slices_v.size() != match_fm_v.size()){
    std::cout << "slice and flashmatch vector length mismatch!" << std::endl;
    return;
  }

  // tpc0
  std::vector<art::Ptr<recob::OpFlash>> flash0_v;
  art::fill_ptr_vector(flash0_v, flash0_h);
  std::vector<art::Ptr<recob::OpFlash>> match_op0 = MatchOpFlash(match_fm_v,flash0_v);

  // tpc1 
  std::vector<art::Ptr<recob::OpFlash>> flash1_v;
  art::fill_ptr_vector(flash1_v, flash1_h);
  std::vector<art::Ptr<recob::OpFlash>> match_op1 = MatchOpFlash(match_fm_v,flash1_v);

  _tree->Fill();
} // end analyze


// define functions 

std::vector<art::Ptr<recob::OpFlash>> sbnd::LightCaloAna::MatchOpFlash(std::vector<art::Ptr<sbn::SimpleFlashMatch>> fm_v,
                                                                       std::vector<art::Ptr<recob::OpFlash>> flash_v){
  std::vector<art::Ptr<recob::OpFlash>> matched_opflashes;
  for (size_t ifm=0; ifm<fm_v.size();ifm++){
    auto fm = fm_v[ifm];
    art::Ptr<recob::OpFlash> nullOpFlash;
    bool found_match = false;
    auto match_time = fm->time;
    for (size_t iop=0; iop<flash_v.size(); iop++){
      auto opflash = flash_v[iop];
      if (abs( opflash->Time() - match_time) < 0.01){
        found_match = true;
        matched_opflashes.push_back(opflash);
        break;
      }
    }
    if (found_match == false) matched_opflashes.push_back(nullOpFlash);
  }
  if (matched_opflashes.size() != fm_v.size()) std::cout << "mismatched opflash and simpleflash vector sizes!" << std::endl;
  return matched_opflashes;
}

DEFINE_ART_MODULE(sbnd::LightCaloAna)

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

  std::vector<std::string> _opflash_producer_v; ///< The OpFlash producers (to be set)
  std::string _slice_producer; ///< The Slice producer (to be set)
  std::string _flashmatch_producer;

  int _run, _subrun, _event;
};


sbnd::LightCaloAna::LightCaloAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  _opflash_producer_v = p.get<std::vector<std::string>>("OpFlashProducers");
  _slice_producer = p.get<std::string>("SliceProducer");
  _flashmatch_producer = p.get<std::string>("FlashMatchProducer");
  // Call appropriate consumes<>() for any products to be retrieved by this module.
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

  std::vector<art::Ptr<recob::Slice>> passed_slices_v; 

  for (size_t n_slice=0; n_slice < slice_v.size(); n_slice++){
    float nu_score = -9999;
    float fm_score = -9999;
    auto slice = slice_v[n_slice];
    std::vector<art::Ptr<recob::PFParticle>> pfp_v = slice_to_pfps.at(n_slice);
    for (size_t n_pfp=0; n_pfp < pfp_v.size(); n_pfp++){
      auto pfp = pfp_v[n_pfp];

      // only select the PRIMARY pfp 
      if(!pfp->IsPrimary() && (abs(pfp->PdgCode()) == 12 || abs(pfp->PdgCode()) == 14|| abs(pfp->PdgCode()) == 16))
        continue;

      // if primary, get nu-score 
      const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> pfpmeta_v = pfp_to_meta.at(pfp->Self());
      const art::Ptr<larpandoraobj::PFParticleMetadata> pfpmeta = pfpmeta_v.front();
      larpandoraobj::PFParticleMetadata::PropertiesMap propmap = pfpmeta->GetPropertiesMap();
      if (propmap.count("NuScore")) nu_score = propmap.at("NuScore");
      else nu_score = -1;
      // select slices that have nuscores > 0.4 
      if (nu_score < 0.4)
        continue;

      // get fm-score 
      std::vector<art::Ptr<sbn::SimpleFlashMatch>> fm_v = pfps_to_sfm.at(pfp.key());
      std::cout << "fm_v size: " << fm_v.size() << std::endl;
      for (size_t n_fm=0; n_fm < fm_v.size(); n_fm++){
        auto fm = fm_v.at(n_fm);
        fm_score = fm->score.total;
        std::cout << "fm_score: " << fm_score << std::endl;
      } // end flashmatch loop
    } // end pfp loop
    if (fm_score < 7)
      passed_slices_v.push_back(slice);
  } // end slice loop
} // end analyze

DEFINE_ART_MODULE(sbnd::LightCaloAna)

////////////////////////////////////////////////////////////////////////
// Class:       SBNDOpT0Finder
// Plugin Type: producer (art v3_05_01)
// File:        SBNDOpT0Finder_module.cc
//
// Generated at Sun Oct  4 17:27:36 2020 by Marco Del Tutto using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcore/Geometry/Geometry.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "sbndcode/OpT0Finder/flashmatch/Base/OpT0FinderTypes.h"
#include "sbndcode/OpT0Finder/flashmatch/Base/FlashMatchManager.h"

#include "TFile.h"
#include "TTree.h"

#include <memory>


class SBNDOpT0Finder;


class SBNDOpT0Finder : public art::EDProducer {
public:
  explicit SBNDOpT0Finder(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SBNDOpT0Finder(SBNDOpT0Finder const&) = delete;
  SBNDOpT0Finder(SBNDOpT0Finder&&) = delete;
  SBNDOpT0Finder& operator=(SBNDOpT0Finder const&) = delete;
  SBNDOpT0Finder& operator=(SBNDOpT0Finder&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  std::vector<flashmatch::QCluster_t> GetLighClusters(art::Event& e);
  float GetNPhotons(const float charge, const art::Ptr<recob::PFParticle> &pfp);

  ::flashmatch::FlashMatchManager _mgr; ///< The flash matching manager
  std::vector<flashmatch::FlashMatch_t> _result; ///< Matching result will be stored here

  std::string _opflash_producer;
  std::string _slice_producer;

  double _flash_trange_start;          ///<
  double _flash_trange_end;            ///<

  float _charge_to_n_photons_track;
  float _charge_to_n_photons_shower;

  TTree* _tree1;
  int _run, _subrun, _event, _matchid, _flashid;
  std::vector<double>               _score, _t0;
  std::vector<double>               _qll_xmin, _tpc_xmin;
  std::vector<double>              _beam_flash_spec;
  std::vector<std::vector<double>> _hypo_flash_spec;
  std::vector<double>              _numc_flash_spec;
  int _fv, _ccnc, _nupdg;
  std::vector<double>    _xfixed_hypo_spec;
  double _xfixed_chi2, _xfixed_ll;

  std::vector<float> _dep_x, _dep_y, _dep_z, _dep_charge, _dep_n_photons;
  std::vector<int> _dep_slice;
};


SBNDOpT0Finder::SBNDOpT0Finder(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  ::art::ServiceHandle<geo::Geometry> geo;

  _opflash_producer = p.get<std::string>("OpFlashProducer");
  _slice_producer = p.get<std::string>("SliceProducer");

  _flash_trange_start = p.get<double>("FlashVetoTimeStart", 0);
  _flash_trange_end = p.get<double>("FlashVetoTimeEnd", 2);

  _charge_to_n_photons_track = p.get<float>("ChargeToNPhotonsTrack");
  _charge_to_n_photons_track = p.get<float>("ChargeToNPhotonsShower");

  _mgr.Configure(p.get<flashmatch::Config_t>("FlashMatchConfig"));

  art::ServiceHandle<art::TFileService> fs;
  _tree1 = fs->make<TTree>("flashmatchtree","");
  _tree1->Branch("run",             &_run,                             "run/I");
  _tree1->Branch("subrun",          &_subrun,                          "subrun/I");
  _tree1->Branch("event",           &_event,                           "event/I");
  // _tree1->Branch("beam_flash_spec", "std::vector<double>",             &_beam_flash_spec);
  // _tree1->Branch("hypo_flash_spec", "std::vector<std::vector<double>>",&_hypo_flash_spec);
  // _tree1->Branch("numc_flash_spec", "std::vector<double>",             &_numc_flash_spec);
  // _tree1->Branch("score",           "std::vector<double>",             &_score);
  // _tree1->Branch("t0",              "std::vector<double>",             &_t0);
  // _tree1->Branch("qll_xmin",        "std::vector<double>",             &_qll_xmin);
  // _tree1->Branch("tpc_xmin",        "std::vector<double>",             &_tpc_xmin);
  // _tree1->Branch("xfixed_hypo_spec","std::vector<double>",             &_xfixed_hypo_spec);
  // _tree1->Branch("fv",              &_fv,                              "fv/I");
  // _tree1->Branch("ccnc",            &_ccnc,                            "ccnc/I");
  // _tree1->Branch("nupdg",           &_nupdg,                           "nupdg/I");
  _tree1->Branch("dep_slice", "std::vector<int>", &_dep_slice);
  _tree1->Branch("dep_x", "std::vector<float>", &_dep_x);
  _tree1->Branch("dep_y", "std::vector<float>", &_dep_y);
  _tree1->Branch("dep_z", "std::vector<float>", &_dep_z);
  _tree1->Branch("dep_charge", "std::vector<float>", &_dep_charge);
  _tree1->Branch("dep_n_photons", "std::vector<float>", &_dep_n_photons);
}

void SBNDOpT0Finder::produce(art::Event& e)
{
  _mgr.Reset();
  _result.clear();
  _mgr.PrintConfig();

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();

  ::art::Handle<std::vector<recob::OpFlash>> flash_h;
  e.getByLabel(_opflash_producer, flash_h);
  if(!flash_h.isValid() || flash_h->empty()) {
    std::cout << "[SBNDOpT0Finder] Don't have good flashes." << std::endl;
    return;
  }

  ::art::ServiceHandle<geo::Geometry> geo;

  int nBeamFlashes = 0;
  std::vector<::flashmatch::Flash_t> beam_flashes;

  for (size_t n = 0; n < flash_h->size(); n++) {

    auto const& flash = (*flash_h)[n];

    std::cout << "[SBNDOpT0Finder] Flash time from " << _opflash_producer << ": " << flash.Time() << std::endl;
    if(flash.Time() < _flash_trange_start || _flash_trange_end < flash.Time()) {
      continue;
    }

    nBeamFlashes++;

    // Construct a Flash_t
    ::flashmatch::Flash_t f;
    f.x = f.x_err = 0;
    f.pe_v.resize(geo->NOpDets());
    f.pe_err_v.resize(geo->NOpDets());
    for (unsigned int i = 0; i < f.pe_v.size(); i++) {
      unsigned int opdet = geo->OpDetFromOpChannel(i);
      f.pe_v[opdet] = flash.PE(i);
      f.pe_err_v[opdet] = sqrt(flash.PE(i));
    }
    f.y = flash.YCenter();
    f.z = flash.ZCenter();
    f.y_err = flash.YWidth();
    f.z_err = flash.ZWidth();
    f.time = flash.Time();
    f.idx = nBeamFlashes-1;
    beam_flashes.resize(nBeamFlashes);
    beam_flashes[nBeamFlashes-1] = f;

  } // flash loop

  // Don't waste other time if there are no flashes in the beam spill
  if (nBeamFlashes == 0) {
    std::cout << "[NeutrinoFlashMatch] Zero beam flashes in this event." << std::endl;
    return;
  }

  // TODO, Only pick one flash for now
  ::flashmatch::Flash_t f = beam_flashes[0];

  _beam_flash_spec.resize(f.pe_v.size());
  _beam_flash_spec = f.pe_v;

  // Emplace flash to Flash Matching Manager
  // _mgr.Emplace(std::move(f));

  // Get all the ligh clusters
  auto light_cluster_v = GetLighClusters(e);
}

std::vector<flashmatch::QCluster_t> SBNDOpT0Finder::GetLighClusters(art::Event& e) {
  // One slice is one QCluster_t.
  // Start from a slice, get all the PFParticles, from there get all the Spacepoints, from
  // there get all the hits on the collection plane.
  // Use the charge on the collection plane to estimate the light, and the 3D spacepoint
  // position to for the 3D location.

  std::vector<flashmatch::QCluster_t> light_cluster_v;

  ::art::Handle<std::vector<recob::Slice>> slice_h;
  e.getByLabel(_slice_producer, slice_h);
  if(!slice_h.isValid() || slice_h->empty()) {
    std::cout << "[SBNDOpT0Finder] Don't have good Slices." << std::endl;
    return light_cluster_v;
  }

  ::art::Handle<std::vector<recob::PFParticle>> pfp_h;
  e.getByLabel(_slice_producer, pfp_h);
  if(!pfp_h.isValid() || pfp_h->empty()) {
    std::cout << "[SBNDOpT0Finder] Don't have good PFParticles." << std::endl;
    return light_cluster_v;
  }

  ::art::Handle<std::vector<recob::SpacePoint>> spacepoint_h;
  e.getByLabel(_slice_producer, spacepoint_h);
  if(!spacepoint_h.isValid() || spacepoint_h->empty()) {
    std::cout << "[SBNDOpT0Finder] Don't have good SpacePoints." << std::endl;
    return light_cluster_v;
  }

  // Construct the vector of Slices
  std::vector<art::Ptr<recob::Slice>> slice_v;
  art::fill_ptr_vector(slice_v, slice_h);

  // Get the associations between slice->pfp->spacepoint->hit
  art::FindManyP<recob::PFParticle> slice_to_pfps (slice_h, e, _slice_producer);
  art::FindManyP<recob::SpacePoint> pfp_to_spacepoints (pfp_h, e, _slice_producer);
  art::FindManyP<recob::Hit> spacepoint_to_hits (spacepoint_h, e, _slice_producer);

  // Loop over the Slices
  for (size_t n_slice = 0; n_slice < slice_h->size(); n_slice++) {
    std::cout << "[SBNDOpT0Finder] At slice " << n_slice << std::endl;

    flashmatch::QCluster_t light_cluster;

    _dep_slice.clear();
    _dep_x.clear();
    _dep_y.clear();
    _dep_z.clear();
    _dep_charge.clear();
    _dep_n_photons.clear();

    // Get the associated PFParticles
    std::vector<art::Ptr<recob::PFParticle>> pfp_v = slice_to_pfps.at(n_slice);

    for (size_t n_pfp = 0; n_pfp < pfp_v.size(); n_pfp++) {
      std::cout << "[SBNDOpT0Finder] At pfp " << n_pfp << std::endl;

      auto pfp = pfp_v[n_pfp];

      // Get the associated SpacePoints
      std::vector<art::Ptr<recob::SpacePoint>> spacepoint_v = pfp_to_spacepoints.at(pfp.key());

      for (size_t n_spacepoint = 0; n_spacepoint < spacepoint_v.size(); n_spacepoint++) {
        // std::cout << "[SBNDOpT0Finder] At spacepoint " << n_spacepoint << std::endl;

        auto spacepoint = spacepoint_v[n_spacepoint];

        // Get the associated hits
        std::vector<art::Ptr<recob::Hit>> hit_v = spacepoint_to_hits.at(spacepoint.key());

        for (size_t n_hit = 0; n_hit < hit_v.size(); n_hit++) {
          // std::cout << "[SBNDOpT0Finder] At hit " << n_hit << std::endl;

          auto hit = hit_v[n_hit];
          if (hit->View() != geo::kZ) {
            continue;
          }
          const auto &position(spacepoint->XYZ());
          const auto charge(hit->Integral());

          light_cluster.emplace_back(position[0],
                                     position[1],
                                     position[2],
                                     GetNPhotons(charge, pfp));

          _dep_slice.push_back(light_cluster_v.size());
          _dep_x.push_back(position[0]);
          _dep_y.push_back(position[1]);
          _dep_z.push_back(position[2]);
          _dep_charge.push_back(charge);
          _dep_n_photons.push_back(GetNPhotons(charge, pfp));
        }
      } // End loop over Spacepoints
    } // End loop over PFParticle

    _tree1->Fill();
    std::cout << "[SBNDOpT0Finder] filled n_slice " << n_slice << std::endl;
    light_cluster_v.emplace_back(light_cluster);

  } // End loop over Slices

  return light_cluster_v;
}

float SBNDOpT0Finder::GetNPhotons(const float charge,
                                  const art::Ptr<recob::PFParticle> &pfp) {
  return charge * (::lar_pandora::LArPandoraHelper::IsTrack(pfp) ? _charge_to_n_photons_track
                                                                 : _charge_to_n_photons_shower);
}


DEFINE_ART_MODULE(SBNDOpT0Finder)














////////////////////////////////////////////////////////////////////////
// Class:       CRTPERecoAnalysis
// Plugin Type: analyzer (Unknown Unknown)
// File:        CRTPERecoAnalysis_module.cc
//
// Generated at Wed Jan 24 12:01:32 2024 by Jiaoyang Li using cetskelgen
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
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "TTree.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "sbnobj/Common/CRT/CRTHit.hh"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "sbndcode/CRT/CRTUtils/CRTBackTracker.h"

class CRTPERecoAnalysis;


class CRTPERecoAnalysis : public art::EDAnalyzer {
public:
  explicit CRTPERecoAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTPERecoAnalysis(CRTPERecoAnalysis const&) = delete;
  CRTPERecoAnalysis(CRTPERecoAnalysis&&) = delete;
  CRTPERecoAnalysis& operator=(CRTPERecoAnalysis const&) = delete;
  CRTPERecoAnalysis& operator=(CRTPERecoAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  // Declare member data here.
  sbnd::CRTBackTracker _crt_back_tracker;
  std::string _crthit_label;
  
  bool _debug;
  bool _data_mode;
  TTree* _tree;

  int _run, _subrun, _event;

  // CRT hit infos. 
  std::vector<std::vector<float>>  _chit_strip_hit_distance_to_readout; ///< CRT hit, distance to readout per strip
  std::vector<std::vector<float>>  _chit_strip_hit_uncorrected_PE; ///< CRT hit, uncorrected PE per strip
  std::vector<std::vector<float>>  _chit_strip_hit_corrected_PE; ///< CRT hit, corrected PE per strip

  // truth info for crt hit backtracking...
  std::vector<int> _chit_backtrack_pdg; ///< CRT hit, truth information of the pdg code 
  std::vector<double> _chit_backtrack_energy; ///< CRT hit, truth information of the particle energy
  std::vector<double> _chit_backtrack_deposited_energy; ///< CRT hit, truth information of the deposited energy for both upstream and downstream
  std::vector<double> _chit_backtrack_purity; ///< CRT hit, truth information of selection purity
};


CRTPERecoAnalysis::CRTPERecoAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  _crthit_label      = p.get<std::string>("CRTHitLabel", "crthit");
  _debug             = p.get<bool>("Debug", false);
  _data_mode         = p.get<bool>("DataMode", false);
  _crt_back_tracker  = p.get<fhicl::ParameterSet>("CRTBackTracker", fhicl::ParameterSet());

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("PERecoAnaTree","CRT PE reco analysis tree");
  _tree->Branch("run",&_run,"run/I");
  _tree->Branch("subrun",&_subrun,"subrun/I");
  _tree->Branch("event",&_event,"event/I");
  _tree->Branch("chit_strip_hit_distance_to_readout",&_chit_strip_hit_distance_to_readout);
  _tree->Branch("chit_strip_hit_uncorrected_PE",&_chit_strip_hit_uncorrected_PE);
  _tree->Branch("chit_strip_hit_corrected_PE",&_chit_strip_hit_corrected_PE);
  if (!_data_mode){
    _tree->Branch("chit_backtrack_pdg",&_chit_backtrack_pdg);
    _tree->Branch("chit_backtrack_energy",&_chit_backtrack_energy);
    _tree->Branch("chit_backtrack_deposited_energy",&_chit_backtrack_deposited_energy);
    _tree->Branch("chit_backtrack_purity",&_chit_backtrack_purity);
  }
}

void CRTPERecoAnalysis::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  _run = e.run();
  _subrun = e.subRun();
  _event = e.event();

  if (!_data_mode) _crt_back_tracker.Initialize(e); // Initialise the backtrack alg. 

  // 
  // Fill the CRT hit info.
  //
  art::Handle<std::vector<sbn::crt::CRTHit> > crt_hit_h;
  e.getByLabel(_crthit_label, crt_hit_h);
  if (!crt_hit_h.isValid()){
    std::cout << "CRTPERecoAnalysis::analyze: invalid CRTHit handle" << std::endl;
    return;
  }
  std::vector<art::Ptr<sbn::crt::CRTHit> > crt_hit_v;
  art::fill_ptr_vector(crt_hit_v, crt_hit_h);
  size_t n_crt_hit = crt_hit_v.size();
  if (_debug) std::cout << "CRTPERecoAnalysis::analyze: number of CRT hits: " << n_crt_hit << std::endl;

  _chit_strip_hit_distance_to_readout.resize(n_crt_hit);
  _chit_strip_hit_corrected_PE.resize(n_crt_hit);
  _chit_strip_hit_uncorrected_PE.resize(n_crt_hit);
  if (!_data_mode){
    _chit_backtrack_pdg.resize(n_crt_hit);
    _chit_backtrack_energy.resize(n_crt_hit);
    _chit_backtrack_deposited_energy.resize(n_crt_hit);
    _chit_backtrack_purity.resize(n_crt_hit);
  }

  for (size_t ihit=0; ihit<n_crt_hit; ++ihit){
    auto crt_hit = crt_hit_v[ihit]; 
    _chit_strip_hit_distance_to_readout[ihit].resize(2);
    _chit_strip_hit_corrected_PE[ihit].resize(2);
    _chit_strip_hit_uncorrected_PE[ihit].resize(2);

    _chit_strip_hit_distance_to_readout[ihit][0] = crt_hit->strip0info[0];
    _chit_strip_hit_distance_to_readout[ihit][1] = crt_hit->strip1info[0];

    _chit_strip_hit_uncorrected_PE[ihit][0] = crt_hit->strip0info[1];
    _chit_strip_hit_uncorrected_PE[ihit][1] = crt_hit->strip1info[1];

    _chit_strip_hit_corrected_PE[ihit][0] = crt_hit->strip0info[2];
    _chit_strip_hit_corrected_PE[ihit][1] = crt_hit->strip1info[2];

    if (!_data_mode){
      const sbnd::CRTBackTracker::TruthMatchMetrics truthMatch = _crt_back_tracker.TruthMatrixFromTotalEnergy(e, crt_hit);
      _chit_backtrack_pdg[ihit] = truthMatch.pdg;
      _chit_backtrack_energy[ihit] = truthMatch.particle_energy;
      _chit_backtrack_deposited_energy[ihit] = truthMatch.depEnergy_total;
      _chit_backtrack_purity[ihit] = truthMatch.purity;
    }
  }

  _tree->Fill();

}


DEFINE_ART_MODULE(CRTPERecoAnalysis)

////////////////////////////////////////////////////////////////////////
// Class:       CRTRecoAnalysis
// Plugin Type: analyzer 
// File:        CRTRecoAnalysis_module.cc
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

class CRTRecoAnalysis;


class CRTRecoAnalysis : public art::EDAnalyzer {
public:
  explicit CRTRecoAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTRecoAnalysis(CRTRecoAnalysis const&) = delete;
  CRTRecoAnalysis(CRTRecoAnalysis&&) = delete;
  CRTRecoAnalysis& operator=(CRTRecoAnalysis const&) = delete;
  CRTRecoAnalysis& operator=(CRTRecoAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:
  // Declare member data here.
  sbnd::CRTBackTracker _crt_back_tracker;
  std::string _crthit_label;
  std::string _febdata_label;
  
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
  std::vector<double> _chit_t1_diff; ///< CRT hit, difference between the hit time and the time of the first particle in the track
  std::vector<double> _chit_x; ///< CRT hit, x position of the hit
  std::vector<double> _chit_y; ///< CRT hit, y position of the hit
  std::vector<double> _chit_z; ///< CRT hit, z position of the hit
  std::vector<double> _chit_t1; ///< CRT hit, time of the hit
  std::vector<double> _chit_true_t1; ///< CRT hit, time of the hit
};


CRTRecoAnalysis::CRTRecoAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  _crthit_label      = p.get<std::string>("CRTHitLabel", "crthit");
  _debug             = p.get<bool>("Debug", false);
  _data_mode         = p.get<bool>("DataMode", false);
  _crt_back_tracker  = p.get<fhicl::ParameterSet>("CRTBackTracker", fhicl::ParameterSet());
  _febdata_label     = p.get<std::string>("FEBDataLabel", "crtsim");

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("PERecoAnaTree","CRT PE reco analysis tree");
  _tree->Branch("run",&_run,"run/I");
  _tree->Branch("subrun",&_subrun,"subrun/I");
  _tree->Branch("event",&_event,"event/I");
  _tree->Branch("chit_strip_hit_distance_to_readout",&_chit_strip_hit_distance_to_readout);
  _tree->Branch("chit_strip_hit_uncorrected_PE",&_chit_strip_hit_uncorrected_PE);
  _tree->Branch("chit_strip_hit_corrected_PE",&_chit_strip_hit_corrected_PE);
  _tree->Branch("chit_t1_diff",&_chit_t1_diff);
  _tree->Branch("chit_x",&_chit_x);
  _tree->Branch("chit_y",&_chit_y);
  _tree->Branch("chit_z",&_chit_z);
  _tree->Branch("chit_t1",&_chit_t1);
  if (!_data_mode){
    _tree->Branch("chit_backtrack_pdg",&_chit_backtrack_pdg);
    _tree->Branch("chit_backtrack_energy",&_chit_backtrack_energy);
    _tree->Branch("chit_backtrack_deposited_energy",&_chit_backtrack_deposited_energy);
    _tree->Branch("chit_backtrack_purity",&_chit_backtrack_purity);
    _tree->Branch("chit_true_t1",&_chit_true_t1);
  }
}

void CRTRecoAnalysis::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  _run    = e.run();
  _subrun = e.subRun();
  _event  = e.event();

  if (!_data_mode) _crt_back_tracker.Initialize(e); // Initialise the backtrack alg. 

  // 
  // Fill the CRT hit info.
  //
  art::Handle<std::vector<sbn::crt::CRTHit> > crt_hit_h;
  e.getByLabel(_crthit_label, crt_hit_h);
  if (!crt_hit_h.isValid()){
    std::cout << "CRTRecoAnalysis::analyze: invalid CRTHit handle" << std::endl;
    return;
  }
  std::vector<art::Ptr<sbn::crt::CRTHit> > crt_hit_v;
  art::fill_ptr_vector(crt_hit_v, crt_hit_h);
  art::FindManyP<sbnd::crt::FEBData> crt_hit_to_feb_data(crt_hit_h, e, _crthit_label);

  //
  // Get the FEB Data
  //
  art::Handle<std::vector<sbnd::crt::FEBData>> feb_data_h;
  art::FindManyP<sim::AuxDetIDE, sbnd::crt::FEBTruthInfo> *febdata_to_ides;
  if (!_data_mode){
    e.getByLabel(_febdata_label, feb_data_h);
    if(!feb_data_h.isValid()){
      std::cout << "FEBData product " << _febdata_label << " not found..." << std::endl;
      throw std::exception();
    }
    febdata_to_ides = new art::FindManyP<sim::AuxDetIDE, sbnd::crt::FEBTruthInfo>(feb_data_h, e, _febdata_label);
  }

  size_t n_crt_hit = crt_hit_v.size();
  if (_debug) std::cout << "CRTRecoAnalysis::analyze: number of CRT hits: " << n_crt_hit << std::endl;

  _chit_strip_hit_distance_to_readout.resize(n_crt_hit);
  _chit_strip_hit_corrected_PE.resize(n_crt_hit);
  _chit_strip_hit_uncorrected_PE.resize(n_crt_hit);
  _chit_t1_diff.resize(n_crt_hit);
  _chit_x.resize(n_crt_hit);
  _chit_y.resize(n_crt_hit);
  _chit_z.resize(n_crt_hit);
  _chit_t1.resize(n_crt_hit);
  if (!_data_mode){
    _chit_backtrack_pdg.resize(n_crt_hit);
    _chit_backtrack_energy.resize(n_crt_hit);
    _chit_backtrack_deposited_energy.resize(n_crt_hit);
    _chit_backtrack_purity.resize(n_crt_hit);
    _chit_true_t1.resize(n_crt_hit);
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
    
    _chit_t1_diff[ihit] = crt_hit->ts0_ns_corr;
    _chit_x[ihit] = crt_hit->x_pos;
    _chit_y[ihit] = crt_hit->y_pos;
    _chit_z[ihit] = crt_hit->z_pos;
    _chit_t1[ihit] = crt_hit->ts1_ns;

    // CRTHit matches to AuxDetIDE.  **To-do** need to find the correct channel for CRT Hits. 
    auto feb_datas = crt_hit_to_feb_data.at(crt_hit.key()); // Use the Assn to find the AuxDetIDE from FEBData. 
    if(feb_datas.size() != 2) std::cout << "ERROR: CRTHit associated to " << feb_datas.size() << " FEBDatas" << std::endl;

    if(!_data_mode) {
      _chit_true_t1[ihit] = 0;
      size_t n_ides = 0;
      for (auto feb_data : feb_datas) {
        auto ide_v = febdata_to_ides->at(feb_data.key());
        for (auto ide : ide_v) {
          _chit_true_t1[ihit] += 0.5 * (ide->entryT + ide->exitT);
          n_ides++;
        }
      }
      _chit_true_t1[ihit] /= n_ides;
      if (_debug) std::cout << "CRTRecoAnalysis::analyze: true hit time: " << _chit_true_t1[ihit] << "; reco hit time: " << _chit_t1[ihit] - 1.7e6 << std::endl;
    }

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


DEFINE_ART_MODULE(CRTRecoAnalysis)

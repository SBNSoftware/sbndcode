////////////////////////////////////////////////////////////////////////
// Class:       SelectionAnaTruth
// Plugin Type: analyzer (Unknown Unknown)
// File:        SelectionAnaTruth_module.cc
//
// Generated at Mon Jan 15 10:39:03 2024 by Jiaoyang Li using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/FileBlock.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "TTree.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/Simulation/AuxDetHit.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTBackTracker.h"
#include <numeric>
#include <algorithm>  // for sort
#include <functional> // for greater

class SelectionAnaTruth;


class SelectionAnaTruth : public art::EDAnalyzer {
public:
  explicit SelectionAnaTruth(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SelectionAnaTruth(SelectionAnaTruth const&) = delete;
  SelectionAnaTruth(SelectionAnaTruth&&) = delete;
  SelectionAnaTruth& operator=(SelectionAnaTruth const&) = delete;
  SelectionAnaTruth& operator=(SelectionAnaTruth&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginSubRun(art::SubRun const& sr) override;
  void respondToOpenInputFile(const art::FileBlock& fb) override;
  void calculateMeanStd(std::vector<double> vec, double &mean, double &std);

private:

  // Declare member data here.
  std::string _auxDetHit_label;
  std::string _mctruth_label;

  bool _debug;
  bool _save_input_file_name;
  int  _interactionMode; // 0 -- dark neutrino 140 MeV; 1 -- dark neutrino 400 MeV; 2 -- dirt; 3 -- cosmic; 
  double _merge_radius; // merge radius for CRT hits.
  
  TTree* _tree;
  int _run, _subrun, _event;
  std::string _file_name;
  float _weight;                                 ///< Signal: the weight which store as the vertex of dark neutrino.
  
  int _n_hits_upstream, _n_hits_downstream;    ///< Number of CRT hits in the upstream/downstream CRTs

  double _timing_difference_betweenHits;       ///< Timing difference between the upstream and downstream CRT hits
  double _distance_betweenHits;                ///< Distance between the upstream and downstream CRT hits

  std::vector<double> _deposited_energy_downstream_vec;   ///< Deposited energy of the CRT hits
  double _deposited_energy_downstream_mean;               ///< Mean of the deposited energy of the CRT hits
  double _deposited_energy_downstream_std;                ///< Standard deviation of the deposited energy of the CRT hits

  double _deposited_energy_downstream_diff;         ///< Difference of the deposited energy between the upstream and downstream CRT hits
  double _deposited_energy_downstream_ratio;              ///< Ratio of the deposited energy between the upstream and downstream CRT hits
  double _deposited_energy_downstream_frac_diff;    ///< Fractional difference of the deposited energy between the upstream and downstream CRT hits

  // POT info. 
  std::string _pot_label;
  TTree* _sr_tree;
  int _sr_run, _sr_subrun;
  double _sr_begintime, _sr_endtime;
  double _sr_pot; ///< Number of POTs per subrun
  double _sr_spills; ///< Number of Spills per subrun
};


SelectionAnaTruth::SelectionAnaTruth(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  _auxDetHit_label = p.get<std::string>("AuxDetHitLabel", "largeant:LArG4DetectorServicevolAuxDetSensitiveCRTStripBERN");
  _debug = p.get<bool>("Debug", false);
  _pot_label = p.get<std::string>("POTLabel", "generator");
  _save_input_file_name = p.get<bool>("SaveInputFileName", true);
  _interactionMode = p.get<int>("InteractionMode");
  _mctruth_label = p.get<std::string>("MCTruthLabel", "generator");
  _merge_radius = p.get<double>("MergeRadius", 4.25);

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("tree", "CRT Selection Analysis Truth");
  _tree->Branch("run", &_run, "run/I");
  _tree->Branch("subrun", &_subrun, "subrun/I");
  _tree->Branch("event", &_event, "event/I");
  _tree->Branch("InteractionMode", &_interactionMode, "InteractionMode/I");
  if (_save_input_file_name) _tree->Branch("file_name", &_file_name);
  _tree->Branch("weight", &_weight, "weight/F");
  _tree->Branch("n_hits_upstream", &_n_hits_upstream, "n_hits_upstream/I");
  _tree->Branch("n_hits_downstream", &_n_hits_downstream, "n_hits_downstream/I");
  _tree->Branch("timing_difference_betweenHits", &_timing_difference_betweenHits, "timing_difference_betweenHits/D");
  _tree->Branch("distance_betweenHits", &_distance_betweenHits, "distance_betweenHits/D");
  
  _tree->Branch("deposited_energy_downstream_vec", &_deposited_energy_downstream_vec);
  _tree->Branch("deposited_energy_downstream_mean", &_deposited_energy_downstream_mean, "deposited_energy_downstream_mean/D");
  _tree->Branch("deposited_energy_downstream_std", &_deposited_energy_downstream_std, "deposited_energy_downstream_std/D");
  _tree->Branch("deposited_energy_downstream_diff", &_deposited_energy_downstream_diff, "deposited_energy_downstream_diff/D");
  _tree->Branch("deposited_energy_downstream_ratio", &_deposited_energy_downstream_ratio, "deposited_energy_downstream_ratio/D");
  _tree->Branch("deposited_energy_downstream_frac_diff", &_deposited_energy_downstream_frac_diff, "_deposited_energy_downstream_frac_diff/D");

  // POT tree. 
  _sr_tree = tfs->make<TTree>("sr_tree", "SubRun POT");
  _sr_tree->Branch("run", &_sr_run, "run/I");
  _sr_tree->Branch("subrun", &_sr_subrun, "subrun/I");
  _sr_tree->Branch("begintime", &_sr_begintime, "begintime/D");
  _sr_tree->Branch("endtime", &_sr_endtime, "endtime/D");
  _sr_tree->Branch("pot", &_sr_pot, "pot/D");
}

void SelectionAnaTruth::analyze(art::Event const& e)
{
  // clear all declared vectors
  _deposited_energy_downstream_vec.clear();

  _run = e.run();
  _subrun = e.subRun();
  _event = e.event();
  if(_debug) {
    std::cout << "This is file: " << _file_name << std::endl;
    std::cout << "This is run:subrun:event " << _run << ":" << _subrun << ":" << _event << std::endl;
  }

  // get the weight info. 
  if (_interactionMode==0 || _interactionMode==1) {
    art::Handle<std::vector<simb::MCTruth>> mctruth_handle;
    e.getByLabel(_mctruth_label, mctruth_handle);
    if(!mctruth_handle.isValid()){
      std::cout << "MCTruth product " << _mctruth_label << " not found..." << std::endl;
      throw std::exception();
    }
    std::vector<art::Ptr<simb::MCTruth>> mct_v;
    art::fill_ptr_vector(mct_v, mctruth_handle);
    for (size_t imct=0; imct< mct_v.size(); imct++){
      auto mct = mct_v[imct];
      for (int ipar = 0; ipar < mct->NParticles(); ipar++) {
        if (ipar==2) {
          auto particle = mct->GetParticle(ipar);
          _weight = particle.Vx();
        }
      }
    }
  }else{
      _weight = 1.;
  }

  //
  // Get the AuxDetHit information. 
  //
  art::Handle<std::vector<sim::AuxDetHit>> auxDetHitHandle;
  e.getByLabel(_auxDetHit_label, auxDetHitHandle);
  if(!auxDetHitHandle.isValid()){
    std::cout << "AuxDetHit product " << _auxDetHit_label << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<sim::AuxDetHit>> auxDetHit_v;
  art::fill_ptr_vector(auxDetHit_v, auxDetHitHandle);

  // Fill in the map for trackID to various info.
  std::map<int, std::vector<double>> trackID_to_t_map;
  std::map<int, std::vector<double>> trackID_to_hitx_map;
  std::map<int, std::vector<double>> trackID_to_hity_map;
  std::map<int, std::vector<double>> trackID_to_hitz_map;
  std::map<int, std::vector<double>> trackID_to_depositedE_map;
  for (size_t iadh = 0; iadh < auxDetHit_v.size(); iadh++) {
    auto const& auxdethit = auxDetHit_v[iadh];
    int trackID = auxdethit->GetTrackID();
    double t =  0.5 * (auxdethit->GetEntryT() + auxdethit->GetExitT());
    double hitx = 0.5 * (auxdethit->GetEntryX() + auxdethit->GetExitX());
    double hity = 0.5 * (auxdethit->GetEntryY() + auxdethit->GetExitY());
    double hitz = 0.5 * (auxdethit->GetEntryZ() + auxdethit->GetExitZ());
    double depositedE = auxdethit->GetEnergyDeposited();

    trackID_to_t_map[trackID].push_back(t);
    trackID_to_hitx_map[trackID].push_back(hitx);
    trackID_to_hity_map[trackID].push_back(hity);
    trackID_to_hitz_map[trackID].push_back(hitz);
    trackID_to_depositedE_map[trackID].push_back(depositedE);
  }

  // Loop over the map to calculate variables. 
  std::vector<double> hit_t_vec_downstream, hit_x_vec_downstream, hit_y_vec_downstream, hit_depositedE_vec_downstream;
  _n_hits_upstream = 0, _n_hits_downstream = 0;
  for (const auto& [trackid, t_vec] : trackID_to_t_map) {
    auto hitx_vec = trackID_to_hitx_map[trackid];
    auto hity_vec = trackID_to_hity_map[trackid];
    auto hitz_vec = trackID_to_hitz_map[trackid]; 
    auto depositedE_vec = trackID_to_depositedE_map[trackid]; 

    double hit_x_upstream(0), hit_y_upstream(0);
    double depositedE_upstream(0), hit_t_uptream(0);
    double hit_x_downstream(0), hit_y_downstream(0);
    double depositedE_downstream(0), hit_t_downstream(0);
    int total_hit_this_track_upstream(0), total_hit_this_track_downstream(0);
    bool IsNewHitDownstream(false);
    // loop through hitz_vec to fill in the information above. 
    for (size_t ihit=0; ihit<hitz_vec.size(); ihit++){
      if (depositedE_vec[ihit]<5e-2) continue; // remove low deposited energy particle. 
      // upstream
      if (hitz_vec[ihit]<0.) { 
          hit_x_upstream += hitx_vec[ihit];
          hit_y_upstream += hity_vec[ihit];  
          hit_t_uptream  += t_vec[ihit];
          depositedE_upstream += depositedE_vec[ihit];  
          total_hit_this_track_upstream++;
          if (total_hit_this_track_upstream==2) { // X-Y coincidence and only fill in once. 
              _n_hits_upstream++; 
          }
      }
      // downstream
      if (hitz_vec[ihit]>0.) { 
        bool isThisHitTooCloseToPreviousHit(false);
        for (size_t iii=0; iii<hit_x_vec_downstream.size(); iii++){
            double distance = std::hypot(hitx_vec[ihit]-hit_x_vec_downstream[iii], hity_vec[ihit]-hit_y_vec_downstream[iii]);
            if (distance < _merge_radius){
                if (_debug) std::cout<<"Hit distance from previous track ID: "<<distance<<std::endl;
                isThisHitTooCloseToPreviousHit = true;
            }
        }
        if (!isThisHitTooCloseToPreviousHit){ // if this hit is not too close to previous hit, then fill in.
            hit_x_downstream += hitx_vec[ihit];
            hit_y_downstream += hity_vec[ihit];  
            hit_t_downstream += t_vec[ihit];
            depositedE_downstream += depositedE_vec[ihit];

            total_hit_this_track_downstream++;
            
            if(total_hit_this_track_downstream==2){ // X-Y coincidence and only fill in once. 
                _n_hits_downstream++;  
                IsNewHitDownstream = true;
            }
        } 
      }
    }
    if (IsNewHitDownstream){
      hit_t_vec_downstream.push_back(hit_t_downstream/total_hit_this_track_downstream);
      hit_x_vec_downstream.push_back(hit_x_downstream/total_hit_this_track_downstream);
      hit_y_vec_downstream.push_back(hit_y_downstream/total_hit_this_track_downstream);
      _deposited_energy_downstream_vec.push_back(depositedE_downstream/total_hit_this_track_downstream);
    }
  }
  if (hit_t_vec_downstream.size()<2){
    _timing_difference_betweenHits = -1; 
    _distance_betweenHits = -1;
    _deposited_energy_downstream_diff = -1;
    _deposited_energy_downstream_ratio = -1;
    _deposited_energy_downstream_frac_diff = -1;
  }
  // calulate the timing difference, distance between hits in the downstream. 
  for (size_t ihit=0; ihit<hit_t_vec_downstream.size(); ihit++){
    for (size_t jhit=ihit+1; jhit<hit_t_vec_downstream.size(); jhit++){
      _timing_difference_betweenHits = std::abs(hit_t_vec_downstream[ihit] - hit_t_vec_downstream[jhit]);

      _distance_betweenHits = std::hypot(hit_x_vec_downstream[ihit]-hit_x_vec_downstream[jhit], hit_y_vec_downstream[ihit]-hit_y_vec_downstream[jhit]);

      _deposited_energy_downstream_diff = std::abs(_deposited_energy_downstream_vec[ihit] - _deposited_energy_downstream_vec[jhit]);
      if (_deposited_energy_downstream_vec[ihit] < _deposited_energy_downstream_vec[jhit]){
        _deposited_energy_downstream_ratio = _deposited_energy_downstream_vec[ihit] / _deposited_energy_downstream_vec[jhit];
      }else{
        _deposited_energy_downstream_ratio = _deposited_energy_downstream_vec[jhit] / _deposited_energy_downstream_vec[ihit];
      }
      _deposited_energy_downstream_frac_diff = _deposited_energy_downstream_diff / (_deposited_energy_downstream_vec[ihit]+_deposited_energy_downstream_vec[jhit]);
      
      if (_debug) std::cout<<"timing difference: "<<_timing_difference_betweenHits<<", distance between hits: "<<_distance_betweenHits<<", deposited energy difference: "<<_deposited_energy_downstream_diff<<", deposited energy ratio: "<<_deposited_energy_downstream_ratio<<", deposited energy fractional difference: "<<_deposited_energy_downstream_frac_diff<<std::endl;
    }
  }

  // calculate the mean and std of the deposited energy of the downstream hits.
  calculateMeanStd(_deposited_energy_downstream_vec, _deposited_energy_downstream_mean, _deposited_energy_downstream_std);

  _tree->Fill();
}

void SelectionAnaTruth::beginSubRun(art::SubRun const& sr) 
{
  _sr_run       = sr.run();
  _sr_subrun    = sr.subRun();
  _sr_begintime = sr.beginTime().value();
  _sr_endtime   = sr.endTime().value();

  art::Handle<sumdata::POTSummary> pot_handle;
  sr.getByLabel(_pot_label, pot_handle);

  if (pot_handle.isValid()) {
    _sr_pot = pot_handle->totpot;
    _sr_spills = pot_handle->totspills;
  } else {
    _sr_pot = 0.;
    _sr_spills = 0.;
  }

  if (_debug) std::cout << "POT for this subrun: " << _sr_pot << " (" << _sr_spills << " spills)" << std::endl;
  _sr_tree->Fill();
}

void SelectionAnaTruth::respondToOpenInputFile(const art::FileBlock& fb)
{
  _file_name = fb.fileName();
}

void SelectionAnaTruth::calculateMeanStd(std::vector<double> vec, double &mean, double &std){
  if (vec.size()>=1){
    mean = std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
    double accum = 0.0;
    std::for_each(vec.begin(), vec.end(), [&](const double d) {
        accum += (d - mean) * (d - mean);
    });
    if (vec.size()>1) std = sqrt(accum / (vec.size()-1));
    else std = -1.;
  }else{
    mean = -1;
    std = -1;
  }
}


DEFINE_ART_MODULE(SelectionAnaTruth)

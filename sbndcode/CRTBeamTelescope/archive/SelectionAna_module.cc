////////////////////////////////////////////////////////////////////////
// Class:       SelectionAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        SelectionAna_module.cc
//
// Generated at Tue Nov 14 06:04:47 2023 by Jiaoyang Li using cetskelgen
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
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTBackTracker.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nurandom/RandomUtils/NuRandomService.h"

#include "TTree.h"
#include "TRandom3.h"
#include <numeric>
#include <algorithm>  // for sort
#include <functional> // for greater

class SelectionAna;

class SelectionAna : public art::EDAnalyzer {
public:
  explicit SelectionAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SelectionAna(SelectionAna const&) = delete;
  SelectionAna(SelectionAna&&) = delete;
  SelectionAna& operator=(SelectionAna const&) = delete;
  SelectionAna& operator=(SelectionAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginSubRun(art::SubRun const& sr) override;
  void respondToOpenInputFile(const art::FileBlock& fb) override;
  void calculateMeanStd(std::vector<double> vec, double &mean, double &std);
  bool canFormSquare(std::vector<double> distances);
  void calculateLenWidthArea(std::vector<std::array<double, 2>> position_x_y, double &length, double &width, double &area);
  void StoreNeutrinoInfo(art::Event const & e, std::string _geant_producer, int geant_track_id, int &nu_pdg, int &ccnc, int &nu_mode, int &nu_interaction_type, double &nu_energy);
  double calculateScalingFactor(double energy, double mass, double boosted_decay_length, double unboosted_decay_length_old, double unboosted_decay_length_new);
  int FindNeutrinoSources(art::Event const & e, std::string _geant_producer, int geant_track_id, int _interactionMode);

private:
  sbnd::CRTBackTracker _crt_back_tracker;

  TRandom3 *fRandomGenerator;
  // Declare member data here.
  std::string _crthit_label;
  std::string _pot_label;
  std::string _mctruth_label;
  std::string _g4_label;

  bool _debug;
  bool _save_input_file_name;
  bool _data_mode; 
  bool _smear_timing; 
  bool _smear_position;
  bool _smear_deposited_energy;
  bool _cosmic_mode;

  int  _interactionMode; // 0 -- signal: dark neutrino; 1 -- background dirt/cosmic/combined (dirt + cosmic)

  TTree* _tree;
  int _run, _subrun, _event;
  std::string _file_name;
  double _weight;                                  ///< Signal: the weight which store as the vertex of dark neutrino.

  double _mct_darkNeutrino_e;                      ///< MCtruth: dark neutrino energy
  double _mct_darkNeutrino_boosted_decay_length;   ///< MCtruth: dark neutrino boosted decay length
  double _mct_darkNeurino_mass;                    ///< MCtruth: dark neutrino mass

  int _n_chits_upstream, _n_chits_downstream;      ///< Number of CRT hits in the upstream/downstream CRTs
  
  // Variables related to geometry for CRT hits.
  int _isSquare;                                    ///< Can CRT hits form a square?
  double _biggest_distance_between_hits;            ///< Mean of the distance between diagonal hits. 
  double _square_side_length;                       ///< Side length of the Square.
  double _square_side_width;                        ///< Side width of the Square.
  double _square_area;                              ///< Area of the Square.

  // Variables realted to deposited energy. 
  std::vector<double> _chit_depositedE;                 ///< CRT hit PEs
  double _chit_depositedE_mean;                         ///< Mean of the CRT hit PE for one event. 
  double _chit_depositedE_sample_std;                   ///< Sample standard deviation of the CRT hit PE for one event.

  double _chit_depositedE_mean_between_hits_biggest_distance;     ///< Mean of the PE diff between any two hits for one event.
  double _chit_depositedE_diff_between_hits_biggest_distance;     ///< Mean of the PE diff between any two hits for one event.
  double _chit_depositedE_ratio_between_hits_biggest_distance;     ///< Mean of the PE ratio between any two hits for one event.
  double _chit_depositedE_frac_diff_between_hits_biggest_distance;     ///< Mean of the PE fractional diff between any two hits for one event.

  // Variables related to t1.
  std::vector<double> _chit_t1;                             ///< CRT hit t1
  std::vector<double> _chit_t1_diff_biggest_distance;       ///< CRT hit t1 diff
  double _t1_diff_between_hits_biggest_distance;       ///< Mean of the t1 between any diagnoal hits for one event.

  // truth info for crt_hit
  std::vector<int>    _chit_backtrack_pdg;              ///< CRT hit, backtracking truth information of the pdg code 
  std::vector<double> _chit_backtrack_energy;           ///< CRT hit, backtracking truth information of the particle energy
  std::vector<double> _chit_backtrack_deposited_energy; ///< CRT hit, backtracking truth information of the deposited energy 
  std::vector<double> _chit_backtrack_purity;           ///< CRT hit, backtracking truth information of selection purity
  std::vector<double> _chit_backtrack_trackID;          ///< CRT hit, backtracking truth information of track id
  std::vector<int>    _chit_backtrack_origin;           ///< CRT hit, 0 -- signal; 1 -- dirt; 2 -- cosmic;

  // truth info for neutrino
  std::vector<int> _chit_backtrack_nu_pdg;
  std::vector<int> _chit_backtrack_nu_ccnc;
  std::vector<int> _chit_backtrack_nu_mode;
  std::vector<int> _chit_backtrack_nu_interaction_type;
  std::vector<double> _chit_backtrack_nu_energy;

  // for smearing to calculate the systematic uncertainty from reconstruction effect. 
  // smear the handles.
  double sigma_pos = 1.88; // sqrt(2)*1.88 cm
  double sigma_timing = 1.88; // ns
  double sigma_depositedE = 1.02; // MeV

  // boundaries for the unphysical values.
  double pos_x_min = -100.419, pos_x_max = 169.296;
  double pos_y_min = -382.778, pos_y_max = -113.064;
  double timing_min = 1.7e6, timing_max = 1.7e6+1600;
  double depositedE_min = 0;

  // POT info. 
  TTree* _sr_tree;
  int _sr_run, _sr_subrun;
  double _sr_begintime, _sr_endtime;
  double _sr_pot; ///< Number of POTs per subrun
  double _sr_spills; ///< Number of Spills per subrun

  // new unboosted decay length
  double _mct_darkNeutrino_new_unboosted_decay_length;
  double _mct_darkNeutrino_unboosted_decay_length=306.53;// cm,
  ULong_t _seednumber;
};


SelectionAna::SelectionAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  _mctruth_label        = p.get<std::string>("MCTruthLabel", "generator");
  _crthit_label         = p.get<std::string>("CRTHitLabel", "crthit");
  _pot_label            = p.get<std::string>("POTLabel", "generator");
  _g4_label             = p.get<std::string>("G4Label", "largeant");
  _crt_back_tracker     = p.get<fhicl::ParameterSet>("CRTBackTracker", fhicl::ParameterSet());
  _interactionMode      = p.get<int>("InteractionMode", 0);
  _debug                = p.get<bool>("Debug", false);
  _smear_timing         = p.get<bool>("SmearTiming", false);
  _smear_position       = p.get<bool>("SmearPosition", false);
  _smear_deposited_energy = p.get<bool>("SmearDepositedEnergy", false);
  _save_input_file_name = p.get<bool>("SaveInputFileName", true);
  _data_mode            = p.get<bool>("DataMode", false);
  _mct_darkNeutrino_new_unboosted_decay_length = p.get<double>("DarkNeutrinoNewUnBoostedDecayLength", 306.53);
  _seednumber           = p.get<ULong_t>("SeedNumber", 14527747);
  _cosmic_mode          = p.get<bool>("CosmicMode", false);

  art::ServiceHandle<art::TFileService> fs;

  _tree = fs->make<TTree>("selection_tree","");
  _tree->Branch("run", &_run, "run/I");
  _tree->Branch("subrun", &_subrun, "subrun/I");
  _tree->Branch("event", &_event, "event/I");
  if (_save_input_file_name) _tree->Branch("file_name", &_file_name);
  if (!_data_mode) {
    _tree->Branch("InteractionMode", &_interactionMode, "InteractionMode/I");
    _tree->Branch("weight", &_weight, "weight/D");
    _tree->Branch("mct_darkNeutrino_e", &_mct_darkNeutrino_e, "mct_darkNeutrino_e/D");
    _tree->Branch("mct_darkNeutrino_boosted_decay_length", &_mct_darkNeutrino_boosted_decay_length, "mct_darkNeutrino_boosted_decay_length/D");
    _tree->Branch("mct_darkNeutrino_mass", &_mct_darkNeurino_mass, "mct_darkNeutrino_mass/D");
  }

  _tree->Branch("n_chits_upstream", &_n_chits_upstream, "n_chits_upstream/I");
  _tree->Branch("n_chits_downstream", &_n_chits_downstream, "n_chits_downstream/I");
  _tree->Branch("isSquare", &_isSquare, "isSquare/I");
  _tree->Branch("square_side_length", &_square_side_length, "square_side_length/D");
  _tree->Branch("square_side_width", &_square_side_width, "square_side_width/D");
  _tree->Branch("square_area", &_square_area, "square_area/D");

  _tree->Branch("biggest_distance_between_hits", &_biggest_distance_between_hits, "biggest_distance_between_hits/D");

  _tree->Branch("chit_depositedE", &_chit_depositedE);
  _tree->Branch("chit_depositedE_mean", &_chit_depositedE_mean, "chit_depositedE_mean/D");
  _tree->Branch("chit_depositedE_sample_std", &_chit_depositedE_sample_std, "chit_depositedE_sample_std/D");

  _tree->Branch("chit_depositedE_mean_between_hits_biggest_distance", &_chit_depositedE_mean_between_hits_biggest_distance, "chit_depositedE_mean_between_hits_biggest_distance/D");
  _tree->Branch("chit_depositedE_ratio_between_hits_biggest_distance", &_chit_depositedE_ratio_between_hits_biggest_distance, "chit_depositedE_ratio_between_hits_biggest_distance/D");
  _tree->Branch("chit_depositedE_frac_diff_between_hits_biggest_distance", &_chit_depositedE_frac_diff_between_hits_biggest_distance, "chit_depositedE_frac_diff_between_hits_biggest_distance/D");

  _tree->Branch("chit_t1", &_chit_t1);
  _tree->Branch("chit_t1_diff_biggest_distance", &_chit_t1_diff_biggest_distance);
  _tree->Branch("t1_diff_between_hits_biggest_distance", &_t1_diff_between_hits_biggest_distance, "t1_diff_between_hits_biggest_distance/D");

  if (!_data_mode){
    _tree->Branch("chit_backtrack_pdg", &_chit_backtrack_pdg);
    _tree->Branch("chit_backtrack_energy", &_chit_backtrack_energy);
    _tree->Branch("chit_backtrack_deposited_energy", &_chit_backtrack_deposited_energy);
    _tree->Branch("chit_backtrack_purity", &_chit_backtrack_purity);
    _tree->Branch("chit_backtrack_trackID", &_chit_backtrack_trackID);
    _tree->Branch("chit_backtrack_origin", &_chit_backtrack_origin);
    _tree->Branch("chit_backtrack_nu_pdg", &_chit_backtrack_nu_pdg);
    _tree->Branch("chit_backtrack_nu_ccnc", &_chit_backtrack_nu_ccnc);
    _tree->Branch("chit_backtrack_nu_mode", &_chit_backtrack_nu_mode);
    _tree->Branch("chit_backtrack_nu_interaction_type", &_chit_backtrack_nu_interaction_type);
    _tree->Branch("chit_backtrack_nu_energy", &_chit_backtrack_nu_energy);
  }

  // POT tree. 
  _sr_tree = fs->make<TTree>("pottree","");
  _sr_tree->Branch("run", &_sr_run, "run/I");
  _sr_tree->Branch("subrun", &_sr_subrun, "subrun/I");
  _sr_tree->Branch("begintime", &_sr_begintime, "begintime/D");
  _sr_tree->Branch("endtime", &_sr_endtime, "endtime/D");
  _sr_tree->Branch("pot", &_sr_pot, "pot/D");
}

void SelectionAna::analyze(art::Event const& e)
{ 
  fRandomGenerator = new TRandom3();     // Initialize the random generator.
  fRandomGenerator->SetSeed(_seednumber); // set seed for the random generator.

  double smear_timing_offset = fRandomGenerator->Gaus(0, sigma_timing);
  double smear_position_offset = fRandomGenerator->Gaus(0, sigma_pos);
  double smear_deposited_energy_offset = fRandomGenerator->Gaus(0, sigma_depositedE);

  if (_smear_timing || _smear_deposited_energy || _smear_position)
    std::cout<<"Random seed: "<<fRandomGenerator->TRandom::GetSeed()<<", _seednumber: "<<_seednumber<<", smear timing: "<<smear_timing_offset<<", smear position: "<<smear_position_offset<<", smear deposited energy: "<<smear_deposited_energy_offset<<std::endl;

  // clear all declared vectors;
  _chit_depositedE.clear();  _chit_t1.clear(); _chit_t1_diff_biggest_distance.clear();
  _chit_backtrack_pdg.clear(); _chit_backtrack_energy.clear(); 
  _chit_backtrack_deposited_energy.clear(); _chit_backtrack_purity.clear(); 
  _chit_backtrack_trackID.clear(); _chit_backtrack_origin.clear();  
  _chit_backtrack_nu_pdg.clear(); _chit_backtrack_nu_ccnc.clear(); 
  _chit_backtrack_nu_mode.clear(); _chit_backtrack_nu_interaction_type.clear(); _chit_backtrack_nu_energy.clear();

  if (!_data_mode) _crt_back_tracker.Initialize(e); // Initialise the backtrack alg. 

  // Get the run, subrun, event information.
  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();

  if(_debug) {
    std::cout << "This is file: " << _file_name << std::endl;
    std::cout << "This is run:subrun:event " << _run << ":" << _subrun << ":" << _event << std::endl;
  }

  // Get the mc truth information for the weight for signal. 
  if(!_data_mode){
    if (_interactionMode==0) {
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
            double _weight_prescale = particle.Vx();
            _mct_darkNeutrino_e = particle.E();
            _mct_darkNeutrino_boosted_decay_length = particle.Vy();
            _mct_darkNeurino_mass = particle.Mass();

            // calculate the scaling factor for the signal.
            double scaling_factor = calculateScalingFactor(_mct_darkNeutrino_e, _mct_darkNeurino_mass, _mct_darkNeutrino_boosted_decay_length, _mct_darkNeutrino_unboosted_decay_length, _mct_darkNeutrino_new_unboosted_decay_length);
            
            // calculate the weight for the signal.
            _weight = _weight_prescale*scaling_factor;
          }
        }
      }
    }else{
        _weight = 1.;
        _mct_darkNeutrino_e = -1.;
        _mct_darkNeutrino_boosted_decay_length = -1.;
        _mct_darkNeurino_mass = -1.;
    }
  }

  // calibration constants to convert from PE to deposited energy.
  std::vector<double> calibration_constant; calibration_constant.resize(2); /*calibration_constant[0], upstream, calibration_constant[1] downstream*/
  calibration_constant[0] = 1.68/(2475.66/40.); // 1.75
  calibration_constant[1] = 1.68/(2239.51/40.);
  if (!_data_mode) {
      calibration_constant[0] = 1.68/(3139.33/40.);
      calibration_constant[1] = 1.68/(3134.69/40.);
  }

  //
  // Get the CRT Hits
  //
  art::Handle<std::vector<sbn::crt::CRTHit>> crt_hit_handle;
  e.getByLabel(_crthit_label, crt_hit_handle);
  if(!crt_hit_handle.isValid()){
    std::cout << "CRTHit product " << _crthit_label << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<sbn::crt::CRTHit>> crt_hit_v;
  art::fill_ptr_vector(crt_hit_v, crt_hit_handle);

  // select within the beam window.
  std::vector<double> beam_time_window = {1.7e6, 1600.}; // 1.7ms, 1600ns
  if (_cosmic_mode) beam_time_window[0] = 1.7e6 - 3000.; // shift the time in the cosmic window. 

  // information to store for the CRT downstream hits.
  std::vector<double> chit_x_downstream, chit_y_downstream, chit_t1_downstream, chit_depoistedE_downstream, chit_corrADC_strip_diff_downstream, chit_t1_strip_diff_downstream;
  chit_x_downstream.clear(); chit_y_downstream.clear(); chit_t1_downstream.clear(); chit_depoistedE_downstream.clear(); chit_corrADC_strip_diff_downstream.clear(); chit_t1_strip_diff_downstream.clear();

  // set the number of hits in the upstream and downstream to 0.
  _n_chits_upstream = 0; _n_chits_downstream = 0;

  // loop over the hits.
  for (size_t ihit = 0; ihit < crt_hit_v.size(); ihit++) {
    
    auto hit = crt_hit_v[ihit];

    double t1 = hit->ts1_ns; // get the time of the hit.
    
    // smear the timing.
    if (_smear_timing)  {
      t1 = t1 + smear_timing_offset;

      if (_debug) std::cout<<"Before hit ts1: "<<hit->ts1_ns<<", smear_timing_offset: "<<smear_timing_offset<<", frac diff: "<<(smear_timing_offset)/hit->ts1_ns<<std::endl;
    }
    // shift t1 inside the beam gate from 1.7e6 for data.
    if (_data_mode) t1 = t1-4420.+1.7e6;

    // select within the beam window.
    if (t1<beam_time_window[0] || t1>beam_time_window[0]+beam_time_window[1]) continue;
    if (_debug) std::cout<<"hit: "<<ihit<<", t1: "<<t1<<", beam_time_window: "<<beam_time_window[0]<<", "<<beam_time_window[1]<<std::endl;
    
    _chit_t1.push_back(t1); // store the t1.

    const std::array<uint16_t,4> corr_adcs = hit->corr_adcs;
    if (hit->tagger == "volTaggerNorth_0") { // downstream
      _n_chits_downstream++; 

      if (_smear_deposited_energy) {
        double chit_depositedE_smeared = hit->peshit*calibration_constant[1] + smear_deposited_energy_offset;
        _chit_depositedE.push_back(chit_depositedE_smeared);
        chit_depoistedE_downstream.push_back(chit_depositedE_smeared);

        if(_debug) std::cout<<"Before hit depositedE: "<<hit->peshit*calibration_constant[1]<<", smear_deposited_energy_offset: "<<smear_deposited_energy_offset<<", frac diff: "<<(smear_deposited_energy_offset)/hit->peshit*calibration_constant[1]<<", sigma_depositedE: "<<sigma_depositedE<<std::endl;
      }else{
        _chit_depositedE.push_back(hit->peshit*calibration_constant[1]);
        chit_depoistedE_downstream.push_back(hit->peshit*calibration_constant[1]);
      }

      chit_x_downstream.push_back(hit->x_pos); chit_y_downstream.push_back(hit->y_pos); 
      chit_t1_downstream.push_back(t1);

      double corrADC_strip_diff = std::abs((corr_adcs[0]+corr_adcs[1])/2. - (corr_adcs[2]+corr_adcs[3])/2.);
      chit_corrADC_strip_diff_downstream.push_back(corrADC_strip_diff);
      chit_t1_strip_diff_downstream.push_back(hit->ts0_ns_corr);

    }else {  // upstream
      _n_chits_upstream++;
      if (_smear_deposited_energy) _chit_depositedE.push_back(hit->peshit*calibration_constant[0] + smear_deposited_energy_offset);
      else _chit_depositedE.push_back(hit->peshit*calibration_constant[0]);
    }

    if(!_data_mode){
      const sbnd::CRTBackTracker::TruthMatchMetrics truthMatch = _crt_back_tracker.TruthMatrixFromTotalEnergy(e, hit);
      _chit_backtrack_pdg.push_back(truthMatch.pdg);
      _chit_backtrack_energy.push_back(truthMatch.particle_energy);
      _chit_backtrack_deposited_energy.push_back(truthMatch.depEnergy_total);
      _chit_backtrack_purity.push_back(truthMatch.purity);
      _chit_backtrack_trackID.push_back(truthMatch.trackid);
      _chit_backtrack_origin.push_back(FindNeutrinoSources(e, _g4_label, truthMatch.trackid, _interactionMode));

      if (_debug) std::cout<<"backtrak pdg: "<<truthMatch.pdg<<", energy: "<<truthMatch.particle_energy<<", depEnergy_total: "<<truthMatch.depEnergy_total<<", purity: "<<truthMatch.purity<<", trackid: "<<truthMatch.trackid<<", origin: "<<FindNeutrinoSources(e, _g4_label, truthMatch.trackid, _interactionMode)<<std::endl;

      // store the neutrino info for the CRT hit.
      int nu_pdg = -1, ccnc = -1, nu_mode = -1, nu_interaction_type = -1; double nu_energy = -1;
      StoreNeutrinoInfo(e, _g4_label, truthMatch.trackid, nu_pdg, ccnc, nu_mode, nu_interaction_type, nu_energy);
      _chit_backtrack_nu_pdg.push_back(nu_pdg);
      _chit_backtrack_nu_ccnc.push_back(ccnc);
      _chit_backtrack_nu_mode.push_back(nu_mode);
      _chit_backtrack_nu_interaction_type.push_back(nu_interaction_type);
      _chit_backtrack_nu_energy.push_back(nu_energy);
    }
  }

  // Calculate the distance between any two hits.
  std::map<double, std::vector<double>>  distance_to_t1_diff_map; distance_to_t1_diff_map.clear();
  std::map<double, std::vector<std::vector<double>>> distance_to_chit_depositedE_mean_map; distance_to_chit_depositedE_mean_map.clear();
  std::map<double, std::vector<double>> distance_to_chit_depositedE_diff_map; distance_to_chit_depositedE_diff_map.clear();
  std::map<double, std::vector<double>> distance_to_chit_depositedE_ratio_map; distance_to_chit_depositedE_ratio_map.clear();
  std::map<double, std::vector<double>> distance_to_chit_depositedE_frac_diff_map; distance_to_chit_depositedE_frac_diff_map.clear();
  std::vector<std::array<double, 2>> position_x_y; position_x_y.clear();
  std::vector<double> distance_between_hits_downstream; distance_between_hits_downstream.clear();
  std::map<double, double>  distance_to_smeared_distance_map; distance_to_smeared_distance_map.clear();

  if (_debug){
    if (chit_x_downstream.size()!=chit_y_downstream.size() || chit_x_downstream.size()!=chit_t1_downstream.size() || chit_x_downstream.size()!=chit_depoistedE_downstream.size()|| chit_x_downstream.size()!=chit_corrADC_strip_diff_downstream.size()|| chit_x_downstream.size()!=chit_t1_strip_diff_downstream.size()){
      std::cout<<"The size of the downstream chits are not the same!"<<std::endl;
    }
  }
  if (chit_x_downstream.size()>=2){
    for (size_t i=0; i<chit_x_downstream.size(); i++){
      for (size_t j=i+1; j<chit_x_downstream.size(); j++){

        double distance = std::hypot(chit_x_downstream[i]-chit_x_downstream[j], chit_y_downstream[i]-chit_y_downstream[j]);
        if (_smear_position) {
          double distance_smeared = std::hypot(chit_x_downstream[i]-chit_x_downstream[j]+smear_position_offset*std::sqrt(2), chit_y_downstream[i]-chit_y_downstream[j]+smear_position_offset*std::sqrt(2));
          distance_to_smeared_distance_map[distance] = distance_smeared;
        }

        double t1_diff_between_hits = std::abs(chit_t1_downstream[i]-chit_t1_downstream[j]);
        if (_smear_timing) {
          t1_diff_between_hits = std::abs((chit_t1_downstream[i]-chit_t1_downstream[j]) + smear_timing_offset*std::sqrt(2));
        }
        
        // deposited E mean. 
        double chit_depositedE_mean_between_hits = (chit_depoistedE_downstream[i]+chit_depoistedE_downstream[j])/2;

        // deposited E diff.
        double chit_depositedE_diff_between_hits = std::abs(chit_depoistedE_downstream[i]-chit_depoistedE_downstream[j]);
        if (_smear_deposited_energy) {
          chit_depositedE_diff_between_hits = std::abs((chit_depoistedE_downstream[i]-chit_depoistedE_downstream[j]) + smear_deposited_energy_offset*std::sqrt(2));
        }

        // deposited E ratio.
        double chit_depositedE_ratio_between_hits = -1; 
        if (chit_depoistedE_downstream[i]<chit_depoistedE_downstream[j]) chit_depositedE_ratio_between_hits = std::abs(chit_depoistedE_downstream[i]/chit_depoistedE_downstream[j]);
        else chit_depositedE_ratio_between_hits = std::abs(chit_depoistedE_downstream[j]/chit_depoistedE_downstream[i]);

        // deposited E fractional diff.
        double chit_depositedE_frac_diff_between_hits = std::abs((chit_depoistedE_downstream[i]-chit_depoistedE_downstream[j]))/((chit_depoistedE_downstream[i]+chit_depoistedE_downstream[j])/2);
        if (_smear_deposited_energy) {
          chit_depositedE_frac_diff_between_hits = std::abs((chit_depoistedE_downstream[i]-chit_depoistedE_downstream[j]) + smear_deposited_energy_offset*std::sqrt(2))/((chit_depoistedE_downstream[i]+chit_depoistedE_downstream[j])/2);   
        }     
        // t1 strip diff 
        double t1_strip_diff_mean = (chit_t1_strip_diff_downstream[i]+chit_t1_strip_diff_downstream[j])/2;

        // corrADC strip diff
        double corrADC_strip_diff_mean = (chit_corrADC_strip_diff_downstream[i]+chit_corrADC_strip_diff_downstream[j])/2;

        distance_between_hits_downstream.push_back(distance);

        distance_to_t1_diff_map[distance].push_back(t1_diff_between_hits);
        distance_to_chit_depositedE_mean_map[distance].push_back({chit_depositedE_mean_between_hits, t1_strip_diff_mean, corrADC_strip_diff_mean});
        distance_to_chit_depositedE_diff_map[distance].push_back(chit_depositedE_diff_between_hits);
        distance_to_chit_depositedE_ratio_map[distance].push_back(chit_depositedE_ratio_between_hits);
        distance_to_chit_depositedE_frac_diff_map[distance].push_back(chit_depositedE_frac_diff_between_hits);
      }
      std::array<double, 2> position = {chit_x_downstream[i], chit_y_downstream[i]};
      position_x_y.push_back(position);
    }
  }

  if (canFormSquare(distance_between_hits_downstream)) {
    _isSquare = 1; // true
    calculateLenWidthArea(position_x_y, _square_side_length, _square_side_width, _square_area);
  }else {
    _isSquare = 0; // false
    _square_side_length = -1.;
    _square_side_width = -1.;
    _square_area = -1.;
  }

  // Calculate the mean and std of the distance between the biggest/smallest distance between any two hits.
  // initialize the variables.
  _biggest_distance_between_hits = -1.;  
  _t1_diff_between_hits_biggest_distance = -1.; 

  _chit_depositedE_mean_between_hits_biggest_distance = -1.;
  _chit_depositedE_diff_between_hits_biggest_distance = -1.;  
  _chit_depositedE_ratio_between_hits_biggest_distance = -1.;  
  _chit_depositedE_frac_diff_between_hits_biggest_distance = -1.;  

  size_t counter=0; //counter_smallest_distance=0, counter_biggest_distance=0;
  for (auto & [distance, chit_depositedE_mean_vec] : distance_to_chit_depositedE_mean_map){
    if (counter==distance_to_chit_depositedE_mean_map.size()-1){ // the maximum distance should belongs to the diagonal hits.  
      // read in the maps.
      auto t1_diff_vec = distance_to_t1_diff_map[distance];
      _chit_t1_diff_biggest_distance = t1_diff_vec;

      auto chit_depositedE_diff_vec = distance_to_chit_depositedE_diff_map[distance];
      auto chit_depositedE_ratio_vec = distance_to_chit_depositedE_ratio_map[distance];
      auto chit_depositedE_frac_diff_vec = distance_to_chit_depositedE_frac_diff_map[distance];

      if (chit_depositedE_mean_vec.size()==2){ // two diagonal pairs. 
        
        _biggest_distance_between_hits = distance;

        if (_smear_position) {
          double distance_smeared = distance_to_smeared_distance_map[distance];
          double scaling_factor = distance_smeared/distance;
          _square_side_length = _square_side_length*scaling_factor;
          _square_side_width = _square_side_width*scaling_factor;
          _square_area = _square_area*scaling_factor*scaling_factor;
        }

        if (_debug) std::cout<<"Before distance: "<<distance<<", smear_position_offset: "<<smear_position_offset<<", frac diff: "<<smear_position_offset/distance<<std::endl;

        // select the pair based on the matrix: t1 strip diff + corrADC strip diff.
        int index = 0; 
        double t1_strip_diff_min = chit_depositedE_mean_vec[0][1];
        if (chit_depositedE_mean_vec[1][1]<t1_strip_diff_min) { index = 1;}
        else if (chit_depositedE_mean_vec[1][1]==t1_strip_diff_min){
          if (chit_depositedE_mean_vec[1][2]<chit_depositedE_mean_vec[0][2]) { index = 1;}
        }
        _t1_diff_between_hits_biggest_distance = t1_diff_vec[index];
        _chit_depositedE_mean_between_hits_biggest_distance = chit_depositedE_mean_vec[index][0];
        _chit_depositedE_diff_between_hits_biggest_distance = chit_depositedE_diff_vec[index];
        _chit_depositedE_ratio_between_hits_biggest_distance = chit_depositedE_ratio_vec[index];
        _chit_depositedE_frac_diff_between_hits_biggest_distance = chit_depositedE_frac_diff_vec[index];
      }
    }
    counter++;
  }
  
  // Calculate the mean and std of the CRT hit PE.
  calculateMeanStd(_chit_depositedE, _chit_depositedE_mean, _chit_depositedE_sample_std);
  _tree->Fill();
}

void SelectionAna::beginSubRun(art::SubRun const& sr) 
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
  std::cout << "POT for this subrun: " << _sr_pot << " (" << _sr_spills << " spills)" << std::endl;
  _sr_tree->Fill();
}

void SelectionAna::respondToOpenInputFile(const art::FileBlock& fb)
{
  _file_name = fb.fileName();
}

void SelectionAna::calculateMeanStd(std::vector<double> vec, double &mean, double &std){
  if (vec.size()>=1){
    mean = std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
    double accum = 0.0;
    std::for_each(vec.begin(), vec.end(), [&](const double d) {
        accum += (d - mean) * (d - mean);
    });
    if (vec.size()>1) std = sqrt(accum / (vec.size()-1));
    else std = -1.;
  }else{
    mean = -1.;
    std = -1.;
  }
}


bool SelectionAna::canFormSquare(std::vector<double> distances) {
  if (distances.size() != 6) {
    return false; // Need exactly 4 points to form a square
  }

  // Sort the distances in non-decreasing order
  std::sort(distances.begin(), distances.end());
   
  // Check if the distances satisfy square properties
  return distances[0] > 0 &&                   // Nonzero distance
    distances[0] == distances[1] &&            // Opposite sides equal
    distances[2] == distances[3] &&
    distances[4] == distances[5] &&            // Diagonals equal
    std::abs(distances[4]-std::sqrt(distances[0]*distances[0] + distances[2]*distances[2]))<0.01;
}

void SelectionAna::calculateLenWidthArea(std::vector<std::array<double, 2>> position_x_y, double &length, double &width, double &area){
  if (position_x_y.size() != 4) { length = -1.; width=-1.; area=-1.; }

  // initialize the variables.
  length = 0.; width = 0.; area = 0.;
  // find the value of length (defined in y direction) and width (defined in x direction)
  for (size_t i=0; i<position_x_y.size(); i++){
    for (size_t j=i+1; j<position_x_y.size(); j++){
        length += std::abs(position_x_y[i][0] - position_x_y[j][0]);
        width += std::abs(position_x_y[i][1] - position_x_y[j][1]);
    }
  }
  length = length/4;
  width = width/4;

  area = length*width;
}


void SelectionAna::StoreNeutrinoInfo(art::Event const & e, std::string _geant_producer, int geant_track_id, int &nu_pdg, int &ccnc, int &nu_mode, int &nu_interaction_type, double &nu_energy){

  lar_pandora::MCTruthToMCParticles truthToParticles;
  lar_pandora::MCParticlesToMCTruth particlesToTruth;

  lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geant_producer, truthToParticles, particlesToTruth);
  bool notFound = true;
  for (auto iter : particlesToTruth) {
    if (iter.first->TrackId() == geant_track_id) {
      simb::Origin_t NuOrigin = iter.second->Origin();
      if (NuOrigin == simb::kBeamNeutrino) { // check if the neutrino is from the beam.
        auto neutrino = iter.second->GetNeutrino();
        auto nu = neutrino.Nu();
        
        nu_pdg = nu.PdgCode();
        ccnc = neutrino.CCNC();
        nu_mode = neutrino.Mode();
        nu_interaction_type = neutrino.InteractionType();
        nu_energy = nu.Momentum().E();

        if (_debug) std::cout<<"Neutrino pdg: "<<nu_pdg<<", nu_mode: "<<nu_mode<<", nu_interaction_type: "<<nu_interaction_type<<", nu_energy: "<<nu_energy<<", outgoing lepton: "<<neutrino.Lepton().PdgCode()<<", HitNuc: "<<neutrino.HitNuc()<<", Target: "<<neutrino.Target()<<", vertex: ("<<neutrino.Nu().Position().X()<<", "<<neutrino.Nu().Position().Y()<<", "<<neutrino.Nu().Position().Z()<<")"<<std::endl;

        notFound = false;
      }
    }
  } 

  if (notFound) {
    nu_pdg = -1;
    ccnc = -1;
    nu_mode = -1;
    nu_interaction_type = -1;
    nu_energy = -1;
  }
}
int SelectionAna::FindNeutrinoSources(art::Event const & e, std::string _geant_producer, int geant_track_id, int _interactionMode) {
  if (_interactionMode==0) return 0; // signal

  lar_pandora::MCTruthToMCParticles truthToParticles;
  lar_pandora::MCParticlesToMCTruth particlesToTruth;

  lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geant_producer, truthToParticles, particlesToTruth);

  for (auto iter : particlesToTruth) {
    if (iter.first->TrackId() == geant_track_id) {
      simb::Origin_t NuOrigin = iter.second->Origin();
      if (NuOrigin == simb::kBeamNeutrino) {
        return 1;
      }else if (NuOrigin == simb::kCosmicRay) {
        return 2;
      }
    }
  }  
  return std::numeric_limits<int>::max();
}

double SelectionAna::calculateScalingFactor(double energy, double mass, double boosted_decay_length, double unboosted_decay_length_old, double unboosted_decay_length_new){
  double gamma = energy/mass;
  double beta = std::sqrt(1-1/(gamma*gamma));
  // exponential function
  double scaling_factor = std::exp(-boosted_decay_length/(gamma*beta*unboosted_decay_length_new))/std::exp(-boosted_decay_length/(gamma*beta*unboosted_decay_length_old));

  return scaling_factor; 
}


DEFINE_ART_MODULE(SelectionAna)
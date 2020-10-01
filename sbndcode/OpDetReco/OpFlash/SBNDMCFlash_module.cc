////////////////////////////////////////////////////////////////////////
// Class:       SBNDMCFlash
// Plugin Type: producer (art v3_04_00)
// File:        SBNDMCFlash_module.cc
//
// Generated at Fri Feb 28 19:56:16 2020 by Marco Del Tutto using cetskelgen
// from cetlib version v3_09_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
// #include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

#include "sbndcode/OpDetReco/OpFlash/FlashFinder/FlashFinderFMWKInterface.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

#include "TVector3.h"
#include "TTree.h"
#include "TH1D.h"
// #include "TRandom3.h"

#include <memory>

class SBNDMCFlash;


class SBNDMCFlash : public art::EDProducer {
public:
  explicit SBNDMCFlash(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SBNDMCFlash(SBNDMCFlash const&) = delete;
  SBNDMCFlash(SBNDMCFlash&&) = delete;
  SBNDMCFlash& operator=(SBNDMCFlash const&) = delete;
  SBNDMCFlash& operator=(SBNDMCFlash&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  std::string _mctruth_label;
  std::string _trigger_label;
  std::string _simphot_label;
  std::string _simphot_insta;
  std::vector<std::string> _pd_to_use;
  int _tpc;
  float _qe_direct, _qe_refl;
  float _window_length, _pre_window;
  bool _debug;

  opdet::sbndPDMapAlg _pds_map;

  // TRandom3 _random;

  int _run, _subrun, _event;
  int _pe_total = 0;
  int _pe_scintillation = 0;
  int _pe_cherenkov = 0;
  std::vector<double> _cherenkov_time_v;
  std::vector<int> _cherenkov_pmt_v;
  std::vector<double> _scintillation_time_v;
  std::vector<int> _scintillation_pmt_v;
  TTree* _tree1;

  void GetFlashLocation(std::vector<double>, double&, double&, double&, double&);

};


SBNDMCFlash::SBNDMCFlash(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  _mctruth_label = p.get<std::string>("MCTruthProduct", "generator");
  _trigger_label = p.get<std::string>("TriggerProduct", "triggersim");
  _simphot_label = p.get<std::string>("SimPhotProduct", "largeant");
  _simphot_insta = p.get<std::string>("SimPhotProductInstance", "");
  _pd_to_use     = p.get<std::vector<std::string>>("PD", _pd_to_use);
  _tpc           = p.get<int>("TPC", -1);
  _qe_direct     = p.get<float>("QEDirect", 0.03);
  _qe_refl       = p.get<float>("QERefl", 0.03);
  _window_length = p.get<float>("WindowLength", 8); // us
  _pre_window    = p.get<float>("PreWindow", 0.1); // us
  _debug         = p.get<bool>("DebugMode",      false);

  art::ServiceHandle<art::TFileService> fs;
  _tree1 = fs->make<TTree>("tree","");
  _tree1->Branch("run",                  &_run,                "run/I");
  _tree1->Branch("subrun",               &_subrun,             "subrun/I");
  _tree1->Branch("event",                &_event,              "event/I");
  _tree1->Branch("pe_total",             &_pe_total,           "pe_total/I");
  _tree1->Branch("pe_scintillation",     &_pe_scintillation,   "pe_scintillation/I");
  _tree1->Branch("pe_cherenkov",         &_pe_cherenkov,       "pe_cherenkov/I");
  _tree1->Branch("cherenkov_time_v",     "std::vector<double>", &_cherenkov_time_v);
  _tree1->Branch("cherenkov_pmt_v",      "std::vector<int>",    &_cherenkov_pmt_v);
  _tree1->Branch("scintillation_time_v", "std::vector<double>", &_scintillation_time_v);
  _tree1->Branch("scintillation_pmt_v",  "std::vector<int>",    &_scintillation_pmt_v);

  produces< std::vector<recob::OpFlash> >();
}

void SBNDMCFlash::produce(art::Event& e)
{

  ::art::ServiceHandle<geo::Geometry> geo;
  auto const *lar_prop = lar::providerFrom<detinfo::LArPropertiesService>();
  auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto const det_prop = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clock_data);

  _qe_direct = _qe_direct / lar_prop->ScintPreScale();
  _qe_refl   = _qe_refl / lar_prop->ScintPreScale();

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();

  // produce OpFlash data-product to be filled within module
  std::unique_ptr< std::vector<recob::OpFlash> > opflashes(new std::vector<recob::OpFlash>);

  if (e.isRealData()) {
    e.put(std::move(opflashes));
    return;
  }

  // art::Handle<std::vector<raw::Trigger> > evt_trigger_h;
  // e.getByLabel(_trigger_label,evt_trigger_h);

  // if( !evt_trigger_h.isValid() || evt_trigger_h->empty() ) {
  //   std::cout << "[NeutrinoMCFlash] Trigger product is not valid or empty." << std::endl;
  //   e.put(std::move(opflashes));
  //   return;
  // }

  art::Handle<std::vector<simb::MCTruth> > evt_mctruth_h;
  e.getByLabel(_mctruth_label,evt_mctruth_h);

  if( !evt_mctruth_h.isValid() || evt_mctruth_h->empty() ) {
    std::cerr << "[NeutrinoMCFlash] MCTruth product is not valid or empty." << std::endl;
    e.put(std::move(opflashes));
    return;
  }

  // art::Handle<std::vector<sim::SimPhotons> > evt_simphot_h;
  // e.getByLabel(_simphot_label,evt_simphot_h);

  // if( !evt_simphot_h.isValid() || evt_simphot_h->empty() ) {
  //   std::cerr << "[NeutrinoMCFlash] SimPhotons product is not valid or empty." << std::endl;
  //   e.put(std::move(opflashes));
  //   return;
  // }


  // if(evt_simphot_h->size() != geo->NOpDets()) {
  //   std::cout << "[NeutrinoMCFlash] Unexpected # of channels in simphotons!" << std::endl;
  //   e.put(std::move(opflashes));
  //   return;
  // }

  // art::Handle<std::vector<sim::SimPhotonsLite> > evt_simphot_h;
  // e.getByLabel(_simphot_label, _simphot_insta, evt_simphot_h);

  // if( !evt_simphot_h.isValid() || evt_simphot_h->empty() ) {
  //   std::cerr << "[NeutrinoMCFlash] SimPhotons product is not valid or empty." << std::endl;
  //   e.put(std::move(opflashes));
  //   return;
  // }

  std::vector<art::Handle<std::vector<sim::SimPhotonsLite> >> evt_simphot_hs;
  e.getManyByType(evt_simphot_hs);
  if( evt_simphot_hs.size() == 0 ) {
    std::cerr << "[NeutrinoMCFlash] SimPhotonsLite product is not valid or empty." << std::endl;
    e.put(std::move(opflashes));
    return;
  }


  // opdet=>opchannel mapping
  std::vector<size_t> opdet2opch(geo->NOpDets(),0);
  for(size_t opch=0; opch<opdet2opch.size(); ++opch){
    opdet2opch[geo->OpDetFromOpChannel(opch)] = opch;
  }

  // Get the OpChannel of the PD to use
  std::vector<int> opch_to_use = lightana::PDNamesToList(_pd_to_use);

  // auto const & evt_trigger = (*evt_trigger_h)[0];
  // auto const trig_time = evt_trigger.TriggerTime();
  auto const trig_time = clock_data.TriggerOffsetTPC();

  if (_debug) std::cout << "trig_time: " << trig_time << std::endl;
  if (_debug) std::cout << "G4ToElecTime(0): " << clock_data.G4ToElecTime(0) << std::endl;
  if (_debug) std::cout << "G4ToElecTime(1000): " << clock_data.G4ToElecTime(1000) << std::endl;
  if (_debug) std::cout << "Number of OpDets: " << geo->NOpDets() << std::endl;

  double nuTime = -1.e9;
  if (_debug) std::cout << "We have " << evt_mctruth_h->size() << " mctruth events." << std::endl;
  for (size_t n = 0; n < evt_mctruth_h->size(); n++) {

    simb::MCTruth const& evt_mctruth = (*evt_mctruth_h)[n];
    if (_debug) std::cout << "Origin: " << evt_mctruth.Origin() << std::endl;
    if (evt_mctruth.Origin() != 1 ) continue;
    if (_debug) std::cout << "We have " << evt_mctruth.NParticles() << " particles." << std::endl;
    for (int p = 0; p < evt_mctruth.NParticles(); p++) {

      simb::MCParticle const& par = evt_mctruth.GetParticle(p);
      //if (par.PdgCode() != 14) continue;
      if (_debug){
        std::cout << "Particle pdg: " << par.PdgCode() << std::endl;
        std::cout << "new Particle time: " << par.T() << std::endl;
        std::cout << "new    converted: " << clock_data.G4ToElecTime(par.T()) - trig_time << std::endl;
        std::cout << std::endl;
      }
      if (std::abs(par.PdgCode()) == 14 || std::abs(par.PdgCode()) == 12) {
        nuTime = par.T();
      }
    }
  }

  if (nuTime == -1.e9) {
    std::cout << "[NeutrinoMCFlash] No neutrino found." << std::endl;
    e.put(std::move(opflashes));
    return;
  }

  std::cout << "[NeutrinoMCFlash] Neutrino G4 interaction time: "  << nuTime << std::endl;

  std::vector<std::vector<double> > pmt_v(1,std::vector<double>(geo->NOpDets(),0));
  _pe_total = _pe_cherenkov = _pe_scintillation = 0;
  _cherenkov_time_v.clear();
  _cherenkov_pmt_v.clear();
  _scintillation_time_v.clear();
  _scintillation_pmt_v.clear();

  // float sampling = time_service->OpticalClock().Frequency() / 1000.0; // GHz
  float start_window = clock_data.TriggerOffsetTPC(); // us
  if (_debug) std::cout << "start_window " << start_window << std::endl;


  // auto const& Photons = *(evt_simphot_h);
  // for (sim::SimPhotonsLite const& photons: Photons) {
  //   unsigned int const nPhotons = std::accumulate(
  //     photons.DetectedPhotons.begin(), photons.DetectedPhotons.end(),
  //     0U, [](auto sum, auto const& entry){ return sum + entry.second; }
  //     );
  //   std::cout << "Channel=" << photons.OpChannel << " has " << nPhotons << " photons (format: [tick] photons):" << std::endl;
  //   for (auto const& pair: photons.DetectedPhotons) {
  //     std::cout << "\t [" << pair.first << "] " << pair.second << std::endl;
  //   }
  // }
  std::cout << "We have " << evt_simphot_hs.size() << " SimPhoton collections" << std::endl;

  float nuTime_elec = clock_data.G4ToElecTime(nuTime) - trig_time;

  for (const art::Handle<std::vector<sim::SimPhotonsLite>> &evt_simphot_h: evt_simphot_hs) {

    bool reflected = (evt_simphot_h.provenance()->productInstanceName() == "Reflected");

    // float qe = _qe_direct;
    // if (reflected) qe = _qe_refl;

    for(sim::SimPhotonsLite const& photons: *(evt_simphot_h)) {

      int opch = photons.OpChannel;


      for(auto const& pair: photons.DetectedPhotons) {

        float photon_time = pair.first - start_window;

        float photon_time_elec = clock_data.G4ToElecTime(photon_time) - trig_time;

        if (_debug && !reflected) {
          std::cout << "Photon time: " << photon_time_elec
                    << " - Neutrino time: " << nuTime_elec << std::endl;
        }


        if (photon_time_elec > nuTime_elec + _window_length) continue;
        if (photon_time_elec > nuTime_elec - _pre_window){

          auto iter = std::find(opch_to_use.begin(), opch_to_use.end(), opch);
          if (iter == opch_to_use.end()) continue;

          auto const& pt = geo->OpDetGeoFromOpChannel(opch).GetCenter();

          if(pt.X() < 0 && _tpc == 1) continue;
          if(pt.X() > 0 && _tpc == 0) continue;

          pmt_v[0][opch] += pair.second;
          _pe_total++;

          // for (int i = 0; i < pair.second; i++) {
          //   float r = _random.Uniform(1.);
          //   if(r < qe) {
          //     pmt_v[0][opch] += 1;
          //     _pe_total ++;
          //   }
          // }
        }
      }
    }
  }

  double Ycenter, Zcenter, Ywidth, Zwidth;
  GetFlashLocation(pmt_v[0], Ycenter, Zcenter, Ywidth, Zwidth);

  recob::OpFlash flash(nuTime_elec,                                // time w.r.t. trigger
                       0,                                          // time width
                       nuTime_elec,                                // flash time in elec clock
                       0.,                                         // frame (?)
                       pmt_v[0],                                   // pe per pmt
                       0, 0, 1,                                    // this are just default values
                       Ycenter, Ywidth, Zcenter, Zwidth);          // flash location

  std::cout << "[NeutrinoMCFlash] MC Flash Time: "  << flash.Time() << std::endl;
  std::cout << "[NeutrinoMCFlash] MC Flash PE:   "  << flash.TotalPE() << std::endl;
  // for (size_t i = 0; i < pmt_v[0].size(); i++) {
  //   std::cout << "ch " << i << " => " << pmt_v[0][i] << std::endl;
  // }

  opflashes->emplace_back(std::move(flash));

  e.put(std::move(opflashes));

  _tree1->Fill();
}

void SBNDMCFlash::GetFlashLocation(std::vector<double> pePerOpChannel,
                                       double& Ycenter,
                                       double& Zcenter,
                                       double& Ywidth,
                                       double& Zwidth)
{

  // Reset variables
  Ycenter = Zcenter = 0.;
  Ywidth  = Zwidth  = -999.;
  double totalPE = 0.;
  double sumy = 0., sumz = 0., sumy2 = 0., sumz2 = 0.;

  for (unsigned int opch = 0; opch < pePerOpChannel.size(); opch++) {

    // Get physical detector location for this opChannel
    double PMTxyz[3];
    ::art::ServiceHandle<geo::Geometry> geo;
    geo->OpDetGeoFromOpChannel(opch).GetCenter(PMTxyz);

    // Add up the position, weighting with PEs
    sumy    += pePerOpChannel[opch]*PMTxyz[1];
    sumy2   += pePerOpChannel[opch]*PMTxyz[1]*PMTxyz[1];
    sumz    += pePerOpChannel[opch]*PMTxyz[2];
    sumz2   += pePerOpChannel[opch]*PMTxyz[2]*PMTxyz[2];

    totalPE += pePerOpChannel[opch];
  }

  Ycenter = sumy/totalPE;
  Zcenter = sumz/totalPE;

  // This is just sqrt(<x^2> - <x>^2)
  if ( (sumy2*totalPE - sumy*sumy) > 0. )
    Ywidth = std::sqrt(sumy2*totalPE - sumy*sumy)/totalPE;

  if ( (sumz2*totalPE - sumz*sumz) > 0. )
    Zwidth = std::sqrt(sumz2*totalPE - sumz*sumz)/totalPE;
}

DEFINE_ART_MODULE(SBNDMCFlash)

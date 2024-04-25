////////////////////////////////////////////////////////////////////////
// Class:       CRTPrintTruth
// Plugin Type: analyzer (Unknown Unknown)
// File:        CRTPrintTruth_module.cc
//
// Generated at Wed Feb 21 10:38:50 2024 by Jiaoyang Li using cetskelgen
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
#include "lardataobj/Simulation/AuxDetHit.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbndcode/CRT/CRTUtils/CRTBackTracker.h"

class CRTPrintTruth;


class CRTPrintTruth : public art::EDAnalyzer {
public:
  explicit CRTPrintTruth(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTPrintTruth(CRTPrintTruth const&) = delete;
  CRTPrintTruth(CRTPrintTruth&&) = delete;
  CRTPrintTruth& operator=(CRTPrintTruth const&) = delete;
  CRTPrintTruth& operator=(CRTPrintTruth&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  art::Ptr<simb::MCParticle> getParticleByID(
          std::vector<art::Ptr<simb::MCParticle>> mclistLARG4,
          int TrackId);

private:

  // Declare member data here.
  sbnd::CRTBackTracker _crt_back_tracker;
  std::string _crthit_label;
  std::string _g4_label;
  std::string _auxdethit_label;

};


CRTPrintTruth::CRTPrintTruth(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  _g4_label          = p.get<std::string>("G4Label", "largeant");
  _crthit_label      = p.get<std::string>("CRTHitLabel", "crthit");
  _auxdethit_label   = p.get<std::string>("AuxDetHitLabel", "largeant:LArG4DetectorServicevolAuxDetSensitiveCRTStripBERN");
  _crt_back_tracker  = p.get<fhicl::ParameterSet>("CRTBackTracker", fhicl::ParameterSet());
}

void CRTPrintTruth::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  _crt_back_tracker.Initialize(e); // Initialise the backtrack alg. 

  // print out the event, run, and subrun numbers
  std::cout<<"PrintTruth This is run: "<<e.id().run()<<", SubRun: "<<e.id().subRun()<<", event: "<<e.id().event()<<std::endl;

  //
  // Get the MCParticles from G4
  //
  art::Handle<std::vector<simb::MCParticle>> mcp_h;
  std::vector<art::Ptr<simb::MCParticle>> mcp_v;

  e.getByLabel(_g4_label, mcp_h);
  if(!mcp_h.isValid()){
    std::cout << "MCTruth product " << _g4_label << " not found..." << std::endl;
    throw std::exception();
  }
  art::fill_ptr_vector(mcp_v, mcp_h);


  //
  // Get the AuxDetHits from G4
  //
  art::Handle<std::vector<sim::AuxDetHit>> adh_h; // geant4
  std::vector<art::Ptr<sim::AuxDetHit>> adh_v;
  e.getByLabel(_auxdethit_label, adh_h);
  if(!adh_h.isValid()){
    std::cout << "AuxDetHit product " << _auxdethit_label << " not found..." << std::endl;
    throw std::exception();
  }
  art::fill_ptr_vector(adh_v, adh_h);

  //
  // Get the CRT Hits
  //
  art::Handle<std::vector<sbn::crt::CRTHit>> crt_hit_h;
  e.getByLabel(_crthit_label, crt_hit_h);
  if(!crt_hit_h.isValid()){
    std::cout << "CRTHit product " << _crthit_label << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<sbn::crt::CRTHit>> crt_hit_v;
  art::fill_ptr_vector(crt_hit_v, crt_hit_h);


  // PRINT OUT ALL MCPARTICLES
  std::map<int, double> mcptrackID_to_startx_map;
  std::map<int, double> mcptrackID_to_starty_map;
  std::map<int, double> mcptrackID_to_startz_map;
  //std::map<int, double> mcptrackID_to_energy_map;
  std::map<int, int> mcptrackID_to_pdg_map;  
  std::map<int, int> mcptrackID_to_motherID_map;  
  std::map<int, double> mcptrackID_to_momentum_map;
  std::map<int, std::string> mcptrackID_to_start_process_map;
  std::map<int, std::string> mcptrackID_to_end_process_map;
  for (auto const& mcp : mcp_v) {
    if (mcp->StatusCode() != 1) continue; 
    //std::cout<<"MCParticle: "<<mcp->TrackId()<<", pdg: "<<mcp->PdgCode()<<", E: "<<mcp->E()<<", Process: "<<mcp->Process()<<", mother: "<<mcp->Mother()/*<<", has "<<mcp->NumberDaughters()<<" daughters. "*/<<std::endl;
    mcptrackID_to_pdg_map[mcp->TrackId()] = mcp->PdgCode();
    mcptrackID_to_motherID_map[mcp->TrackId()] = mcp->Mother();
    mcptrackID_to_startx_map[mcp->TrackId()] = mcp->Vx();
    mcptrackID_to_starty_map[mcp->TrackId()] = mcp->Vy();
    mcptrackID_to_startz_map[mcp->TrackId()] = mcp->Vz();
    //mcptrackID_to_energy_map[mcp->TrackId()] = mcp->E();
    mcptrackID_to_momentum_map[mcp->TrackId()] = mcp->P()*1000.;
    mcptrackID_to_start_process_map[mcp->TrackId()] = mcp->Process();
    mcptrackID_to_end_process_map[mcp->TrackId()] = mcp->EndProcess();
  }

  int counter = 0;
  for (auto const& adh : adh_v) {
    std::cout<<"PrintTruth AuxDetHit "<<counter<<": pdg: "<<mcptrackID_to_pdg_map[adh->GetTrackID()]<<", P: "<<mcptrackID_to_momentum_map[adh->GetTrackID()]<<" MeV, Process: "<<mcptrackID_to_start_process_map[adh->GetTrackID()]<<", hit postion: ("<<0.5*(adh->GetEntryX()+adh->GetExitX())<<", "<<0.5*(adh->GetEntryY()+adh->GetExitY())<<", "<<0.5*(adh->GetEntryZ()+adh->GetExitZ())<<"); track id: "<<adh->GetTrackID();

    std::string process = mcptrackID_to_start_process_map[adh->GetTrackID()];
    int thistrackid = adh->GetTrackID();
    while (process!="primary" && mcptrackID_to_pdg_map[thistrackid]!=0) {
      int motherid = mcptrackID_to_motherID_map[thistrackid];
      process = mcptrackID_to_start_process_map[motherid];
      thistrackid = motherid;
      int mother_pdg = mcptrackID_to_pdg_map[motherid];
      std::cout<<", which has mother pdg of: "<<mother_pdg<<", process: "<<process;
    }
    std::cout<<std::endl;
    counter++;
  }

  size_t n_hits = crt_hit_v.size();
  for (size_t i = 0; i < n_hits; i++) {
    auto hit = crt_hit_v[i];
    const sbnd::CRTBackTracker::TruthMatchMetrics truthMatch = _crt_back_tracker.TruthMatrixFromTotalEnergy(e, hit);

    std::cout<<"PrintTruth CRTHit "<<i/*<<", time: "<<hit->ts1_ns-1.7e6*/<<"; hit postion: ("<<hit->x_pos<<", "<<hit->y_pos<<", "<<hit->z_pos<<"); pdg: "<<truthMatch.pdg<<"; deposited energy: "<<truthMatch.depEnergy_total*1000.<<" MeV; purity: "<<truthMatch.purity<<"; trackID: "<<truthMatch.trackid<<std::endl;
    
    /*<<", mother: "<<getParticleByID(mcp_v, truthMatch.trackid)->Mother()<<", "<<getParticleByID(mcp_v, getParticleByID(mcp_v, truthMatch.trackid)->Mother())->PdgCode()*/

    /*lar_pandora::MCTruthToMCParticles truthToParticles;
    lar_pandora::MCParticlesToMCTruth particlesToTruth;

    lar_pandora::LArPandoraHelper::CollectMCParticles(e, _g4_label, truthToParticles, particlesToTruth);
    for (auto iter : particlesToTruth) {
      if (iter.first->TrackId() == truthMatch.trackid) {
        art::Ptr<simb::MCTruth> mcNeutrino = iter.second;
        int mcdaughter = mcNeutrino->NParticles();
        std::cout<<"Neutrino sources: "<<mcNeutrino->Origin()<<std::endl;
        
        for (int idaughter=0; idaughter<mcdaughter; ++idaughter) {
          const simb::MCParticle mcpart = mcNeutrino->GetParticle(idaughter);
          if (mcpart.StatusCode() != 1) continue; // Exclude particles that are not propagated
          
          std::cout<<mcpart.TrackId()<<", pdg: "<<mcpart.PdgCode()<<", E: "<<mcpart.E()<<", T: "<<mcpart.T()<<", P: "<<mcpart.P()<<", Process: "<<mcpart.Process()<<", has "<<mcpart.NumberDaughters()<<" daughters. "<<std::endl;

          // acess with the track id
          art::Ptr<simb::MCParticle> mcpart2 = getParticleByID(mcp_v, mcpart.TrackId());

          std::cout<<"2 "<<mcpart2->TrackId()<<", pdg: "<<mcpart2->PdgCode()<<", E: "<<mcpart2->E()<<", T: "<<mcpart2->T()<<", P: "<<mcpart2->P()<<", Process: "<<mcpart2->Process()<<", has "<<mcpart2->NumberDaughters()<<" daughters. "<<std::endl;
          
          // loop through all daughter particles
          for (int igranddaughter=0; igranddaughter<mcpart.NumberDaughters(); ++igranddaughter) {
            int mcgranddaughter_trackid = mcpart.Daughter(igranddaughter);
            std::cout<<mcgranddaughter_trackid<<", mcgranddaughter_trackid"<<std::endl;
            auto mcgranddaughter = getParticleByID(mcp_v, mcgranddaughter_trackid);
            if (mcgranddaughter->StatusCode() != 1) continue; // Exclude particles that are not propagated
            // print out the daughter particles
            std::cout<<"  "<<igranddaughter<<", pdg: "<<mcgranddaughter->PdgCode()<<", E: "<<mcgranddaughter->E()<<", T: "<<mcgranddaughter->T()<<", P: "<<mcgranddaughter->P()<<", Process: "<<mcgranddaughter->Process()<<std::endl;
          }
        }
      }
    }*/  
  }
  std::cout<<"PrintTruth Hit =================================================================="<<std::endl;
}

art::Ptr<simb::MCParticle> CRTPrintTruth::getParticleByID(
          std::vector<art::Ptr<simb::MCParticle>> mclistLARG4,
          int TrackId) 
{
  for(auto particle : mclistLARG4) {
    if ( particle->TrackId() == TrackId){
      return particle;
    }
  }
  return art::Ptr<simb::MCParticle>();
}


DEFINE_ART_MODULE(CRTPrintTruth)

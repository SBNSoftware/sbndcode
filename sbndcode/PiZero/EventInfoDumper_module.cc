////////////////////////////////////////////////////////////////////////
// Class:       EventInfoDumper
// Plugin Type: analyzer
// File:        EventInfoDumper_module.cc
// Author:      Henry Lay (h.lay@lancaster.ac.uk)
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
#include "canvas/Persistency/Common/FindOneP.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"

namespace sbnd {
  class EventInfoDumper;
}

class sbnd::EventInfoDumper : public art::EDAnalyzer {
public:
  explicit EventInfoDumper(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  EventInfoDumper(EventInfoDumper const&) = delete;
  EventInfoDumper(EventInfoDumper&&) = delete;
  EventInfoDumper& operator=(EventInfoDumper const&) = delete;
  EventInfoDumper& operator=(EventInfoDumper&&) = delete;

  // Required functions.
  void analyze(const art::Event &e) override;

  bool SignalEvent(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles);

  void DumpEvent(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles);

  bool VolumeCheck(const TVector3 &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);

private:

  art::ServiceHandle<cheat::ParticleInventoryService> particleInv;
  art::ServiceHandle<cheat::BackTrackerService>       backTracker;

  std::string fMCParticleModuleLabel;

  bool fAVOnly;

  int  _run;
  int  _subrun;
  int  _event;
  bool _signal;
};

sbnd::EventInfoDumper::EventInfoDumper(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fMCParticleModuleLabel(p.get<std::string>("MCParticleModuleLabel", "largeant"))
  , fAVOnly(p.get<bool>("AVOnly", true))
  {}

void sbnd::EventInfoDumper::analyze(const art::Event &e)
{
  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  =  e.id().event();

  std::cout << "This is event " << _run << "-" << _subrun << "-" << _event << std::endl;

  // Get MCTruths
  std::vector<art::Handle<std::vector<simb::MCTruth>>> MCTruthHandles = e.getMany<std::vector<simb::MCTruth>>();
  _signal = SignalEvent(e, MCTruthHandles);

  if(_signal)
    std::cout << "Is signal event" << std::endl;
  else
    std::cout << "Is not signal event" << std::endl;

  DumpEvent(e, MCTruthHandles);
  
  std::cout << std::endl;
}

bool sbnd::EventInfoDumper::SignalEvent(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles)
{
  std::vector<std::pair<bool, double>> events;

  for(auto const& MCTruthHandle : MCTruthHandles)
    {
      std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
      art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

      art::FindManyP<simb::MCParticle> MCTruthToMCParticles(MCTruthHandle, e, fMCParticleModuleLabel);

      for(auto const& mct : MCTruthVec)
        {
          if(mct->Origin() != 1)
            continue;

          const simb::MCNeutrino mcn = mct->GetNeutrino();
          const simb::MCParticle nu  = mcn.Nu();

          const bool nc = mcn.CCNC() == 1;
          const bool fv = VolumeCheck(nu.Position().Vect(), 20., 5., 10., 50.);

          unsigned pizeros = 0;

          for(int i = 0; i < mct->NParticles(); ++i)
            {
              const auto mcp = mct->GetParticle(i);

              if(mcp.PdgCode() == 111 && mcp.StatusCode() != 1)
                ++pizeros;
            }

          const bool pizero = pizeros > 0;

          const std::vector<art::Ptr<simb::MCParticle>> MCParticleVec = MCTruthToMCParticles.at(mct.key());
          double total_en = 0.;

          for(auto const& mcp : MCParticleVec)
            {
              std::vector<const sim::IDE*> ides = backTracker->TrackIdToSimIDEs_Ps(mcp->TrackId());

              for(auto const& ide : ides)
                total_en += ide->energy / 1000.;
            }

          events.push_back({nc && fv && pizero, total_en});
        }
    }

  if(events.size() == 0)
    return false;

  std::sort(events.begin(), events.end(),
            [](const auto &a, const auto &b)
            { return a.second > b.second; });

  return events.at(0).first;
}

void sbnd::EventInfoDumper::DumpEvent(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles)
{
  for(auto const& MCTruthHandle : MCTruthHandles)
    {
      std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
      art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

      art::FindManyP<simb::MCParticle> MCTruthToMCParticles(MCTruthHandle, e, fMCParticleModuleLabel);

      for(auto const& mct : MCTruthVec)
        {
          if(mct->Origin() != 1)
            continue;

          const simb::MCNeutrino mcn = mct->GetNeutrino();
          const simb::MCParticle nu  = mcn.Nu();

          if(fAVOnly && !VolumeCheck(nu.Position().Vect()))
            continue;

          const bool nc = mcn.CCNC() == 1;

          std::cout << '\n'
                    << "===== Neutrino Event =====\n"
                    << "PDG: " << nu.PdgCode() << '\n'
                    << "CCNC: " << (nc ? "NC" : "CC") << '\n'
                    << "E: " << nu.E() << '\n'
                    << "Vtx: (" << nu.Vx() << ", " << nu.Vy() << ", " << nu.Vz() << ") cm\n";

          const std::vector<art::Ptr<simb::MCParticle>> MCParticleVec = MCTruthToMCParticles.at(mct.key());

          int i = 0;

          for(auto const& mcp : MCParticleVec)
            {
              if(mcp->StatusCode() != 1 || mcp->Mother() != nu.TrackId() + 10000000)
                continue;

              std::cout << "\tParticle " << i << '\n'
                        << "\t\tPDG: " << mcp->PdgCode() << '\n'
                        << "\t\tStatusCode: " << mcp->StatusCode() << '\n'
                        << "\t\tTrackID: " << mcp->TrackId() << '\n'
                        << "\t\tMother: " << mcp->Mother() << '\n'
                        << "\t\tNumDaughters: " << mcp->NumberDaughters() << '\n'
                        << "\t\tLength: " << (mcp->EndPosition().Vect() - mcp->Position().Vect()).Mag() << '\n'
                        << "\t\tE: " << mcp->E() << '\n'
                        << "\t\tp: (" << mcp->Px() << ", " << mcp->Py() << ", " << mcp->Pz() << ") cm\n";

              if(mcp->PdgCode() == 111 && mcp->NumberDaughters() == 2)
                {
                  const simb::MCParticle* gamma1 = particleInv->TrackIdToParticle_P(mcp->Daughter(0));
                  const simb::MCParticle* gamma2 = particleInv->TrackIdToParticle_P(mcp->Daughter(1));

                  const TVector3 gamma1_mom = gamma1->Momentum().Vect();
                  const TVector3 gamma2_mom = gamma2->Momentum().Vect();
                  const double open_angle   = TMath::RadToDeg() * gamma1_mom.Angle(gamma2_mom);

                  std::cout << "\t\t\tPDG: " << gamma1->PdgCode() << '\n'
                            << "\t\t\tStatusCode: " << gamma1->StatusCode() << '\n'
                            << "\t\t\tE: " << gamma1->E() << '\n'
                            << "\t\t\tp: (" << gamma1->Px() << ", " << gamma1->Py() << ", " << gamma1->Pz() << ") cm\n"
                            << "\t\t\tLength: " << (gamma1->EndPosition().Vect() - gamma1->Position().Vect()).Mag() << '\n'
                            << "\t\t\tPDG: " << gamma2->PdgCode() << '\n'
                            << "\t\t\tStatusCode: " << gamma2->StatusCode() << '\n'
                            << "\t\t\tE: " << gamma2->E() << '\n'
                            << "\t\t\tp: (" << gamma2->Px() << ", " << gamma2->Py() << ", " << gamma2->Pz() << ") cm\n"
                            << "\t\t\tLength: " << (gamma2->EndPosition().Vect() - gamma2->Position().Vect()).Mag() << '\n'
                            << "\t\tOpen Angle: " << open_angle << '\n';
                }
              ++i;
            }
          std::cout << std::endl;
        }
    }
}

bool sbnd::EventInfoDumper::VolumeCheck(const TVector3 &pos, const double &walls, const double &cath, const double &front, const double &back)
{
  const bool xedges = pos.X() < (200. - walls) && pos.X() > (-200. + walls);
  const bool yedges = pos.Y() < (200. - walls) && pos.Y() > (-200. + walls);
  const bool zedges = pos.Z() < (500. - back)  && pos.Z() > (0. + front);
  const bool caths  = pos.X() > cath || pos.X() < -cath;

  return xedges && yedges && zedges && caths;
}

DEFINE_ART_MODULE(sbnd::EventInfoDumper)

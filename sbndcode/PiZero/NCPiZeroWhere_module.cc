////////////////////////////////////////////////////////////////////////
// Class:       NCPiZeroWhere
// Plugin Type: analyzer
// File:        NCPiZeroWhere_module.cc
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
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

#include "NCPiZeroStructs.h"

namespace sbnd {
  class NCPiZeroWhere;
}

class sbnd::NCPiZeroWhere : public art::EDAnalyzer {
public:
  explicit NCPiZeroWhere(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NCPiZeroWhere(NCPiZeroWhere const&) = delete;
  NCPiZeroWhere(NCPiZeroWhere&&) = delete;
  NCPiZeroWhere& operator=(NCPiZeroWhere const&) = delete;
  NCPiZeroWhere& operator=(NCPiZeroWhere&&) = delete;

  // Required functions.
  void analyze(const art::Event &e) override;

  void AnalyseNeutrinos(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles);

  void AnalyseMCTruth(const art::Event &e, const art::Ptr<simb::MCTruth> &mct);

  bool VolumeCheck(const geo::Point_t &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);
  bool VolumeCheck(const TVector3 &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);

private:

  art::InputTag fMCParticleModuleLabel;
};

sbnd::NCPiZeroWhere::NCPiZeroWhere(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fMCParticleModuleLabel(p.get<art::InputTag>("MCParticleModuleLabel"))
  {
  }

void sbnd::NCPiZeroWhere::analyze(const art::Event &e)
{
  // Get MCTruths
  std::vector<art::Handle<std::vector<simb::MCTruth>>> MCTruthHandles = e.getMany<std::vector<simb::MCTruth>>();

  AnalyseNeutrinos(e, MCTruthHandles);
}

void sbnd::NCPiZeroWhere::AnalyseNeutrinos(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles)
{
  for(auto const& MCTruthHandle : MCTruthHandles)
    {
      std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
      art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

      for(auto const& mct : MCTruthVec)
        {
          if(mct->Origin() != 1)
            continue;

          AnalyseMCTruth(e, mct);
        }
    }
}

void sbnd::NCPiZeroWhere::AnalyseMCTruth(const art::Event &e, const art::Ptr<simb::MCTruth> &mct)
{
  const simb::MCNeutrino mcn = mct->GetNeutrino();
  const simb::MCParticle nu  = mcn.Nu();

  const bool nc = mcn.CCNC() == 1;
  const bool fv = VolumeCheck(nu.Position().Vect(), 20., 5., 10., 50.);

  art::FindManyP<simb::MCParticle> MCTruthToMCParticles( { mct }, e, fMCParticleModuleLabel);
  const std::vector<art::Ptr<simb::MCParticle>> MCParticleVec = MCTruthToMCParticles.at(0);

  int protons = 0, charged_pions = 0, neutral_pions = 0;

  for(auto const& mcp : MCParticleVec)
    {
      if(mcp->Process() == "primary" && mcp->StatusCode() == 1)
        {
          switch(abs(mcp->PdgCode()))
            {
            case 2212:
              if(mcp->P() > .4)
                ++protons;
              break;
            case 211:
              if(mcp->P() > .15)
                ++charged_pions;
              break;
            case 111:
              if(mcp->NumberDaughters() == 2)
                ++neutral_pions;
              break;
            default:
              break;
            }
        }
    }

  const bool pizero = neutral_pions == 1;

  if(nc && fv && pizero)
    {
      const bool av = VolumeCheck(nu.Position().Vect());

      bool pizero2 = false;

      for(int i = 0; i < mct->NParticles(); ++i)
        {
          const simb::MCParticle mcp = mct->GetParticle(i);

          if(mcp.PdgCode() == 111)
            pizero2 = true;
        }

      if(!(nc && av && pizero2))
        {
          std::cout << "WHAAAAAAAAAAAATTTTTT" << std::endl;
          std::cout << e.run() << "-" << e.subRun() << " " << e.event() << std::endl;
          throw std::exception();
        }
      else
        std::cout << "Good signal" << std::endl;
    }
}

bool sbnd::NCPiZeroWhere::VolumeCheck(const geo::Point_t &pos, const double &walls, const double &cath, const double &front, const double &back)
{
  const TVector3 posVec(pos.X(), pos.Y(), pos.Z());
  return VolumeCheck(posVec, walls, cath, front, back);
}

bool sbnd::NCPiZeroWhere::VolumeCheck(const TVector3 &pos, const double &walls, const double &cath, const double &front, const double &back)
{
  const bool xedges = pos.X() < (200. - walls) && pos.X() > (-200. + walls);
  const bool yedges = pos.Y() < (200. - walls) && pos.Y() > (-200. + walls);
  const bool zedges = pos.Z() < (500. - back)  && pos.Z() > (0. + front);
  const bool caths  = pos.X() > cath || pos.X() < -cath;

  return xedges && yedges && zedges && caths;
}

DEFINE_ART_MODULE(sbnd::NCPiZeroWhere)

////////////////////////////////////////////////////////////////////////
// Class:       NCPiZeroFilter
// Plugin Type: filter (Unknown Unknown)
// File:        NCPiZeroFilter_module.cc
//
// Generated at Fri Sep 16 13:35:13 2022 by Henry Lay using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCTruth.h"

#include <memory>

namespace sbnd {
  class NCPiZeroFilter;
}


class sbnd::NCPiZeroFilter : public art::EDFilter {
public:
  explicit NCPiZeroFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NCPiZeroFilter(NCPiZeroFilter const&) = delete;
  NCPiZeroFilter(NCPiZeroFilter&&) = delete;
  NCPiZeroFilter& operator=(NCPiZeroFilter const&) = delete;
  NCPiZeroFilter& operator=(NCPiZeroFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

private:

  // Declare member data here.
  bool fVerbose;

  std::string fNuGenModuleLabel;

  float fXEdgeCut, fYEdgeCut, fZFrontCut, fZBackCut, fCathodeXCut;

  float fXOutEdge, fYEdge, fZFrontEdge, fZBackEdge, fXCathodeEdge;
};


sbnd::NCPiZeroFilter::NCPiZeroFilter(fhicl::ParameterSet const& pset)
  : EDFilter{pset},
  fVerbose (pset.get<bool>("Verbose", false)),
  fNuGenModuleLabel (pset.get<std::string>("NuGenModuleLabel")),
  fXEdgeCut (pset.get<float>("XEdgeCut",15.f)),
  fYEdgeCut (pset.get<float>("YEdgeCut",15.f)),
  fZFrontCut (pset.get<float>("ZFrontCut",30.f)),
  fZBackCut (pset.get<float>("ZBackCut",65.f)),
  fCathodeXCut (pset.get<float>("CathodeXCut",1.5f))
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

bool sbnd::NCPiZeroFilter::filter(art::Event& e)
{
  art::Handle<std::vector<simb::MCTruth> > handleNeutrinos;
  e.getByLabel(fNuGenModuleLabel, handleNeutrinos);

  if(fVerbose) std::cout << "===============================================================\n"
			 << "Run: " << e.run() << " SubRun: " << e.subRun() << " Event: "
			 << e.event() << '\n'
			 << "===============================================================\n"
			 << '\n'
			 << "Contains " << handleNeutrinos->size() << " neutrinos\n"
			 << std::endl;
  
  std::vector<std::pair<int, double> > neutrino_list;

  for(unsigned int nu_i = 0; nu_i < handleNeutrinos->size(); ++nu_i)
    {
      const art::Ptr<simb::MCTruth> truthNeutrino(handleNeutrinos, nu_i);
      if(truthNeutrino.isNull()) continue;

      const double nu_en = truthNeutrino->GetNeutrino().Nu().E();
      
      neutrino_list.push_back(std::pair<int, double>(nu_i, nu_en));
    }

  std::sort(neutrino_list.begin(), neutrino_list.end(), [](const std::pair<int, double> &x,
							   const std::pair<int, double> &y)
	    { return x.second > y.second; });

  const art::Ptr<simb::MCTruth> leadNeutrino(handleNeutrinos, neutrino_list[0].first);
  const simb::MCNeutrino nu = leadNeutrino->GetNeutrino();
  unsigned int charged_leptons(0), pizeros(0), pions(0);
  
  for(int daughter_i = 0; daughter_i < leadNeutrino->NParticles(); ++daughter_i)
    {
      const simb::MCParticle daughter = leadNeutrino->GetParticle(daughter_i);
      const int pdg                   = std::abs(daughter.PdgCode());
      
      if(pdg == 11 || pdg == 13 || pdg == 15)
	++charged_leptons;
      else if(pdg == 111)
	++pizeros;
      else if(pdg == 211)
	++pions;
    }
  
  if(fVerbose) std::cout << "Leading Neutrino\n"
			 << "\tCCNC:            " << nu.CCNC() << '\n'
			 << "\tCharged Leptons: " << charged_leptons << '\n'
			 << "\tPiZeros:         " << pizeros << '\n'
			 << "\tCharged Pions:   " << pions << '\n';
  
  if(nu.CCNC() != 1)
    {
      if(fVerbose) std::cout << "\tPasses?          No\n" << std::endl;
      return false;
    }
  else if(charged_leptons != 0)
    {
      std::cout << "---------------------------------------------------------------\n"
		<< "ERROR: NC events contains " << charged_leptons << " charged leptons?"
		<< "---------------------------------------------------------------\n" 
		<< std::endl;
      if(fVerbose) std::cout << "\tPasses?          No\n" << std::endl;
      return false;
    }
  else if(pizeros == 1 && pions == 0) 
    {
      if(fVerbose) std::cout << "\tPasses?          Yes\n" << std::endl;
      return true;
    }
  else
    {
      if(fVerbose) std::cout << "\tPasses?          No\n" << std::endl;
      return false;
    }
}

DEFINE_ART_MODULE(sbnd::NCPiZeroFilter)

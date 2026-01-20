////////////////////////////////////////////////////////////////////////
// Class:       DeuteronFilter
// Plugin Type: filter (Unknown Unknown)
// File:        DeuteronFilter_module.cc
//
// Generated Mon Jan 5th 2026 by Anna Beever
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

#include <memory>

// Additional framework includes
#include "art_root_io/TFileService.h"

// Additional LArSoft includes

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

namespace nuclearFragments {
  class DeuteronFilter;
}


class nuclearFragments::DeuteronFilter : public art::EDFilter {
public:
  explicit DeuteronFilter(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DeuteronFilter(DeuteronFilter const&) = delete;
  DeuteronFilter(DeuteronFilter&&) = delete;
  DeuteronFilter& operator=(DeuteronFilter const&) = delete;
  DeuteronFilter& operator=(DeuteronFilter&&) = delete;

  // Required functions.
  bool filter(art::Event& e) override;

  bool VolumeCheck(const TVector3 &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);

private:

  // Define input labels
  std::string fNuGenLabel;

};


nuclearFragments::DeuteronFilter::DeuteronFilter(fhicl::ParameterSet const& p)
  : EDFilter{p},
  // More initializers here.
  fNuGenLabel(p.get<std::string>("NuGenLabel"))
{
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

bool nuclearFragments::DeuteronFilter::filter(art::Event& e)
{ 

  double minP = 0.4959;

  //Truth level
  art::ValidHandle<std::vector<simb::MCTruth>> truNuHandle = e.getValidHandle<std::vector<simb::MCTruth>>(fNuGenLabel);
  std::vector<art::Ptr<simb::MCTruth>> truNuVector;

  if(truNuHandle.isValid())
  {
    art::fill_ptr_vector(truNuVector,truNuHandle);
  }

  for(const art::Ptr<simb::MCTruth> &truNu : truNuVector)
  {
    if(truNu.isNull()){
      continue;
    }

    const simb::MCNeutrino mcn = truNu->GetNeutrino();
    const simb::MCParticle nu  = mcn.Nu();

    bool signal_deuteron = false;

    for(int i = 0; i < truNu->NParticles(); ++i)
    {
      
      const simb::MCParticle mcp = truNu->GetParticle(i);

      if(mcp.PdgCode() == 1000010020 && mcp.P() > minP)
      {
        signal_deuteron = true;
      }
    
    }

    if(signal_deuteron)
      return true;
    }

  return false;

}

bool nuclearFragments::DeuteronFilter::VolumeCheck(const TVector3 &pos, const double &walls, const double &cath, const double &front, const double &back)
{
  const bool xedges = pos.X() < (200. - walls) && pos.X() > (-200. + walls);
  const bool yedges = pos.Y() < (200. - walls) && pos.Y() > (-200. + walls);
  const bool zedges = pos.Z() < (500. - back)  && pos.Z() > (0. + front);
  const bool caths  = pos.X() > cath || pos.X() < -cath;

  return xedges && yedges && zedges && caths;
}

DEFINE_ART_MODULE(nuclearFragments::DeuteronFilter)

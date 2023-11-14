#ifndef SECONDSHOWERFINDERALG_H_SEEN
#define SECONDSHOWERFINDERALG_H_SEEN

///////////////////////////////////////////////
// SecondShowerFinderAlg.h
//
// Alg intended to find the subleading photon
// missed from the reconstruction of neutral
// pion candidates.
//
// Author: Henry Lay (h.lay@lancaster.ac.uk)
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"

class SecondShowerFinderAlg
{
 public:
  SecondShowerFinderAlg();

  SecondShowerFinderAlg(fhicl::ParameterSet const& p);

  bool FindSecondShower(const art::Ptr<recob::Shower> &shower, const std::vector<art::Ptr<recob::Hit>> &hits);

};

#endif

////////////////////////////////////////////////////////////////////////
// Class:       NCPiZeroFilter
// Plugin Type: filter (Unknown Unknown)
// File:        NCPiZeroFilter_module.cc
//
// Generated at Mon Nov 27 04:33:34 2023 by Henry Lay using cetskelgen
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
#include "canvas/Persistency/Common/Ptr.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include <memory>

namespace sbnd {
  class NCPiZeroFilter;
}

class sbnd::NCPiZeroFilter : public art::EDFilter {
public:
  explicit NCPiZeroFilter(fhicl::ParameterSet const& p);

  NCPiZeroFilter(NCPiZeroFilter const&) = delete;
  NCPiZeroFilter(NCPiZeroFilter&&) = delete;
  NCPiZeroFilter& operator=(NCPiZeroFilter const&) = delete;
  NCPiZeroFilter& operator=(NCPiZeroFilter&&) = delete;

  bool filter(art::Event& e) override;

  bool VolumeCheck(const TVector3 &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);

private:

  art::InputTag fMCTruthModuleLabel;

};


sbnd::NCPiZeroFilter::NCPiZeroFilter(fhicl::ParameterSet const& p)
  : EDFilter{p}
  , fMCTruthModuleLabel(p.get<art::InputTag>("MCTruthModuleLabel"))
  {
  }

bool sbnd::NCPiZeroFilter::filter(art::Event& e)
{
  art::Handle<std::vector<simb::MCTruth>> mcTruthHandle;
  e.getByLabel(fMCTruthModuleLabel, mcTruthHandle);
  if(!mcTruthHandle.isValid()){
    std::cout << "MCTruth product " << fMCTruthModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  std::vector<art::Ptr<simb::MCTruth>> mcTruthVec;
  art::fill_ptr_vector(mcTruthVec, mcTruthHandle);

  for(const art::Ptr<simb::MCTruth> &mct : mcTruthVec)
    {
      if(mct->Origin() != 1)
        continue;

      const simb::MCNeutrino mcn = mct->GetNeutrino();
      const simb::MCParticle nu  = mcn.Nu();

      const bool nc = mcn.CCNC() == 1;
      const bool av = VolumeCheck(nu.Position().Vect());

      bool pizero = false;

      for(int i = 0; i < mct->NParticles(); ++i)
        {
          const simb::MCParticle mcp = mct->GetParticle(i);

          if(mcp.PdgCode() == 111)
            pizero = true;
        }

      if(nc && av && pizero)
        return true;
    }

  return false;
}

bool sbnd::NCPiZeroFilter::VolumeCheck(const TVector3 &pos, const double &walls, const double &cath, const double &front, const double &back)
{
  const bool xedges = pos.X() < (200. - walls) && pos.X() > (-200. + walls);
  const bool yedges = pos.Y() < (200. - walls) && pos.Y() > (-200. + walls);
  const bool zedges = pos.Z() < (500. - back)  && pos.Z() > (0. + front);
  const bool caths  = pos.X() > cath || pos.X() < -cath;

  return xedges && yedges && zedges && caths;
}

DEFINE_ART_MODULE(sbnd::NCPiZeroFilter)

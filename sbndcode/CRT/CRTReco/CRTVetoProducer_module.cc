////////////////////////////////////////////////////////////////////////
// Class:       CRTVetoProducer
// Plugin Type: producer
// File:        CRTVetoProducer_module.cc
//
// Author:      Alex Antonakis (aantonak@ucsb.edu)
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Math/GenVector/PositionVector2D.h"
#include "Math/GenVector/DisplacementVector2D.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"
#include "sbnobj/SBND/CRT/CRTVeto.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"


#include "TMath.h"

#include <memory>


namespace sbnd::crt {
  class CRTVetoProducer;
}


class sbnd::crt::CRTVetoProducer : public art::EDProducer {
public:
  explicit CRTVetoProducer(fhicl::ParameterSet const& p);

  CRTVetoProducer(CRTVetoProducer const&) = delete;
  CRTVetoProducer(CRTVetoProducer&&) = delete;
  CRTVetoProducer& operator=(CRTVetoProducer const&) = delete;
  CRTVetoProducer& operator=(CRTVetoProducer&&) = delete;

  void produce(art::Event& e) override;

  // Functions to check if a point is in a specific Tagger
  bool inSouth(const double &x, const double &y, const double &z);
  bool inNorth(const double &x, const double &y, const double &z);
  bool inEast(const double &x, const double &y, const double &z);
  bool inWest(const double &x, const double &y, const double &z);
  bool inTopLow(const double &x, const double &y, const double &z);
  bool inTopHigh(const double &x, const double &y, const double &z);

private:

  CRTGeoAlg        fCRTGeoAlg;
  double           fWindowStart;
  double           fWindowEnd;
  std::string      fCRTSpacePointModuleLabel;
};


sbnd::crt::CRTVetoProducer::CRTVetoProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fCRTGeoAlg(p.get<fhicl::ParameterSet>("CRTGeoAlg"))
  , fWindowStart(p.get<double>("WindowStart"))
  , fWindowEnd(p.get<double>("WindowEnd"))
  , fCRTSpacePointModuleLabel(p.get<std::string>("CRTSpacePointModuleLabel"))
{
  produces<CRTVeto>();
 // produces<art::Assns<CRTSpacePoint, CRTVeto>>();
}

void sbnd::crt::CRTVetoProducer::produce(art::Event& e)
{
  //auto crtveto            = std::make_unique<CRTVeto>();
  //auto vetoSpacePointAssn = std::make_unique<art::Assns<CRTSpacePoint, CRTVeto>>();

  // Get the Space Points for this event
  art::Handle<std::vector<CRTSpacePoint>> CRTSpacePointHandle;
  e.getByLabel(fCRTSpacePointModuleLabel, CRTSpacePointHandle);
  if(!CRTSpacePointHandle.isValid()){
    std::cout << "CRTSpacePoint product " << fCRTSpacePointModuleLabel << " not found..." << std::endl;
    throw std::exception();
  } 
  std::vector<art::Ptr<CRTSpacePoint>> CRTSpacePointVec;
  art::fill_ptr_vector(CRTSpacePointVec, CRTSpacePointHandle);

  int S = 0;
  int N = 0;
  int E = 0;
  int W = 0;
  int TL = 0;
  int TH = 0;

  // loop over the space points
  for (unsigned int i = 0; i < CRTSpacePointVec.size(); i++){

    double x = CRTSpacePointVec[i]->X();
    double y = CRTSpacePointVec[i]->Y();
    double z = CRTSpacePointVec[i]->Z();
    double t = CRTSpacePointVec[i]->Ts0();
    // check if in the beam window --> Make configureable
    if ((t > fWindowEnd) || (t < fWindowStart)) {
      continue;
    }
    // We have a valid space point --> make assn
    //util::CreateAssn(*this, e, *crtveto, CRTSpacePointVec[i], *vetoSpacePointAssn);
    
    if (inSouth(x, y, z)) {
      S += 1;
    }
    else if (inNorth(x, y, z)) {
      N += 1;
    }
    else if (inEast(x, y, z)) {
      E += 1;
    }
    else if (inWest(x, y, z)) {
      W += 1;
    }
    else if (inTopLow(x, y, z)) {
      TL += 1;
    }
    else if (inTopHigh(x, y, z)) {
      TH += 1;
    }
    else {
      std::cout << "Not a valid SpacePoint --> Outside Geometry !!!" << std::endl;
    }

  } // end loop over space points

  int Nwalls_all = S + N + E + W;
  int Nwalls_noNorth = S + E + W;
  int Ntop = TL + TH;

  bool v0 = false;
  bool v1 = false;
  bool v2 = false;
  bool v3 = false;

  // V0
  if ((Nwalls_all + Ntop) > 0) {
    v0 = true;  
  }
  // V1
  if ( (Nwalls_all > 0) || ((TL > 0) && (TH > 0)) ) {
    v1 = true;
  }
  // V2
  if ((Nwalls_noNorth + Ntop) > 0) {
    v2 = true;  
  }
  // V3
  if ( (Nwalls_noNorth > 0) || ((TL > 0) && (TH > 0)) ) {
    v3 = true;
  }

  //crtveto =  CRTVeto(v0, v1, v2, v3); 
  auto crtveto            = std::make_unique<CRTVeto>(v0, v1, v2, v3);
 
  e.put(std::move(crtveto));
  //e.put(std::move(vetoSpacePointAssn));

} // end of produce

bool sbnd::crt::CRTVetoProducer::inSouth(const double &x, const double &y, const double &z)
{
  if ((z > -1000) && (z < -179) && (y > -370) && (y < 400)) {
    return true;
  }
  return false;
}
bool sbnd::crt::CRTVetoProducer::inNorth(const double &x, const double &y, const double &z)
{
  // check that this isn't affected by x from east and west --> TODO
  if ((z > 746) && (y > -370) && (y < 400)) {
    return true;
  }
  return false;
}

bool sbnd::crt::CRTVetoProducer::inEast(const double &x, const double &y, const double &z)
{
  //east_sel = "x < -380 and z < 770 and  -370 < y < 400"
  if ((x < -380) && (z < 770) && (y > -370) && (y < 400)) {
    return true;
  }
  return false;
}

bool sbnd::crt::CRTVetoProducer::inWest(const double &x, const double &y, const double &z)
{
  //west_sel = "x < -380 and z < 770 and  -370 < y < 400"
  if ((x > 380) && (z < 770) && (y > -370) && (y < 400)) {
    return true;
  }
  return false;
}

bool sbnd::crt::CRTVetoProducer::inTopLow(const double &x, const double &y, const double &z)
{
  if ((y > 400) && (y < 700)) {
    return true;
  }
  return false;
}

bool sbnd::crt::CRTVetoProducer::inTopHigh(const double &x, const double &y, const double &z)
{
  if (y > 700) {
    return true;
  }
  return false;
}



DEFINE_ART_MODULE(sbnd::crt::CRTVetoProducer)

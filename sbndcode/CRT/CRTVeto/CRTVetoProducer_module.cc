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

#include "sbnobj/SBND/CRT/CRTEnums.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"
#include "sbnobj/SBND/CRT/CRTVeto.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "lardata/Utilities/AssociationUtil.h"

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
  bool inSouth(const CRTTagger t);
  bool inNorth(const CRTTagger t);
  bool inEast(const CRTTagger t);
  bool inWest(const CRTTagger t);
  bool inTopLow(const CRTTagger t);
  bool inTopHigh(const CRTTagger t);

private:

  CRTGeoAlg        fCRTGeoAlg;
  double           fWindowStart;
  double           fWindowEnd;
  std::string      fCRTClusterModuleLabel;
  std::string      fCRTSpacePointModuleLabel;
  bool		   fIsData;

  std::vector<int> fVetoSpacePointIndices;
  std::vector<art::Ptr<CRTSpacePoint>> fVetoSpacePoints;
};


sbnd::crt::CRTVetoProducer::CRTVetoProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fCRTGeoAlg(p.get<fhicl::ParameterSet>("CRTGeoAlg"))
  , fWindowStart(p.get<double>("WindowStart"))
  , fWindowEnd(p.get<double>("WindowEnd"))
  , fCRTClusterModuleLabel(p.get<std::string>("CRTClusterModuleLabel"))
  , fCRTSpacePointModuleLabel(p.get<std::string>("CRTSpacePointModuleLabel"))
  , fIsData(p.get<bool>("IsData"))
{
  produces<std::vector<CRTVeto>>();
  produces<art::Assns<CRTVeto, CRTSpacePoint>>();
}

void sbnd::crt::CRTVetoProducer::produce(art::Event& e)
{

  fVetoSpacePoints.clear();
  fVetoSpacePointIndices.clear();
 
  if (fIsData) {
    std::cout << "Filling CRT Veto for Data" << std::endl;

  }

 
  // Get the Clusters for this event
  art::Handle<std::vector<CRTCluster>> CRTClusterHandle;
  e.getByLabel(fCRTClusterModuleLabel, CRTClusterHandle);
  if(!CRTClusterHandle.isValid()){
    std::cout << "CRTCluster product " << fCRTClusterModuleLabel << " not found..." << std::endl;
    throw std::exception();
  } 
  std::vector<art::Ptr<CRTCluster>> CRTClusterVec;
  art::fill_ptr_vector(CRTClusterVec, CRTClusterHandle);

  // Get the Space Points for this event
  art::Handle<std::vector<CRTSpacePoint>> CRTSpacePointHandle;
  e.getByLabel(fCRTSpacePointModuleLabel, CRTSpacePointHandle);
  if(!CRTSpacePointHandle.isValid()){
    std::cout << "CRTSpacePoint product " << fCRTSpacePointModuleLabel << " not found..." << std::endl;
    throw std::exception();
  } 
  std::vector<art::Ptr<CRTSpacePoint>> CRTSpacePointVec;
  art::fill_ptr_vector(CRTSpacePointVec, CRTSpacePointHandle);

  art::FindManyP<CRTSpacePoint> cluster_sps(CRTClusterVec, e, fCRTSpacePointModuleLabel);
  
  // Get indices of Veto Space Points
  /*
  for (unsigned int i = 0; i < CRTSpacePointVec.size(); ++i) {

    double t = CRTSpacePointVec.at(i)->Ts0()/1000; // Convert to us
    if ((t > fWindowEnd) || (t < fWindowStart)) {
      continue;
    }
    fVetoSpacePointIndices.push_back(i);

  }
  */

  // initialize counters for veto hits in different taggers
  int S = 0;
  int N = 0;
  int E = 0;
  int W = 0;
  int TL = 0;
  int TH = 0;

  // loop over clusters in this event
  // Note: there are usually more clusters than space points
  for(auto const &cluster : CRTClusterVec) {
  
    std::vector<art::Ptr<CRTSpacePoint>> clusterSpacePoints = cluster_sps.at(cluster.key());
  
    // loop over the space points in this cluster --> should be one per cluster
    for (auto const &sp : clusterSpacePoints) {
     
      // Get the time and tagger for this space point
      double t = sp->Ts0()/1000; // Convert to us
      CRTTagger tag = cluster->Tagger();  

      // check if in the beam window --> Make configureable
      if ((t > fWindowEnd) || (t < fWindowStart)) {
        continue;
      }
      std::cout << std::endl;
      std::cout << "Found Space Point in time with the beam !!!" << std::endl;
      std::cout << std::endl;
    
      fVetoSpacePoints.push_back(sp);
      fVetoSpacePointIndices.push_back(sp.key());

      if (inSouth(tag)) {
        S += 1;
      }
      else if (inNorth(tag)) {
        N += 1;
      }
      else if (inEast(tag)) {
        E += 1;
      }
      else if (inWest(tag)) {
        W += 1;
      }
      else if (inTopLow(tag)) {
        TL += 1;
      }
      else if (inTopHigh(tag)) {
        TH += 1;
      }
      else {
        std::cout << "Not a valid SpacePoint --> Outside Geometry !!!" << std::endl;
      }
    } // end loop over space points
  } // end of loop over clusters
 
  int Nwalls_all = S + N + E + W;
  int Nwalls_noNorth = S + E + W;
  int Ntop = TL + TH;

  bool v0 = false;
  bool v1 = false;
  bool v2 = false;
  bool v3 = false;
  bool v4 = false;

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
  if (S > 0)  {
    v4 = true;
  }

  //crtveto =  CRTVeto(v0, v1, v2, v3); 
  //
  //
  // Make Products for the Event
  auto crtveto            = std::make_unique<CRTVeto>(v0, v1, v2, v3, v4);
  auto vetoVec            = std::make_unique<std::vector<CRTVeto>>();
  vetoVec->push_back(*crtveto);

  auto vetoSpacePointAssn = std::make_unique<art::Assns<CRTVeto, CRTSpacePoint>>();

  //auto const vetoPtrMaker = art::PtrMaker<CRTVeto>(e);
  //auto const vetoPtr = vetoPtrMaker(0); 
  //auto const CRTSpacePointPtrMaker = art::PtrMaker<CRTSpacePoint>(e, CRTSpacePointHandle.id());

  if (fVetoSpacePoints.size() != fVetoSpacePointIndices.size()) {
    std::cout << std::endl;
    std::cout << "Problem: Number of Space Points does not match!" << std::endl;
    std::cout << std::endl;
  }

  util::CreateAssn(*this, e, *vetoVec, fVetoSpacePoints, *vetoSpacePointAssn);
  
  /*
  for (unsigned int i = 0; i < fVetoSpacePoints.size(); ++i) {
    std::cout << "Veto SP index " << fVetoSpacePointIndices.at(i) << std::endl;
    //auto const CRTSpacePointPtr = CRTSpacePointPtrMaker(fVetoSpacePointIndices.at(i));
    //vetoSpacePointAssn->addSingle(vetoPtr, CRTSpacePointPtr);
    //vetoSpacePointAssn->addSingle(crtveto, CRTSpacePointPtr);
    util::CreateAssn(*this, e, *vetoVec, fVetoSpacePoints.at(i), *vetoSpacePointAssn);
  }
  */

  //e.put(std::move(crtveto));
  e.put(std::move(vetoVec));
  e.put(std::move(vetoSpacePointAssn));

} // end of produce

bool sbnd::crt::CRTVetoProducer::inSouth(const CRTTagger t)
{
  if (t == kSouthTagger) {
    return true;
  }
  return false;
}
bool sbnd::crt::CRTVetoProducer::inNorth(const CRTTagger t)
{
  if (t == kNorthTagger) {
    return true;
  }
  return false;
}

bool sbnd::crt::CRTVetoProducer::inEast(const CRTTagger t)
{
  //east_sel = "x < -380 and z < 770 and  -370 < y < 400"
  if (t == kEastTagger) {
    return true;
  }
  return false;
}

bool sbnd::crt::CRTVetoProducer::inWest(const CRTTagger t)
{
  //west_sel = "x < -380 and z < 770 and  -370 < y < 400"
  if (t == kWestTagger) {
    return true;
  }
  return false;
}

bool sbnd::crt::CRTVetoProducer::inTopLow(const CRTTagger t)
{
  if (t == kTopLowTagger) {
    return true;
  }
  return false;
}

bool sbnd::crt::CRTVetoProducer::inTopHigh(const CRTTagger t)
{
  if (t == kTopHighTagger) {
    return true;
  }
  return false;
}



DEFINE_ART_MODULE(sbnd::crt::CRTVetoProducer)

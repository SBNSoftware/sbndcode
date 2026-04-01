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

#include "artdaq-core/Data/RawEvent.hh"

#include "Math/GenVector/PositionVector2D.h"
#include "Math/GenVector/DisplacementVector2D.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "sbndcode/Timing/SBNDRawTimingObj.h"

#include "sbnobj/SBND/Timing/DAQTimestamp.hh"
#include "sbnobj/SBND/CRT/CRTEnums.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"
#include "sbnobj/SBND/CRT/CRTVeto.hh"
#include "sbnobj/SBND/Timing/FrameShiftInfo.hh"

#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "art_root_io/TFileService.h"
#include "TNtuple.h"

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

  double           fWindowStart;
  double           fWindowEnd;
  std::string      fCRTClusterModuleLabel;
  std::string      fCRTSpacePointModuleLabel;
  std::string      fCRTTimingReferenceInfoLabel;
  bool		   fIsData;
  bool 		   fDebug; 
 
  std::string      fFrameShiftModuleLabel;

  unsigned int fEventID;
  unsigned int fRun;
  unsigned int fSubRun;

};


sbnd::crt::CRTVetoProducer::CRTVetoProducer(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fWindowStart(p.get<double>("WindowStart"))
  , fWindowEnd(p.get<double>("WindowEnd"))
  , fCRTClusterModuleLabel(p.get<std::string>("CRTClusterModuleLabel"))
  , fCRTSpacePointModuleLabel(p.get<std::string>("CRTSpacePointModuleLabel"))
  , fCRTTimingReferenceInfoLabel(p.get<std::string>("CRTTimingReferenceInfoLabel"))
  , fIsData(p.get<bool>("IsData"))
  , fDebug(p.get<bool>("Debug", false))
  , fFrameShiftModuleLabel(p.get<std::string>("FrameShiftModuleLabel", ""))
{
  produces<std::vector<CRTVeto>>();
  produces<art::Assns<CRTVeto, CRTSpacePoint>>();
}

void sbnd::crt::CRTVetoProducer::produce(art::Event& e)
{
  // Declare Space Point Vector. Will fill this with Space Points that we veto on
  std::vector<art::Ptr<CRTSpacePoint>> fVetoSpacePoints;

  fEventID = e.id().event();
  fRun = e.run();
  fSubRun = e.subRun();

  if (fDebug) {
    std::cout << std::endl;
    std::cout << "Filling CRTVeto for Run " << fRun << " SubRun " << fSubRun << " Event " << fEventID << std::endl;
    std::cout << std::endl;
  }

  // Data needs to be handled differently than MC  
  if (fIsData) {
    
    // Check that FrameShift module was applied successfully
    art::Handle<sbnd::timing::FrameShiftInfo> frameShiftHandle;
    e.getByLabel(fFrameShiftModuleLabel, frameShiftHandle);
    
    if(!frameShiftHandle.isValid())
      throw std::runtime_error("Frame Shift Info object is invalid, check data quality");

    art::Handle<raw::TimingReferenceInfo> TimingReferenceHandle;
    e.getByLabel(fCRTTimingReferenceInfoLabel, TimingReferenceHandle);
    if(TimingReferenceHandle.isValid())
    {
      uint16_t TType = TimingReferenceHandle->timingType;
      // Check if the timing is valid
      if (TType != 0 && TType != 1) {
        
        if (fDebug) { 
          std::cout << std::endl;
          std::cout << "TType is not valid for this event --> Don't Flag Event!" << std::endl;
          std::cout << std::endl;
        }    
        auto crtveto            = std::make_unique<CRTVeto>(false, false, false, false, false);
        auto vetoVec            = std::make_unique<std::vector<CRTVeto>>();
        vetoVec->push_back(*crtveto);
        auto vetoSpacePointAssn = std::make_unique<art::Assns<CRTVeto, CRTSpacePoint>>();
        e.put(std::move(vetoVec));
        e.put(std::move(vetoSpacePointAssn));
        return;

      }

    }

  } // end of check for Data

  // Get the Clusters for this event
  art::Handle<std::vector<CRTCluster>> CRTClusterHandle;
  e.getByLabel(fCRTClusterModuleLabel, CRTClusterHandle);
  if(!CRTClusterHandle.isValid()){
    if (fDebug)
      std::cout << "CRTCluster product " << fCRTClusterModuleLabel << " not found..." << std::endl;
    throw std::exception();
  } 
  std::vector<art::Ptr<CRTCluster>> CRTClusterVec;
  art::fill_ptr_vector(CRTClusterVec, CRTClusterHandle);

  art::FindManyP<CRTSpacePoint> cluster_sps(CRTClusterVec, e, fCRTSpacePointModuleLabel);
  
  // initialize counters for veto hits in different taggers
  int S = 0;
  int N = 0;
  int E = 0;
  int W = 0;
  int TL = 0;
  int TH = 0;

  // loop over clusters in this event
  // Note: there are usually more clusters than space points
  //
  if (fDebug) {
    std::cout << std::endl;
    std::cout << "Start Loop over Clusters ..." << std::endl;
    std::cout << std::endl;
  }
  for(auto const &cluster : CRTClusterVec) {
  
    std::vector<art::Ptr<CRTSpacePoint>> clusterSpacePoints = cluster_sps.at(cluster.key());
  
    // loop over the space points in this cluster --> should be one per cluster
    for (auto const &sp : clusterSpacePoints) {

      // Get the time and tagger for this space point
      double t = sp->Ts0(); // Starts in ns

      t /= 1000; // convert to us

      CRTTagger tag = cluster->Tagger();  

      // check if in the beam window --> Make configureable
      if ((t > fWindowEnd) || (t < fWindowStart)) {
        continue;
      }
   
      // We Found a Space Point in our Veto Window !!!

      fVetoSpacePoints.push_back(sp);
      //fVetoSpacePointIndices.push_back(sp.key());

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
        if (fDebug)
          std::cout << "Not a valid SpacePoint --> Outside Geometry !!!" << std::endl;
      }
    } // end loop over space points
  } // end of loop over clusters
  
  if (fDebug) {
    std::cout << "Finished Loop over Clusters --> Produce CRTVeto products" << std::endl;  
    std::cout << std::endl; 
  }

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
  // V4
  if (S > 0)  {
    v4 = true;
  }

  // Make Products for the Event
  // See sbnobj/SBND/CRT/CRTVeto.hh for description of V0, V1, V2, V3, V4
  auto crtveto            = std::make_unique<CRTVeto>(v0, v1, v2, v3, v4);
  auto vetoVec            = std::make_unique<std::vector<CRTVeto>>();
  vetoVec->push_back(*crtveto);

  auto vetoSpacePointAssn = std::make_unique<art::Assns<CRTVeto, CRTSpacePoint>>();

  util::CreateAssn(*this, e, *vetoVec, fVetoSpacePoints, *vetoSpacePointAssn);
  
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
  if (t == kEastTagger) {
    return true;
  }
  return false;
}

bool sbnd::crt::CRTVetoProducer::inWest(const CRTTagger t)
{
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

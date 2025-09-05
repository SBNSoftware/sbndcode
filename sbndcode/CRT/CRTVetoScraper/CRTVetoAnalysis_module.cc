////////////////////////////////////////////////////////////////////////
// Class:       CRTVetoAnalysis
// Plugin Type: analyzer
// File:        CRTVetoAnalysis_module.cc
// Author:      Alex Antonakis (aantonakis@ucsb.edu)
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

#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/ChannelMaps/CRT/CRTChannelMapService.h"
#include "sbndcode/CRT/CRTBackTracker/CRTBackTrackerAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"
#include "sbndcode/Decoders/PTB/sbndptb.h"
#include "sbndcode/Timing/SBNDRawTimingObj.h"
#include "sbnobj/SBND/CRT/CRTVeto.hh"
#include "TNtuple.h"

namespace sbnd::crt {
  class CRTVetoAnalysis;
}

class sbnd::crt::CRTVetoAnalysis : public art::EDAnalyzer {
public:
  explicit CRTVetoAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTVetoAnalysis(CRTVetoAnalysis const&) = delete;
  CRTVetoAnalysis(CRTVetoAnalysis&&) = delete;
  CRTVetoAnalysis& operator=(CRTVetoAnalysis const&) = delete;
  CRTVetoAnalysis& operator=(CRTVetoAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void beginJob() override;

private:

  unsigned int fEventID;
  unsigned int fRun;
  unsigned int fSubRun;

  TNtuple *fVeto_tree;
};

sbnd::crt::CRTVetoAnalysis::CRTVetoAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} 
{

}

void sbnd::crt::CRTVetoAnalysis::analyze(art::Event const& e)
{

  fRun = e.id().run();
  fSubRun = e.id().subRun();
  fEventID =  e.id().event();

  std::cout << "Run " << fRun << " SubRun " << fSubRun << " Event " << fEventID << std::endl;
  std::cout << std::endl;  

  art::Handle<std::vector<CRTVeto>> CRTVetoHandle;
  std::vector<art::Ptr<CRTVeto>> CRTVetoVec;
  if (e.getByLabel("crtveto", CRTVetoHandle))
    art::fill_ptr_vector(CRTVetoVec, CRTVetoHandle);

  if (CRTVetoVec.empty()) {
    std::cout << "empty CRTVeto in this event :(" << std::endl;
    return;
  }
  fVeto_tree->Fill(fRun, fSubRun, fEventID, CRTVetoVec[0]->V0(), CRTVetoVec[0]->V1(), CRTVetoVec[0]->V2(), CRTVetoVec[0]->V3(), CRTVetoVec[0]->V4());
  
  if (CRTVetoVec[0]->V0()) {
    std::cout << "Flagged an Event with V0" << std::endl;
  }

}

void sbnd::crt::CRTVetoAnalysis::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;

  fVeto_tree = tfs->make<TNtuple>("veto_tree", "veto_tree", "run:subrun:evt:v0:v1:v2:v3:v4"); 

}

DEFINE_ART_MODULE(sbnd::crt::CRTVetoAnalysis)



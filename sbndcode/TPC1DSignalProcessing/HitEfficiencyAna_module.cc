// HitEfficiencyAna_module.cc
//
//Shamelessly stolen from ICARUS by A. Scarff
//
// A basic "skeleton" to read in art::Event records from a file,
// access their information, and do something with them.
// See
// <https://cdcvs.fnal.gov/redmine/projects/larsoftsvn/wiki/Using_the_Framework>
// for a description of the ART classes used here.
// Almost everything you see in the code below may have to be changed
// by you to suit your task. The example task is to make histograms
// and n-tuples related to dE/dx of particle tracks in the detector.
// As you try to understand why things are done a certain way in this
// example ("What's all this stuff about 'auto const&'?"), it will help
// to read ADDITIONAL_NOTES.txt in the same directory as this file.
#ifndef HitEfficiencyAna_module
#define HitEfficiencyAna_module
// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
//#include "cetlib/search_path.h"
#include "cetlib/cpu_timer.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "sbndcode/TPC1DSignalProcessing/tools/IHitEfficiencyHistogramTool.h"
#include "TTree.h"
// C++ Includes
#include <map>
#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
namespace HitEfficiencyAna
{
  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition
  class HitEfficiencyAna : public art::EDAnalyzer
  {
  public:
    // Standard constructor and destructor for an ART module.
    explicit HitEfficiencyAna(fhicl::ParameterSet const& pset);
    virtual ~HitEfficiencyAna();
    // This method is called once, at the start of the job. In this
    // example, it will define the histograms and n-tuples we'll write.
    void beginJob();
    void endJob();
    // This method is called once, at the start of each run. It's a
    // good place to read databases or files that may have
    // run-dependent information.
    void beginRun(const art::Run& run);
    // This method reads in any parameters from the .fcl files. This
    // method is called 'reconfigure' because it might be called in the
    // middle of a job; e.g., if the user changes parameter values in an
    // interactive event display.
    void reconfigure(fhicl::ParameterSet const& pset);
    // The analysis routine, called once per event.
    void analyze (const art::Event& evt);
  private:
    // The variables that will go into the n-tuple.
    int fEvent;
    int fRun;
    int fSubRun;
    int fNumEvents;
    
    // Introduce a temporary event counter to handle the current issue where all input
    // files start with the same run/subrun/event sequence
    int fNumCalls;
    
    // Keep track of the hit histogramming tools here
    std::vector<std::unique_ptr<IHitEfficiencyHistogramTool>> fHitHistogramToolVec;
    
    // And our tuple
    mutable TTree* fTree;
    // Other variables that will be shared between different methods.
    const geo::GeometryCore*           fGeometry;       // pointer to Geometry service
    //const lariov::DetPedestalProvider& fPedestalRetrievalAlg; ///< Keep track of an instance to the pedestal retrieval alg
  }; // class HitEfficiencyAna
  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation
  //-----------------------------------------------------------------------
  // Constructor
  HitEfficiencyAna::HitEfficiencyAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet),
      fNumCalls(0)
      //      fPedestalRetrievalAlg(*lar::providerFrom<lariov::DetPedestalService>())
  {
    fGeometry = lar::providerFrom<geo::Geometry>();
    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
  }
  //-----------------------------------------------------------------------
  // Destructor
  HitEfficiencyAna::~HitEfficiencyAna()
  {}
  //-----------------------------------------------------------------------
  void HitEfficiencyAna::beginJob()
  {
    // Get the detector length, to determine the maximum bin edge of one
    // of the histograms.
    //double detectorLength = fGeometry->DetLength();
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us.
    art::ServiceHandle<art::TFileService> tfs;
    std::cout << "beginjob running" << std::endl;
    // First task is to create a TTree..
    fTree = tfs->makeAndRegister<TTree>("HitEffic_t","Hit Efficiency Tuple");
    fTree->Branch("Run",      &fRun,       "Run/I");
    fTree->Branch("SubRun",   &fSubRun,    "SubRun/I");
    fTree->Branch("Event",    &fEvent,     "Event/I");
    // The arguments to 'make<whatever>' are the same as those passed
    // to the 'whatever' constructor, provided 'whatever' is a ROOT
    // class that TFileService recognizes.
    for (auto& hitHistTool : fHitHistogramToolVec)
      {
        hitHistTool->initializeHists(tfs, "HitEffic");
        hitHistTool->initializeTuple(fTree);
      }
    // zero out the event counter
    fNumEvents = 0;
  }
  //-----------------------------------------------------------------------
  void HitEfficiencyAna::beginRun(const art::Run& /*run*/)
  {
    // How to convert from number of electrons to GeV.  The ultimate
    // source of this conversion factor is
    // ${LARSIM_DIR}/include/SimpleTypesAndConstants/PhysicalConstants.h.
    //    art::ServiceHandle<sim::LArG4Parameters> larParameters;
    //    fElectronsToGeV = 1./larParameters->GeVToElectrons();
  }
  //-----------------------------------------------------------------------
  void HitEfficiencyAna::reconfigure(fhicl::ParameterSet const& p)
  {
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    // Implement the tools for handling the responses
    const std::vector<fhicl::ParameterSet>& hitHistogramToolVec = p.get<std::vector<fhicl::ParameterSet>>("HitEfficiencyHistogramToolList");
    
    for(auto& hitHistogramTool : hitHistogramToolVec)
      fHitHistogramToolVec.push_back(art::make_tool<IHitEfficiencyHistogramTool>(hitHistogramTool));
    std::cout << "reconfigure running" << std::endl;
    return;
  }
  //-----------------------------------------------------------------------
  void HitEfficiencyAna::analyze(const art::Event& event)
  {
    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event();
    fRun    = ++fNumCalls; //event.run();
    fSubRun = event.subRun();
    fNumEvents++;
    
    // Make a pass through all hits to make contrasting plots
    std::cout << "-- Run: " << fRun << ", SubRun: " << fSubRun << ", Event: " << fEvent << " -------" << std::endl;
    std::cout << "HELLO!" << std::endl;

    for(auto& hitHistTool : fHitHistogramToolVec) hitHistTool->fillHistograms(event);
    
    fTree->Fill();
    return;
  }
  void HitEfficiencyAna::endJob()
  {
    // Make a call to normalize histograms if so desired
    for(auto& hitHistTool : fHitHistogramToolVec) hitHistTool->endJob(fNumEvents);
    std::cout << "endjob running" << std::endl;
    return;
  }
  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see HitEfficiencyAna.fcl for more information.
  DEFINE_ART_MODULE(HitEfficiencyAna)
} // namespace HitEfficiencyAna
#endif // HitEfficiencyAna_module

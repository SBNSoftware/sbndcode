////////////////////////////////////////////////////////////////////////
// Class:       OpDetAnalyzer
// Plugin Type: analyzer (Unknown Unknown)
// File:        OpDetAnalyzer_module.cc
//
// Generated at Tue Oct 19 04:32:01 2021 by Patrick Green using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include <TTree.h>

#include "larcore/Geometry/Geometry.h"

#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/sim.h"

#include "sbncode/OpDet/PDMapAlg.h"

#include <iostream>

namespace opdet {
  class OpDetAnalyzer;
}


class opdet::OpDetAnalyzer : public art::EDAnalyzer {
public:
  explicit OpDetAnalyzer(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OpDetAnalyzer(OpDetAnalyzer const&) = delete;
  OpDetAnalyzer(OpDetAnalyzer&&) = delete;
  OpDetAnalyzer& operator=(OpDetAnalyzer const&) = delete;
  OpDetAnalyzer& operator=(OpDetAnalyzer&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  std::vector<std::string> fInputModule;      // Input tag for OpDet collection

  // PDS map
  std::unique_ptr<opdet::PDMapAlg> fPDSMapPtr;
  std::vector<std::string> fPMTMapLabel, fArapucaMapLabel;

  std::vector<int> fNOpDetAll;
  std::vector<int> fNOpDetDirect;
  std::vector<int> fNOpDetReflected;

  TTree *fAllPhotonsTree;
  TTree *fTheOpDetTree;
  TTree *fTheEventTree;

  int fEventID;
  int fOpChannel;
  float fWavelength;
  float fTime;

  int fCountOpDetAll;
  int fCountOpDetDirect;
  int fCountOpDetReflected;
  float fOpDetX, fOpDetY, fOpDetZ;
  bool fisPMT; 

  int fCountEventAll;
  int fCountEventDirect;
  int fCountEventReflected;

  int fVerbosity; // debugging

  /// Value used when a typical ultraviolet light wavelength is needed.
  static constexpr double kVUVWavelength = 128.0; // nm

  /// Value used when a typical visible light wavelength is needed.
  static constexpr double kVisibleWavelength = 450.0; // nm

};


opdet::OpDetAnalyzer::OpDetAnalyzer(fhicl::ParameterSet const& pset)
  : EDAnalyzer{pset},
    fPDSMapPtr{art::make_tool<opdet::PDMapAlg>(pset.get<fhicl::ParameterSet>("PDSMapTool"))}
{
  
  fVerbosity = pset.get<int>("Verbosity");
  
  try
  {
     fInputModule = pset.get<std::vector<std::string>>("InputModule");
  }
  catch(...)
  {
     fInputModule.push_back(pset.get<std::string>("InputModule"));
  }
}

void opdet::OpDetAnalyzer::beginJob()
{
  // create output trees
  art::ServiceHandle<art::TFileService> tfs;
  
  // all photons tree
  fAllPhotonsTree = tfs->make<TTree>("AllPhotons", "AllPhotons");
  fAllPhotonsTree->Branch("EventID", &fEventID, "EventID/I");
  fAllPhotonsTree->Branch("Wavelength", &fWavelength, "Wavelength/F");
  fAllPhotonsTree->Branch("OpChannel", &fOpChannel, "OpChannel/I");
  fAllPhotonsTree->Branch("Time", &fTime, "Time/F");
  
  // photons per opdet tree
  fTheOpDetTree = tfs->make<TTree>("PhotonsPerOpDet","PhotonsPerOpDet");
  fTheOpDetTree->Branch("EventID", &fEventID, "EventID/I");
  fTheOpDetTree->Branch("OpChannel", &fOpChannel, "OpChannel/I");
  fTheOpDetTree->Branch("OpDetX", &fOpDetX, "OpDetX/F");
  fTheOpDetTree->Branch("OpDetY", &fOpDetY, "OpDetY/F");
  fTheOpDetTree->Branch("OpDetZ", &fOpDetZ, "OpDetZ/F");
  fTheOpDetTree->Branch("isPMT", &fisPMT);
  fTheOpDetTree->Branch("CountAll", &fCountOpDetAll, "CountAll/I");
  fTheOpDetTree->Branch("CountDirect", &fCountOpDetDirect, "CountDirect/I");
  fTheOpDetTree->Branch("CountReflected", &fCountOpDetReflected, "CountReflected/I");

  // photons per event tree
  fTheEventTree  = tfs->make<TTree>("PhotonsPerEvent","PhotonsPerEvent");
  fTheEventTree->Branch("EventID", &fEventID, "EventID/I");
  fTheEventTree->Branch("CountAll", &fCountEventAll, "CountAll/I");
  fTheEventTree->Branch("CountDirect", &fCountEventDirect, "CountDirect/I");
  fTheEventTree->Branch("CountReflected", &fCountEventReflected, "CountReflected/I");

}

void opdet::OpDetAnalyzer::analyze(art::Event const& evt)
{

  // geometry information
  art::ServiceHandle<geo::Geometry> geo;

  // get event number
  fEventID = evt.id().event();
  
  // adapted from standard SimPhotonCounter code

  // Get *ALL* SimPhotonsCollection from Event
  auto photon_handles = evt.getMany<std::vector<sim::SimPhotonsLite>>();
  if (photon_handles.size() == 0)
    throw art::Exception(art::errors::ProductNotFound)<<"sim SimPhotons retrieved and you requested them.";

  // reset counters
  // event
  fCountEventAll=0;
  fCountEventDirect=0;
  fCountEventReflected=0;
  // opdet
  fCountOpDetAll=0;
  fCountOpDetDirect=0;
  fCountOpDetReflected=0;

  // reset vectors
  fNOpDetAll.clear(); fNOpDetAll.resize(geo->NOpChannels(), 0);
  fNOpDetDirect.clear(); fNOpDetDirect.resize(geo->NOpChannels(), 0);
  fNOpDetReflected.clear(); fNOpDetReflected.resize(geo->NOpChannels(), 0);

  // Get SimPhotonsLite from Event
  for(auto const& mod : fInputModule){

    // Loop over direct/reflected photons
    for (auto const& ph_handle: photon_handles) {
      // Do some checking before we proceed
      if (!ph_handle.isValid()) continue;
      if (ph_handle.provenance()->moduleLabel() != mod) continue;   //not the most efficient way of doing this, but preserves the logic of the module. Andrzej

      bool Reflected = (ph_handle.provenance()->productInstanceName() == "Reflected");
      
      if(fVerbosity > 0) std::cout<<"Found OpDet hit collection of size "<< (*ph_handle).size()<<std::endl;

      if((*ph_handle).size()>0)
      {

        for ( auto const& photon : (*ph_handle) )
        {
          //Get data from HitCollection entry
          fOpChannel=photon.OpChannel;
          std::map<int, int> PhotonsMap = photon.DetectedPhotons;

          // do not save if PD is not sensitive to this light
          std::string pd_type=fPDSMapPtr->pdType(fOpChannel);
          if (Reflected && pd_type=="xarapuca_vuv") continue;
          if(!Reflected && (pd_type=="xarapuca_vis" || pd_type=="pmt_uncoated")) continue;

          for(auto it = PhotonsMap.begin(); it!= PhotonsMap.end(); it++)
          {

            // Calculate wavelength in nm
            if (Reflected) fWavelength = kVisibleWavelength;
            else fWavelength= kVUVWavelength;

            // Get arrival time from phot
            fTime= it->first;

            for(int i = 0; i < it->second ; i++)
            {
              // Increment per OpDet counters and fill all photon tree
              fCountOpDetAll++;
              fAllPhotonsTree->Fill();

              // all
              fNOpDetAll[fOpChannel]++;
  
              // direct light
              if (!Reflected){
                fNOpDetDirect[fOpChannel]++;
              }
              // reflected light
              else if (Reflected) {
                fNOpDetReflected[fOpChannel]++;
              }            
            }
          }
        } // eng loop over photons        
      }
    } // end loop over photon handles
  } // end loop over input modules


  // fill information in perOpDet and perEvent trees
  for (unsigned int OpDet = 0; OpDet < geo->NOpChannels(); OpDet++) {
    
    // perOpDet tree
    fOpChannel = OpDet;
    geo::Point_t OpDetCenter = geo->OpDetGeoFromOpDet(OpDet).GetCenter();
    fOpDetX = OpDetCenter.X();
    fOpDetY = OpDetCenter.Y();
    fOpDetZ = OpDetCenter.Z();
    fisPMT =  geo->OpDetGeoFromOpDet(OpDet).isSphere();
    fCountOpDetAll = fNOpDetAll[OpDet];
    fCountOpDetDirect = fNOpDetDirect[OpDet];
    fCountOpDetReflected = fNOpDetReflected[OpDet];
    fTheOpDetTree->Fill();

    // increment counters for perEvent tree
    fCountEventAll+=fCountOpDetAll;
    fCountEventDirect+=fCountOpDetDirect;
    fCountEventReflected+=fCountOpDetReflected;

     if(fVerbosity >2) std::cout<<"OpDetResponseInterface PerOpDet : Event "<<fEventID<<" OpDet " << fOpChannel << " All " << fCountOpDetAll << " Dir " <<fCountOpDetDirect << ", Refl " << fCountOpDetReflected<<std::endl;
  }

  // perEvent tree
  fTheEventTree->Fill();

  if(fVerbosity >1) std::cout<<"OpDetResponseInterface PerEvent : Event "<<fEventID<<" All " << fCountEventAll << " Dir " <<fCountEventDirect << ", Refl " << fCountEventReflected <<std::endl;

}

void opdet::OpDetAnalyzer::endJob()
{
  // Implementation of optional member function here.
}


DEFINE_ART_MODULE(opdet::OpDetAnalyzer)

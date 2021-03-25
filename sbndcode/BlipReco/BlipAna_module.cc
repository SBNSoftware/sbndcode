////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////

#ifndef BLIPANA_H
#define BLIPANA_H

//##########################
//### Framework includes ###
//##########################
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//########################
//### LArSoft includes ###
//########################
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
//#include "nusimdata/SimulationBase/GTruth.h"
//#include "nusimdata/SimulationBase/MCTruth.h"
//#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
//#include "lardataobj/RawData/BeamInfo.h"
//#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
//#include "lardataobj/AnalysisBase/CosmicTag.h"
//#include "lardataobj/AnalysisBase/FlashMatch.h"
//#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

//########################
//### SBNDCode includes###
//########################
#include "sbndcode/RecoUtils/RecoUtils.h"

//########################
//### C++ includes     ###
//########################
#include <cstring> // std::memcpy()
#include <vector>
#include <map>
#include <utility>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional> // std::mem_fn()
#include <typeinfo>
#include <cmath>

//########################
//### ROOT includes    ###
//########################
#include "TTree.h"
#include "TTimeStamp.h"

// Set global constants and max array sizes
const int kNplanes  = 3;     //number of wire planes
const int kMaxHits  = 25000; //maximum number of hits;

namespace{  
  template <typename ITER, typename TYPE> 
    inline void FillWith(ITER from, ITER to, TYPE value)
    { std::fill(from, to, value); }
  template <typename ITER, typename TYPE>
    inline void FillWith(ITER from, size_t n, TYPE value)
    { std::fill(from, from + n, value); }
  template <typename CONT, typename V>
    inline void FillWith(CONT& data, const V& value)
    { FillWith(std::begin(data), std::end(data), value); }
}

namespace sbnd { 

  //====================================================
  //  Data storage structure
  //====================================================
  class BlipAnaTreeDataStruct 
  {
    public:
    // === Function for resetting data ===
    void Clear();
    // === Hit information ===
    int	  no_hits;                  //number of hits
    int   hit_id[kMaxHits];         //hit ID (in case only subset of hits are saved)
    int	  hit_tpc[kMaxHits];        //tpc number
    int	  hit_plane[kMaxHits];      //plane number
    int	  hit_wire[kMaxHits];       //wire number
    int	  hit_channel[kMaxHits];    //channel ID
    float	hit_peakT[kMaxHits];      //peak time (tick)
    float hit_width[kMaxHits];      //shape RMS
    float	hit_ph[kMaxHits];         //amplitude
    float	hit_area[kMaxHits];     //charge (area) in ADC units
    int	  hit_trkid[kMaxHits];      //is this hit associated with a reco track?
    int	  hit_mcid[kMaxHits];       //TrackID of leading MCParticle that created the hit
    float hit_mcfrac[kMaxHits];       //fraction of hit energy from leading MCParticle
    float hit_mcenergy[kMaxHits];     //true energy
    float hit_mccharge[kMaxHits];    //true number of electrons (drift-attenuated)
    float hit_charge[kMaxHits];      //reconstructed number of electrons
  };
  
  
  //====================================================
  //  BlipAna class definition
  //====================================================
  class BlipAna : public art::EDAnalyzer 
  { 
    public:
    explicit BlipAna(fhicl::ParameterSet const& pset);
    virtual ~BlipAna();
    void analyze(const art::Event& evt);
//    void reconfigure(fhicl::ParameterSet const &p);
    
    private:
    BlipAnaTreeDataStruct* fData;
    TTree* fTree;
    std::string fSomeLabel;
  };

}

//====================================================
//  Method for resetting data structure
//====================================================
void sbnd::BlipAnaTreeDataStruct::Clear() 
{
  no_hits = 0;
  FillWith(hit_tpc,     -999);
  FillWith(hit_id,      -999);
  FillWith(hit_plane,   -999);
  FillWith(hit_wire,    -999);
  FillWith(hit_channel, -999);
  FillWith(hit_peakT,   -999);
  FillWith(hit_width,   -999);
  FillWith(hit_ph,      -999);
  FillWith(hit_area,    -999);
  FillWith(hit_trkid,   -999);
  FillWith(hit_mcid,    -999);
  FillWith(hit_mcfrac,  -999);
  FillWith(hit_mcenergy,-999);
  FillWith(hit_mccharge,-999);
  FillWith(hit_charge,  -999);
}



//====================================================
//  BlipAna constructor
//====================================================
sbnd::BlipAna::BlipAna(fhicl::ParameterSet const& pset) : 
  EDAnalyzer(pset),
  fData(nullptr),
  fSomeLabel(pset.get<std::string>("SomeLabel"))
{
  if(!fData) fData = new BlipAnaTreeDataStruct();
  fData->Clear();
}
// Destructor
sbnd::BlipAna::~BlipAna(){}


//====================================================
//  Main event-by-event analysis
//====================================================
void sbnd::BlipAna::analyze(const art::Event& evt)
{
  std::cout<<"Processing event "<<evt.id().event()<<"\n";
  fData->Clear();
  //fTree->Fill();
}

DEFINE_ART_MODULE(sbnd::BlipAna)

#endif

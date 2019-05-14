#ifndef CRTBACKTRACKER_H_SEEN
#define CRTBACKTRACKER_H_SEEN


///////////////////////////////////////////////
// CRTBackTracker.h
//
// Quick and dirty backtracker for SBND CRT
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"

#include "lardataobj/Simulation/AuxDetSimChannel.h"

#include "sbndcode/CRT/CRTProducts/CRTData.hh"
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"

// c++
#include <vector>


namespace sbnd{

  class CRTBackTracker {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> CRTDataLabel {
        Name("CRTDataLabel")
      };
      fhicl::Atom<art::InputTag> CRTHitLabel {
        Name("CRTHitLabel")
      };
      fhicl::Atom<art::InputTag> CRTTrackLabel {
        Name("CRTTrackLabel")
      };

    };

    CRTBackTracker(const Config& config);

    CRTBackTracker(const fhicl::ParameterSet& pset) :
      CRTBackTracker(fhicl::Table<Config>(pset, {})()) {}

    CRTBackTracker();

    ~CRTBackTracker();

    void reconfigure(const Config& config);

    // Check that two CRT data products are the same
    bool DataCompare(const crt::CRTData& data1, const crt::CRTData& data2);

    // Check that two CRT hits are the same
    bool HitCompare(const crt::CRTHit& hit1, const crt::CRTHit& hit2);

    // Check that two CRT tracks are the same
    bool TrackCompare(const crt::CRTTrack& hit1, const crt::CRTTrack& track2);

    // Get all the true particle IDs that contributed to the CRT data product
    std::vector<int> AllTrueIds(const art::Event& event, const crt::CRTData& data);

    // Get all the true particle IDs that contributed to the CRT hit
    std::vector<int> AllTrueIds(const art::Event& event, const crt::CRTHit& hit);

    // Get all the true particle IDs that contributed to the CRT track
    std::vector<int> AllTrueIds(const art::Event& event, const crt::CRTTrack& track);

    // Get the true particle ID that contributed the most energy to the CRT data product
    int TrueIdFromTotalEnergy(const art::Event& event, const crt::CRTData& data);

    // Get the true particle ID that contributed the most energy to the CRT hit
    int TrueIdFromTotalEnergy(const art::Event& event, const crt::CRTHit& hit);

    // Get the true particle ID that contributed the most energy to the CRT track
    int TrueIdFromTotalEnergy(const art::Event& event, const crt::CRTTrack& track);

  private:

    art::InputTag fCRTDataLabel;
    art::InputTag fCRTHitLabel;
    art::InputTag fCRTTrackLabel;

  };

}

#endif

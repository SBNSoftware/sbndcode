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
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// Utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"

#include "lardataobj/Simulation/AuxDetSimChannel.h"

#include "sbnobj/SBND/CRT/CRTData.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"

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
      fhicl::Atom<bool> RollupUnsavedIds {
        Name("RollupUnsavedIds")
      };

    };

    CRTBackTracker(const Config& config);

    CRTBackTracker(const fhicl::ParameterSet& pset) :
      CRTBackTracker(fhicl::Table<Config>(pset, {})()) {}

    CRTBackTracker();

    ~CRTBackTracker();

    void reconfigure(const Config& config);

    // Initialize to speed things up
    void Initialize(const art::Event& event);

    // Check that two CRT data products are the same
    bool DataCompare(const sbnd::crt::CRTData& data1, const sbnd::crt::CRTData& data2);

    // Check that two CRT hits are the same
    bool HitCompare(const sbn::crt::CRTHit& hit1, const sbn::crt::CRTHit& hit2);

    // Check that two CRT tracks are the same
    bool TrackCompare(const sbn::crt::CRTTrack& track1, const sbn::crt::CRTTrack& track2);

    // Get all the true particle IDs that contributed to the CRT data product
    std::vector<int> AllTrueIds(const art::Event& event, const sbnd::crt::CRTData& data);

    // Get all the true particle IDs that contributed to the CRT hit
    std::vector<int> AllTrueIds(const art::Event& event, const sbn::crt::CRTHit& hit);

    // Get all the true particle IDs that contributed to the CRT track
    std::vector<int> AllTrueIds(const art::Event& event, const sbn::crt::CRTTrack& track);

    // Get the true particle ID that contributed the most energy to the CRT data product
    int TrueIdFromTotalEnergy(const art::Event& event, const sbnd::crt::CRTData& data);
    // Faster function - needs Initialize() to be called first
    int TrueIdFromDataId(const art::Event& event, int data_i);

    // Get the true particle ID that contributed the most energy to the CRT hit
    int TrueIdFromTotalEnergy(const art::Event& event, const sbn::crt::CRTHit& hit);
    // Faster function - needs Initialize() to be called first
    int TrueIdFromHitId(const art::Event& event, int hit_i);

    // Get the true particle ID that contributed the most energy to the CRT track
    int TrueIdFromTotalEnergy(const art::Event& event, const sbn::crt::CRTTrack& track);
    // Faster function - needs Initialize() to be called first
    int TrueIdFromTrackId(const art::Event& event, int track_i);

  private:

    art::InputTag fCRTDataLabel;
    art::InputTag fCRTHitLabel;
    art::InputTag fCRTTrackLabel;

    bool fRollupUnsavedIds;

    std::map<int, std::map<int, double>> fDataTrueIds;
    std::map<int, std::map<int, double>> fHitTrueIds;
    std::map<int, std::map<int, double>> fTrackTrueIds;

  };

}

#endif

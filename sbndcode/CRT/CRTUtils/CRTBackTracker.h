#ifndef CRTBACKTRACKER_H_SEEN
#define CRTBACKTRACKER_H_SEEN


///////////////////////////////////////////////////////////////////////////
// CRTBackTracker.h
//
// Quick and dirty backtracker for SBND CRT
// T Brooks (tbrooks@fnal.gov), November 2018
//
// Modified by Jiaoyang Li/李 娇瑒 (jiaoyang.li@ed.ac.uk), May 2023
///////////////////////////////////////////////////////////////////////////

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

// sbnobj
#include "sbnobj/SBND/CRT/CRTData.hh"
#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/FEBTruthInfo.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"

// larsim
#include "larsim/MCCheater/ParticleInventoryService.h"

// lardataobj
#include "lardataobj/Simulation/AuxDetSimChannel.h" 
#include "lardataobj/Simulation/ParticleAncestryMap.h"
// c++
#include <vector>


namespace sbnd{

  class CRTBackTracker {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<art::InputTag> FEBDataLabel {
        Name("FEBDataLabel")
      };
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
      fhicl::Atom<art::InputTag> SimModuleLabel {
        Name("SimModuleLabel")
      };

    };

    struct TruthMatchMetrics {
      int           trackid;
      int           pdg;
      double        completeness;
      double        purity;
      double        hitcompleteness;
      double        hitpurity;
      double        particle_energy;
      double        depEnergy_total;

      TruthMatchMetrics(int _trackid, double _completeness, double _purity,
                        double _hitcompleteness, double _hitpurity, int _pdg,
                        double _particle_energy, double _depEnergy_total)
      {
        trackid                = _trackid;
        completeness           = _completeness;
        purity                 = _purity;
        hitcompleteness        = _hitcompleteness;
        hitpurity              = _hitpurity;
        pdg                    = _pdg;
        particle_energy        = _particle_energy;
        depEnergy_total        = _depEnergy_total;
      }
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
    bool FEBDataCompare(const sbnd::crt::FEBData& data1, const sbnd::crt::FEBData& data2);

    // Check that two CRT hits are the same
    bool HitCompare(const sbn::crt::CRTHit& hit1, const sbn::crt::CRTHit& hit2);

    // Check that two CRT tracks are the same
    bool TrackCompare(const sbn::crt::CRTTrack& track1, const sbn::crt::CRTTrack& track2);

    // Get all the true particle IDs that contributed to the CRT data product
    std::vector<int> AllTrueIds(const art::Event& event, const sbnd::crt::FEBData& data);

    // Get all the true particle IDs that contributed to the CRT hit
    std::vector<int> AllTrueIds(const art::Event& event, const sbn::crt::CRTHit& hit);

    // Get all the true particle IDs that contributed to the CRT track
    std::vector<int> AllTrueIds(const art::Event& event, const art::Ptr<sbn::crt::CRTTrack> &track);

    // Get the true particle ID that contributed the most energy to the CRT data product
    int TrueIdFromTotalEnergy(const art::Event& event, const sbnd::crt::FEBData& data);
    // Faster function - needs Initialize() to be called first
    int TrueIdFromFEBDataId(const art::Event& event, int data_i);

    // Get the true particle ID that contributed the most energy to the CRT hit
    int TrueIdFromTotalEnergy(const art::Event& event, const sbn::crt::CRTHit& hit);
    // Faster function - needs Initialize() to be called first
    int TrueIdFromHitId(const art::Event& event, int hit_i);

    // Get the true particle ID that contributed the most energy to the CRT track
    int TrueIdFromTotalEnergy(const art::Event& event, const sbn::crt::CRTTrack& track);

    // Faster function - needs Initialize() to be called first
    int TrueIdFromTrackId(const art::Event& event, int track_i);



    // new functions: 
    TruthMatchMetrics TruthMatrixFromTotalEnergy(const art::Event& event, const art::Ptr<sbn::crt::CRTTrack> &crtTrack);//const sbn::crt::CRTTrack& track);
    TruthMatchMetrics TruthMatrixFromTotalEnergy(const art::Event& event, const art::Ptr<sbn::crt::CRTHit> &crtHit);//const sbn::crt::CRTTrack& track);
    void TrueParticlePDGEnergyTime(const int trackID, int &pdg, double &energy, double &time);
    int  RollUpID(const int &id);

    
  private:

    art::ServiceHandle<cheat::ParticleInventoryService> particleInv;

    art::InputTag fFEBDataLabel;
    art::InputTag fCRTDataLabel;
    art::InputTag fCRTHitLabel;
    art::InputTag fCRTTrackLabel;

    art::InputTag fSimModuleLabel;

    bool fRollupUnsavedIds;

    std::map<int, std::map<int, double>> fCRTDataTrueIds;
    std::map<int, std::map<int, double>> fFEBDataTrueIds;
    std::map<int, std::map<int, double>> fHitTrueIds;
    std::map<int, std::map<int, double>> fTrackTrueIds;
    std::map<int, int>                   fTrackIDMotherMap;


  };

}

#endif
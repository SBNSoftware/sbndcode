#ifndef CRTBACKTRACKERALG_H_SEEN
#define CRTBACKTRACKERALG_H_SEEN

///////////////////////////////////////////////
// CRTBackTrackerAlg.h
//
// Truth Matching Utilities for CRT analysis
// Henry Lay (h.lay@lancaster.ac.uk)
// November 2022
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// Utility libraries
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

// sbnobj
#include "sbnobj/SBND/CRT/FEBData.hh"
#include "sbnobj/SBND/CRT/FEBTruthInfo.hh"
#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"

// sbndcode
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

// larsim
#include "larsim/MCCheater/ParticleInventoryService.h"

// larcorealg
#include "larcorealg/CoreUtils/enumerate.h"

// lardataobj
#include "lardataobj/Simulation/ParticleAncestryMap.h"

namespace sbnd::crt {
  
  class CRTBackTrackerAlg {
  public:
    
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      
      fhicl::Atom<art::InputTag> SimModuleLabel {
        Name("SimModuleLabel"),
          };
      fhicl::Atom<art::InputTag> SimDepositModuleLabel {
        Name("SimDepositModuleLabel"),
          };
      fhicl::Atom<art::InputTag> FEBDataModuleLabel {
        Name("FEBDataModuleLabel"),
          };
      fhicl::Atom<art::InputTag> StripHitModuleLabel {
        Name("StripHitModuleLabel"),
          };
      fhicl::Atom<art::InputTag> ClusterModuleLabel {
        Name("ClusterModuleLabel"),
          };
      fhicl::Atom<art::InputTag> SpacePointModuleLabel {
        Name("SpacePointModuleLabel"),
          };
    };

    struct TrueDeposit {
      int       trackid;
      int       pdg;
      CRTTagger tagger;
      double    x;
      double    y;
      double    z;
      double    energy;
      double    time;
      double    coreX;
      double    coreY;
      double    coreZ;
      double    coreEnergy;
      double    coreTime;
      
      TrueDeposit(int _trackid = -999999, int _pdg = -999999, CRTTagger _tagger = kUndefinedTagger, 
                  double _x = -999999., double _y = -999999., double _z = -999999., 
                  double _energy = -999999., double _time = -999999.,
                  double _coreX = -999999., double _coreY = -999999., double _coreZ = -999999.,
                  double _coreEnergy = -999999., double _coreTime = -999999.)
      {
        trackid    = _trackid;
        pdg        = _pdg;
        tagger     = _tagger;
        x          = _x;
        y          = _y;
        z          = _z;
        energy     = _energy;
        time       = _time;
        coreX      = _coreX;
        coreY      = _coreY;
        coreZ      = _coreZ;
        coreEnergy = _coreEnergy;
        coreTime   = _coreTime;
      }
    };

    struct TruthMatchMetrics {
      int         trackid;
      double      completeness;
      double      purity;
      double      hitcompleteness;
      double      hitpurity;
      TrueDeposit deposit;

      TruthMatchMetrics(int _trackid, double _completeness, double _purity,
                        double _hitcompleteness, double _hitpurity, TrueDeposit _deposit)
      {
        trackid         = _trackid;
        completeness    = _completeness;
        purity          = _purity;
        hitcompleteness = _hitcompleteness;
        hitpurity       = _hitpurity;
        deposit         = _deposit;
      }
    };

    CRTBackTrackerAlg(const Config& config);
    
  CRTBackTrackerAlg(const fhicl::ParameterSet& pset) :
    CRTBackTrackerAlg(fhicl::Table<Config>(pset, {})()) {}
    
    CRTBackTrackerAlg();

    ~CRTBackTrackerAlg();

    void reconfigure(const Config& config);

    void SetupMaps(const art::Event &event);

    int RollUpID(const int &id);

    void RunRecoStatusChecks(const art::Event &event);
    
    std::map<std::pair<int, CRTTagger>, bool> GetRecoStatusMap();

    TrueDeposit GetTrueDeposit(std::pair<int, CRTTagger> category);

    TruthMatchMetrics TruthMatching(const art::Event &event, const art::Ptr<CRTStripHit> &stripHit);

    TruthMatchMetrics TruthMatching(const art::Event &event, const art::Ptr<CRTCluster> &cluster);

  private:
    
    CRTGeoAlg fCRTGeoAlg;
    art::ServiceHandle<cheat::ParticleInventoryService> particleInv;

    art::InputTag fSimModuleLabel;
    art::InputTag fSimDepositModuleLabel;
    art::InputTag fFEBDataModuleLabel;
    art::InputTag fStripHitModuleLabel;
    art::InputTag fClusterModuleLabel;
    art::InputTag fSpacePointModuleLabel;

    std::map<std::pair<int, CRTTagger>, TrueDeposit> fTrueDepositsMap;
    std::map<std::pair<int, CRTTagger>, bool>        fTrackIDRecoMap;

    std::map<std::pair<int, CRTTagger>, double> fMCPIDEsEnergyMap;
    std::map<std::pair<int, CRTTagger>, int>    fMCPStripHitsMap;

    std::map<int, int> fTrackIDMotherMap;
    std::map<int, int> fStripHitMCPMap;
    
  };
}

#endif

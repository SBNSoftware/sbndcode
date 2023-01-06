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

// sbndcode
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

#include "larsim/MCCheater/ParticleInventoryService.h"

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
    };

    struct TruthMatchMetrics {
      int    trackid;
      double completeness;
      double purity;
      double hitcompleteness;
      double hitpurity;
      double x;
      double y;
      double z;
      double energy;
      double time;

      TruthMatchMetrics(int _trackid, double _completeness, double _purity,
                        double _hitcompleteness, double _hitpurity,
                        double _x, double _y, double _z,
                        double _energy, double _time)
      {
        trackid         = _trackid;
        completeness    = _completeness;
        purity          = _purity;
        hitcompleteness = _hitcompleteness;
        hitpurity       = _hitpurity;
        x               = _x;
        y               = _y;
        z               = _z;
        energy          = _energy;
        time            = _time;
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

    TruthMatchMetrics TruthMatching(const art::Event &event, const art::Ptr<CRTStripHit> &stripHit);

    TruthMatchMetrics TruthMatching(const art::Event &event, const art::Ptr<CRTCluster> &cluster);

  private:
    
    CRTGeoAlg fCRTGeoAlg;
    art::ServiceHandle<cheat::ParticleInventoryService> particleInv;

    art::InputTag fSimModuleLabel;
    art::InputTag fSimDepositModuleLabel;
    art::InputTag fStripHitModuleLabel;
    art::InputTag fFEBDataModuleLabel;
    art::InputTag fClusterModuleLabel;

    std::map<std::pair<int, CRTTagger>, double> fMCPIDEsEnergyMap;
    std::map<std::pair<int, CRTTagger>, int>    fMCPStripHitsMap;

    std::map<int, int> fTrackIDMotherMap;
    std::map<int, int> fStripHitMCPMap;
  };
}

#endif

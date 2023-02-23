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
#include "canvas/Persistency/Common/FindOneP.h"

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
#include "sbnobj/SBND/CRT/CRTTrack.hh"

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

      fhicl::Atom<art::InputTag> TrackModuleLabel {
        Name("TrackModuleLabel"),
          };
    };

    struct TrueDeposit {
      int       trackid;
      int       pdg;
      CRTTagger tagger;
      double    energy;
      double    time;
      double    x;
      double    y;
      double    z;
      double    coreEnergy;
      double    coreTime;
      double    coreX;
      double    coreY;
      double    coreZ;

      TrueDeposit(int _trackid = -999999, int _pdg = -999999, CRTTagger _tagger = kUndefinedTagger, 
                  double _energy = -999999., double _time = -999999.,
                  double _x = -999999., double _y = -999999., double _z = -999999., 
                  double _coreEnergy = -999999., double _coreTime = -999999.,
                  double _coreX = -999999., double _coreY = -999999., double _coreZ = -999999.)
      {
        trackid    = _trackid;
        pdg        = _pdg;
        tagger     = _tagger;
        energy     = _energy;
        time       = _time;
        x          = _x;
        y          = _y;
        z          = _z;
        coreEnergy = _coreEnergy;
        coreTime   = _coreTime;
        coreX      = _coreX;
        coreY      = _coreY;
        coreZ      = _coreZ;
      }
    };

    struct TrueTrackInfo {
      int trackid;
      int pdg;
      double startx;
      double starty;
      double startz;
      double endx;
      double endy;
      double endz;
      double energy;
      double tof;
      bool found_ends;
      TrueDeposit deposit1;
      TrueDeposit deposit2;

      TrueTrackInfo()
      {
        trackid    = -std::numeric_limits<int>::max();
        pdg        = -std::numeric_limits<int>::max();
        startx     = -std::numeric_limits<double>::max();
        starty     = -std::numeric_limits<double>::max();
        startz     = -std::numeric_limits<double>::max();
        endx       = -std::numeric_limits<double>::max();
        endy       = -std::numeric_limits<double>::max();
        endz       = -std::numeric_limits<double>::max();
        energy     = -std::numeric_limits<double>::max();
        tof        = -std::numeric_limits<double>::max();
        found_ends = false;
        deposit1   = TrueDeposit();
        deposit2   = TrueDeposit();
      }

      TrueTrackInfo(int _trackid, int _pdg, geo::Point_t _start, geo::Point_t _end,
                    double _energy, double _tof, bool _found_ends, TrueDeposit &_deposit1,
                    TrueDeposit &_deposit2)
      {
        trackid    = _trackid;
        pdg        = _pdg;
        startx     = _start.X();
        starty     = _start.Y();
        startz     = _start.Z();
        endx       = _end.X();
        endy       = _end.Y();
        endz       = _end.Z();
        energy     = _energy;
        tof        = _tof;
        found_ends = _found_ends;
        deposit1   = _deposit1;
        deposit2   = _deposit2;
      }
    };

    struct TruthMatchMetrics {
      int           trackid;
      double        completeness;
      double        purity;
      double        hitcompleteness;
      double        hitpurity;
      TrueDeposit   deposit;
      TrueTrackInfo trackinfo;

      TruthMatchMetrics(int _trackid, double _completeness, double _purity,
                        double _hitcompleteness, double _hitpurity, TrueDeposit _deposit,
                        TrueTrackInfo _trackinfo = TrueTrackInfo())
      {
        trackid         = _trackid;
        completeness    = _completeness;
        purity          = _purity;
        hitcompleteness = _hitcompleteness;
        hitpurity       = _hitpurity;
        deposit         = _deposit;
        trackinfo       = _trackinfo;
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
    
    std::map<std::pair<int, CRTTagger>, bool> GetSpacePointRecoStatusMap();

    TrueDeposit GetTrueDeposit(std::pair<int, CRTTagger> category);

    TruthMatchMetrics TruthMatching(const art::Event &event, const art::Ptr<CRTStripHit> &stripHit);

    TruthMatchMetrics TruthMatching(const art::Event &event, const art::Ptr<CRTCluster> &cluster);

    TruthMatchMetrics TruthMatching(const art::Event &event, const art::Ptr<CRTTrack> &track);

    std::pair<double, geo::Point_t> LineTaggerIntersectionPoint(const geo::Point_t &start, 
                                                                const geo::Vector_t &dir, 
                                                                const CRTTagger &tagger);

  private:
    
    CRTGeoAlg fCRTGeoAlg;
    art::ServiceHandle<cheat::ParticleInventoryService> particleInv;

    art::InputTag fSimModuleLabel;
    art::InputTag fSimDepositModuleLabel;
    art::InputTag fFEBDataModuleLabel;
    art::InputTag fStripHitModuleLabel;
    art::InputTag fClusterModuleLabel;
    art::InputTag fSpacePointModuleLabel;
    art::InputTag fTrackModuleLabel;

    std::map<std::pair<int, CRTTagger>, TrueDeposit> fTrueDepositsPerTaggerMap;
    std::map<int, TrueDeposit>                       fTrueDepositsMap;
    std::map<int, TrueTrackInfo>                     fTrueTrackInfosMap;
    std::map<std::pair<int, CRTTagger>, bool>        fTrackIDSpacePointRecoMap;
    std::map<std::pair<int, CRTTagger>, double>      fMCPIDEsEnergyPerTaggerMap;
    std::map<int, double>                            fMCPIDEsEnergyMap;
    std::map<std::pair<int, CRTTagger>, int>         fMCPStripHitsMap;
    std::map<int, int>                               fTrackIDMotherMap;
    std::map<int, int>                               fStripHitMCPMap;
    
  };
}

#endif

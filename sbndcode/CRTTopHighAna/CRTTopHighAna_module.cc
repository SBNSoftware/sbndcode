////////////////////////////////////////////////////////////////////////
// Class:       CRTTopHighAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        CRTTopHighAna_module.cc
//
// Generated at Thu May 25 09:40:08 2023 by Henry Lay using cetskelgen
// from  version .
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
#include "TTree.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "sbnobj/SBND/CRT/CRTData.hh"
#include "sbnobj/SBND/CRT/FEBData.hh"
#include "lardataobj/AnalysisBase/T0.h"
#include "sbnobj/SBND/CRT/FEBTruthInfo.hh"
#include <chrono>
#include "lardataobj/Simulation/ParticleAncestryMap.h"

class CRTTopHighAna;


class CRTTopHighAna : public art::EDAnalyzer {
public:
  explicit CRTTopHighAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTTopHighAna(CRTTopHighAna const&) = delete;
  CRTTopHighAna(CRTTopHighAna&&) = delete;
  CRTTopHighAna& operator=(CRTTopHighAna const&) = delete;
  CRTTopHighAna& operator=(CRTTopHighAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void ResetVars();

  bool IntersectsTagger(const std::string name, const art::Ptr<simb::MCParticle> &particle);

  bool IntersectsModules(const std::string name, const art::Ptr<simb::MCParticle> &particle);
  
  bool IntersectsStrips(const std::string name, const art::Ptr<simb::MCParticle> &particle);
  
  bool PassesThroughVolume(const double rmin[3], const double rmax[3], const art::Ptr<simb::MCParticle> &particle);

  bool IsInsideVolume(const double rmin[3], const double rmax[3], const geo::Point_t &point);

  bool IsInsideVolume(const double rmin[3], const double rmax[3], const TVector3 &point);

  double RelativeDistanceToIntersect(const double plane, const double dir, const double start);

  bool PassesThroughVolume(const double rmin[3], const double rmax[3], const geo::Vector_t dir, const geo::Point_t start);

  bool IsGoodMatch(const art::Event &e, const int trackid, const art::Ptr<sbn::crt::CRTHit> &crthit);

  bool IsGoodMatch(const art::Event &e, const int trackid, const art::Ptr<sbn::crt::CRTTrack> &crttrack);

private:

  sbnd::CRTGeoAlg fCRTGeoAlg;

  std::set<int> fTrackRecoSet, fTrackRecoLongSet, fHitMatchSet, fGoodHitMatchSet, fTrackMatchSet, fGoodTrackMatchSet;
  std::map<int, int> fHitMatchMap;
  std::map<int, std::set<int>> fTrackMatchMap;
  art::Handle<sim::ParticleAncestryMap> fDroppedTrackIDMapHandle;

  std::string fMCTruthLabel, fMCParticleLabel, fPFPLabel, fTPCTrackLabel, fHitMatchLabel,
    fTrackMatchLabel, fCRTHitLabel, fCRTTrackLabel, fFEBDataLabel, fCRTDataLabel;

  TTree* fTree;
  
  int pdg;
  float energy, angleToVertical, startAngleToVertical;
  bool primary;

  bool startContained, endContained, intersectsTPC, enters, exits, entersOrExits,
    startCRTContained, endCRTContained, intersectsCRT, entersCRT, exitsCRT, entersOrExitsCRT,
    intersectsBottom, intersectsBottomStrips, intersectsEast, intersectsEastStrips,
    intersectsWest, intersectsWestStrips, intersectsNorth, intersectsNorthStrips, intersectsSouth, intersectsSouthStrips,
    intersectsTopLow, intersectsTopLowStrips, intersectsTopHigh, intersectsTopHighStrips, intersectsAny, intersectsAnyStrips,
    intersectsAnyNoBottom, intersectsAnyStripsNoBottom, intersectsExtendedTopLow, intersectsRaisedTopLow, intersectsExtendedRaisedTopLow,
    intersectsExtendedRaisedTopLowFullNorth, intersectsExtendedRaisedTopLowFullSouth,
    hasTPCTrack, hasTPCLongTrack, hasHitMatch, hasTrackMatch, hasGoodHitMatch, hasGoodTrackMatch;

  int hitMatchTagger, trackMatchTagger1, trackMatchTagger2, trackMatchTagger3;
};


CRTTopHighAna::CRTTopHighAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fMCTruthLabel    (p.get<std::string>("MCTruthLabel", "generator"))
  , fMCParticleLabel (p.get<std::string>("MCParticleLabel", "largeant"))
  , fPFPLabel        (p.get<std::string>("PFPLabel", "pandoraSCE"))
  , fTPCTrackLabel   (p.get<std::string>("TPCTrackLabel", "pandoraSCETrack"))
  , fHitMatchLabel   (p.get<std::string>("HitMatchLabel", "crthitt0SCE"))
  , fTrackMatchLabel (p.get<std::string>("TrackMatchLabel", "crttrackt0SCE"))
  , fCRTHitLabel     (p.get<std::string>("CRTHitLabel", "crthit"))
  , fCRTTrackLabel   (p.get<std::string>("CRTTrackLabel", "crttrack"))
  , fFEBDataLabel    (p.get<std::string>("FEBDataLabel", "crtsim"))
  , fCRTDataLabel    (p.get<std::string>("CRTDataLabel", "crt"))
  {
    art::ServiceHandle<art::TFileService> fs;

    fTree = fs->make<TTree>("tree", "");

    fTree->Branch("pdg", &pdg);
    fTree->Branch("energy", &energy);
    fTree->Branch("primary", &primary);
    fTree->Branch("angleToVertical", &angleToVertical);
    fTree->Branch("startAngleToVertical", &startAngleToVertical);

    fTree->Branch("startContained", &startContained);
    fTree->Branch("endContained", &endContained);
    fTree->Branch("intersectsTPC", &intersectsTPC);
    fTree->Branch("enters", &enters);
    fTree->Branch("exits", &exits);
    fTree->Branch("entersOrExits", &entersOrExits);

    fTree->Branch("startCRTContained", &startCRTContained);
    fTree->Branch("endCRTContained", &endCRTContained);
    fTree->Branch("intersectsCRT", &intersectsCRT);
    fTree->Branch("entersCRT", &entersCRT);
    fTree->Branch("exitsCRT", &exitsCRT);
    fTree->Branch("entersOrExitsCRT", &entersOrExitsCRT);

    fTree->Branch("intersectsBottom", &intersectsBottom);
    fTree->Branch("intersectsBottomStrips", &intersectsBottomStrips);
    fTree->Branch("intersectsEast", &intersectsEast);
    fTree->Branch("intersectsEastStrips", &intersectsEastStrips);
    fTree->Branch("intersectsWest", &intersectsWest);
    fTree->Branch("intersectsWestStrips", &intersectsWestStrips);
    fTree->Branch("intersectsNorth", &intersectsNorth);
    fTree->Branch("intersectsNorthStrips", &intersectsNorthStrips);
    fTree->Branch("intersectsSouth", &intersectsSouth);
    fTree->Branch("intersectsSouthStrips", &intersectsSouthStrips);
    fTree->Branch("intersectsTopLow", &intersectsTopLow);
    fTree->Branch("intersectsTopLowStrips", &intersectsTopLowStrips);
    fTree->Branch("intersectsTopHigh", &intersectsTopHigh);
    fTree->Branch("intersectsTopHighStrips", &intersectsTopHighStrips);
    fTree->Branch("intersectsAny", &intersectsAny);
    fTree->Branch("intersectsAnyNoBottom", &intersectsAnyNoBottom);
    fTree->Branch("intersectsAnyStrips", &intersectsAnyStrips);
    fTree->Branch("intersectsAnyStripsNoBottom", &intersectsAnyStripsNoBottom);

    fTree->Branch("intersectsExtendedTopLow", &intersectsExtendedTopLow);
    fTree->Branch("intersectsRaisedTopLow", &intersectsRaisedTopLow);
    fTree->Branch("intersectsExtendedRaisedTopLow", &intersectsExtendedRaisedTopLow);
    fTree->Branch("intersectsExtendedRaisedTopLowFullNorth", &intersectsExtendedRaisedTopLowFullNorth);
    fTree->Branch("intersectsExtendedRaisedTopLowFullSouth", &intersectsExtendedRaisedTopLowFullSouth);

    fTree->Branch("hasTPCTrack", &hasTPCTrack);
    fTree->Branch("hasTPCLongTrack", &hasTPCLongTrack);
    fTree->Branch("hasHitMatch", &hasHitMatch);
    fTree->Branch("hasTrackMatch", &hasTrackMatch);
    fTree->Branch("hasGoodHitMatch", &hasGoodHitMatch);
    fTree->Branch("hasGoodTrackMatch", &hasGoodTrackMatch);
    fTree->Branch("hitMatchTagger", &hitMatchTagger);
    fTree->Branch("trackMatchTagger1", &trackMatchTagger1);
    fTree->Branch("trackMatchTagger2", &trackMatchTagger2);
    fTree->Branch("trackMatchTagger3", &trackMatchTagger3);
  }

void CRTTopHighAna::analyze(art::Event const& e)
{
  fTrackRecoSet.clear(); fTrackRecoLongSet.clear(); fHitMatchSet.clear(); fGoodHitMatchSet.clear(); fTrackMatchSet.clear();
  fGoodTrackMatchSet.clear(); fHitMatchMap.clear(); fTrackMatchMap.clear();

  fDroppedTrackIDMapHandle.clear();

  art::Handle<std::vector<simb::MCTruth>> MCTruthHandle;
  e.getByLabel(fMCTruthLabel, MCTruthHandle);
  if(!MCTruthHandle.isValid()){
    std::cout << "MCTruth product " << fMCTruthLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
  art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

  art::FindManyP<simb::MCParticle> MCTruthToMCParticles(MCTruthHandle, e, fMCParticleLabel);

  art::Handle<std::vector<recob::PFParticle>> PFPHandle;
  e.getByLabel(fPFPLabel, PFPHandle);
  if(!PFPHandle.isValid()){
    std::cout << "PFP product " << fPFPLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<recob::PFParticle>> PFPVec;
  art::fill_ptr_vector(PFPVec, PFPHandle);

  art::Handle<std::vector<recob::Track>> TPCTrackHandle;
  e.getByLabel(fTPCTrackLabel, TPCTrackHandle);
  if(!TPCTrackHandle.isValid()){
    std::cout << "TPCTrack product " << fTPCTrackLabel << " not found..." << std::endl;
    throw std::exception();
  }

  art::Handle<std::vector<anab::T0>> HitMatchHandle;
  e.getByLabel(fHitMatchLabel, HitMatchHandle);
  if(!HitMatchHandle.isValid()){
    std::cout << "HitMatch product " << fHitMatchLabel << " not found..." << std::endl;
    throw std::exception();
  }

  art::Handle<std::vector<anab::T0>> TrackMatchHandle;
  e.getByLabel(fTrackMatchLabel, TrackMatchHandle);
  if(!TrackMatchHandle.isValid()){
    std::cout << "TrackMatch product " << fTrackMatchLabel << " not found..." << std::endl;
    throw std::exception();
  }

  art::Handle<std::vector<sbn::crt::CRTTrack>> CRTTrackHandle;
  e.getByLabel(fCRTTrackLabel, CRTTrackHandle);
  if(!CRTTrackHandle.isValid()){
    std::cout << "CRTTrack product " << fCRTTrackLabel << " not found..." << std::endl;
    throw std::exception();
  }

  e.getByLabel(fMCParticleLabel, fDroppedTrackIDMapHandle);
  if(!fDroppedTrackIDMapHandle.isValid()){
    std::cout << "fDroppedTrackIDMap product " << fMCParticleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  
  art::FindOneP<recob::Track> pfpToTrack(PFPHandle, e, fTPCTrackLabel);
  art::FindManyP<recob::Hit> tracksToHits(TPCTrackHandle, e, fTPCTrackLabel);
  art::FindOneP<anab::T0> trackToHitMatch(TPCTrackHandle, e, fHitMatchLabel);
  art::FindOneP<anab::T0> trackToTrackMatch(TPCTrackHandle, e, fTrackMatchLabel);
  art::FindOneP<sbn::crt::CRTHit> matchToHit(HitMatchHandle, e, fHitMatchLabel);
  art::FindOneP<sbn::crt::CRTTrack> matchToTrack(TrackMatchHandle, e, fTrackMatchLabel);
  art::FindManyP<sbn::crt::CRTHit> crtTracksToHits(CRTTrackHandle, e, fCRTTrackLabel);

  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  for(auto const& pfp : PFPVec)
    {
      if(pfp->PdgCode() != 13)
        continue;

      const art::Ptr<recob::Track> track = pfpToTrack.at(pfp.key());

      if(track.isNull())
        {
          std::cout << "Bad track pointer" << std::endl;
          //      throw std::exception();
          continue;
        }

      const std::vector<art::Ptr<recob::Hit>> hits = tracksToHits.at(track.key());

      int trackid = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, hits, true);
      int rolledtrackid = fDroppedTrackIDMapHandle->GetAncestor(trackid);
      if(fDroppedTrackIDMapHandle->Exists(rolledtrackid))
        trackid = rolledtrackid;
      fTrackRecoSet.insert(trackid);

      if(track->Length() > 5)
        fTrackRecoLongSet.insert(trackid);

      const art::Ptr<anab::T0> hitmatch   = trackToHitMatch.at(track.key());
      const art::Ptr<anab::T0> trackmatch = trackToTrackMatch.at(track.key());

      if(hitmatch.isNonnull())
        {
          const art::Ptr<sbn::crt::CRTHit> crthit = matchToHit.at(hitmatch.key());

          if(fHitMatchSet.count(trackid) != 0)
            std::cout << "Already got a hit match" << std::endl;

          fHitMatchSet.insert(trackid);
          fHitMatchMap[trackid] = crthit->plane;

          if(IsGoodMatch(e, trackid, crthit))
            fGoodHitMatchSet.insert(trackid);
        }

      if(trackmatch.isNonnull())
        {
          const art::Ptr<sbn::crt::CRTTrack> crttrack           = matchToTrack.at(trackmatch.key());
          const std::vector<art::Ptr<sbn::crt::CRTHit>> crthits = crtTracksToHits.at(crttrack.key());

          if(fTrackMatchSet.count(trackid) != 0)
            std::cout << "Already got a track match" << std::endl;

          fTrackMatchSet.insert(trackid);
          for(auto const& crthit : crthits)
            fTrackMatchMap[trackid].insert(crthit->plane);

          if(IsGoodMatch(e, trackid, crttrack))
            fGoodTrackMatchSet.insert(trackid);
        }
    }

  for(auto const& truth : MCTruthVec)
    {
      const std::vector<art::Ptr<simb::MCParticle>> particles = MCTruthToMCParticles.at(truth.key());

      for(auto const& particle : particles)
        {
          ResetVars();

          pdg     = particle->PdgCode();
          energy  = particle->E();
          primary = particle->Mother() == 0;

          const TVector3 start = particle->Position().Vect();
          const TVector3 end   = particle->EndPosition().Vect();
          angleToVertical      = TMath::RadToDeg() * (end - start).Angle(TVector3(0, -1, 0));

          const TVector3 startdir = particle->Momentum().Vect();
          startAngleToVertical = TMath::RadToDeg() * startdir.Angle(TVector3(0, -1, 0));

          const double tpcmin[3] = { -200, -200, 0 };
          const double tpcmax[3] = { 200, 200, 500 };
          startContained = IsInsideVolume(tpcmin, tpcmax, start);
          endContained   = IsInsideVolume(tpcmin, tpcmax, end);
          intersectsTPC  = PassesThroughVolume(tpcmin, tpcmax, particle);

          enters = intersectsTPC && !startContained;
          exits  = intersectsTPC && !endContained;

          entersOrExits = intersectsTPC && !(startContained && endContained);

          const std::vector<double> limits = fCRTGeoAlg.CRTLimitsStrips();
          const double crtmin[3] = {limits[0], limits[1], limits[2]};
          const double crtmax[3] = {limits[3], limits[4], limits[5]};
          startCRTContained = IsInsideVolume(crtmin, crtmax, start);
          endCRTContained   = IsInsideVolume(crtmin, crtmax, end);
          intersectsCRT     = PassesThroughVolume(crtmin, crtmax, particle);

          entersCRT = intersectsCRT && !startCRTContained;
          exitsCRT  = intersectsCRT && !endCRTContained;

          entersOrExitsCRT = intersectsCRT && !(startCRTContained && endCRTContained);

          intersectsBottom        = IntersectsTagger("volTaggerBot_0", particle);
          intersectsBottomStrips  = IntersectsStrips("volTaggerBot_0", particle);
          intersectsEast          = IntersectsTagger("volTaggerEast_0", particle);
          intersectsEastStrips    = IntersectsStrips("volTaggerEast_0", particle);
          intersectsWest          = IntersectsTagger("volTaggerWest_0", particle);
          intersectsWestStrips    = IntersectsStrips("volTaggerWest_0", particle);
          intersectsNorth         = IntersectsTagger("volTaggerNorth_0", particle);
          intersectsNorthStrips   = IntersectsStrips("volTaggerNorth_0", particle);
          intersectsSouth         = IntersectsTagger("volTaggerSouth_0", particle);
          intersectsSouthStrips   = IntersectsStrips("volTaggerSouth_0", particle);
          intersectsTopLow        = IntersectsTagger("volTaggerTopLow_0", particle);
          intersectsTopLowStrips  = IntersectsStrips("volTaggerTopLow_0", particle);
          intersectsTopHigh       = IntersectsTagger("volTaggerTopHigh_0", particle);
          intersectsTopHighStrips = IntersectsStrips("volTaggerTopHigh_0", particle);

          intersectsAnyNoBottom = intersectsEast || intersectsWest || intersectsNorth
            || intersectsSouth || intersectsTopLow || intersectsTopHigh;

          intersectsAnyStripsNoBottom = intersectsEastStrips || intersectsWestStrips || intersectsNorthStrips
            || intersectsSouthStrips || intersectsTopLowStrips || intersectsTopHighStrips;

          intersectsAny = intersectsAnyNoBottom || intersectsBottom;

          intersectsAnyStrips = intersectsAnyStripsNoBottom || intersectsBottomStrips;

          const sbnd::CRTTaggerGeo toplow = fCRTGeoAlg.GetTagger("volTaggerTopLow_0");

          const double toplowextmin[3] = {toplow.minX - 45.5, toplow.minY, toplow.minZ - 74.5};
          const double toplowextmax[3] = {toplow.maxX + 45.5, toplow.maxY, toplow.maxZ + 74.5};
          intersectsExtendedTopLow = PassesThroughVolume(toplowextmin, toplowextmax, particle);

          const double toplowraismin[3] = {toplow.minX, toplow.minY + 42, toplow.minZ};
          const double toplowraismax[3] = {toplow.maxX, toplow.maxY + 42, toplow.maxZ};
          intersectsRaisedTopLow = PassesThroughVolume(toplowraismin, toplowraismax, particle);

          const double toplowextraismin[3] = {toplow.minX - 45.5, toplow.minY + 42, toplow.minZ - 74.5};
          const double toplowextraismax[3] = {toplow.maxX + 45.5, toplow.maxY + 42, toplow.maxZ + 74.5};
          intersectsExtendedRaisedTopLow = PassesThroughVolume(toplowextraismin, toplowextraismax, particle);

          const double toplowextraisfullnorthmin[3] = {toplow.minX - 45.5, toplow.minY + 42, toplow.minZ};
          const double toplowextraisfullnorthmax[3] = {toplow.maxX + 45.5, toplow.maxY + 42, toplow.maxZ + 149};
          intersectsExtendedRaisedTopLowFullNorth = PassesThroughVolume(toplowextraisfullnorthmin, toplowextraisfullnorthmax, particle);

          const double toplowextraisfullsouthmin[3] = {toplow.minX - 45.5, toplow.minY + 42, toplow.minZ - 149};
          const double toplowextraisfullsouthmax[3] = {toplow.maxX + 45.5, toplow.maxY + 42, toplow.maxZ};
          intersectsExtendedRaisedTopLowFullSouth = PassesThroughVolume(toplowextraisfullsouthmin, toplowextraisfullsouthmax, particle);

          int trackid = particle->TrackId();
          int rolledtrackid = fDroppedTrackIDMapHandle->GetAncestor(trackid);
          if(fDroppedTrackIDMapHandle->Exists(rolledtrackid))
            trackid = rolledtrackid;

          hasTPCTrack       = fTrackRecoSet.count(trackid) != 0;
          hasTPCLongTrack   = fTrackRecoLongSet.count(trackid) != 0;
          hasHitMatch       = fHitMatchSet.count(trackid) != 0;
          hasTrackMatch     = fTrackMatchSet.count(trackid) != 0;
          hasGoodHitMatch   = fGoodHitMatchSet.count(trackid) != 0;
          hasGoodTrackMatch = fGoodTrackMatchSet.count(trackid) != 0;

          if(hasHitMatch)
            hitMatchTagger = fHitMatchMap.at(trackid);

          if(hasTrackMatch)
            {
              auto const& taggers = fTrackMatchMap.at(trackid);
              std::vector<int> taggersVec(taggers.begin(), taggers.end());

              if(taggers.size() != 2 && taggers.size() != 3)
                {
                  std::cout << "Number of taggers used for track = " << taggers.size() << std::endl;
                  std::cout << "(";
                  for(auto const& tagger : taggers)
                    std::cout << tagger << ", ";
                  std::cout << ")" << std::endl;
                  //              throw std::exception();
                }
              else if(taggers.size() == 2)
                {
                  trackMatchTagger1 = taggersVec[0];
                  trackMatchTagger2 = taggersVec[1];
                }
              else if(taggers.size() > 2)
                {
                  trackMatchTagger1 = taggersVec[0];
                  trackMatchTagger2 = taggersVec[1];
                  trackMatchTagger3 = taggersVec[2];
                }
            }

          fTree->Fill();
        }
    }
}

void CRTTopHighAna::ResetVars()
{
  pdg = -1; energy = -1.; angleToVertical = -400.;

  intersectsTPC = false; startContained = false; endContained     = false;
  enters        = false; exits          = false; entersOrExits    = false;
  entersCRT     = false; exitsCRT       = false; entersOrExitsCRT = false;

  startCRTContained = false; endCRTContained = false; entersCRT = false; exitsCRT = false; entersOrExitsCRT = false;

  intersectsBottom      = false; intersectsBottomStrips      = false; intersectsEast   = false; intersectsEastStrips   = false;
  intersectsWest        = false; intersectsWestStrips        = false; intersectsNorth  = false; intersectsNorthStrips  = false;
  intersectsSouth       = false; intersectsSouthStrips       = false; intersectsTopLow = false; intersectsTopLowStrips = false;
  intersectsTopHigh     = false; intersectsTopHighStrips     = false; intersectsAny    = false; intersectsAnyStrips    = false;
  intersectsAnyNoBottom = false; intersectsAnyStripsNoBottom = false;

  intersectsExtendedTopLow                = false; intersectsRaisedTopLow                  = false;
  intersectsExtendedRaisedTopLow          = false; intersectsExtendedRaisedTopLowFullNorth = false;
  intersectsExtendedRaisedTopLowFullSouth = false; 

  hasTPCTrack = false; hasTPCLongTrack = false; hasHitMatch = false; hasTrackMatch = false;

  hitMatchTagger = -1; trackMatchTagger1 = -1; trackMatchTagger2 = -1; trackMatchTagger3 = -1;
}

bool CRTTopHighAna::IntersectsTagger(const std::string name, const art::Ptr<simb::MCParticle> &particle)
{
  const sbnd::CRTTaggerGeo tagger = fCRTGeoAlg.GetTagger(name);
  
  const double rmin[3] = {tagger.minX, tagger.minY, tagger.minZ};
  const double rmax[3] = {tagger.maxX, tagger.maxY, tagger.maxZ};

  return PassesThroughVolume(rmin, rmax, particle);
}

bool CRTTopHighAna::IntersectsModules(const std::string name, const art::Ptr<simb::MCParticle> &particle)
{
  const sbnd::CRTTaggerGeo tagger = fCRTGeoAlg.GetTagger(name);
  
  float minx = std::numeric_limits<float>::max(), miny = std::numeric_limits<float>::max(), minz = std::numeric_limits<float>::max();
  float maxx = -std::numeric_limits<float>::max(), maxy = -std::numeric_limits<float>::max(), maxz = -std::numeric_limits<float>::max();
  
  for(auto const& [modname, module] : tagger.modules)
    {
      const double rmin[3] = {module.minX, module.minY, module.minZ};
      const double rmax[3] = {module.maxX, module.maxY, module.maxZ};

      if(minx > rmin[0]) minx = rmin[0];
      if(miny > rmin[1]) miny = rmin[1];
      if(minz > rmin[2]) minz = rmin[2];
      if(maxx < rmax[0]) maxx = rmax[0];
      if(maxy < rmax[1]) maxy = rmax[1];
      if(maxz < rmax[2]) maxz = rmax[2];

      if(PassesThroughVolume(rmin, rmax, particle))
        return true;
    }

  return false;
}

bool CRTTopHighAna::IntersectsStrips(const std::string name, const art::Ptr<simb::MCParticle> &particle)
{
  const sbnd::CRTTaggerGeo tagger = fCRTGeoAlg.GetTagger(name);
  
  float minx = std::numeric_limits<float>::max(), miny = std::numeric_limits<float>::max(), minz = std::numeric_limits<float>::max();
  float maxx = -std::numeric_limits<float>::max(), maxy = -std::numeric_limits<float>::max(), maxz = -std::numeric_limits<float>::max();
  
  for(auto const& [modname, module] : tagger.modules)
    {
      for(auto const& [stripname, strip] : module.strips)
        {
          const double rmin[3] = {strip.minX, strip.minY, strip.minZ};
          const double rmax[3] = {strip.maxX, strip.maxY, strip.maxZ};

          if(minx > rmin[0]) minx = rmin[0];
          if(miny > rmin[1]) miny = rmin[1];
          if(minz > rmin[2]) minz = rmin[2];
          if(maxx < rmax[0]) maxx = rmax[0];
          if(maxy < rmax[1]) maxy = rmax[1];
          if(maxz < rmax[2]) maxz = rmax[2];

          if(PassesThroughVolume(rmin, rmax, particle))
            return true;
        }
    }

  return false;
}

bool CRTTopHighAna::PassesThroughVolume(const double rmin[3], const double rmax[3], const art::Ptr<simb::MCParticle> &particle)
{
  for(unsigned i = 0; i < particle->NumberTrajectoryPoints(); i++)
    {
      geo::Point_t point {particle->Vx(i), particle->Vy(i), particle->Vz(i)};
      if(IsInsideVolume(rmin, rmax, point))
        return true;

      if(i != 0)
        {
          geo::Point_t prevPoint {particle->Vx(i - 1), particle->Vy(i - 1), particle->Vz(i - 1)};
          geo::Vector_t dir = point - prevPoint;

          if(PassesThroughVolume(rmin, rmax, dir, prevPoint))
            return true;
        }
    }
  return false;
}

bool CRTTopHighAna::IsInsideVolume(const double rmin[3], const double rmax[3], const geo::Point_t &point)
{
  return (point.X() >= rmin[0] && point.X() <= rmax[0] &&
          point.Y() >= rmin[1] && point.Y() <= rmax[1] &&
          point.Z() >= rmin[2] && point.Z() <= rmax[2]);
}

bool CRTTopHighAna::IsInsideVolume(const double rmin[3], const double rmax[3], const TVector3 &point)
{
  return (point.X() >= rmin[0] && point.X() <= rmax[0] &&
          point.Y() >= rmin[1] && point.Y() <= rmax[1] &&
          point.Z() >= rmin[2] && point.Z() <= rmax[2]);
}

double CRTTopHighAna::RelativeDistanceToIntersect(const double plane, const double dir, const double start)
{
  return (plane - start) / dir;
}

bool CRTTopHighAna::PassesThroughVolume(const double rmin[3], const double rmax[3], const geo::Vector_t dir, const geo::Point_t start)
{
  const double minkx = RelativeDistanceToIntersect(rmin[0], dir.X(), start.X());
  const double minky = RelativeDistanceToIntersect(rmin[1], dir.Y(), start.Y());
  const double minkz = RelativeDistanceToIntersect(rmin[2], dir.Z(), start.Z());

  const double maxkx = RelativeDistanceToIntersect(rmax[0], dir.X(), start.X());
  const double maxky = RelativeDistanceToIntersect(rmax[1], dir.Y(), start.Y());
  const double maxkz = RelativeDistanceToIntersect(rmax[2], dir.Z(), start.Z());

  const geo::Point_t intersectminx = start + minkx * dir;
  const geo::Point_t intersectminy = start + minky * dir;
  const geo::Point_t intersectminz = start + minkz * dir;

  const geo::Point_t intersectmaxx = start + maxkx * dir;
  const geo::Point_t intersectmaxy = start + maxky * dir;
  const geo::Point_t intersectmaxz = start + maxkz * dir;

  return (IsInsideVolume(rmin, rmax, intersectminx) && minkx > 0 && minkx < 1) ||
    (IsInsideVolume(rmin, rmax, intersectminy) && minky > 0 && minky < 1) ||
    (IsInsideVolume(rmin, rmax, intersectminz) && minkz > 0 && minkz < 1) ||
    (IsInsideVolume(rmin, rmax, intersectmaxx) && maxkx > 0 && maxkx < 1) ||
    (IsInsideVolume(rmin, rmax, intersectmaxy) && maxky > 0 && maxky < 1) ||
    (IsInsideVolume(rmin, rmax, intersectmaxz) && maxkz > 0 && maxkz < 1);
}

bool CRTTopHighAna::IsGoodMatch(const art::Event &e, const int trackid, const art::Ptr<sbn::crt::CRTHit> &crthit)
{
  if(crthit.isNull())
    return false;

  art::Handle<std::vector<sbn::crt::CRTHit>> CRTHitHandle;
  e.getByLabel(fCRTHitLabel, CRTHitHandle);

  art::Handle<std::vector<sbnd::crt::CRTData>> CRTDataHandle;
  e.getByLabel(fCRTDataLabel, CRTDataHandle);

  art::Handle<std::vector<sbnd::crt::FEBData>> FEBDataHandle;
  e.getByLabel(fFEBDataLabel, FEBDataHandle);

  art::FindManyP<sbnd::crt::CRTData> crtHitToCRTData(CRTHitHandle, e, fCRTHitLabel);
  art::FindOneP<sbnd::crt::FEBData> crtDataToFEBData(CRTDataHandle, e, fCRTDataLabel);
  art::FindManyP<sim::AuxDetIDE, sbnd::crt::FEBTruthInfo> febDataToIDEs(FEBDataHandle, e, fFEBDataLabel);

  std::map<int, double> idemap;

  const std::vector<art::Ptr<sbnd::crt::CRTData>> crtdatas = crtHitToCRTData.at(crthit.key());

  for(auto const& crtdata : crtdatas)
    {
      const art::Ptr<sbnd::crt::FEBData> febdata = crtDataToFEBData.at(crtdata.key());
      const std::vector<art::Ptr<sim::AuxDetIDE>> ides = febDataToIDEs.at(febdata.key());

      unsigned i = 0;
      for(auto const& ide : ides)
        {
          const sbnd::crt::FEBTruthInfo *febtruthinfo = febDataToIDEs.data(febdata.key())[i];
          if((uint) febtruthinfo->GetChannel() == (crtdata->Channel() % 32))
            {
              int idetrackid = ide->trackID;
              int rolledidetrackid = fDroppedTrackIDMapHandle->GetAncestor(idetrackid);
              if(fDroppedTrackIDMapHandle->Exists(rolledidetrackid))
                idetrackid = rolledidetrackid;
      
              if(idemap.find(idetrackid) == idemap.end())
                idemap[idetrackid] = 0.;

              idemap[idetrackid] += ide->energyDeposited;
            }
        }
    }

  if(idemap.size() == 0)
    return false;

  std::vector<std::pair<int, double>> idevec(idemap.begin(), idemap.end());

  std::sort(idevec.begin(), idevec.end(),
            [](auto const& a, auto const& b)
            { return a.second > b.second; });

  return idevec.front().first == trackid;
}

bool CRTTopHighAna::IsGoodMatch(const art::Event &e, const int trackid, const art::Ptr<sbn::crt::CRTTrack> &crttrack)
{
  if(crttrack.isNull())
    return false;

  art::Handle<std::vector<sbn::crt::CRTTrack>> CRTTrackHandle;
  e.getByLabel(fCRTTrackLabel, CRTTrackHandle);

  art::FindManyP<sbn::crt::CRTHit> crtTracksToCRTHits(CRTTrackHandle, e, fCRTTrackLabel);
  const std::vector<art::Ptr<sbn::crt::CRTHit>> crthits = crtTracksToCRTHits.at(crttrack.key());

  if(crthits.size() == 0)
    return false;

  bool goodmatches = IsGoodMatch(e, trackid, crthits.front());

  for(unsigned i = 1; i < crthits.size(); ++i)
    goodmatches &= IsGoodMatch(e, trackid, crthits[i]);

  return goodmatches;
}

DEFINE_ART_MODULE(CRTTopHighAna)

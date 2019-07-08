#include "CosmicIdAlg.h"

namespace sbnd{

CosmicIdAlg::CosmicIdAlg(const Config& config){

  this->reconfigure(config);

}


CosmicIdAlg::CosmicIdAlg(){

}


CosmicIdAlg::~CosmicIdAlg(){

}


void CosmicIdAlg::reconfigure(const Config& config){

  fTpcTrackModuleLabel = config.TpcTrackModuleLabel();
  fPandoraLabel        = config.PandoraLabel();
  fCrtHitModuleLabel   = config.CrtHitModuleLabel();
  fCrtTrackModuleLabel = config.CrtTrackModuleLabel();
  fCaloModuleLabel     = config.CaloModuleLabel();

  fApplyFiducialCut  = config.ApplyFiducialCut();
  fApplyStoppingCut  = config.ApplyStoppingCut();
  fApplyGeometryCut  = config.ApplyGeometryCut();
  fApplyCpaCrossCut  = config.ApplyCpaCrossCut();
  fApplyApaCrossCut  = config.ApplyApaCrossCut();
  fApplyCrtTrackCut  = config.ApplyCrtTrackCut();
  fApplyCrtHitCut    = config.ApplyCrtHitCut();
  fApplyPandoraT0Cut = config.ApplyPandoraT0Cut();

  fUseTrackAngleVeto    = config.UseTrackAngleVeto();
  fMinSecondTrackLength = config.MinSecondTrackLength();
  fMinVertexDistance    = config.MinVertexDistance();
  fMinMergeAngle        = config.MinMergeAngle();

  fvTag = config.FVTagAlg();
  spTag = config.SPTagAlg();
  ccTag = config.CCTagAlg();
  acTag = config.ACTagAlg();
  chTag = config.CHTagAlg();
  ctTag = config.CTTagAlg();
  ptTag = config.PTTagAlg();

  return;
}

void CosmicIdAlg::SetCuts(bool FV, bool SP, bool Geo, bool CC, bool AC, bool CT, bool CH, bool PT){

  fApplyFiducialCut = FV;
  fApplyStoppingCut = SP;
  fApplyGeometryCut = Geo;
  fApplyCpaCrossCut = CC;
  fApplyApaCrossCut = AC;
  fApplyCrtTrackCut = CT;
  fApplyCrtHitCut = CH;
  fApplyPandoraT0Cut = PT;

}

bool CosmicIdAlg::CosmicId(recob::Track track, const art::Event& event, std::vector<double> t0Tpc0, std::vector<double> t0Tpc1){

  auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTpcTrackModuleLabel);
  art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTpcTrackModuleLabel);
  art::FindManyP<anab::Calorimetry> findManyCalo(tpcTrackHandle, event, fCaloModuleLabel);

  std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(track.ID());

  if(fApplyPandoraT0Cut){
    if(ptTag.PandoraT0CosmicId(track, event)) return true;
  }    

  if(fApplyFiducialCut){
    if(fvTag.FiducialVolumeCosmicId(track)) return true;
  }

  if(fApplyStoppingCut){
    std::vector<art::Ptr<anab::Calorimetry>> calos = findManyCalo.at(track.ID());
    if(spTag.StoppingParticleCosmicId(track, calos)) return true;
  }

  if(fApplyGeometryCut){
    bool tpc0Flash = CosmicIdUtils::BeamFlash(t0Tpc0, 4.);
    bool tpc1Flash = CosmicIdUtils::BeamFlash(t0Tpc1, 4.);
    if(geoTag.GeometryCosmicId(track, hits, tpc0Flash, tpc1Flash)) return true;
  }

  if(fApplyCpaCrossCut){
    std::vector<recob::Track> tracks;
    for(auto const& tpcTrack : (*tpcTrackHandle)){
      tracks.push_back(tpcTrack);
    }

    if(ccTag.CpaCrossCosmicId(track, tracks, findManyHits)) return true;
  }

  if(fApplyApaCrossCut){
    if(acTag.ApaCrossCosmicId(track, hits, t0Tpc0, t0Tpc1)) return true;
  }

  if(fApplyCrtTrackCut){
    auto crtTrackHandle = event.getValidHandle<std::vector<crt::CRTTrack>>(fCrtTrackModuleLabel);
    std::vector<crt::CRTTrack> crtTracks = (*crtTrackHandle);

    if(ctTag.CrtTrackCosmicId(track, crtTracks, event)) return true;
  }

  if(fApplyCrtHitCut){
    auto crtHitHandle = event.getValidHandle<std::vector<crt::CRTHit>>(fCrtHitModuleLabel);
    std::vector<crt::CRTHit> crtHits;
    for(auto const& crtHit : (*crtHitHandle)){
      if(crtHit.tagger != "volTaggerTopHigh_0"){
        crtHits.push_back(crtHit);
      }
    }

    if(chTag.CrtHitCosmicId(track, crtHits, event)) return true;
  }

  return false;

}

bool CosmicIdAlg::CosmicId(recob::PFParticle pfparticle, std::map< size_t, art::Ptr<recob::PFParticle> > pfParticleMap, const art::Event& event, std::vector<double> t0Tpc0, std::vector<double> t0Tpc1){

  art::Handle< std::vector<recob::PFParticle> > pfParticleHandle;
  event.getByLabel(fPandoraLabel, pfParticleHandle);
  art::FindManyP< recob::Track > pfPartToTrackAssoc(pfParticleHandle, event, fTpcTrackModuleLabel);

  auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTpcTrackModuleLabel);
  art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTpcTrackModuleLabel);
  art::FindManyP<anab::Calorimetry> findManyCalo(tpcTrackHandle, event, fCaloModuleLabel);

  std::vector<recob::Track> nuTracks;
  for (const size_t daughterId : pfparticle.Daughters()){
  
    // Get tracks associated with daughter
    art::Ptr<recob::PFParticle> pParticle = pfParticleMap.at(daughterId);
    const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pParticle.key()));
    if(associatedTracks.size() != 1) continue;

    recob::Track track = *associatedTracks.front();
    nuTracks.push_back(track);
    
    
  }

  if(fApplyPandoraT0Cut){
    if(ptTag.PandoraT0CosmicId(pfparticle, pfParticleMap, event)) return true;
  }

  if(nuTracks.size() == 0) return false;

  std::sort(nuTracks.begin(), nuTracks.end(), [](auto& left, auto& right){
              return left.Length() > right.Length();});

  recob::Track track = nuTracks[0];
  std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(track.ID());

  if(fApplyFiducialCut){
    if(fvTag.FiducialVolumeCosmicId(track)) return true;
  }

  if(fApplyGeometryCut){
    bool tpc0Flash = CosmicIdUtils::BeamFlash(t0Tpc0, 4.);
    bool tpc1Flash = CosmicIdUtils::BeamFlash(t0Tpc1, 4.);
    if(geoTag.GeometryCosmicId(track, hits, tpc0Flash, tpc1Flash)) return true;
  }

  if(fApplyCrtTrackCut){
    auto crtTrackHandle = event.getValidHandle<std::vector<crt::CRTTrack>>(fCrtTrackModuleLabel);
    std::vector<crt::CRTTrack> crtTracks = (*crtTrackHandle);

    if(ctTag.CrtTrackCosmicId(track, crtTracks, event)) return true;
  }

  if(fApplyCpaCrossCut){
    std::vector<recob::Track> tracks;
    for(auto const& tpcTrack : (*tpcTrackHandle)){
      tracks.push_back(tpcTrack);
    }

    if(ccTag.CpaCrossCosmicId(track, tracks, findManyHits)) return true;
  }

  std::vector<std::pair<recob::Track, double>> secondaryTracks;
  if(fUseTrackAngleVeto && nuTracks.size() > 1){
    TVector3 start = track.Vertex<TVector3>();
    TVector3 end = track.End<TVector3>();

    // Loop over the secondary tracks
    // Find smallest angle between primary track and any secondary tracks above a certain length
    for(size_t i = 1; i < nuTracks.size(); i++){
      recob::Track track2 = nuTracks[i];
      // Only consider secondary tracks longer than 5 cm (try to exclude michel electrons)
      if(track2.Length() < fMinSecondTrackLength) continue;
      TVector3 start2 = track2.Vertex<TVector3>();
      TVector3 end2 = track2.End<TVector3>();
      // Do they share the same vertex? (no delta rays)
      if((start-start2).Mag() < fMinVertexDistance){ 
        double angle = (end - start).Angle(end2 - start2);
        secondaryTracks.push_back(std::make_pair(track2, angle));
      }
    }
  }

  if(secondaryTracks.size() > 0){
    // Sort tracks by smallest angle
    std::sort(secondaryTracks.begin(), secondaryTracks.end(), [](auto& left, auto& right){
              return left.second < right.second;});
    // If secondary track angle is compatible with split track (near 180) then try to merge
    if(secondaryTracks[0].second > fMinMergeAngle){
      recob::Track track2 = secondaryTracks[0].first;

      if(fApplyFiducialCut){
        if(!fvTag.InFiducial(track.End()) && !fvTag.InFiducial(track2.End())) return true;
      }

      if(fApplyStoppingCut){
        // Apply stopping cut to the longest track
        std::vector<art::Ptr<anab::Calorimetry>> calos = findManyCalo.at(track.ID());
        if(spTag.StoppingParticleCosmicId(track, calos)) return true;
        // Apply stopping cut assuming the tracks are split
        std::vector<art::Ptr<anab::Calorimetry>> calos2 = findManyCalo.at(track2.ID());
        if(spTag.StoppingParticleCosmicId(track, track2, calos, calos2)) return true;
      }

      if(fApplyApaCrossCut){
        // Apply apa crossing cut to the longest track
        if(acTag.ApaCrossCosmicId(track, hits, t0Tpc0, t0Tpc1)) return true;
        // Also apply to secondary track FIXME need to check primary track doesn't go out of bounds
        std::vector<art::Ptr<recob::Hit>> hits2 = findManyHits.at(track2.ID());
        if(acTag.ApaCrossCosmicId(track2, hits2, t0Tpc0, t0Tpc1)) return true;
      }

      if(fApplyCrtHitCut){
        // Apply crt hit match cut to both tracks
        auto crtHitHandle = event.getValidHandle<std::vector<crt::CRTHit>>(fCrtHitModuleLabel);
        std::vector<crt::CRTHit> crtHits;
        for(auto const& crtHit : (*crtHitHandle)){
          if(crtHit.tagger != "volTaggerTopHigh_0"){
            crtHits.push_back(crtHit);
          }
        }
        if(chTag.CrtHitCosmicId(track, crtHits, event)) return true;
        if(chTag.CrtHitCosmicId(track2, crtHits, event)) return true;
      }
    }
    // Don't apply other cuts if angle between tracks is consistent with neutrino interaction
  }

  // If there's only one track from vertex then apply cuts as normal
  if(secondaryTracks.size() == 0){

    if(fApplyStoppingCut){
      std::vector<art::Ptr<anab::Calorimetry>> calos = findManyCalo.at(track.ID());
      if(spTag.StoppingParticleCosmicId(track, calos)) return true;
    }

    if(fApplyApaCrossCut){
      if(acTag.ApaCrossCosmicId(track, hits, t0Tpc0, t0Tpc1)) return true;
    }

    if(fApplyCrtHitCut){
      auto crtHitHandle = event.getValidHandle<std::vector<crt::CRTHit>>(fCrtHitModuleLabel);
      std::vector<crt::CRTHit> crtHits;
      for(auto const& crtHit : (*crtHitHandle)){
        if(crtHit.tagger != "volTaggerTopHigh_0"){
          crtHits.push_back(crtHit);
        }
      }

      if(chTag.CrtHitCosmicId(track, crtHits, event)) return true;
    }
  }

  return false;

}

}

#include "PandoraNuScoreCosmicIdAlg.h"

namespace sbnd{

  PandoraNuScoreCosmicIdAlg::PandoraNuScoreCosmicIdAlg(const Config& config){

    this->reconfigure(config);

  }


  PandoraNuScoreCosmicIdAlg::PandoraNuScoreCosmicIdAlg(){

  }


  PandoraNuScoreCosmicIdAlg::~PandoraNuScoreCosmicIdAlg(){

  }


  void PandoraNuScoreCosmicIdAlg::reconfigure(const Config& config){

    fPandoraLabel        = config.PandoraLabel();
    fTpcTrackModuleLabel = config.TpcTrackModuleLabel();
    fNuScoreCut          = config.NuScoreCut();

    return;
  }

  // Finds any t0s associated with track by pandora, tags if outside beam
  bool PandoraNuScoreCosmicIdAlg::PandoraNuScoreCosmicId(recob::Track track, const art::Event& event){

    // Get the pfps and associations
    art::Handle< std::vector<recob::PFParticle> > pfParticleHandle;
    event.getByLabel(fPandoraLabel, pfParticleHandle);
    art::FindManyP<recob::Track> pfPartToTrackAssoc(pfParticleHandle, event, fTpcTrackModuleLabel);
    art::FindManyP<larpandoraobj::PFParticleMetadata> PFPMetaDataAssoc(pfParticleHandle, event, fPandoraLabel);

    // Loop over all the pfps
    for(auto const &pfp : (*pfParticleHandle)){
      // Get the associated track if there is one
      const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pfp.Self()));
      if(associatedTracks.size() != 1) continue;
      recob::Track trk = *associatedTracks.front();
      if(trk.ID() != track.ID()) continue;

      recob::PFParticle PFPNeutrino = GetPFPNeutrino(pfp, (*pfParticleHandle));

      float pfpNuScore = GetPandoraNuScore(PFPNeutrino, PFPMetaDataAssoc);
      if (pfpNuScore < fNuScoreCut){
        return true;
      } else {
        return false;
      }
    }
    return false;

  }

  // Finds any t0s associated with pfparticle by pandora, tags if outside beam
  bool PandoraNuScoreCosmicIdAlg::PandoraNuScoreCosmicId(recob::PFParticle pfparticle,
      std::map< size_t, art::Ptr<recob::PFParticle> > pfParticleMap, const art::Event& event){

    // Get pfp associations to t0s
    art::Handle< std::vector<recob::PFParticle> > pfParticleHandle;
    event.getByLabel(fPandoraLabel, pfParticleHandle);

    art::FindManyP<larpandoraobj::PFParticleMetadata> PFPMetaDataAssoc(pfParticleHandle, event, fPandoraLabel);

    recob::PFParticle PFPNeutrino = GetPFPNeutrino(pfparticle, pfParticleMap);

    float pfpNuScore = GetPandoraNuScore(PFPNeutrino, PFPMetaDataAssoc);

    if (pfpNuScore < fNuScoreCut){
      return true;
    }
    return false;
  }


  recob::PFParticle PandoraNuScoreCosmicIdAlg::GetPFPNeutrino(recob::PFParticle pfparticle, std::map< size_t, art::Ptr<recob::PFParticle> >& pfParticleMap){

    if ((pfparticle.PdgCode()==12) ||(pfparticle.PdgCode()==14)){
      return pfparticle;
    } else {
      art::Ptr<recob::PFParticle> parentPFP= pfParticleMap.at(pfparticle.Parent());
      return GetPFPNeutrino(*parentPFP, pfParticleMap);
    }
  }

  recob::PFParticle PandoraNuScoreCosmicIdAlg::GetPFPNeutrino(recob::PFParticle pfparticle,
      const std::vector<recob::PFParticle>& pfpVec){

    if ((pfparticle.PdgCode()==12) ||(pfparticle.PdgCode()==14)){
      return pfparticle;
    } else {
      size_t parentID = pfparticle.Parent();
      auto parentPFPIter = std::find_if(pfpVec.begin(), pfpVec.end(),
          [&](const auto& pfp){return pfp.Self()==parentID;});

      if (parentPFPIter==pfpVec.end()){
        return pfparticle;
      } else {
        return GetPFPNeutrino(*parentPFPIter, pfpVec);
      }
    }
  }

  float PandoraNuScoreCosmicIdAlg::GetPandoraNuScore(recob::PFParticle pfparticle,
      art::FindManyP<larpandoraobj::PFParticleMetadata> PFPMetaDataAssoc){

    const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata> > pfpMetaVec =
      PFPMetaDataAssoc.at(pfparticle.Self());

    if (pfpMetaVec.size() !=1){
      std::cout<<"Cannot get PFPMetadata"<<std::endl;
      return 99999;
    }

    art::Ptr<larpandoraobj::PFParticleMetadata> pfpMeta = pfpMetaVec.front();

    larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap = pfpMeta->GetPropertiesMap();
    auto propertiesMapIter = propertiesMap.find("NuScore");
    if (propertiesMapIter == propertiesMap.end()){
      std::cout<<"Cannot get PFP Nu Score in Metadata"<<std::endl;
      std::cout<<"PFP pdg: "<<pfparticle.PdgCode()<<std::endl;
      return 99999;
    }

    return propertiesMapIter->second;
  }

}

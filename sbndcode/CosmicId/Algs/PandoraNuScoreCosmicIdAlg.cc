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
    for(auto const pfp : (*pfParticleHandle)){
      // Get the associated track if there is one
      const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pfp.Self()));
      if(associatedTracks.size() != 1) continue;
      recob::Track trk = *associatedTracks.front();
      if(trk.ID() != track.ID()) continue;

      recob::PFParticle PFPNeutrino = GetPFPNeutrino(pfp, (*pfParticleHandle));

      const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata> > pfpMetaVec = PFPMetaDataAssoc.at(PFPNeutrino.Self());
      if (pfpMetaVec.size() !=1){
        std::cout<<"Cannot get PFPMetadata"<<std::endl;
        return false;
      }

      art::Ptr<larpandoraobj::PFParticleMetadata> pfpMeta = pfpMetaVec.front();

      larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap = pfpMeta->GetPropertiesMap();
      auto propertiesMapIter = propertiesMap.find("NuScore");
      if (propertiesMapIter == propertiesMap.end())
        return false;
      float pfpNuScore = propertiesMapIter->second;

      if (pfpNuScore < fNuScoreCut){
        return true;
      }

    }
    return false;

  }

  // Finds any t0s associated with pfparticle by pandora, tags if outside beam
  bool PandoraNuScoreCosmicIdAlg::PandoraNuScoreCosmicId(recob::PFParticle pfparticle, std::map< size_t, art::Ptr<recob::PFParticle> > pfParticleMap, const art::Event& event){

    // Get pfp associations to t0s
    art::Handle< std::vector<recob::PFParticle> > pfParticleHandle;
    event.getByLabel(fPandoraLabel, pfParticleHandle);

    art::FindManyP<larpandoraobj::PFParticleMetadata> PFPMetaDataAssoc(pfParticleHandle, event, fPandoraLabel);

    recob::PFParticle PFPNeutrino = GetPFPNeutrino(pfparticle, pfParticleMap);

    const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata> > pfpMetaVec = PFPMetaDataAssoc.at(PFPNeutrino.Self());
    if (pfpMetaVec.size() !=1){
      std::cout<<"Cannot get PFPMetadata"<<std::endl;
      return false;
    }

    art::Ptr<larpandoraobj::PFParticleMetadata> pfpMeta = pfpMetaVec.front();

    larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap = pfpMeta->GetPropertiesMap();
    auto propertiesMapIter = propertiesMap.find("NuScore");
    if (propertiesMapIter == propertiesMap.end())
      return false;
    float pfpNuScore = propertiesMapIter->second;

    if (pfpNuScore < fNuScoreCut){
      return true;
    }
    return false;
  }


  recob::PFParticle PandoraNuScoreCosmicIdAlg::GetPFPNeutrino(recob::PFParticle pfparticle, std::map< size_t, art::Ptr<recob::PFParticle> >& pfParticleMap){

    if ((pfparticle.PdgCode()==12) ||(pfparticle.PdgCode()==14)){
      return pfparticle;
    } else {
      art::Ptr<recob::PFParticle> parentPFP = pfParticleMap.at(pfparticle.Parent());
      return GetPFPNeutrino(*parentPFP, pfParticleMap);
    }
  }

  recob::PFParticle PandoraNuScoreCosmicIdAlg::GetPFPNeutrino(recob::PFParticle pfparticle, const std::vector<recob::PFParticle>& pfpVec){

    if ((pfparticle.PdgCode()==12) ||(pfparticle.PdgCode()==14)){
      return pfparticle;
    } else {
      size_t parentID = pfparticle.Parent();
      recob::PFParticle parentPFP;
      for (auto const& pfpIter: pfpVec){
        if (pfpIter.Self() == parentID){
          parentPFP = pfpIter;
          break;
        }
      }
      return GetPFPNeutrino(parentPFP, pfpVec);
    }
  }


}

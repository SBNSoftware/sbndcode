#include "CRTBackTrackerAlg.h"

namespace sbnd::crt {
  
  CRTBackTrackerAlg::CRTBackTrackerAlg(const Config& config)
  {
    this->reconfigure(config);
  }
  
  CRTBackTrackerAlg::CRTBackTrackerAlg(){}
  
  CRTBackTrackerAlg::~CRTBackTrackerAlg(){}

  void CRTBackTrackerAlg::reconfigure(const Config& config)
  {
    fSimModuleLabel = config.SimModuleLabel();
    fSimDepositModuleLabel = config.SimDepositModuleLabel();
    fFEBDataModuleLabel = config.FEBDataModuleLabel();
    fStripHitModuleLabel = config.StripHitModuleLabel();
    fClusterModuleLabel = config.ClusterModuleLabel();

    return;
  }

  void CRTBackTrackerAlg::SetupMaps(const art::Event &event)
  {
    fMCPIDEsEnergyMap.clear();
    fMCPStripHitsMap.clear();

    fTrackIDMotherMap.clear();
    fStripHitMCPMap.clear();

    auto droppedTrackIdMaps = event.getMany<std::map<int, std::set<int>>>();

    for(auto const& droppedTrackIdMap : droppedTrackIdMaps)
      {
        for(auto const& [mother, ids] : *droppedTrackIdMap)
          {
            for(auto const& id : ids)
              fTrackIDMotherMap[id] = mother;
          }
      }

    art::Handle<std::vector<sim::AuxDetIDE>> ideHandle;
    event.getByLabel(fFEBDataModuleLabel, ideHandle);
    std::vector<art::Ptr<sim::AuxDetIDE>> ideVec;
    art::fill_ptr_vector(ideVec, ideHandle);

    std::set<std::pair<int, CRTTagger>> depositCategories;

    for(auto const ide : ideVec)
      {
        const double x = (ide->entryX + ide->exitX) / 2.;
        const double y = (ide->entryY + ide->exitY) / 2.;
        const double z = (ide->entryZ + ide->exitZ) / 2.;
        const CRTTagger tagger = fCRTGeoAlg.WhichTagger(x, y, z);

        depositCategories.insert({RollUpID(ide->trackID), tagger});
      }
        
    for(auto const category : depositCategories)
      fMCPIDEsEnergyMap[category] = 0.;
      
    for(auto const ide : ideVec)
      {
        const double x = (ide->entryX + ide->exitX) / 2.;
        const double y = (ide->entryY + ide->exitY) / 2.;
        const double z = (ide->entryZ + ide->exitZ) / 2.;
        const CRTTagger tagger = fCRTGeoAlg.WhichTagger(x, y, z);

        fMCPIDEsEnergyMap[{RollUpID(ide->trackID), tagger}] += ide->energyDeposited;
      }

    art::Handle<std::vector<CRTStripHit>> stripHitHandle;
    event.getByLabel(fStripHitModuleLabel, stripHitHandle);
    std::vector<art::Ptr<CRTStripHit>> stripHitVec;
    art::fill_ptr_vector(stripHitVec, stripHitHandle);

    for(auto const stripHit : stripHitVec)
      {
        const CRTTagger tagger = fCRTGeoAlg.ChannelToTaggerEnum(stripHit->Channel());
        TruthMatchMetrics truthMatch = TruthMatching(event, stripHit);

        fStripHitMCPMap[stripHit.key()] = truthMatch.trackid;

        if(fMCPStripHitsMap.find({truthMatch.trackid, tagger}) == fMCPStripHitsMap.end())
          fMCPStripHitsMap[{truthMatch.trackid, tagger}] = 0;

        ++fMCPStripHitsMap[{truthMatch.trackid, tagger}];
      }
  }

  int CRTBackTrackerAlg::RollUpID(const int &id)
  {
    if(fTrackIDMotherMap.find(id) != fTrackIDMotherMap.end())
      return fTrackIDMotherMap[id];

    return id;
  }

  CRTBackTrackerAlg::TruthMatchMetrics CRTBackTrackerAlg::TruthMatching(const art::Event &event, const art::Ptr<CRTStripHit> &stripHit)
  {  
    art::Handle<std::vector<FEBData>> febDataHandle;
    event.getByLabel(fFEBDataModuleLabel, febDataHandle);

    art::Handle<std::vector<CRTStripHit>> stripHitHandle;
    event.getByLabel(fStripHitModuleLabel, stripHitHandle);
    std::vector<art::Ptr<CRTStripHit>> stripHitVec;
    art::fill_ptr_vector(stripHitVec, stripHitHandle);

    art::FindManyP<sim::AuxDetIDE, FEBTruthInfo> febDataToIDEs(febDataHandle, event, fFEBDataModuleLabel);
    art::FindManyP<FEBData> stripHitToFEBData(stripHitHandle, event, fStripHitModuleLabel);
    const CRTTagger tagger = fCRTGeoAlg.ChannelToTaggerEnum(stripHit->Channel());

    std::map<int, double> idToEnergyMap;
    double totalEnergy = 0.;

    auto const febData = stripHitToFEBData.at(stripHit.key());
    assert(febData.size() == 1);
    auto const assnIDEVec = febDataToIDEs.at(febData.at(0).key());
    for(unsigned i = 0; i < assnIDEVec.size(); ++i)
      {
        const art::Ptr<sim::AuxDetIDE> ide = assnIDEVec[i];
        const FEBTruthInfo *febTruthInfo = febDataToIDEs.data(febData.at(0).key())[i];
        if((uint) febTruthInfo->GetChannel() == (stripHit->Channel() % 32))
          {
            idToEnergyMap[RollUpID(ide->trackID)] += ide->energyDeposited;
            totalEnergy                           += ide->energyDeposited;
          }
      }

    double bestPur = 0., comp = 0.;
    int trackid = -99999;

    std::cout << "New Strip Hit" << std::endl;

    for(auto const [id, en] : idToEnergyMap)
      {
        const simb::MCParticle* particle = particleInv->TrackIdToParticle_P(id);
        const int pdg = particle == NULL ? -1 : particle->PdgCode();
        const std::string endprocess = particle == NULL ? "" : particle->EndProcess();
        std::cout << "\tTrackID: " << id << " En: " << en << " PDG: " << pdg << " EndProc: " <<  endprocess << std::endl;
        double pur = en / totalEnergy;
        if(pur > bestPur)
          {
            trackid = id;
            bestPur = pur;
            comp    = en / fMCPIDEsEnergyMap[{id, tagger}];
          }
      }
    std::cout << "\tPur: " << bestPur << " Comp: " << comp << std::endl;
    std::cout << std::endl;
    return TruthMatchMetrics(trackid, comp, bestPur);
  }

  CRTBackTrackerAlg::TruthMatchMetrics CRTBackTrackerAlg::TruthMatching(const art::Event &event, const art::Ptr<CRTCluster> &cluster)
  {
    art::Handle<std::vector<FEBData>> febDataHandle;
    event.getByLabel(fFEBDataModuleLabel, febDataHandle);

    art::Handle<std::vector<CRTStripHit>> stripHitHandle;
    event.getByLabel(fStripHitModuleLabel, stripHitHandle);

    art::Handle<std::vector<CRTCluster>> clusterHandle;
    event.getByLabel(fClusterModuleLabel, clusterHandle);
    std::vector<art::Ptr<CRTCluster>> clusterVec;
    art::fill_ptr_vector(clusterVec, clusterHandle);

    art::FindManyP<sim::AuxDetIDE, FEBTruthInfo> febDataToIDEs(febDataHandle, event, fFEBDataModuleLabel);
    art::FindManyP<FEBData> stripHitToFEBData(stripHitHandle, event, fStripHitModuleLabel);
    art::FindManyP<CRTStripHit> clusterToStripHits(clusterHandle, event, fClusterModuleLabel);

    std::map<int, double> idToEnergyMap;
    double totalEnergy = 0.;
    std::map<int, int> idToNHitsMap;

    auto const assnStripHitVec = clusterToStripHits.at(cluster.key());

    for(auto const stripHit : assnStripHitVec)
      {
        auto const febData = stripHitToFEBData.at(stripHit.key());
        assert(febData.size() == 1);
        auto const assnIDEVec = febDataToIDEs.at(febData.at(0).key());
        for(unsigned i = 0; i < assnIDEVec.size(); ++i)
          {
            const art::Ptr<sim::AuxDetIDE> ide = assnIDEVec[i];
            const FEBTruthInfo *febTruthInfo = febDataToIDEs.data(febData.at(0).key())[i];
            if((uint) febTruthInfo->GetChannel() == (stripHit->Channel() % 32))
              {
                idToEnergyMap[RollUpID(ide->trackID)] += ide->energyDeposited;
                totalEnergy                           += ide->energyDeposited;
              }
          }

        if(idToNHitsMap.find(fStripHitMCPMap[stripHit.key()]) == idToNHitsMap.end())
           idToNHitsMap[fStripHitMCPMap[stripHit.key()]] = 0;
        
        ++idToNHitsMap[fStripHitMCPMap[stripHit.key()]];
      }

    double bestPur = 0., comp = 0.;
    int trackid = -99999;

    for(auto const [id, en] : idToEnergyMap)
      {
        double pur = en / totalEnergy;
        if(pur > bestPur)
          {
            trackid = id;
            bestPur = pur;
            comp    = en / fMCPIDEsEnergyMap[{id, cluster->Tagger()}];
          }
      }

    double hitComp = idToNHitsMap[trackid] / (double) fMCPStripHitsMap[{trackid, cluster->Tagger()}];
    double hitPur  = idToNHitsMap[trackid] / (double) assnStripHitVec.size();

    return TruthMatchMetrics(trackid, comp, bestPur, hitComp, hitPur);
  }
}

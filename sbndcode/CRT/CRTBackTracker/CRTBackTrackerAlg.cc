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
    fSpacePointModuleLabel = config.SpacePointModuleLabel();
    fTrackModuleLabel = config.TrackModuleLabel();

    return;
  }

  void CRTBackTrackerAlg::SetupMaps(const art::Event &event)
  {
    fTrueDepositsPerTaggerMap.clear();
    fTrueDepositsMap.clear();
    fTrueTrackInfosMap.clear();
    fTrackIDSpacePointRecoMap.clear();
    fTrackIDTrackRecoMap.clear();
    fMCPIDEsEnergyPerTaggerMap.clear();
    fMCPIDEsEnergyMap.clear();
    fMCPStripHitsMap.clear();
    fTrackIDMotherMap.clear();
    fStripHitMCPMap.clear();

    art::Handle<std::vector<sim::ParticleAncestryMap>> droppedTrackIDMapVecHandle;
    event.getByLabel(fSimModuleLabel, droppedTrackIDMapVecHandle);

    if(droppedTrackIDMapVecHandle.isValid())
      {
	for(auto const& droppedTrackIdMap : *droppedTrackIDMapVecHandle)
	  {
	    for(auto const& [mother, ids] : droppedTrackIdMap.GetMap())
	      {
		for(auto const& id : ids)
		  fTrackIDMotherMap[id] = mother;
	      }
	  }
      }

    art::Handle<sim::ParticleAncestryMap> droppedTrackIDMapHandle;
    event.getByLabel(fSimModuleLabel, droppedTrackIDMapHandle);

    if(droppedTrackIDMapHandle.isValid())
      {
	auto const& droppedTrackIdMap = droppedTrackIDMapHandle->GetMap();

	for(auto const& [mother, ids] : droppedTrackIdMap)
	  {
	    for(auto const& id : ids)
	      fTrackIDMotherMap[id] = mother;
	  }
      }

    art::Handle<std::vector<sim::AuxDetIDE>> ideHandle;
    event.getByLabel(fFEBDataModuleLabel, ideHandle);
    std::vector<art::Ptr<sim::AuxDetIDE>> ideVec;
    art::fill_ptr_vector(ideVec, ideHandle);

    std::set<Category> depositCategories;
    std::set<int>      depositTrackIDs;

    std::map<Category, double> categoryToEnergyMap, categoryToXMap, categoryToYMap,
      categoryToZMap, categoryToTimeMap;
    std::map<Category, uint> categoryToNIDEsMap;

    std::map<int, double> idToEnergyMap, idToXMap, idToYMap,
      idToZMap, idToTimeMap;
    std::map<int, uint> idToNIDEsMap;

    std::map<Category, double> coreCategoryToEnergyMap, coreCategoryToXMap, coreCategoryToYMap,
      coreCategoryToZMap, coreCategoryToTimeMap;
    std::map<Category, uint> coreCategoryToNIDEsMap;

    for(auto const ide : ideVec)
      {
        const double x = (ide->entryX + ide->exitX) / 2.;
        const double y = (ide->entryY + ide->exitY) / 2.;
        const double z = (ide->entryZ + ide->exitZ) / 2.;
        const CRTTagger tagger = fCRTGeoAlg.WhichTagger(x, y, z);

        const int rollUpID = RollUpID(ide->trackID);

        Category category(rollUpID, tagger);
        depositCategories.insert(category);
        depositTrackIDs.insert(rollUpID);

        fTrackIDSpacePointRecoMap[category] = false;
        fTrackIDTrackRecoMap[rollUpID]      = { false, false };

        categoryToEnergyMap[category] += ide->energyDeposited;
        categoryToXMap[category]      += x;
        categoryToYMap[category]      += y;
        categoryToZMap[category]      += z;
        categoryToTimeMap[category]   += (ide->entryT + ide->exitT) / 2.;

        idToEnergyMap[rollUpID] += ide->energyDeposited;
        idToXMap[rollUpID]      += x;
        idToYMap[rollUpID]      += y;
        idToZMap[rollUpID]      += z;
        idToTimeMap[rollUpID]   += (ide->entryT + ide->exitT) / 2.;

        Category coreCategory(ide->trackID, tagger);
        coreCategoryToEnergyMap[coreCategory] += ide->energyDeposited;
        coreCategoryToXMap[coreCategory]      += x;
        coreCategoryToYMap[coreCategory]      += y;
        coreCategoryToZMap[coreCategory]      += z;
        coreCategoryToTimeMap[coreCategory]   += (ide->entryT + ide->exitT) / 2.;

        if(categoryToNIDEsMap.find(category) == categoryToNIDEsMap.end())
          categoryToNIDEsMap[category] = 0;

        ++categoryToNIDEsMap[category];

        if(idToNIDEsMap.find(rollUpID) == idToNIDEsMap.end())
          idToNIDEsMap[rollUpID] = 0;

        ++idToNIDEsMap[rollUpID];

        if(coreCategoryToNIDEsMap.find(coreCategory) == coreCategoryToNIDEsMap.end())
          coreCategoryToNIDEsMap[coreCategory] = 0;

        ++coreCategoryToNIDEsMap[coreCategory];
      }

    for(auto const category : depositCategories)
      {
        fMCPIDEsEnergyPerTaggerMap[category] = 0.;
        
        const double x      = categoryToXMap[category] / categoryToNIDEsMap[category];
        const double y      = categoryToYMap[category] / categoryToNIDEsMap[category];
        const double z      = categoryToZMap[category] / categoryToNIDEsMap[category];
        const double time   = categoryToTimeMap[category] / categoryToNIDEsMap[category];
        const double energy = categoryToEnergyMap[category];

        const double coreX      = coreCategoryToXMap[category] / coreCategoryToNIDEsMap[category];
        const double coreY      = coreCategoryToYMap[category] / coreCategoryToNIDEsMap[category];
        const double coreZ      = coreCategoryToZMap[category] / coreCategoryToNIDEsMap[category];
        const double coreTime   = coreCategoryToTimeMap[category] / coreCategoryToNIDEsMap[category];
        const double coreEnergy = coreCategoryToEnergyMap[category];

        int pdg;
        double particle_energy, particle_time;
        TrueParticlePDGEnergyTime(category.trackid, pdg, particle_energy, particle_time);

        fTrueDepositsPerTaggerMap[category] = TrueDeposit(category.trackid, pdg, category.tagger,
                                                          energy, time, x, y, z, true,
                                                          coreEnergy, coreTime, coreX, coreY, coreZ);
      }

    for(auto const trackID : depositTrackIDs)
      {
        fMCPIDEsEnergyMap[trackID] = 0.;
        
        const double x      = idToXMap[trackID] / idToNIDEsMap[trackID];
        const double y      = idToYMap[trackID] / idToNIDEsMap[trackID];
        const double z      = idToZMap[trackID] / idToNIDEsMap[trackID];
        const double time   = idToTimeMap[trackID] / idToNIDEsMap[trackID];
        const double energy = idToEnergyMap[trackID];

        int pdg;
        double particle_energy, particle_time;
        TrueParticlePDGEnergyTime(trackID, pdg, particle_energy, particle_time);

        struct SortTagger {
          double    time;
          double    energy;
          CRTTagger tagger;
        };

        std::vector<SortTagger> taggers;

        for(auto const& [category, deposit] : fTrueDepositsPerTaggerMap)
          {
            if(category.trackid == trackID)
              taggers.push_back({deposit.time, deposit.energy, category.tagger});
          }

        std::sort(taggers.begin(), taggers.end(),
                  [](const SortTagger &a, const SortTagger &b)
                  { return a.time < b.time; });

        Category category1, category2;

        if(taggers.size() == 1)
          {         
            category1 = {trackID, taggers[0].tagger};
            category2 = {trackID, taggers[0].tagger};
          }
        else if(taggers.size() == 2)
          {
            category1 = {trackID, taggers[0].tagger};
            category2 = {trackID, taggers[1].tagger};
          }
        else if(taggers.size() == 3)
          {
            category1 = {trackID, taggers[0].tagger};
            category2 = {trackID, taggers[2].tagger};
          }
        else if(taggers.size() > 3)
          {
            std::sort(taggers.begin(), taggers.end(),
                      [](const auto &a, const auto &b)
                      { return a.energy < b.energy; });

            const unsigned n = taggers.size();

            if(taggers[n-1].time < taggers[n-2].time)
              {
                category1 = {trackID, taggers[n-1].tagger};
                category2 = {trackID, taggers[n-2].tagger};
              }
            else
              {
                category1 = {trackID, taggers[n-2].tagger};
                category2 = {trackID, taggers[n-1].tagger};
              }
          }

        fTrueTrackInfosMap[trackID] = TrueTrackInfo(trackID, pdg, particle_energy,
                                                    fTrueDepositsPerTaggerMap[category1],
                                                    fTrueDepositsPerTaggerMap[category2]);

        fTrueDepositsMap[trackID] = TrueDeposit(trackID, pdg, kUndefinedTagger,
                                                energy, time, x, y, z, taggers.size() > 1);
      }

    for(auto const ide : ideVec)
      {
        const double x = (ide->entryX + ide->exitX) / 2.;
        const double y = (ide->entryY + ide->exitY) / 2.;
        const double z = (ide->entryZ + ide->exitZ) / 2.;
        const CRTTagger tagger = fCRTGeoAlg.WhichTagger(x, y, z);

        const int rollUpID = RollUpID(ide->trackID);
        Category category(rollUpID, tagger);

        fMCPIDEsEnergyPerTaggerMap[category] += ide->energyDeposited;
        fMCPIDEsEnergyMap[rollUpID]          += ide->energyDeposited;
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

        Category category(truthMatch.trackid, tagger);
        if(fMCPStripHitsMap.find(category) == fMCPStripHitsMap.end())
          fMCPStripHitsMap[category] = 0;

        ++fMCPStripHitsMap[{truthMatch.trackid, tagger}];
      }
  }

  int CRTBackTrackerAlg::RollUpID(const int &id)
  {
    if(fTrackIDMotherMap.find(id) != fTrackIDMotherMap.end())
      return fTrackIDMotherMap[id];

    return id;
  }

  void CRTBackTrackerAlg::RunSpacePointRecoStatusChecks(const art::Event &event)
  {
    art::Handle<std::vector<CRTSpacePoint>> spacePointHandle;
    event.getByLabel(fSpacePointModuleLabel, spacePointHandle);

    art::FindManyP<CRTCluster> spacePointsToClusters(spacePointHandle, event, fSpacePointModuleLabel);

    for(unsigned i = 0; i < spacePointHandle->size(); ++i)
      {
        const art::Ptr<CRTSpacePoint> spacePoint(spacePointHandle, i);
        const art::Ptr<CRTCluster> cluster = spacePointsToClusters.at(spacePoint.key())[0];

        TruthMatchMetrics truthMatch = TruthMatching(event, cluster);

        Category category(truthMatch.trackid, cluster->Tagger());
        fTrackIDSpacePointRecoMap[category] = true;
      }
  }

  void CRTBackTrackerAlg::RunTrackRecoStatusChecks(const art::Event &event)
  {
    art::Handle<std::vector<CRTTrack>> trackHandle;
    event.getByLabel(fTrackModuleLabel, trackHandle);

    for(unsigned i = 0; i < trackHandle->size(); ++i)
      {
        const art::Ptr<CRTTrack> track(trackHandle, i);

        TruthMatchMetrics truthMatch = TruthMatching(event, track);

        fTrackIDTrackRecoMap[truthMatch.trackid] = { true, track->Triple()};
      }
  }

  std::map<CRTBackTrackerAlg::Category, bool> CRTBackTrackerAlg::GetSpacePointRecoStatusMap()
  {
    return fTrackIDSpacePointRecoMap;
  }

  std::map<int, std::pair<bool, bool>> CRTBackTrackerAlg::GetTrackRecoStatusMap()
  {
    return fTrackIDTrackRecoMap;
  }

  CRTBackTrackerAlg::TrueDeposit CRTBackTrackerAlg::GetTrueDeposit(Category category)
  {
    return fTrueDepositsPerTaggerMap[category];
  }

  CRTBackTrackerAlg::TrueDeposit CRTBackTrackerAlg::GetTrueDeposit(int trackid)
  {
    return fTrueDepositsMap[trackid];
  }

  CRTBackTrackerAlg::TruthMatchMetrics CRTBackTrackerAlg::TruthMatching(const art::Event &event, const art::Ptr<CRTStripHit> &stripHit)
  {  
    art::Handle<std::vector<FEBData>> febDataHandle;
    event.getByLabel(fFEBDataModuleLabel, febDataHandle);

    art::Handle<std::vector<CRTStripHit>> stripHitHandle;
    event.getByLabel(fStripHitModuleLabel, stripHitHandle);

    art::FindManyP<sim::AuxDetIDE, FEBTruthInfo> febDataToIDEs(febDataHandle, event, fFEBDataModuleLabel);
    art::FindOneP<FEBData> stripHitToFEBData(stripHitHandle, event, fStripHitModuleLabel);
    const CRTTagger tagger = fCRTGeoAlg.ChannelToTaggerEnum(stripHit->Channel());

    auto const febData = stripHitToFEBData.at(stripHit.key());
    auto const assnIDEVec = febDataToIDEs.at(febData.key());

    std::map<int, double> idToEnergyMap;
    double totalEnergy = 0., x = 0., y = 0., z = 0., t = 0.;
    uint nides = 0;

    for(unsigned i = 0; i < assnIDEVec.size(); ++i)
      {
        const art::Ptr<sim::AuxDetIDE> ide = assnIDEVec[i];
        const FEBTruthInfo *febTruthInfo = febDataToIDEs.data(febData.key())[i];
        if((uint) febTruthInfo->GetChannel() == (stripHit->Channel() % 32))
          {
            idToEnergyMap[RollUpID(ide->trackID)] += ide->energyDeposited;
            totalEnergy                           += ide->energyDeposited;

            x                                     += (ide->entryX + ide->exitX) / 2.;
            y                                     += (ide->entryY + ide->exitY) / 2.;
            z                                     += (ide->entryZ + ide->exitZ) / 2.;
            t                                     += (ide->entryT + ide->exitT) / 2.;

            ++nides;
          }
      }

    x /= nides;
    y /= nides;
    z /= nides;
    t /= nides;

    double bestPur = 0., comp = 0.;
    int trackid = -99999;

    for(auto const [id, en] : idToEnergyMap)
      {
        double pur = en / totalEnergy;
        if(pur > bestPur)
          {
            Category category(id, tagger);

            trackid = id;
            bestPur = pur;
            comp    = en / fMCPIDEsEnergyPerTaggerMap[category];
          }
      }

    TrueDeposit deposit(trackid, -999999, tagger, totalEnergy, t, x, y, z, true);

    return TruthMatchMetrics(trackid, comp, bestPur, 1., 1., deposit);
  }

  CRTBackTrackerAlg::TruthMatchMetrics CRTBackTrackerAlg::TruthMatching(const art::Event &event, const art::Ptr<CRTCluster> &cluster)
  {
    art::Handle<std::vector<FEBData>> febDataHandle;
    event.getByLabel(fFEBDataModuleLabel, febDataHandle);

    art::Handle<std::vector<CRTStripHit>> stripHitHandle;
    event.getByLabel(fStripHitModuleLabel, stripHitHandle);

    art::Handle<std::vector<CRTCluster>> clusterHandle;
    event.getByLabel(fClusterModuleLabel, clusterHandle);

    art::FindManyP<sim::AuxDetIDE, FEBTruthInfo> febDataToIDEs(febDataHandle, event, fFEBDataModuleLabel);
    art::FindManyP<FEBData> stripHitToFEBData(stripHitHandle, event, fStripHitModuleLabel);
    art::FindManyP<CRTStripHit> clusterToStripHits(clusterHandle, event, fClusterModuleLabel);

    std::map<int, double> idToEnergyMap;
    double totalEnergy = 0.;
    std::map<int, uint> idToNHitsMap;

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
            Category category(id, cluster->Tagger());

            trackid = id;
            bestPur = pur;
            comp    = en / fMCPIDEsEnergyPerTaggerMap[category];
          }
      }

    Category category(trackid, cluster->Tagger());

    double hitComp = idToNHitsMap[trackid] / (double) fMCPStripHitsMap[category];
    double hitPur  = idToNHitsMap[trackid] / (double) assnStripHitVec.size();

    return TruthMatchMetrics(trackid, comp, bestPur, hitComp, hitPur, 
                             fTrueDepositsPerTaggerMap[category]);
  }

  CRTBackTrackerAlg::TruthMatchMetrics CRTBackTrackerAlg::TruthMatching(const art::Event &event, const art::Ptr<CRTTrack> &track)
  {
    art::Handle<std::vector<FEBData>> febDataHandle;
    event.getByLabel(fFEBDataModuleLabel, febDataHandle);

    art::Handle<std::vector<CRTStripHit>> stripHitHandle;
    event.getByLabel(fStripHitModuleLabel, stripHitHandle);

    art::Handle<std::vector<CRTCluster>> clusterHandle;
    event.getByLabel(fClusterModuleLabel, clusterHandle);

    art::Handle<std::vector<CRTSpacePoint>> spacePointHandle;
    event.getByLabel(fSpacePointModuleLabel, spacePointHandle);

    art::Handle<std::vector<CRTTrack>> trackHandle;
    event.getByLabel(fTrackModuleLabel, trackHandle);

    art::FindManyP<sim::AuxDetIDE, FEBTruthInfo> febDataToIDEs(febDataHandle, event, fFEBDataModuleLabel);
    art::FindManyP<FEBData> stripHitToFEBData(stripHitHandle, event, fStripHitModuleLabel);
    art::FindManyP<CRTStripHit> clusterToStripHits(clusterHandle, event, fClusterModuleLabel);
    art::FindOneP<CRTCluster> spacePointToCluster(spacePointHandle, event, fSpacePointModuleLabel);
    art::FindManyP<CRTSpacePoint> trackToSpacePoints(trackHandle, event, fTrackModuleLabel);

    std::map<int, double> idToEnergyMap;
    double totalEnergy = 0.;

    auto const spacePointVec = trackToSpacePoints.at(track.key());

    for(auto const spacePoint : spacePointVec)
      {
        auto const cluster     = spacePointToCluster.at(spacePoint.key());
        auto const stripHitVec = clusterToStripHits.at(cluster.key());
        
        for(auto const stripHit : stripHitVec)
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
          }
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
            comp    = en / fMCPIDEsEnergyMap[id];
          }
      }

    return TruthMatchMetrics(trackid, comp, bestPur, 1., 1., fTrueDepositsMap[trackid],
                             fTrueTrackInfosMap[trackid]);
  }

  std::pair<double, geo::Point_t> CRTBackTrackerAlg::LineTaggerIntersectionPoint(const geo::Point_t &start, 
                                                                                 const geo::Vector_t &dir, 
                                                                                 const CRTTagger &tagger)
  {
    const CoordSet constrainedPlane = CRTCommonUtils::GetTaggerDefinedCoordinate(tagger);
    const CRTTaggerGeo taggerGeo    = fCRTGeoAlg.GetTagger(CRTCommonUtils::GetTaggerName(tagger));
    double k;

    switch(constrainedPlane)
      {
      case kX:
        {
          const double x = (taggerGeo.maxX + taggerGeo.minX) / 2.;
          k =  (x - start.X()) / dir.X();
        }
        break;
      case kY:
        {
          const double y = (taggerGeo.maxY + taggerGeo.minY) / 2.;
          k = (y - start.Y()) / dir.Y();
        }
        break;
      case kZ:
        {
          const double z = (taggerGeo.maxZ + taggerGeo.minZ) / 2.;
          k = (z - start.Z()) / dir.Z();
        }
        break;
      default:
        std::cout << "Tagger not defined in one plane" << std::endl;
        k = std::numeric_limits<double>::max();
        break;
      }
    
    if(!fCRTGeoAlg.IsPointInsideCRTLimits(start + k * dir))
      return {999999., {999999., 999999., 999999.}};
    
    return {k, start + k * dir};
  }

  void CRTBackTrackerAlg::TrueParticlePDGEnergyTime(const int trackID, int &pdg, double &energy, double &time)
  {
    const simb::MCParticle* particle = particleInv->TrackIdToParticle_P(trackID);

    pdg    = particle == NULL ? -std::numeric_limits<int>::max()    : particle->PdgCode();
    energy = particle == NULL ? -std::numeric_limits<double>::max() : particle->E();
    time   = particle == NULL ? -std::numeric_limits<double>::max() : particle->T();
  }
}

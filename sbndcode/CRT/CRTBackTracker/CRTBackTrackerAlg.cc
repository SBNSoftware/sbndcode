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
    fMCPIDEsEnergyPerTaggerMap.clear();
    fMCPIDEsEnergyMap.clear();
    fMCPStripHitsMap.clear();
    fTrackIDMotherMap.clear();
    fStripHitMCPMap.clear();

    art::Handle<std::vector<sim::ParticleAncestryMap>> droppedTrackIDMapHandle;
    event.getByLabel(fSimModuleLabel, droppedTrackIDMapHandle);

    for(auto const& droppedTrackIdMap : *droppedTrackIDMapHandle)
      {
        for(auto const& [mother, ids] : droppedTrackIdMap.GetMap())
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
    std::set<int>                       depositTrackIDs;

    std::map<std::pair<int, CRTTagger>, double> categoryToEnergyMap, categoryToXMap, categoryToYMap, 
      categoryToZMap, categoryToTimeMap;
    std::map<std::pair<int, CRTTagger>, uint> categoryToNIDEsMap;

    std::map<int, double> idToEnergyMap, idToXMap, idToYMap, 
      idToZMap, idToTimeMap;
    std::map<int, uint> idToNIDEsMap;

    std::map<std::pair<int, CRTTagger>, double> coreCategoryToEnergyMap, coreCategoryToXMap, coreCategoryToYMap, 
      coreCategoryToZMap, coreCategoryToTimeMap;
    std::map<std::pair<int, CRTTagger>, uint> coreCategoryToNIDEsMap;

    for(auto const ide : ideVec)
      {
        const double x = (ide->entryX + ide->exitX) / 2.;
        const double y = (ide->entryY + ide->exitY) / 2.;
        const double z = (ide->entryZ + ide->exitZ) / 2.;
        const CRTTagger tagger = fCRTGeoAlg.WhichTagger(x, y, z);

        const int rollUpID = RollUpID(ide->trackID);

        depositCategories.insert({rollUpID, tagger});
        depositTrackIDs.insert(rollUpID);
        
        fTrackIDSpacePointRecoMap[{rollUpID, tagger}] = false;

        categoryToEnergyMap[{rollUpID, tagger}] += ide->energyDeposited;
        categoryToXMap[{rollUpID, tagger}]      += x;
        categoryToYMap[{rollUpID, tagger}]      += y;
        categoryToZMap[{rollUpID, tagger}]      += z;
        categoryToTimeMap[{rollUpID, tagger}]   += (ide->entryT + ide->exitT) / 2.;

        idToEnergyMap[rollUpID] += ide->energyDeposited;
        idToXMap[rollUpID]      += x;
        idToYMap[rollUpID]      += y;
        idToZMap[rollUpID]      += z;
        idToTimeMap[rollUpID]   += (ide->entryT + ide->exitT) / 2.;

        coreCategoryToEnergyMap[{ide->trackID, tagger}] += ide->energyDeposited;
        coreCategoryToXMap[{ide->trackID, tagger}]      += x;
        coreCategoryToYMap[{ide->trackID, tagger}]      += y;
        coreCategoryToZMap[{ide->trackID, tagger}]      += z;
        coreCategoryToTimeMap[{ide->trackID, tagger}]   += (ide->entryT + ide->exitT) / 2.;

        if(categoryToNIDEsMap.find({rollUpID, tagger}) == categoryToNIDEsMap.end())
          categoryToNIDEsMap[{rollUpID, tagger}] = 0;

        ++categoryToNIDEsMap[{rollUpID, tagger}];

        if(idToNIDEsMap.find(rollUpID) == idToNIDEsMap.end())
          idToNIDEsMap[rollUpID] = 0;

        ++idToNIDEsMap[rollUpID];

        if(coreCategoryToNIDEsMap.find({ide->trackID, tagger}) == coreCategoryToNIDEsMap.end())
          coreCategoryToNIDEsMap[{ide->trackID, tagger}] = 0;

        ++coreCategoryToNIDEsMap[{ide->trackID, tagger}];
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

        const simb::MCParticle* particle = particleInv->TrackIdToParticle_P(category.first);
        const int pdg = particle == NULL ? -999999 : particle->PdgCode();

        fTrueDepositsPerTaggerMap[category] = TrueDeposit(category.first, pdg, category.second,
                                                          energy, time, x, y, z,
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

        const simb::MCParticle* particle = particleInv->TrackIdToParticle_P(trackID);
        const int pdg = particle == NULL ? -999999 : particle->PdgCode();

        bool found_start = false, found_end = false;
        geo::Point_t start, end;
        double start_time = 0., end_time = 0.;
        unsigned start_i = 0;

        for(unsigned i = 0; i < particle->NumberTrajectoryPoints(); ++i)
          {
            const TVector3 pos = particle->Position(i).Vect();
            const geo::Point_t point(pos.X(), pos.Y(), pos.Z());

            if(!found_start && fCRTGeoAlg.IsPointInsideCRTLimits(point))
              {
                start       = point;
                start_time  = particle->T(i);
                start_i     = i;
                found_start = true;
              }
            else if(fCRTGeoAlg.IsPointInsideCRTLimits(point))
              {
                end       = point;
                end_time  = particle->T(i);
                found_end = true;
              }
          }

        if(!found_start || !found_end)
          {
            const TVector3 pos0 = particle->Position().Vect();
            const geo::Point_t point0(pos0.X(), pos0.Y(), pos0.Z());
            const TVector3 posN = particle->EndPosition().Vect();
            const geo::Point_t pointN(posN.X(), posN.Y(), posN.Z());
            
            const double kmax       = (pointN - point0).R();
            const geo::Vector_t dir = (pointN - point0).Unit();

            const double time0      = particle->T();
            const double timeN      = particle->EndT();
            const double time_diff  = timeN - time0;

            std::vector<std::pair<double, geo::Point_t>> intersects, back_intersects;

            for(unsigned tag_i = 0; tag_i < fCRTGeoAlg.NumTaggers(); ++tag_i)
              {
                CRTTagger tag = static_cast<CRTTagger>(tag_i);

                intersects.push_back(LineTaggerIntersectionPoint(point0, dir, tag));
                back_intersects.push_back(LineTaggerIntersectionPoint(pointN, -dir, tag));
              }

            std::sort(intersects.begin(), intersects.end(),
                      [](const auto &a, const auto &b) { 
                        if(a.first < 0 && b.first > 0 ) return false;
                        else if(a.first > 0 && b.first < 0) return true;
                        else
                          return a.first < b.first;
                      });

            std::sort(back_intersects.begin(), back_intersects.end(),
                      [](const auto &a, const auto &b) { 
                        if(a.first < 0 && b.first > 0 ) return false;
                        else if(a.first > 0 && b.first < 0) return true;
                        else
                          return a.first < b.first;
                      });

            for(auto const inte : intersects)
              {
                std::cout << inte.first << "\t" << inte.second << std::endl;
              }
            std::cout << "Chose: " << intersects.begin()->first << "\t" << intersects.begin()->second << std::endl;
            std::cout << std::endl;
            start = intersects.begin()->second;
            end   = back_intersects.begin()->second;
            
            start_time = time0 + (intersects.begin()->first / kmax) * time_diff;
            end_time   = timeN - (back_intersects.begin()->first / kmax) * time_diff;
          }

        fTrueTrackInfosMap[trackID] = TrueTrackInfo(trackID, pdg, start, end, particle->E(start_i),
                                                    end_time - start_time);

        fTrueDepositsMap[trackID] = TrueDeposit(trackID, pdg, kUndefinedTagger,
                                                energy, time, x, y, z);
      }

    for(auto const ide : ideVec)
      {
        const double x = (ide->entryX + ide->exitX) / 2.;
        const double y = (ide->entryY + ide->exitY) / 2.;
        const double z = (ide->entryZ + ide->exitZ) / 2.;
        const CRTTagger tagger = fCRTGeoAlg.WhichTagger(x, y, z);

        fMCPIDEsEnergyPerTaggerMap[{RollUpID(ide->trackID), tagger}] += ide->energyDeposited;
        fMCPIDEsEnergyMap[RollUpID(ide->trackID)]                    += ide->energyDeposited;   
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

  void CRTBackTrackerAlg::RunRecoStatusChecks(const art::Event &event)
  {
    art::Handle<std::vector<CRTSpacePoint>> spacePointHandle;
    event.getByLabel(fSpacePointModuleLabel, spacePointHandle);

    art::FindManyP<CRTCluster> spacePointsToClusters(spacePointHandle, event, fSpacePointModuleLabel);

    for(unsigned i = 0; i < spacePointHandle->size(); ++i)
      {
        const art::Ptr<CRTSpacePoint> spacepoint(spacePointHandle, i);
        const art::Ptr<CRTCluster> cluster = spacePointsToClusters.at(spacepoint.key())[0];

        TruthMatchMetrics truthmatch = TruthMatching(event, cluster);

        fTrackIDSpacePointRecoMap[{truthmatch.trackid, cluster->Tagger()}] = true;
      } 
  }

  std::map<std::pair<int, CRTTagger>, bool> CRTBackTrackerAlg::GetSpacePointRecoStatusMap()
  {
    return fTrackIDSpacePointRecoMap;
  }

  CRTBackTrackerAlg::TrueDeposit CRTBackTrackerAlg::GetTrueDeposit(std::pair<int, CRTTagger> category)
  {
    return fTrueDepositsPerTaggerMap[category];
  }

  CRTBackTrackerAlg::TruthMatchMetrics CRTBackTrackerAlg::TruthMatching(const art::Event &event, const art::Ptr<CRTStripHit> &stripHit)
  {  
    art::Handle<std::vector<FEBData>> febDataHandle;
    event.getByLabel(fFEBDataModuleLabel, febDataHandle);

    art::Handle<std::vector<CRTStripHit>> stripHitHandle;
    event.getByLabel(fStripHitModuleLabel, stripHitHandle);

    art::FindManyP<sim::AuxDetIDE, FEBTruthInfo> febDataToIDEs(febDataHandle, event, fFEBDataModuleLabel);
    art::FindManyP<FEBData> stripHitToFEBData(stripHitHandle, event, fStripHitModuleLabel);
    const CRTTagger tagger = fCRTGeoAlg.ChannelToTaggerEnum(stripHit->Channel());

    auto const febData = stripHitToFEBData.at(stripHit.key());
    auto const assnIDEVec = febDataToIDEs.at(febData.at(0).key());

    std::map<int, double> idToEnergyMap;
    double totalEnergy = 0., x = 0., y = 0., z = 0., t = 0.;
    uint nides = 0;

    for(unsigned i = 0; i < assnIDEVec.size(); ++i)
      {
        const art::Ptr<sim::AuxDetIDE> ide = assnIDEVec[i];
        const FEBTruthInfo *febTruthInfo = febDataToIDEs.data(febData.at(0).key())[i];
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
            trackid = id;
            bestPur = pur;
            comp    = en / fMCPIDEsEnergyPerTaggerMap[{id, tagger}];
          }
      }

    TrueDeposit deposit(-999999, -999999, tagger, totalEnergy, t, x, y, z);

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
            trackid = id;
            bestPur = pur;
            comp    = en / fMCPIDEsEnergyPerTaggerMap[{id, cluster->Tagger()}];
          }
      }

    double hitComp = idToNHitsMap[trackid] / (double) fMCPStripHitsMap[{trackid, cluster->Tagger()}];
    double hitPur  = idToNHitsMap[trackid] / (double) assnStripHitVec.size();

    return TruthMatchMetrics(trackid, comp, bestPur, hitComp, hitPur, 
                             fTrueDepositsPerTaggerMap[{trackid, cluster->Tagger()}]);
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
}

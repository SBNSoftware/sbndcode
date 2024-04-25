#include "CRTBackTracker.h"

namespace sbnd{

CRTBackTracker::CRTBackTracker(const Config& config){

  this->reconfigure(config);

}


CRTBackTracker::CRTBackTracker(){

}


CRTBackTracker::~CRTBackTracker(){

}

void CRTBackTracker::reconfigure(const Config& config){

  fFEBDataLabel         = config.FEBDataLabel();
  fCRTDataLabel         = config.CRTDataLabel();
  fCRTHitLabel          = config.CRTHitLabel();
  fCRTTrackLabel        = config.CRTTrackLabel();
  fSimModuleLabel       = config.SimModuleLabel();

  fRollupUnsavedIds     = config.RollupUnsavedIds();
  fSearchDistanceRadius = config.SearchDistanceRadius();

  return;

}

void CRTBackTracker::Initialize(const art::Event& event){

  // Clear those data structures!
  fCRTDataTrueIds.clear();
  fFEBDataTrueIds.clear();
  fHitTrueIds.clear();
  fTrackTrueIds.clear();
  fTrackIDMotherMap.clear();

  // new modifications ..
  art::Handle<std::vector<sim::ParticleAncestryMap>> droppedTrackIDMapVecHandle;
  event.getByLabel(fSimModuleLabel, droppedTrackIDMapVecHandle);

  if(droppedTrackIDMapVecHandle.isValid()){
	  for(auto const& droppedTrackIdMap : *droppedTrackIDMapVecHandle){
	    for(auto const& [mother, ids] : droppedTrackIdMap.GetMap()){
        for(auto const& id : ids) fTrackIDMotherMap[id] = mother;
	    }
	  }
  }

  art::Handle<sim::ParticleAncestryMap> droppedTrackIDMapHandle;
  event.getByLabel(fSimModuleLabel, droppedTrackIDMapHandle);

  if(droppedTrackIDMapHandle.isValid()){
	  auto const& droppedTrackIdMap = droppedTrackIDMapHandle->GetMap();
	  for(auto const& [mother, ids] : droppedTrackIdMap){
	    for(auto const& id : ids) fTrackIDMotherMap[id] = mother;
	  }
  }
  // end of new modifications ..



  // Get a true track_id map with FEBData
  art::Handle< std::vector<sbnd::crt::FEBData>> febDataHandle; // Get a handle to the FEB data in the event
  std::vector<art::Ptr<sbnd::crt::FEBData> > febDataList;
  if (event.getByLabel(fFEBDataLabel, febDataHandle))
    art::fill_ptr_vector(febDataList, febDataHandle);
  
  art::FindManyP<sim::AuxDetIDE> findManyIdes(febDataHandle, event, fFEBDataLabel);

  std::map<art::Ptr<sbnd::crt::FEBData>, int> dataPtrMap;

  for(size_t data_i = 0; data_i < febDataList.size(); data_i++){

    dataPtrMap[febDataList[data_i]] = data_i;

    // Get all the true IDs from all the IDEs in the FEBData
    std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(data_i);
    for(size_t i = 0; i < ides.size(); i++){
      int id = ides[i]->trackID;
      if(fRollupUnsavedIds) id = std::abs(id);
      fFEBDataTrueIds[data_i][id] += ides[i]->energyDeposited; // std::map<int, std::map<int, double>> fFEBDataTrueIds
    }

  }

  // Get a true track_id map with CRTHit
  art::Handle< std::vector<sbn::crt::CRTHit>> crtHitHandle;
  std::vector<art::Ptr<sbn::crt::CRTHit> > crtHitList;
  if (event.getByLabel(fCRTHitLabel, crtHitHandle))
    art::fill_ptr_vector(crtHitList, crtHitHandle);

  art::FindManyP<sbnd::crt::FEBData> findManyData(crtHitHandle, event, fCRTHitLabel);

  std::map<art::Ptr<sbn::crt::CRTHit>, int> hitPtrMap;

  for(size_t hit_i = 0; hit_i < crtHitList.size(); hit_i++){

    hitPtrMap[crtHitList[hit_i]] = hit_i;

    std::vector<art::Ptr<sbnd::crt::FEBData>> data = findManyData.at(hit_i);
    for(size_t data_i = 0; data_i < data.size(); data_i++){

      int dataID = dataPtrMap[data[data_i]];

      for(auto const& di : fFEBDataTrueIds[dataID]){
        fHitTrueIds[hit_i][di.first] += di.second;
      }

    }
  }

  // Get a true track_id map with CRTTrack
  art::Handle< std::vector<sbn::crt::CRTTrack>> crtTrackHandle;
  std::vector<art::Ptr<sbn::crt::CRTTrack> > crtTrackList;
  if (event.getByLabel(fCRTTrackLabel, crtTrackHandle))
    art::fill_ptr_vector(crtTrackList, crtTrackHandle);

  art::FindManyP<sbn::crt::CRTHit> findManyHits(crtTrackHandle, event, fCRTTrackLabel);

  for(size_t track_i = 0; track_i < crtTrackList.size(); track_i++){

    std::vector<art::Ptr<sbn::crt::CRTHit>> hits = findManyHits.at(track_i);
    for(size_t hit_i = 0; hit_i < hits.size(); hit_i++){

      int hitID = hitPtrMap[hits[hit_i]];

      for(auto const& hi : fHitTrueIds[hitID]){
        fTrackTrueIds[track_i][hi.first] += hi.second;
      }

    }
  }
}

// Check that two CRT data products are the same
bool CRTBackTracker::FEBDataCompare(const sbnd::crt::FEBData& data1, const sbnd::crt::FEBData& data2){

  if(data1.Coinc() != data2.Coinc()) return false;
  if(data1.Mac5() != data2.Mac5())   return false;
  if(data1.Ts0() != data2.Ts0())     return false;
  if(data1.Ts1() != data2.Ts1())     return false;
  if(data1.ADC() != data2.ADC())     return false;

  return true;
}

// Check that two CRT hits are the same
bool CRTBackTracker::HitCompare(const sbn::crt::CRTHit& hit1, const sbn::crt::CRTHit& hit2){

  if(hit1.ts1_ns != hit2.ts1_ns) return false;
  if(hit1.plane != hit2.plane) return false;
  if(hit1.x_pos != hit2.x_pos) return false;
  if(hit1.y_pos != hit2.y_pos) return false;
  if(hit1.z_pos != hit2.z_pos) return false;
  if(hit1.x_err != hit2.x_err) return false;
  if(hit1.y_err != hit2.y_err) return false;
  if(hit1.z_err != hit2.z_err) return false;
  if(hit1.tagger != hit2.tagger) return false;

  return true;

}

// Check that two CRT tracks are the same
bool CRTBackTracker::TrackCompare(const sbn::crt::CRTTrack& track1, const sbn::crt::CRTTrack& track2){

  if(track1.ts1_ns != track2.ts1_ns) return false;

  if(track1.plane1 != track2.plane1) return false;
  if(track1.x1_pos != track2.x1_pos) return false;
  if(track1.y1_pos != track2.y1_pos) return false;
  if(track1.z1_pos != track2.z1_pos) return false;
  if(track1.x1_err != track2.x1_err) return false;
  if(track1.y1_err != track2.y1_err) return false;
  if(track1.z1_err != track2.z1_err) return false;

  if(track1.plane2 != track2.plane2) return false;
  if(track1.x2_pos != track2.x2_pos) return false;
  if(track1.y2_pos != track2.y2_pos) return false;
  if(track1.z2_pos != track2.z2_pos) return false;
  if(track1.x2_err != track2.x2_err) return false;
  if(track1.y2_err != track2.y2_err) return false;
  if(track1.z2_err != track2.z2_err) return false;

  return true;

}

// Get all the true particle IDs that contributed to the FEB data product 
// Comments: not all 
std::vector<int> CRTBackTracker::AllTrueIds(const art::Event& event, const sbnd::crt::FEBData& data){

  std::vector<int> ids;

  // Get a handle to the CRT data in the event
  auto febDataHandle = event.getValidHandle<std::vector<sbnd::crt::FEBData>>(fFEBDataLabel);
  
  art::FindManyP<sim::AuxDetIDE> findManyIdes(febDataHandle, event, fFEBDataLabel);

  // Find which one matches the data passed to the function
  int data_i = 0, index = 0;
  for(auto const& crtData : (*febDataHandle)){
    if(FEBDataCompare(crtData, data)) data_i = index;
    index++;
  }

  // Get all the true IDs from all the IDEs in the hit
  std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(data_i);
  for(size_t i = 0; i < ides.size(); i++){
    int id = ides[i]->trackID;
    if(fRollupUnsavedIds) id = std::abs(id);
    ids.push_back(id);
  }

  // Remove any repeated IDs
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

  return ids;

}

// Get all the true particle IDs that contributed to the CRT hit
std::vector<int> CRTBackTracker::AllTrueIds(const art::Event& event, const sbn::crt::CRTHit& hit){

  std::vector<int> ids;

  // Get a handle to the CRT hits in the event
  auto crtHitHandle = event.getValidHandle<std::vector<sbn::crt::CRTHit>>(fCRTHitLabel);
  
  art::FindManyP<sbnd::crt::FEBData> findManyData(crtHitHandle, event, fCRTHitLabel);

  // Find which one matches the hit passed to the function
  int hit_i = 0, index = 0;
  for(auto const& crtHit : (*crtHitHandle)){
    if(HitCompare(crtHit, hit)) hit_i = index;
    index++;
  }

  // Get the crt data associated to that hit and the IDEs associate to the data
  std::vector<art::Ptr<sbnd::crt::FEBData>> data = findManyData.at(hit_i);
  art::FindManyP<sim::AuxDetIDE> findManyIdes(data, event, fFEBDataLabel);

  // Get all the true IDs from all the IDEs in the hit
  for(size_t i = 0; i < data.size(); i++){
    std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(i);
    for(size_t j = 0; j < ides.size(); j++){
      int id = ides[j]->trackID;
      if(fRollupUnsavedIds) id = std::abs(id);
      ids.push_back(id);
    }
  }

  // Remove any repeated IDs
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

  return ids;

}

// Get all the true particle IDs that contributed to the CRT track
std::vector<int> CRTBackTracker::AllTrueIds(const art::Event& event, const art::Ptr<sbn::crt::CRTTrack> &track){

  std::vector<int> ids;

  // Get a handle to the CRT tracks in the event
  auto crtTrackHandle = event.getValidHandle<std::vector<sbn::crt::CRTTrack>>(fCRTTrackLabel);
  
  art::FindManyP<sbn::crt::CRTHit> findManyHits(crtTrackHandle, event, fCRTTrackLabel);

  // Find which one matches the track passed to the function
  /*int track_i = 0, index = 0;
  for(auto const& crtTrack : (*crtTrackHandle)){
    if(TrackCompare(crtTrack, track)) track_i = index;
    index++;
  }*/

  // Get the crt hits associated to that hit and the data associate to the hits
  std::vector<art::Ptr<sbn::crt::CRTHit>> hits = findManyHits.at(track.key());
  //art::Handle<std::vector<sbn::crt::CRTHit>> hits = findManyHits.at(track_i);
  art::FindManyP<sbnd::crt::FEBData> findManyFEBData(hits, event, fCRTHitLabel); 

  // Get all the true IDs from all the IDEs in the track
  for(size_t ihit = 0; ihit < hits.size(); ihit++){
    /*std::vector<art::Ptr<sbnd::crt::FEBData>> FEBdata = findManyFEBData.at(ihit);
    art::FindManyP<sim::AuxDetIDE> findManyIdes(FEBdata, event, fFEBDataLabel); */

    std::vector<art::Ptr<sbnd::crt::FEBData>> FEBdata = findManyFEBData.at(ihit);
    art::FindManyP<sim::AuxDetIDE, sbnd::crt::FEBTruthInfo> febdata_to_ides(FEBdata, event, fFEBDataLabel);

    for(size_t feb_i = 0; feb_i < FEBdata.size(); feb_i++){
      //std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(feb_i); 
      
      auto const feb_data = FEBdata[feb_i];
      mf::LogDebug("CrtBackTrack") << "FEB " << feb_i << " with mac " << feb_data->Mac5() << std::endl;
      
      auto ides = febdata_to_ides.at(feb_data.key());
      mf::LogDebug("CrtBackTrack") << "We have " << ides.size() << " IDEs." << std::endl;
      
      for(size_t ide_i = 0; ide_i < ides.size(); ide_i++){ // Loop through AuxDetIDE. 

        const sbnd::crt::FEBTruthInfo *fti = febdata_to_ides.data(feb_i)[ide_i]; // use .data() to access the metadata/third-layer data
        // to check if channel numbers from FEBTruthInfo match with CRTHit.
        if ( fti->GetChannel()== (int)hits.at(ihit)->channel0 || fti->GetChannel()==(int)hits.at(ihit)->channel1 ){ 
          int id = ides[ide_i]->trackID;
          if(fRollupUnsavedIds) id = std::abs(id);
          ids.push_back(id);
          mf::LogDebug("CrtBackTrack") << "ide_i " << ide_i << " ene " << ides[ide_i]->energyDeposited << " fti->GetChannel() " << fti->GetChannel() << std::endl;
        }
      }
    }
  }

  // Remove any repeated IDs
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

  return ids;
}

// Get the true particle ID that contributed the most energy to the CRT data product
int CRTBackTracker::TrueIdFromTotalEnergy(const art::Event& event, const sbnd::crt::FEBData & data){ 

  std::map<int, double> ids;

  // Get a handle to the CRT data in the event
  auto febDataHandle = event.getValidHandle<std::vector<sbnd::crt::FEBData>>(fFEBDataLabel);
  
  art::FindManyP<sim::AuxDetIDE> findManyIdes(febDataHandle, event, fFEBDataLabel);

  // Find which one matches the data passed to the function
  int data_i = 0, index = 0;
  for(auto const& febData : (*febDataHandle)){
    if(FEBDataCompare(febData, data)) data_i = index;
    index++;
  }

  // Get all the true IDs from all the IDEs in the hit
  std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(data_i);
  for(size_t i = 0; i < ides.size(); i++){
    int id = ides[i]->trackID;
    if(fRollupUnsavedIds) id = std::abs(id);
    ids[id] += ides[i]->energyDeposited;
  }

  // Find the true ID that contributed the most energy
  double maxEnergy = -1;
  int trueId = -99999;
  for(auto &id : ids){
    if(id.second > maxEnergy){
      maxEnergy = id.second;
      trueId = id.first;       // Might want to store the purity and completeness info. 
    }
  }

  return trueId;

}

int CRTBackTracker::TrueIdFromFEBDataId(const art::Event& event, int data_i){

  if(fCRTDataTrueIds.find(data_i) != fCRTDataTrueIds.end()){ 
    double maxEnergy = -1;
    int trueId = -99999;
    for(auto &id : fCRTDataTrueIds[data_i]){
      if(id.second > maxEnergy){
        maxEnergy = id.second;
        trueId = id.first;
      }
    }
    return trueId;
  }

  return -99999;
}


// Get the true particle ID that contributed the most energy to the CRT hit
int CRTBackTracker::TrueIdFromTotalEnergy(const art::Event& event, const sbn::crt::CRTHit& hit){

  std::map<int, double> ids;

  // Get a handle to the CRT hits in the event
  auto crtHitHandle = event.getValidHandle<std::vector<sbn::crt::CRTHit>>(fCRTHitLabel);
  
  art::FindManyP<sbnd::crt::FEBData> findManyData(crtHitHandle, event, fCRTHitLabel);

  // Find which one matches the hit passed to the function
  int hit_i = 0, index = 0;
  for(auto const& crtHit : (*crtHitHandle)){
    if(HitCompare(crtHit, hit)) hit_i = index;
    index++;
  }

  // Get the crt data associated to that hit and the IDEs associate to the data
  std::vector<art::Ptr<sbnd::crt::FEBData>> data = findManyData.at(hit_i);
  art::FindManyP<sim::AuxDetIDE> findManyIdes(data, event, fFEBDataLabel);

  // Get all the true IDs from all the IDEs in the hit
  for(size_t i = 0; i < data.size(); i++){
    std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(i);
    for(size_t j = 0; j < ides.size(); j++){
      int id = ides[j]->trackID;
      if(fRollupUnsavedIds) id = std::abs(id);
      ids[id] += ides[j]->energyDeposited;
    }
  }

  // Find the true ID that contributed the most energy
  double maxEnergy = -1;
  int trueId = -99999;
  for(auto &id : ids){
    if(id.second > maxEnergy){
      maxEnergy = id.second;
      trueId = id.first;
    }
  }

  return trueId;

}

int CRTBackTracker::TrueIdFromHitId(const art::Event& event, int hit_i){

  if(fHitTrueIds.find(hit_i) != fHitTrueIds.end()){ 
    double maxEnergy = -1;
    int trueId = -99999;
    for(auto &id : fHitTrueIds[hit_i]){
      if(id.second > maxEnergy){
        maxEnergy = id.second;
        trueId = id.first;
      }
    }
    return trueId;
  }

  return -99999;
}
// Get the true particle ID that contributed the most energy to the CRT track
int CRTBackTracker::TrueIdFromTotalEnergy(const art::Event& event, const sbn::crt::CRTTrack& track){

  std::map<int, double> ids;

  // Get a handle to the CRT tracks in the event
  auto crtTrackHandle = event.getValidHandle<std::vector<sbn::crt::CRTTrack>>(fCRTTrackLabel);
  
  art::FindManyP<sbn::crt::CRTHit> findManyHits(crtTrackHandle, event, fCRTTrackLabel);

  // Find which one matches the track passed to the function
  int track_i = 0, index = 0;
  for(auto const& crtTrack : (*crtTrackHandle)){
    if(TrackCompare(crtTrack, track)) track_i = index;
    index++;
  }

  // Get the crt hits associated to that hit and the data associate to the hits
  std::vector<art::Ptr<sbn::crt::CRTHit>> hits = findManyHits.at(track_i);
  art::FindManyP<sbnd::crt::CRTData> findManyData(hits, event, fCRTHitLabel);

  // Get all the true IDs from all the IDEs in the track
  for(size_t i = 0; i < hits.size(); i++){
    std::vector<art::Ptr<sbnd::crt::CRTData>> data = findManyData.at(i);
    art::FindManyP<sim::AuxDetIDE> findManyIdes(data, event, fCRTDataLabel);
    for(size_t j = 0; j < data.size(); j++){
      std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(j);
      for(size_t k = 0; k < ides.size(); k++){
        int id = ides[k]->trackID;
        if(fRollupUnsavedIds) id = std::abs(id);
        ids[id] += ides[k]->energyDeposited;
      }
    }
  }

  // Find the true ID that contributed the most energy
  double maxEnergy = -1;
  int trueId = -99999;
  for(auto &id : ids){
    if(id.second > maxEnergy){
      maxEnergy = id.second;
      trueId = id.first;
    }
  }

  return trueId;

}

int CRTBackTracker::TrueIdFromTrackId(const art::Event& event, int track_i){

  if(fTrackTrueIds.find(track_i) != fTrackTrueIds.end()){ 
    double maxEnergy = -1;
    int trueId = -99999;
    for(auto &id : fTrackTrueIds[track_i]){
      if(id.second > maxEnergy){
        maxEnergy = id.second;
        trueId = id.first;
      }
    }

    return trueId;
  }

  return -99999;

}






// ------------------------------------ new additions: ------------------------------------
// Get the true particle ID that contributed the most energy to the CRT hit
CRTBackTracker::TruthMatchMetrics CRTBackTracker::TruthMatrixFromTotalEnergy(const art::Event& event, const art::Ptr<sbn::crt::CRTHit> &crtHit){

  // Get a handle to the CRT hits in the event
  auto crtHitHandle   = event.getValidHandle<std::vector<sbn::crt::CRTHit>>(fCRTHitLabel); 
  auto febDataHandle  = event.getValidHandle<std::vector<sbnd::crt::FEBData>>(fFEBDataLabel); 
  
  art::FindManyP<sbnd::crt::FEBData> findManyFEBData(crtHitHandle, event, fCRTHitLabel);
  art::FindManyP<sim::AuxDetIDE, sbnd::crt::FEBTruthInfo> febdata_to_ides(febDataHandle, event, fFEBDataLabel);

  std::map<int, double> idToEnergyMap;
  std::map<int, double> idToDistanceMap;
  double totalEnergy = 0.;

  // Get the FEB data associated to that hit and the data associate to the hits
  std::vector<art::Ptr<sbnd::crt::FEBData>> FEBdataVec = findManyFEBData.at(crtHit.key());
  for(auto const FEBdata : FEBdataVec){
    std::vector<art::Ptr<sim::AuxDetIDE>> ides = febdata_to_ides.at(FEBdata.key());
    //std::cout<<"For FEBdata.key() == "<<FEBdata.key()<<std::endl;
    for(size_t ide_i = 0; ide_i < ides.size(); ide_i++){
      const sbnd::crt::FEBTruthInfo *fti = febdata_to_ides.data(FEBdata.key())[ide_i]; // use .data() to access the metadata/third-layer data
      if ( fti->GetChannel()== (int)(crtHit->channel0 % 32) || fti->GetChannel()==(int)(crtHit->channel1 % 32) ){ 
        
        int id = ides[ide_i]->trackID;
        if(fRollupUnsavedIds) { id = RollUpID(id);} //id = std::abs(id); 

        // check if the ide position is close enough to the reco crt hit. 
        double ide_distance_from_crthit_2d = std::hypot(crtHit->x_pos-ides[ide_i]->entryX, crtHit->y_pos-ides[ide_i]->entryY);
        if( ide_distance_from_crthit_2d < fSearchDistanceRadius){
          idToEnergyMap[id]       += ides[ide_i]->energyDeposited;
          totalEnergy             += ides[ide_i]->energyDeposited;
        }
      }
    }
  }

  // Find the true ID that contributed the most energy
  double bestPur = 0., comp = 0.;
  int trackid_true = -99999;
  //double particle_energy = 0.; int pdg = -1; double deposited_energy_track = 0.; double time = -1.;
  double particle_energy = -99999.; int pdg = -99999; double deposited_energy_track = -99999.; double time = -99999.;
  for(auto const [id, en] : idToEnergyMap){
    double pur = en / totalEnergy;
    if(pur > bestPur){
      trackid_true = id;
      TrueParticlePDGEnergyTime(trackid_true, pdg, particle_energy, time);
      deposited_energy_track = en;
      bestPur = pur;
      comp    = en / particle_energy; 
    }
  }
  return TruthMatchMetrics(trackid_true, comp, bestPur, 1., 1., pdg, particle_energy, deposited_energy_track);
}


// Get the true particle ID that contributed the most energy to the CRT track
CRTBackTracker::TruthMatchMetrics CRTBackTracker::TruthMatrixFromTotalEnergy(const art::Event& event, const art::Ptr<sbn::crt::CRTTrack> &crtTrack){

  // Get a handle to the CRT tracks in the event
  auto crtTrackHandle = event.getValidHandle<std::vector<sbn::crt::CRTTrack>>(fCRTTrackLabel);
  auto crtHitHandle   = event.getValidHandle<std::vector<sbn::crt::CRTHit>>(fCRTHitLabel); 
  auto febDataHandle  = event.getValidHandle<std::vector<sbnd::crt::FEBData>>(fFEBDataLabel); 
  
  art::FindManyP<sbn::crt::CRTHit> findManyHits(crtTrackHandle, event, fCRTTrackLabel);
  art::FindManyP<sbnd::crt::FEBData> findManyFEBData(crtHitHandle, event, fCRTHitLabel);
  art::FindManyP<sim::AuxDetIDE, sbnd::crt::FEBTruthInfo> febdata_to_ides(febDataHandle, event, fFEBDataLabel);

  std::map<int, double> idToEnergyMap;
  std::map<int, double> idToTimeMap;
  double totalEnergy = 0.;

  // Get the crt hits associated to that hit and the data associate to the hits
  std::vector<art::Ptr<sbn::crt::CRTHit>> crtHitVec = findManyHits.at(crtTrack.key());
  for(auto const crtHit : crtHitVec){
    std::vector<art::Ptr<sbnd::crt::FEBData>> FEBdataVec = findManyFEBData.at(crtHit.key());
    for(auto const FEBdata : FEBdataVec){
      std::vector<art::Ptr<sim::AuxDetIDE>> ides = febdata_to_ides.at(FEBdata.key());
      for(size_t ide_i = 0; ide_i < ides.size(); ide_i++){
        const sbnd::crt::FEBTruthInfo *fti = febdata_to_ides.data(FEBdata.key())[ide_i]; // use .data() to access the metadata/third-layer data
        if ( fti->GetChannel()== (int)(crtHit->channel0 % 32) || fti->GetChannel()==(int)(crtHit->channel1 % 32) ){ 
          int id = ides[ide_i]->trackID;
          if(fRollupUnsavedIds) { id = RollUpID(id);} //id = std::abs(id); 

          double ide_distance_from_crthit_2d = std::hypot(crtHit->x_pos-ides[ide_i]->entryX, crtHit->y_pos-ides[ide_i]->entryY);
          if( ide_distance_from_crthit_2d < fSearchDistanceRadius){
            idToEnergyMap[id]       += ides[ide_i]->energyDeposited;
            totalEnergy             += ides[ide_i]->energyDeposited;
            idToTimeMap[id]         += 0.5 * (ides[ide_i]->entryT + ides[ide_i]->exitT);
          }
        }
      }
    }
  }

  // Find the true ID that contributed the most energy
  double bestPur = 0., comp = 0.;
  int trackid_true = -99999;
  double particle_energy = -99999.; int pdg = -99999; double deposited_energy_track = -99999.; double time = -99999.;

  for(auto const [id, en] : idToEnergyMap){
    double pur = en / totalEnergy;
    //double depositedE_time = idToTimeMap[id];
    if(pur > bestPur){
      trackid_true = id;
      TrueParticlePDGEnergyTime(trackid_true, pdg, particle_energy, time);
      //std::cout<<"mcparticle time: "<<time<<"; energy deposition time: "<<depositedE_time<<std::endl;
      deposited_energy_track = en;
      bestPur = pur;
      comp    = en / particle_energy; 
    }
  }
  return TruthMatchMetrics(trackid_true, comp, bestPur, 1., 1., pdg, particle_energy, deposited_energy_track);
}

int CRTBackTracker::RollUpID(const int &id)
{
  if(fTrackIDMotherMap.find(id) != fTrackIDMotherMap.end())
    return fTrackIDMotherMap[id];

  return id;
}

void CRTBackTracker::TrueParticlePDGEnergyTime(const int trackID, int &pdg, double &energy, double &time)
{
  const simb::MCParticle* particle = particleInv->TrackIdToParticle_P(trackID);

  pdg    = particle == NULL ? -std::numeric_limits<int>::max()    : particle->PdgCode();
  energy = particle == NULL ? -std::numeric_limits<double>::max() : particle->E();
  time   = particle == NULL ? -std::numeric_limits<double>::max() : particle->T();
}

}
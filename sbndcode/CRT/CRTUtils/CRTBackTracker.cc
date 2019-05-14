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

  fCRTDataLabel = config.CRTDataLabel();
  fCRTHitLabel = config.CRTHitLabel();
  fCRTTrackLabel = config.CRTTrackLabel();

  return;

}

// Check that two CRT data products are the same
bool CRTBackTracker::DataCompare(const crt::CRTData& data1, const crt::CRTData& data2){

  if(data1.Channel() != data2.Channel()) return false;
  if(data1.T0() != data2.T0()) return false;
  if(data1.T1() != data2.T1()) return false;
  if(data1.ADC() != data2.ADC()) return false;

  return true;

}

// Check that two CRT hits are the same
bool CRTBackTracker::HitCompare(const crt::CRTHit& hit1, const crt::CRTHit& hit2){

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
bool CRTBackTracker::TrackCompare(const crt::CRTTrack& track1, const crt::CRTTrack& track2){

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

// Get all the true particle IDs that contributed to the CRT data product
std::vector<int> CRTBackTracker::AllTrueIds(const art::Event& event, const crt::CRTData& data){

  std::vector<int> ids;

  // Get a handle to the CRT data in the event
  auto crtDataHandle = event.getValidHandle<std::vector<crt::CRTData>>(fCRTDataLabel);
  
  art::FindManyP<sim::AuxDetIDE> findManyIdes(crtDataHandle, event, fCRTDataLabel);

  // Find which one matches the data passed to the function
  int data_i = 0, index = 0;
  for(auto const& crtData : (*crtDataHandle)){
    if(DataCompare(crtData, data)) data_i = index;
    index++;
  }

  // Get all the true IDs from all the IDEs in the hit
  std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(data_i);
  for(size_t i = 0; i < ides.size(); i++){
    ids.push_back(ides[i]->trackID);
  }

  // Remove any repeated IDs
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

  return ids;

}

// Get all the true particle IDs that contributed to the CRT hit
std::vector<int> CRTBackTracker::AllTrueIds(const art::Event& event, const crt::CRTHit& hit){

  std::vector<int> ids;

  // Get a handle to the CRT hits in the event
  auto crtHitHandle = event.getValidHandle<std::vector<crt::CRTHit>>(fCRTHitLabel);
  
  art::FindManyP<crt::CRTData> findManyData(crtHitHandle, event, fCRTHitLabel);

  // Find which one matches the hit passed to the function
  int hit_i = 0, index = 0;
  for(auto const& crtHit : (*crtHitHandle)){
    if(HitCompare(crtHit, hit)) hit_i = index;
    index++;
  }

  // Get the crt data associated to that hit and the IDEs associate to the data
  std::vector<art::Ptr<crt::CRTData>> data = findManyData.at(hit_i);
  art::FindManyP<sim::AuxDetIDE> findManyIdes(data, event, fCRTDataLabel);

  // Get all the true IDs from all the IDEs in the hit
  for(size_t i = 0; i < data.size(); i++){
    std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(i);
    for(size_t j = 0; j < ides.size(); j++){
      ids.push_back(ides[j]->trackID);
    }
  }

  // Remove any repeated IDs
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

  return ids;

}

// Get all the true particle IDs that contributed to the CRT track
std::vector<int> CRTBackTracker::AllTrueIds(const art::Event& event, const crt::CRTTrack& track){

  std::vector<int> ids;

  // Get a handle to the CRT tracks in the event
  auto crtTrackHandle = event.getValidHandle<std::vector<crt::CRTTrack>>(fCRTTrackLabel);
  
  art::FindManyP<crt::CRTHit> findManyHits(crtTrackHandle, event, fCRTTrackLabel);

  // Find which one matches the track passed to the function
  int track_i = 0, index = 0;
  for(auto const& crtTrack : (*crtTrackHandle)){
    if(TrackCompare(crtTrack, track)) track_i = index;
    index++;
  }

  // Get the crt hits associated to that hit and the data associate to the hits
  std::vector<art::Ptr<crt::CRTHit>> hits = findManyHits.at(track_i);
  art::FindManyP<crt::CRTData> findManyData(hits, event, fCRTHitLabel);

  // Get all the true IDs from all the IDEs in the track
  for(size_t i = 0; i < hits.size(); i++){
    std::vector<art::Ptr<crt::CRTData>> data = findManyData.at(i);
    art::FindManyP<sim::AuxDetIDE> findManyIdes(data, event, fCRTDataLabel);
    for(size_t j = 0; j < data.size(); j++){
      std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(j);
      for(size_t k = 0; k < ides.size(); k++){
        ids.push_back(ides[k]->trackID);
      }
    }
  }

  // Remove any repeated IDs
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

  return ids;

}

// Get the true particle ID that contributed the most energy to the CRT data product
int CRTBackTracker::TrueIdFromTotalEnergy(const art::Event& event, const crt::CRTData& data){

  std::map<int, double> ids;

  // Get a handle to the CRT data in the event
  auto crtDataHandle = event.getValidHandle<std::vector<crt::CRTData>>(fCRTDataLabel);
  
  art::FindManyP<sim::AuxDetIDE> findManyIdes(crtDataHandle, event, fCRTDataLabel);

  // Find which one matches the data passed to the function
  int data_i = 0, index = 0;
  for(auto const& crtData : (*crtDataHandle)){
    if(DataCompare(crtData, data)) data_i = index;
    index++;
  }

  // Get all the true IDs from all the IDEs in the hit
  std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(data_i);
  for(size_t i = 0; i < ides.size(); i++){
    ids[ides[i]->trackID] += ides[i]->energyDeposited;
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

// Get the true particle ID that contributed the most energy to the CRT hit
int CRTBackTracker::TrueIdFromTotalEnergy(const art::Event& event, const crt::CRTHit& hit){

  std::map<int, double> ids;

  // Get a handle to the CRT hits in the event
  auto crtHitHandle = event.getValidHandle<std::vector<crt::CRTHit>>(fCRTHitLabel);
  
  art::FindManyP<crt::CRTData> findManyData(crtHitHandle, event, fCRTHitLabel);

  // Find which one matches the hit passed to the function
  int hit_i = 0, index = 0;
  for(auto const& crtHit : (*crtHitHandle)){
    if(HitCompare(crtHit, hit)) hit_i = index;
    index++;
  }

  // Get the crt data associated to that hit and the IDEs associate to the data
  std::vector<art::Ptr<crt::CRTData>> data = findManyData.at(hit_i);
  art::FindManyP<sim::AuxDetIDE> findManyIdes(data, event, fCRTDataLabel);

  // Get all the true IDs from all the IDEs in the hit
  for(size_t i = 0; i < data.size(); i++){
    std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(i);
    for(size_t j = 0; j < ides.size(); j++){
      ids[ides[j]->trackID] += ides[j]->energyDeposited;
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

// Get the true particle ID that contributed the most energy to the CRT track
int CRTBackTracker::TrueIdFromTotalEnergy(const art::Event& event, const crt::CRTTrack& track){

  std::map<int, double> ids;

  // Get a handle to the CRT tracks in the event
  auto crtTrackHandle = event.getValidHandle<std::vector<crt::CRTTrack>>(fCRTTrackLabel);
  
  art::FindManyP<crt::CRTHit> findManyHits(crtTrackHandle, event, fCRTTrackLabel);

  // Find which one matches the track passed to the function
  int track_i = 0, index = 0;
  for(auto const& crtTrack : (*crtTrackHandle)){
    if(TrackCompare(crtTrack, track)) track_i = index;
    index++;
  }

  // Get the crt hits associated to that hit and the data associate to the hits
  std::vector<art::Ptr<crt::CRTHit>> hits = findManyHits.at(track_i);
  art::FindManyP<crt::CRTData> findManyData(hits, event, fCRTHitLabel);

  // Get all the true IDs from all the IDEs in the track
  for(size_t i = 0; i < hits.size(); i++){
    std::vector<art::Ptr<crt::CRTData>> data = findManyData.at(i);
    art::FindManyP<sim::AuxDetIDE> findManyIdes(data, event, fCRTDataLabel);
    for(size_t j = 0; j < data.size(); j++){
      std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(j);
      for(size_t k = 0; k < ides.size(); k++){
        ids[ides[k]->trackID] += ides[k]->energyDeposited;
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

}

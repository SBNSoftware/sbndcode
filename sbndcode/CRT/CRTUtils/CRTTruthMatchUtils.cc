#include "CRTTruthMatchUtils.h"

// =============================== UTILITY FUNCTIONS ==============================

std::vector<int> CRTTruthMatchUtils::AllTrueIds(art::Handle<std::vector<sbnd::crt::CRTHit>> hitHandle, const art::Event& event, art::InputTag hitLabel, int hit_i){

  std::vector<int> ids;

  art::FindManyP<sbnd::crt::CRTData> findManyData(hitHandle, event, hitLabel);
  std::vector<art::Ptr<sbnd::crt::CRTData>> data = findManyData.at(hit_i);
  for(size_t i = 0; i < data.size(); i++){
    ids.push_back(data[i]->TrackID());
  }
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

  return ids;

}


std::vector<int> CRTTruthMatchUtils::AllTrueIds(art::Handle<std::vector<sbnd::crt::CRTTrack>> trackHandle, const art::Event& event, art::InputTag trackLabel, art::InputTag hitLabel, int track_i){

  std::vector<int> ids;

  art::FindManyP<sbnd::crt::CRTHit> findManyHits(trackHandle, event, trackLabel);
  std::vector<art::Ptr<sbnd::crt::CRTHit>> hits = findManyHits.at(track_i);
  art::FindManyP<sbnd::crt::CRTData> findManyData(hits, event, hitLabel);
  for(size_t i = 0; i < hits.size(); i++){
    std::vector<art::Ptr<sbnd::crt::CRTData>> data = findManyData.at(i);
    for(size_t j = 0; j < data.size(); j++){
      ids.push_back(data[j]->TrackID());
    }
  }
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

  return ids;

}

std::vector<int> CRTTruthMatchUtils::AllTrueIds(art::Handle<std::vector<sbnd::crt::CRTData>> dataHandle, const art::Event& event, art::InputTag dataLabel, int data_i){

  std::vector<int> ids;

  art::FindManyP<sim::AuxDetIDE> findManyIdes(dataHandle, event, dataLabel);
  std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(data_i);
  for(size_t i = 0; i < ides.size(); i++){
    ids.push_back(ides[i]->trackID);
  }
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

  return ids;

}


std::vector<int> CRTTruthMatchUtils::AllTrueIds(art::Handle<std::vector<sbnd::crt::CRTHit>> hitHandle, const art::Event& event, art::InputTag hitLabel, art::InputTag dataLabel, int hit_i){

  std::vector<int> ids;

  art::FindManyP<sbnd::crt::CRTData> findManyData(hitHandle, event, hitLabel);
  std::vector<art::Ptr<sbnd::crt::CRTData>> data = findManyData.at(hit_i);
  art::FindManyP<sim::AuxDetIDE> findManyIdes(data, event, dataLabel);
  for(size_t i = 0; i < data.size(); i++){
    std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(i);
    for(size_t j = 0; j < ides.size(); j++){
      ids.push_back(ides[j]->trackID);
    }
  }
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

  return ids;

}


std::vector<int> CRTTruthMatchUtils::AllTrueIds(art::Handle<std::vector<sbnd::crt::CRTTrack>> trackHandle, const art::Event& event, art::InputTag trackLabel, art::InputTag hitLabel, art::InputTag dataLabel, int track_i){

  std::vector<int> ids;

  art::FindManyP<sbnd::crt::CRTHit> findManyHits(trackHandle, event, trackLabel);
  std::vector<art::Ptr<sbnd::crt::CRTHit>> hits = findManyHits.at(track_i);
  art::FindManyP<sbnd::crt::CRTData> findManyData(hits, event, hitLabel);
  for(size_t i = 0; i < hits.size(); i++){
    std::vector<art::Ptr<sbnd::crt::CRTData>> data = findManyData.at(i);
    art::FindManyP<sim::AuxDetIDE> findManyIdes(data, event, dataLabel);
    for(size_t j = 0; j < data.size(); j++){
      std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(j);
      for(size_t k = 0; k < ides.size(); k++){
        ids.push_back(ides[k]->trackID);
      }
    }
  }
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

  return ids;

}

int CRTTruthMatchUtils::TrueIdFromTotalEnergy(art::Handle<std::vector<sbnd::crt::CRTData>> dataHandle, const art::Event& event, art::InputTag dataLabel, int data_i){

  std::map<int, double> ids;

  art::FindManyP<sim::AuxDetIDE> findManyIdes(dataHandle, event, dataLabel);
  std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(data_i);
  for(size_t i = 0; i < ides.size(); i++){
    ids[ides[i]->trackID] += ides[i]->energyDeposited;
  }

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

int CRTTruthMatchUtils::TrueIdFromTotalEnergy(art::Handle<std::vector<sbnd::crt::CRTHit>> hitHandle, const art::Event& event, art::InputTag hitLabel, art::InputTag dataLabel, int hit_i){

  std::map<int, double> ids;

  art::FindManyP<sbnd::crt::CRTData> findManyData(hitHandle, event, hitLabel);
  std::vector<art::Ptr<sbnd::crt::CRTData>> data = findManyData.at(hit_i);
  art::FindManyP<sim::AuxDetIDE> findManyIdes(data, event, dataLabel);
  for(size_t i = 0; i < data.size(); i++){
    std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(i);
    for(size_t j = 0; j < ides.size(); j++){
      ids[ides[j]->trackID] += ides[j]->energyDeposited;
    }
  }

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


int CRTTruthMatchUtils::TrueIdFromTotalEnergy(art::Handle<std::vector<sbnd::crt::CRTTrack>> trackHandle, const art::Event& event, art::InputTag trackLabel, art::InputTag hitLabel, art::InputTag dataLabel, int track_i){

  std::map<int, double> ids;

  art::FindManyP<sbnd::crt::CRTHit> findManyHits(trackHandle, event, trackLabel);
  std::vector<art::Ptr<sbnd::crt::CRTHit>> hits = findManyHits.at(track_i);
  art::FindManyP<sbnd::crt::CRTData> findManyData(hits, event, hitLabel);
  for(size_t i = 0; i < hits.size(); i++){
    std::vector<art::Ptr<sbnd::crt::CRTData>> data = findManyData.at(i);
    art::FindManyP<sim::AuxDetIDE> findManyIdes(data, event, dataLabel);
    for(size_t j = 0; j < data.size(); j++){
      std::vector<art::Ptr<sim::AuxDetIDE>> ides = findManyIdes.at(j);
      for(size_t k = 0; k < ides.size(); k++){
        ids[ides[k]->trackID] += ides[k]->energyDeposited;
      }
    }
  }

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

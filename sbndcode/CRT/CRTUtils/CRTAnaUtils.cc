#include "CRTAnaUtils.h"

namespace sbnd{

// =============================== UTILITY FUNCTIONS ==============================

std::vector<std::vector<art::Ptr<sbnd::crt::CRTHit>>> CRTAnaUtils::CreateCRTTzeros(std::vector<art::Ptr<sbnd::crt::CRTHit>> crtHits, double fTimeLimit){

  std::vector<std::vector<art::Ptr<sbnd::crt::CRTHit>>> CRTTzeroVect;
  int iflag[2000] = {};

  // Loop over crt hits
  for(size_t i = 0; i < crtHits.size(); i++){
    if(iflag[i] == 0){
      std::vector<art::Ptr<sbnd::crt::CRTHit>> CRTTzero;
      double time_ns_A = crtHits[i]->ts1_ns;
      iflag[i]=1;
      CRTTzero.push_back(crtHits[i]);

      // Sort into a Tzero collection
      // Loop over all the other CRT hits
      for(size_t j = i+1; j < crtHits.size(); j++){
        if(iflag[j] == 0){
          // If ts1_ns - ts1_ns < diff then put them in a vector
          double time_ns_B = crtHits[j]->ts1_ns;
          double diff = std::abs(time_ns_B - time_ns_A) * 1e-3; // [us]
          if(diff < fTimeLimit){
            iflag[j] = 1;
            CRTTzero.push_back(crtHits[j]);
          }
        }
      }
      CRTTzeroVect.push_back(CRTTzero);
    }
  }

  return CRTTzeroVect;

}

std::vector<sbnd::crt::CRTTrack> CRTAnaUtils::CreateCRTTracks(std::vector<std::vector<art::Ptr<sbnd::crt::CRTHit>>> crtTzeros, double fAverageHitDist, bool fUseTopPlane, double fDistanceLimit){

  std::vector<sbnd::crt::CRTTrack> returnTracks;

  sbnd::CRTTrackRecoAlg trackAlg(fAverageHitDist, fDistanceLimit);

  // Loop over tzeros
  for(size_t i = 0; i < crtTzeros.size(); i++){

    //loop over hits for this tzero, sort by tagger
    std::map<std::string, std::vector<art::Ptr<sbnd::crt::CRTHit>>> tagHits;
    for (size_t ah = 0; ah< crtTzeros[i].size(); ++ah){        
      std::string tag = crtTzeros[i][ah]->tagger;       
      tagHits[tag].push_back(crtTzeros[i][ah]);
    } // loop over hits
    
    //loop over planes and calculate average hits
    std::vector<sbnd::crt::CRTHit> hits;
    for (auto &keyVal : tagHits){
      std::string tag = keyVal.first;
      std::vector<sbnd::crt::CRTHit> ahits = trackAlg.AverageHits(tagHits[tag]);
      if(fUseTopPlane && tag == "volTaggerTopHigh_0"){ 
        hits.insert(hits.end(), ahits.begin(), ahits.end());
      }
      else if(tag != "volTaggerTopHigh_0"){ 
        hits.insert(hits.end(), ahits.begin(), ahits.end());
      }
    }

    std::vector<sbnd::crt::CRTTrack> tempTracks = trackAlg.CreateTracks(hits);
    returnTracks.insert(returnTracks.end(), tempTracks.begin(), tempTracks.end());

   }

  return returnTracks;

}

std::vector<sbnd::crt::CRTTrack> CRTAnaUtils::CreateCRTTracks(std::vector<art::Ptr<sbnd::crt::CRTHit>> crtHits, double fTimeLimit, double fAverageHitDist, bool fUseTopPlane, double fDistanceLimit){

    auto crtTzeros = CreateCRTTzeros(crtHits, fTimeLimit);
    auto crtTracks = CreateCRTTracks(crtTzeros, fAverageHitDist, fUseTopPlane, fDistanceLimit);

    return crtTracks;

}

double CRTAnaUtils::T0FromCRTHits(recob::Track tpcTrack, std::vector<sbnd::crt::CRTHit> crtHits, int tpc, double fMinTrackLength, double fTrackDirectionFrac, double fDistanceLimit){

  sbnd::CRTT0MatchAlg t0Alg;
  
  if (tpcTrack.Length() < fMinTrackLength) return -99999; 

  // Calculate direction as an average over directions
  std::pair<TVector3, TVector3> startEndDir = t0Alg.TrackDirectionAverage(tpcTrack, fTrackDirectionFrac);
  TVector3 startDir = startEndDir.first;
  TVector3 endDir = startEndDir.second;

  auto start = tpcTrack.Vertex<TVector3>();
  auto end = tpcTrack.End<TVector3>();

  // ====================== Matching Algorithm ========================== //
  // Get the allowed t0 range
  std::pair<double, double> t0MinMax = t0Alg.TrackT0Range(start.X(), end.X(), tpc);
  std::vector<std::pair<double, double>> t0Candidates;

  // Loop over all the CRT hits
  for(auto &crtHit : crtHits){
    // Check if hit is within the allowed t0 range
    double crtTime = ((double)(int)crtHit.ts1_ns) * 1e-3; // [us]
    if (!(crtTime >= t0MinMax.first - 10. && crtTime <= t0MinMax.second + 10.)) continue;
    TVector3 crtPoint(crtHit.x_pos, crtHit.y_pos, crtHit.z_pos);
  
    // Calculate the distance between the crossing point and the CRT hit
    double startDist = t0Alg.DistOfClosestApproach(start, startDir, crtHit, tpc, crtTime);
    // If the distance is less than some limit record the time
    if (startDist < fDistanceLimit){ 
      t0Candidates.push_back(std::make_pair(startDist, crtTime));
    }
  
    // Calculate the distance between the crossing point and the CRT hit
    double endDist = t0Alg.DistOfClosestApproach(end, endDir, crtHit, tpc, crtTime);
    // If the distance is less than some limit record the time
    if (endDist < fDistanceLimit){ 
      t0Candidates.push_back(std::make_pair(endDist, crtTime));
    }
  
  }

  // Sort the candidates by distance
  std::sort(t0Candidates.begin(), t0Candidates.end(), [](auto& left, auto& right){
            return left.first < right.first;});
  double bestTime = -99999;
  if(t0Candidates.size()>0) {
    bestTime = t0Candidates[0].second;
  }

  return bestTime;

}

std::pair<crt::CRTHit, double> CRTAnaUtils::ClosestCRTHit(recob::Track tpcTrack, std::vector<sbnd::crt::CRTHit> crtHits, int tpc, double fTrackDirectionFrac){

  sbnd::CRTT0MatchAlg t0Alg;
  
  // Calculate direction as an average over directions
  std::pair<TVector3, TVector3> startEndDir = t0Alg.TrackDirectionAverage(tpcTrack, fTrackDirectionFrac);
  TVector3 startDir = startEndDir.first;
  TVector3 endDir = startEndDir.second;

  auto start = tpcTrack.Vertex<TVector3>();
  auto end = tpcTrack.End<TVector3>();

  // ====================== Matching Algorithm ========================== //
  // Get the allowed t0 range
  std::pair<double, double> t0MinMax = t0Alg.TrackT0Range(start.X(), end.X(), tpc);
  std::vector<std::pair<crt::CRTHit, double>> t0Candidates;

  // Loop over all the CRT hits
  for(auto &crtHit : crtHits){
    // Check if hit is within the allowed t0 range
    double crtTime = ((double)(int)crtHit.ts1_ns) * 1e-3; // [us]
    if (!(crtTime >= t0MinMax.first - 10. && crtTime <= t0MinMax.second + 10.)) continue;
    TVector3 crtPoint(crtHit.x_pos, crtHit.y_pos, crtHit.z_pos);
  
    // Calculate the distance between the crossing point and the CRT hit
    double startDist = t0Alg.DistOfClosestApproach(start, startDir, crtHit, tpc, crtTime);
    double endDist = t0Alg.DistOfClosestApproach(end, endDir, crtHit, tpc, crtTime);
    // If the distance is less than some limit record the time
    if ((crtPoint-start).Mag() < (crtPoint-end).Mag()){ 
      t0Candidates.push_back(std::make_pair(crtHit, startDist));
    }
    else{
      t0Candidates.push_back(std::make_pair(crtHit, endDist));
    }
  
  }

  // Sort the candidates by distance
  std::sort(t0Candidates.begin(), t0Candidates.end(), [](auto& left, auto& right){
            return left.second < right.second;});

  if(t0Candidates.size() > 0){
    return t0Candidates[0];
  }
  crt::CRTHit hit;
  return std::make_pair(hit, -99999);

}

double CRTAnaUtils::T0FromCRTTracks(recob::Track tpcTrack, std::vector<sbnd::crt::CRTTrack> crtTracks, int tpc, double fMaxAngleDiff, double fMaxDistance){

  sbnd::CRTTrackMatchAlg trackAlg;

  // Get the length, angle and start and end position of the TPC track
  TVector3 tpcStart = tpcTrack.Vertex<TVector3>();
  TVector3 tpcEnd = tpcTrack.End<TVector3>();
  double tpcTheta = (tpcStart - tpcEnd).Theta();
  double tpcPhi = (tpcStart - tpcEnd).Phi();

  int crtIndex = 0;
  std::vector<sbnd::RecoCRTTrack> recoCrtTracks;

  // Transform CRT tracks to expected TPC reconstructed tracks
  for (auto& crtTrack : crtTracks){

    //Check that crt track crosses tpc volume, if not skip it
    if(!trackAlg.CrossesTPC(crtTrack)){ crtIndex++; continue; }
 
    std::vector<sbnd::RecoCRTTrack> tempTracks = trackAlg.CrtToRecoTrack(crtTrack, crtIndex);
    recoCrtTracks.insert(recoCrtTracks.end(), tempTracks.begin(), tempTracks.end());
 
    crtIndex++;
  }
 
  std::vector<std::pair<double, double>> crtTpcMatchCandidates;
  // Loop over the reco crt tracks
  for (auto const& recoCrtTrack : recoCrtTracks){
 
    // Get the length, angle and start and end position of the TPC track
    TVector3 crtStart = recoCrtTrack.start;
    TVector3 crtEnd = recoCrtTrack.end;
    double crtTheta = (crtStart - crtEnd).Theta();
    double crtPhi = (crtStart - crtEnd).Phi();
 
    // Find the difference with the CRT track
    double dDist1 = (crtStart-tpcStart).Mag();
    double dDist2 = (crtEnd-tpcEnd).Mag();
    if(std::max((crtStart-tpcStart).Mag(), (crtEnd-tpcEnd).Mag()) > 
       std::max((crtStart-tpcEnd).Mag(), (crtEnd-tpcStart).Mag())){
      crtTheta = (crtEnd - crtStart).Theta();
      crtPhi = (crtEnd - crtStart).Phi();
      dDist1 = (crtEnd-tpcStart).Mag();
      dDist2 = (crtStart-tpcEnd).Mag();
    }
    double dTheta = atan2(sin(tpcTheta - crtTheta), cos(tpcTheta - crtTheta));
    double dPhi = atan2(sin(tpcPhi - crtPhi), cos(tpcPhi - crtPhi));
 
    // Do the actual matching
    if(std::abs(dTheta) < fMaxAngleDiff && std::abs(dPhi) < fMaxAngleDiff && 
       tpc == recoCrtTrack.tpc && (dDist1<fMaxDistance||dDist2<fMaxDistance)){
      crtTpcMatchCandidates.push_back(std::make_pair(recoCrtTrack.trueTime, std::abs(dTheta)));
    }
 
  }

  // Choose the track which matches the closest
  double bestTime = -99999;
  if(crtTpcMatchCandidates.size() > 0){
    std::sort(crtTpcMatchCandidates.begin(), crtTpcMatchCandidates.end(), [](auto& left, auto& right){
              return left.second < right.second;});
    bestTime = crtTpcMatchCandidates[0].first;
  }
  
  return bestTime;

}

std::vector<double> CRTAnaUtils::ApaT0sFromCRTHits(std::vector<art::Ptr<crt::CRTHit>> crtHits, double fTimeLimit){

  TPCGeoAlg fTpcGeo;

  std::vector<std::vector<art::Ptr<crt::CRTHit>>> crtT0Ptr = CreateCRTTzeros(crtHits, fTimeLimit);
  std::vector<double> stopT0;
  for(size_t i = 0; i < crtT0Ptr.size(); i++){
    double t0 = 0;
    double npts = 0;
    for(size_t j = 0; j < crtT0Ptr[i].size(); j++){

      if(crtT0Ptr[i][j]->tagger == "volTaggerTopHigh_0") continue;

      if(crtT0Ptr[i][j]->tagger == "volTaggerBot_0") continue;

      if(crtT0Ptr[i][j]->y_pos < fTpcGeo.MinY()) continue;

      if(std::abs(crtT0Ptr[i][j]->x_pos) < fTpcGeo.MaxX()) continue;

      t0 += (double)(int)crtT0Ptr[i][0]->ts1_ns*1e-3; // [us]
      npts++;
    }

    if(t0 != 0) stopT0.push_back(t0/npts);
  }
  return stopT0;
}

std::vector<double> CRTAnaUtils::ApaT0sFromCRTTracks(std::vector<crt::CRTTrack> crtTracks){

  CRTTrackMatchAlg trackAlg;

  std::vector<double> throughT0;

  for(auto const& crtTrack : crtTracks){
    double trackT0 = (double)(int)crtTrack.ts1_ns*1e-3;
    if(trackAlg.CrossesAPA(crtTrack)) throughT0.push_back(trackT0);
  }

  return throughT0;
}

}

#include "CRTAnaUtils.h"

// ==================== FORWARD DECLARATIONS OF INTERNAL FUNCTIONS ===============

struct RecoCRTTrack{
  int crtID;
  int tpc;
  TVector3 start; // [cm]
  TVector3 end; // [cm]
  double trueTime; // [us]
  bool complete;
};

// FOR TRACKING

sbnd::crt::CRTHit FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t, 
                                         std::vector<std::pair<int,float>>> tpesmap, float peshit, double time, 
                                         double x, double ex, double y, double ey, double z, double ez, std::string tagger);

sbnd::crt::CRTTrack FillCrtTrack(sbnd::crt::CRTHit hit1, sbnd::crt::CRTHit hit2, bool complete);

std::vector<sbnd::crt::CRTHit> AverageHits(std::vector<sbnd::crt::CRTHit> hits, double fAverageHitDist);

sbnd::crt::CRTHit DoAverage(std::vector<sbnd::crt::CRTHit> hits);

TVector3 CrossPoint(sbnd::crt::CRTHit hit, TVector3 start, TVector3 diff);

// FOR T0 MATCHING

std::pair<double, double> TrackT0Range(double startX, double endX, int tpc);

double DistOfClosestApproach(TVector3 trackPos, TVector3 trackDir, sbnd::crt::CRTHit crtHit, int tpc, double t0);

std::vector<RecoCRTTrack> CrtToRecoTrack(sbnd::crt::CRTTrack, int id);

std::vector<RecoCRTTrack> CreateRecoCRTTrack(TVector3 start, TVector3 end, double shift, 
                                             int tpc, int id, double time, bool complete);
                                                 
bool CrossesTPC(sbnd::crt::CRTTrack track);

// =============================== UTILITY FUNCTIONS ==============================

std::vector<std::vector<sbnd::crt::CRTHit>> CRTAnaUtils::CreateCRTTzeros(std::vector<sbnd::crt::CRTHit> crtHits, double fTimeLimit){

  std::vector<std::vector<sbnd::crt::CRTHit>> CRTTzeroVect;
  int iflag[2000] = {};

  // Loop over crt hits
  for(size_t i = 0; i < crtHits.size(); i++){
    if(iflag[i] == 0){
      std::vector<sbnd::crt::CRTHit> CRTTzero;
      double time_ns_A = crtHits[i].ts1_ns;
      iflag[i]=1;
      CRTTzero.push_back(crtHits[i]);

      // Sort into a Tzero collection
      // Loop over all the other CRT hits
      for(size_t j = i+1; j < crtHits.size(); j++){
        if(iflag[j] == 0){
          // If ts1_ns - ts1_ns < diff then put them in a vector
          double time_ns_B = crtHits[j].ts1_ns;
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

std::vector<sbnd::crt::CRTTrack> CRTAnaUtils::CreateCRTTracks(std::vector<std::vector<sbnd::crt::CRTHit>> crtTzeros, double fAverageHitDist, bool fUseTopPlane, double fDistanceLimit){

  std::vector<sbnd::crt::CRTTrack> returnTracks;

  // Loop over tzeros
  for(size_t i = 0; i < crtTzeros.size(); i++){

    //loop over hits for this tzero, sort by tagger
    std::map<std::string, std::vector<sbnd::crt::CRTHit>> tagHits;
    for (size_t ah = 0; ah< crtTzeros[i].size(); ++ah){        
      std::string tag = crtTzeros[i][ah].tagger;       
      tagHits[tag].push_back(crtTzeros[i][ah]);
    } // loop over hits
    
    //loop over planes and calculate average hits
    std::vector<sbnd::crt::CRTHit> hits;
    for (auto &keyVal : tagHits){
      std::string tag = keyVal.first;
      std::vector<sbnd::crt::CRTHit> ahits = AverageHits(tagHits[tag], fAverageHitDist);
      if(fUseTopPlane && tag == "volTaggerTopHigh_0"){ 
        hits.insert(hits.end(), ahits.begin(), ahits.end());
      }
      else if(tag != "volTaggerTopHigh_0"){ 
        hits.insert(hits.end(), ahits.begin(), ahits.end());
      }
    }

    //Store list of hit pairs with distance between them
    std::vector<std::pair<std::pair<size_t, size_t>, double>> hitPairDist;
    std::vector<std::pair<size_t, size_t>> usedPairs;
 
    //Calculate the distance between all hits on different planes
    for(size_t i = 0; i < hits.size(); i++){
      sbnd::crt::CRTHit hit1 = hits[i];
 
      for(size_t j = 0; j < hits.size(); j++){
        sbnd::crt::CRTHit hit2 = hits[j];
        std::pair<size_t, size_t> hitPair = std::make_pair(i, j);
        std::pair<size_t, size_t> rhitPair = std::make_pair(j, i);
 
        //Only compare hits on different taggers and don't reuse hits
        if(hit1.tagger!=hit2.tagger && std::find(usedPairs.begin(), usedPairs.end(), rhitPair)==usedPairs.end()){
          //Calculate the distance between hits and store
          TVector3 pos1(hit1.x_pos, hit1.y_pos, hit1.z_pos);
          TVector3 pos2(hit2.x_pos, hit2.y_pos, hit2.z_pos);
          double dist = (pos1 - pos2).Mag();
          usedPairs.push_back(hitPair);
          hitPairDist.push_back(std::make_pair(hitPair, dist));
        }
      }
    }
 
    //Sort map by distance
    std::sort(hitPairDist.begin(), hitPairDist.end(), [](auto& left, auto& right){
              return left.second > right.second;});
 
    //Store potential hit collections + distance along 1D hit
    std::vector<std::pair<std::vector<size_t>, double>> tracks;
    for(size_t i = 0; i < hitPairDist.size(); i++){
      size_t hit_i = hitPairDist[i].first.first;
      size_t hit_j = hitPairDist[i].first.second;
      //Make sure bottom plane hit is always hit_i
      if(hits[hit_j].tagger=="volTaggerBot_0") std::swap(hit_i, hit_j);
      sbnd::crt::CRTHit ihit = hits[hit_i];
      sbnd::crt::CRTHit jhit = hits[hit_j];
 
      //If the bottom plane hit is a 1D hit
      if(ihit.x_err>100. || ihit.z_err>100.){
        double facMax = 1;
        std::vector<size_t> nhitsMax;
        double minDist = 99999;
 
        //Loop over the length of the 1D hit
        for(int i = 0; i<21; i++){
          double fac = (i)/10.;
          std::vector<size_t> nhits;
          double totalDist = 0.;
          TVector3 start(ihit.x_pos-(1.-fac)*ihit.x_err, ihit.y_pos, ihit.z_pos-(1.-fac)*ihit.z_err);
          TVector3 end(jhit.x_pos, jhit.y_pos, jhit.z_pos);
          TVector3 diff = start - end;
 
          //Loop over the rest of the hits
          for(size_t k = 0; k < hits.size(); k++){
            if(k == hit_i || k == hit_j || hits[k].tagger == ihit.tagger || hits[k].tagger == jhit.tagger) continue;
            //Calculate the distance between the track crossing point and the true hit
            TVector3 mid(hits[k].x_pos, hits[k].y_pos, hits[k].z_pos);
            TVector3 cross = CrossPoint(hits[k], start, diff);
            double dist = (cross-mid).Mag();
 
            //If the distance is less than some limit add the hit to the track and record the distance
            if(dist < fDistanceLimit){
              nhits.push_back(k);
              totalDist += dist;
            }
          }
 
          //If the distance down the 1D hit means more hits are included and they are closer to the track record it
          if(nhits.size()>=nhitsMax.size() && totalDist/nhits.size() < minDist){
            nhitsMax = nhits;
            facMax = fac;
            minDist = totalDist/nhits.size();
          }
          nhits.clear();
        }
 
        //Record the track candidate
        std::vector<size_t> trackCand;
        trackCand.push_back(hit_i);
        trackCand.push_back(hit_j);
        trackCand.insert(trackCand.end(), nhitsMax.begin(), nhitsMax.end());
        tracks.push_back(std::make_pair(trackCand, facMax));
      }
 
      //If there is no 1D hit
      else{
        TVector3 start(ihit.x_pos, ihit.y_pos, ihit.z_pos);
        TVector3 end(jhit.x_pos, jhit.y_pos, jhit.z_pos);
        TVector3 diff = start - end;
        std::vector<size_t> trackCand;
        trackCand.push_back(hit_i);
        trackCand.push_back(hit_j);
 
        //Loop over all the other hits
        for(size_t k = 0; k < hits.size(); k++){
          if(k == hit_i || k == hit_j || hits[k].tagger == ihit.tagger || hits[k].tagger == jhit.tagger) continue;
          //Calculate distance to other hits not on the planes of the track hits
          TVector3 mid(hits[k].x_pos, hits[k].y_pos, hits[k].z_pos);
          TVector3 cross = CrossPoint(hits[k], start, diff);
          double dist = (cross-mid).Mag();
 
          //Record any within a certain distance
          if(dist < fDistanceLimit){
            trackCand.push_back(k);
          }
        }
        tracks.push_back(std::make_pair(trackCand, 1));
      }
    }
 
    //Sort track candidates by number of hits
    std::sort(tracks.begin(), tracks.end(), [](auto& left, auto& right){
              return left.first.size() > right.first.size();});
 
    //Record used hits
    std::vector<size_t> usedHits;
    //Loop over candidates
    for(auto& track : tracks){
      size_t hit_i = track.first[0];
      size_t hit_j = track.first[1];
 
      // Make sure the first hit is the top high tagger if there are only two hits
      if(hits[hit_j].tagger=="volTaggerTopHigh_0") std::swap(hit_i, hit_j);
      sbnd::crt::CRTHit ihit = hits[track.first[0]];
      sbnd::crt::CRTHit jhit = hits[track.first[1]];
 
      //Check no hits in track have been used
      bool used = false;
 
      //Loop over hits in track candidate
      for(size_t i = 0; i < track.first.size(); i++){
        //Check if any of the hits have been used
        if(std::find(usedHits.begin(), usedHits.end(), track.first[i]) != usedHits.end()) used=true;
      }
 
      //If any of the hits have already been used skip this track
      if(used) continue;
      ihit.x_pos -= (1.-track.second)*ihit.x_err;
      ihit.z_pos -= (1.-track.second)*ihit.z_err;
      //Create track
      sbnd::crt::CRTTrack crtTrack = FillCrtTrack(ihit, jhit, true);
 
      //If only the top two planes are hit create an incomplete/stopping track
      if(track.first.size()==2 && ihit.tagger == "volTaggerTopHigh_0" && jhit.tagger == "volTaggerTopLow_0"){ 
        crtTrack.complete = false;
      }
 
      returnTracks.push_back(crtTrack);
 
      //Record which hits were used only if the track has more than two hits
      for(size_t i = 0; i < track.first.size(); i++){
        if(track.first.size()>2) usedHits.push_back(track.first[i]);
      }
    }
  }

  return returnTracks;

}

std::vector<sbnd::crt::CRTTrack> CRTAnaUtils::CreateCRTTracks(std::vector<sbnd::crt::CRTHit> crtHits, double fTimeLimit, double fAverageHitDist, bool fUseTopPlane, double fDistanceLimit){

    std::vector<std::vector<sbnd::crt::CRTHit>> crtTzeros = CreateCRTTzeros(crtHits, fTimeLimit);

    std::vector<sbnd::crt::CRTTrack> crtTracks = CreateCRTTracks(crtTzeros, fAverageHitDist, fUseTopPlane, fDistanceLimit);

    return crtTracks;

}

double CRTAnaUtils::T0FromCRTHits(recob::Track tpcTrack, std::vector<sbnd::crt::CRTHit> crtHits, int tpc, double fMinTrackLength, double fTrackDirectionFrac, double fDistanceLimit){
  
  if (tpcTrack.Length() < fMinTrackLength) return -99999; 

  // Calculate direction as an average over directions
  size_t nTrackPoints = tpcTrack.NumberTrajectoryPoints();
  int endPoint = (int)floor(nTrackPoints*fTrackDirectionFrac);
  double xTotStart = 0; double yTotStart = 0; double zTotStart = 0;
  double xTotEnd = 0; double yTotEnd = 0; double zTotEnd = 0;
  for(int i = 0; i < endPoint; i++){
    xTotStart += tpcTrack.DirectionAtPoint(i)[0];
    yTotStart += tpcTrack.DirectionAtPoint(i)[1];
    zTotStart += tpcTrack.DirectionAtPoint(i)[2];
    xTotEnd += tpcTrack.DirectionAtPoint(nTrackPoints - (i+1))[0];
    yTotEnd += tpcTrack.DirectionAtPoint(nTrackPoints - (i+1))[1];
    zTotEnd += tpcTrack.DirectionAtPoint(nTrackPoints - (i+1))[2];
  }
  TVector3 startDir = {-xTotStart/endPoint, -yTotStart/endPoint, -zTotStart/endPoint};
  TVector3 endDir = {xTotEnd/endPoint, yTotEnd/endPoint, zTotEnd/endPoint};

  TVector3 start = tpcTrack.Vertex();
  TVector3 end = tpcTrack.End();

  // ====================== Matching Algorithm ========================== //
  // Get the allowed t0 range
  std::pair<double, double> t0MinMax = TrackT0Range(start.X(), end.X(), tpc);
  std::vector<std::pair<double, double>> t0Candidates;

  // Loop over all the CRT hits
  for(auto &crtHit : crtHits){
    // Check if hit is within the allowed t0 range
    double crtTime = ((double)(int)crtHit.ts1_ns) * 1e-3;
    if (!(crtTime >= t0MinMax.first - 10. && crtTime <= t0MinMax.second + 10.)) continue;
    TVector3 crtPoint(crtHit.x_pos, crtHit.y_pos, crtHit.z_pos);
  
    // Calculate the distance between the crossing point and the CRT hit
    double startDist = DistOfClosestApproach(start, startDir, crtHit, tpc, crtTime);
    // If the distance is less than some limit record the time
    if (startDist < fDistanceLimit){ 
      t0Candidates.push_back(std::make_pair(startDist, crtTime));
    }
  
    // Calculate the distance between the crossing point and the CRT hit
    double endDist = DistOfClosestApproach(end, endDir, crtHit, tpc, crtTime);
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

double CRTAnaUtils::T0FromCRTTracks(recob::Track tpcTrack, std::vector<sbnd::crt::CRTTrack> crtTracks, int tpc, double fMaxAngleDiff, double fMaxDistance){

  // Get the length, angle and start and end position of the TPC track
  TVector3 tpcStart = tpcTrack.Vertex();
  TVector3 tpcEnd = tpcTrack.End();
  double tpcTheta = (tpcStart - tpcEnd).Theta();
  double tpcPhi = (tpcStart - tpcEnd).Phi();

  int crtIndex = 0;
  std::vector<RecoCRTTrack> recoCrtTracks;

  // Transform CRT tracks to expected TPC reconstructed tracks
  for (auto& crtTrack : crtTracks){

    //Check that crt track crosses tpc volume, if not skip it
    if(!CrossesTPC(crtTrack)){ crtIndex++; continue; }
 
    std::vector<RecoCRTTrack> tempTracks = CrtToRecoTrack(crtTrack, crtIndex);
    recoCrtTracks.insert(recoCrtTracks.end(), tempTracks.begin(), tempTracks.end());
 
    crtIndex++;
  }
 
  std::vector<std::pair<int, double>> crtTpcMatchCandidates;
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
  int bestTime = -99999;
  if(crtTpcMatchCandidates.size() > 0){
    std::sort(crtTpcMatchCandidates.begin(), crtTpcMatchCandidates.end(), [](auto& left, auto& right){
              return left.second < right.second;});
    bestTime = crtTpcMatchCandidates[0].first;
  }
  
  return bestTime;

}

// ============================== INTERNAL FUNCTIONS ===========================

// FOR TRACKING

// Function to make creating CRTHits easier
sbnd::crt::CRTHit FillCrtHit(std::vector<uint8_t> tfeb_id, std::map<uint8_t, 
                                         std::vector<std::pair<int,float>>> tpesmap, float peshit, double time, 
                                         double x, double ex, double y, double ey, double z, double ez, std::string tagger)
{

    sbnd::crt::CRTHit crtHit;

    crtHit.feb_id      = tfeb_id;
    crtHit.pesmap      = tpesmap;
    crtHit.peshit      = peshit;
    crtHit.ts0_s_corr  = 0;
    crtHit.ts0_ns      = time * 1e3;
    crtHit.ts0_ns_corr = 0;
    crtHit.ts1_ns      = time * 1e3;
    crtHit.ts0_s       = time * 1e-6;
    crtHit.x_pos       = x;
    crtHit.x_err       = ex;
    crtHit.y_pos       = y; 
    crtHit.y_err       = ey;
    crtHit.z_pos       = z;
    crtHit.z_err       = ez;
    crtHit.tagger      = tagger;

    return crtHit;

} // CRTTrackProducer::FillCrtHit()


// Function to make creating CRTTracks easier
sbnd::crt::CRTTrack FillCrtTrack(sbnd::crt::CRTHit hit1, sbnd::crt::CRTHit hit2, bool complete)
{

  sbnd::crt::CRTTrack newtr;

  newtr.ts0_s         = (hit1.ts0_s + hit2.ts0_s)/2.;
  newtr.ts0_s_err     = (uint32_t)((hit1.ts0_s - hit2.ts0_s)/2.);
  newtr.ts0_ns_h1     = hit1.ts0_ns;
  newtr.ts0_ns_err_h1 = hit1.ts0_ns_corr;
  newtr.ts0_ns_h2     = hit2.ts0_ns;
  newtr.ts0_ns_err_h2 = hit2.ts0_ns_corr;
  newtr.ts0_ns        = (uint32_t)((hit1.ts0_ns + hit2.ts0_ns)/2.);
  newtr.ts0_ns_err    = (uint16_t)(sqrt(hit1.ts0_ns_corr*hit1.ts0_ns_corr + hit2.ts0_ns_corr*hit2.ts0_ns_corr)/2.);
  newtr.ts1_ns        = (int32_t)(((double)(int)hit1.ts1_ns + (double)(int)hit2.ts1_ns)/2.);
  newtr.ts1_ns_err    = (uint16_t)(sqrt(hit1.ts0_ns_corr*hit1.ts0_ns_corr + hit2.ts0_ns_corr*hit2.ts0_ns_corr)/2.);
  newtr.peshit        = hit1.peshit+hit2.peshit;
  newtr.x1_pos        = hit1.x_pos;
  newtr.x1_err        = hit1.x_err;
  newtr.y1_pos        = hit1.y_pos;
  newtr.y1_err        = hit1.y_err;
  newtr.z1_pos        = hit1.z_pos;
  newtr.z1_err        = hit1.z_err;
  newtr.x2_pos        = hit2.x_pos;
  newtr.x2_err        = hit2.x_err;
  newtr.y2_pos        = hit2.y_pos;
  newtr.y2_err        = hit2.y_err;
  newtr.z2_pos        = hit2.z_pos;
  newtr.z2_err        = hit2.z_err;
  float deltax        = hit1.x_pos - hit2.x_pos;
  float deltay        = hit1.y_pos - hit2.y_pos;
  float deltaz        = hit1.z_pos - hit2.z_pos;
  newtr.length        = sqrt(deltax*deltax + deltay*deltay+deltaz*deltaz);
  newtr.thetaxy       = atan2(deltax,deltay);
  newtr.phizy         = atan2(deltaz,deltay);
  newtr.plane1        = hit1.plane;
  newtr.plane2        = hit2.plane;
  newtr.complete      = complete;

  return(newtr);

} // FillCrtTrack()

// Function to average hits within a certain distance of each other
std::vector<sbnd::crt::CRTHit> AverageHits(std::vector<sbnd::crt::CRTHit> hits, double fAverageHitDist)
{

  std::vector<sbnd::crt::CRTHit> returnHits;
  std::vector<sbnd::crt::CRTHit> aveHits;
  std::vector<sbnd::crt::CRTHit> spareHits;

  if (hits.size()>0){
    // loop over size of tx
    bool first = true;
    TVector3 middle(0., 0., 0.);
    for (size_t i = 0; i < hits.size(); i++){
      // Get the position of the hit
      TVector3 pos(hits[i].x_pos, hits[i].y_pos, hits[i].z_pos);
      // If first then set average = hit pos
      if(first){
        middle = pos;
        first = false;
      }
      // If distance from average < limit then add to average
      if((pos-middle).Mag() < fAverageHitDist){
        aveHits.push_back(hits[i]);
      }
      // Else add to another vector
      else{
        spareHits.push_back(hits[i]);
      }
    }

    sbnd::crt::CRTHit aveHit = DoAverage(aveHits);
    returnHits.push_back(aveHit);

    //Do this recursively
    std::vector<sbnd::crt::CRTHit> moreHits = AverageHits(spareHits, fAverageHitDist);
    returnHits.insert(returnHits.end(), moreHits.begin(), moreHits.end());
    return returnHits;
  }
  else {
    return returnHits;
  }

} // CRTTrackProducer::AverageHits()

  
// Take a list of hits and find average parameters
sbnd::crt::CRTHit DoAverage(std::vector<sbnd::crt::CRTHit> hits)
{

  // Initialize values
  std::string tagger = hits[0].tagger;
  double xpos = 0.; 
  double ypos = 0.;
  double zpos = 0.;
  double xmax = -99999; double xmin = 99999;
  double ymax = -99999; double ymin = 99999;
  double zmax = -99999; double zmin = 99999;
  double ts1_ns = 0.;
  int nhits = 0;

  // Loop over hits
  for( auto& hit : hits ){
    // Get the mean x,y,z and times
    xpos += hit.x_pos;
    ypos += hit.y_pos;
    zpos += hit.z_pos;
    ts1_ns += (double)(int)hit.ts1_ns;
    // For the errors get the maximum limits
    if(hit.x_pos + hit.x_err > xmax) xmax = hit.x_pos + hit.x_err;
    if(hit.x_pos - hit.x_err < xmin) xmin = hit.x_pos - hit.x_err;
    if(hit.y_pos + hit.y_err > ymax) ymax = hit.y_pos + hit.y_err;
    if(hit.y_pos - hit.y_err < ymin) ymin = hit.y_pos - hit.y_err;
    if(hit.z_pos + hit.z_err > zmax) zmax = hit.z_pos + hit.z_err;
    if(hit.z_pos - hit.z_err < zmin) zmin = hit.z_pos - hit.z_err;
    // Add all the unique IDs in the vector
    nhits++;
  }

  // Create a hit
  sbnd::crt::CRTHit crtHit = FillCrtHit(hits[0].feb_id, hits[0].pesmap, hits[0].peshit, 
                                  (ts1_ns/nhits)*1e-3, xpos/nhits, (xmax-xmin)/2,
                                  ypos/nhits, (ymax-ymin)/2., zpos/nhits, (zmax-zmin)/2., tagger);

  return crtHit;

} // DoAverage()

// Function to calculate the crossing point of a track and tagger
TVector3 CrossPoint(sbnd::crt::CRTHit hit, TVector3 start, TVector3 diff)
{
  TVector3 cross;
  // Use the error to get the fixed coordinate of a tagger
  if(hit.x_err > 0.39 && hit.x_err < 0.41){
    double xc = hit.x_pos;
    TVector3 crossp(xc, 
                    ((xc - start.X()) / (diff.X()) * diff.Y()) + start.Y(), 
                    ((xc - start.X()) / (diff.X()) * diff.Z()) + start.Z());
    cross = crossp;
  }

  else if(hit.y_err > 0.39 && hit.y_err < 0.41){
    double yc = hit.y_pos;
    TVector3 crossp(((yc - start.Y()) / (diff.Y()) * diff.X()) + start.X(), 
                    yc, 
                    ((yc - start.Y()) / (diff.Y()) * diff.Z()) + start.Z());
    cross = crossp;
  }

  else if(hit.z_err > 0.39 && hit.z_err < 0.41){
    double zc = hit.z_pos;
    TVector3 crossp(((zc - start.Z()) / (diff.Z()) * diff.X()) + start.X(), 
                    ((zc - start.Z()) / (diff.Z()) * diff.Y()) + start.Y(), 
                    zc);
    cross = crossp;
  }

  return cross;

} // CrossPoint()

// FOR T0 MATCHING

// Utility function that determines the possible t0 range of a track
std::pair<double, double> TrackT0Range(double startX, double endX, int tpc){

  geo::GeometryCore const* fGeometryService = lar::providerFrom<geo::Geometry>();
  detinfo::DetectorProperties const* fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
  std::pair<double, double> result;
  double Vd = fDetectorProperties->DriftVelocity();

  // Whole track must be within tpc
  // Find which TPC the track hits are detected in
  if(tpc == 0){
    // Lowest |X| is furthest from APA
    double lowX = std::max(startX, endX);
    // xmin is shift from furthest to 0 (the CPA)
    double xmax = 0 - lowX;
    // Highest |X| is closest to APA
    double highX = std::min(startX, endX);
    // xmax is shift from closest to APA
    double xmin = -(2.0*fGeometryService->DetHalfWidth()+3.) - highX;
    double t0max = -(xmin/Vd);
    double t0min = -(xmax/Vd);
    result = std::make_pair(t0min, t0max);
  }

  else{
    // Lowest |X| is furthest from APA
    double lowX = std::min(startX, endX);
    // xmin is shift from furthest to 0 (the CPA)
    double xmin = 0 - lowX;
    // Highest |X| is closest to APA
    double highX = std::max(startX, endX);
    // xmax is shift from closest to APA
    double xmax = (2.0*fGeometryService->DetHalfWidth()+3.) - highX;
    double t0min = xmin/Vd;
    double t0max = xmax/Vd;
    result = std::make_pair(t0min, t0max);
  }

  return result;

} // TrackT0Range()

double DistOfClosestApproach(TVector3 trackPos, TVector3 trackDir, sbnd::crt::CRTHit crtHit, int tpc, double t0){

  detinfo::DetectorProperties const* fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
  double minDist = 99999;

  // Convert the t0 into an x shift
  double shift = t0 * fDetectorProperties->DriftVelocity();
  // Apply the shift depending on which TPC the track is in
  if (tpc == 1) trackPos[0] += shift;
  if (tpc == 0) trackPos[0] -= shift;

  TVector3 endPos = trackPos + trackDir;
  double denominator = (endPos - trackPos).Mag();
  // 1D hits should only have a lot of variance in one direction
  if(crtHit.x_err > 50.){
    // Loop over size of hit to find the min dist
    for(int i = 0; i < 20.; i++){
      double xpos = crtHit.x_pos + ((i+1.)/10. - 1.)*crtHit.x_err;
      TVector3 crtPoint(xpos, crtHit.y_pos, crtHit.z_pos);
      double numerator = ((crtPoint - trackPos).Cross(crtPoint-endPos)).Mag();
      double dca = numerator/denominator;
      if(dca < minDist) minDist = dca;
    }
  }
  else if(crtHit.y_err > 50.){
    // Loop over size of hit to find the min dist
    for(int i = 0; i < 20.; i++){
      double ypos = crtHit.y_pos + ((i+1.)/10. - 1.)*crtHit.y_err;
      TVector3 crtPoint(crtHit.x_pos, ypos, crtHit.z_pos);
      double numerator = ((crtPoint - trackPos).Cross(crtPoint-endPos)).Mag();
      double dca = numerator/denominator;
      if(dca < minDist) minDist = dca;
    }
  }
  else if(crtHit.y_err > 50.){
    // Loop over size of hit to find the min dist
    for(int i = 0; i < 20.; i++){
      double zpos = crtHit.z_pos + ((i+1.)/10. - 1.)*crtHit.z_err;
      TVector3 crtPoint(crtHit.x_pos, crtHit.y_pos, zpos);
      double numerator = ((crtPoint - trackPos).Cross(crtPoint-endPos)).Mag();
      double dca = numerator/denominator;
      if(dca < minDist) minDist = dca;
    }
  }
  else{
    TVector3 crtPoint(crtHit.x_pos, crtHit.y_pos, crtHit.z_pos);
    double numerator = ((crtPoint - trackPos).Cross(crtPoint-endPos)).Mag();
    double dca = numerator/denominator;
    if(dca < minDist) minDist = dca;
  }

  return minDist;

} // DistToOfClosestApproach()

// Function to transform a CRTTrack into an expected reconstructed track
std::vector<RecoCRTTrack> CrtToRecoTrack(sbnd::crt::CRTTrack track, int id){

  detinfo::DetectorProperties const* fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
  std::vector<RecoCRTTrack> recoCrtTracks;
  // Get the time of the track
  double crtTime = ((double)(int)track.ts1_ns) * 1e-3;
  // Convert time into a x shift
  double xShift = crtTime * fDetectorProperties->DriftVelocity();

  // Shift track, remembering to take into account the tpc, if the track crosses the cpa and 
  //the size of the readout window
  TVector3 crtStart(track.x1_pos, track.y1_pos, track.z1_pos);
  TVector3 crtEnd(track.x2_pos, track.y2_pos, track.z2_pos);
  if(crtStart.Y() < crtEnd.Y()) std::swap(crtStart, crtEnd);
  if(!track.complete){
    // Find point where track crosses bottom plane
    TVector3 diff = (crtEnd - crtStart).Unit();
    TVector3 newEnd = crtStart + 5000*diff; 
    crtEnd = newEnd;
  }

  TVector3 cpaCrossStart, cpaCrossEnd;
  // Calculate the expected reconstructed length, angle and start and end position
  std::vector<int> crtTpc;
  if(crtStart.X() < 0. && crtEnd.X() < 0.){

    // Track in TPC 0
    std::vector<RecoCRTTrack> tempTracks = CreateRecoCRTTrack(crtStart, crtEnd, xShift, 0, 
                                                              id, crtTime, track.complete); 
    recoCrtTracks.insert(recoCrtTracks.end(), tempTracks.begin(), tempTracks.end());
  }
  else if(crtStart.X() > 0. && crtEnd.X() > 0.){

    // Track in TPC 1
    std::vector<RecoCRTTrack> tempTracks = CreateRecoCRTTrack(crtStart, crtEnd, xShift, 1, 
                                                              id, crtTime, track.complete); 
    recoCrtTracks.insert(recoCrtTracks.end(), tempTracks.begin(), tempTracks.end());
  }
  else {
    // Track in both TPCs and will be split
    TVector3 direction = crtStart - crtEnd;
    double step = (0. - crtStart.X())/direction.X(); 
    TVector3 cpaCross(0., crtStart.Y() + step*direction.Y(), crtStart.Z() + step*direction.Z());
    cpaCrossStart = cpaCross;
    cpaCrossEnd = cpaCross;

    if(crtStart.X() < 0.){ 
      std::vector<RecoCRTTrack> tempTracks0 = CreateRecoCRTTrack(crtStart, cpaCrossStart, xShift, 0, 
                                                                 id, crtTime, track.complete); 
      recoCrtTracks.insert(recoCrtTracks.end(), tempTracks0.begin(), tempTracks0.end());

      std::vector<RecoCRTTrack> tempTracks1 = CreateRecoCRTTrack(crtEnd, cpaCrossEnd, xShift, 1, 
                                                                 id, crtTime, track.complete); 
      recoCrtTracks.insert(recoCrtTracks.end(), tempTracks1.begin(), tempTracks1.end());
    }
    else {
      std::vector<RecoCRTTrack> tempTracks0 = CreateRecoCRTTrack(crtEnd, cpaCrossEnd, xShift, 0, 
                                                                 id, crtTime, track.complete); 
      recoCrtTracks.insert(recoCrtTracks.end(), tempTracks0.begin(), tempTracks0.end());

      std::vector<RecoCRTTrack> tempTracks1 = CreateRecoCRTTrack(crtStart, cpaCrossStart, xShift, 1, 
                                                                 id, crtTime, track.complete); 
      recoCrtTracks.insert(recoCrtTracks.end(), tempTracks1.begin(), tempTracks1.end());
    }

  }

  return recoCrtTracks;

} // CRTToRecoTrack()


// Function to shift CRTTrack in X and work out how much is reconstructed
std::vector<RecoCRTTrack> CreateRecoCRTTrack(TVector3 start, TVector3 end, double shift, 
                                             int tpc, int id, double time, bool complete){

  geo::GeometryCore const* fGeometryService = lar::providerFrom<geo::Geometry>();
  detinfo::DetectorClocks const* fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>(); 
  detinfo::DetectorProperties const* fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
  std::vector<RecoCRTTrack> recoCrtTracks;

  // Get the true entry and exit points in the TPC
  double xmin = -2.0 * fGeometryService->DetHalfWidth();
  double xmax = 2.0 * fGeometryService->DetHalfWidth();
  double ymin = -fGeometryService->DetHalfHeight();
  double ymax = fGeometryService->DetHalfHeight();
  double zmin = 0.;
  double zmax = fGeometryService->DetLength();

  // Get track info
  TVector3 diff = end - start;
  TVector3 startTPC, endTPC;
  bool first = true;

  // Loop over trajectory points
  int npts = 1000;
  for (int traj_i = 0; traj_i <= npts; traj_i++){
    TVector3 trajPoint = start + ((traj_i)/((double)npts))*diff;
    // Check if point is within reconstructable volume
    if (trajPoint[0] >= xmin && trajPoint[0] <= xmax && trajPoint[1] >= ymin && trajPoint[1] <= ymax 
        && trajPoint[2] >= zmin && trajPoint[2] <= zmax && first){
      first = false;
      startTPC = trajPoint;
    }
    if (trajPoint[0] >= xmin && trajPoint[0] <= xmax && trajPoint[1] >= ymin && trajPoint[1] <= ymax 
        && trajPoint[2] >= zmin && trajPoint[2] <= zmax){
      endTPC = trajPoint;
    }
  }

  // Don't shift if not inside TPC
  if(startTPC == endTPC){
    return recoCrtTracks;
  }

  // Shift in x depending on TPC
  if(tpc == 1){
    // Track in TPC 1
    startTPC[0] -= shift;
    endTPC[0] -= shift;
  }
  else if(tpc == 0){
    // Track in TPC 0
    startTPC[0] += shift;
    endTPC[0] += shift;
  }
  
  double readoutWindow  = fDetectorClocks->TPCTick2Time((double)fDetectorProperties->ReadOutWindowSize());
  double driftTimeTicks = (2.*fGeometryService->DetHalfWidth())/fDetectorProperties->DriftVelocity();
  double deltaX = (readoutWindow - driftTimeTicks) * fDetectorProperties->DriftVelocity();

  if(tpc == 0) xmax = deltaX;
  if(tpc == 1) xmin = -deltaX;

  // Get track info
  TVector3 diffTPC = endTPC - startTPC;
  TVector3 startCut, endCut;
  first = true;
  // Loop over trajectory points
  for (int traj_i = 0; traj_i <= npts; traj_i++){
    TVector3 trajPoint = startTPC + ((traj_i)/((double)npts))*diffTPC;
    // Check if point is within reconstructable volume
    if (trajPoint[0] >= xmin && trajPoint[0] <= xmax && trajPoint[1] >= ymin && trajPoint[1] <= ymax 
        && trajPoint[2] >= zmin && trajPoint[2] <= zmax && first){
      first = false;
      startCut = trajPoint;
    }
    if (trajPoint[0] >= xmin && trajPoint[0] <= xmax && trajPoint[1] >= ymin && trajPoint[1] <= ymax 
        && trajPoint[2] >= zmin && trajPoint[2] <= zmax){
      endCut = trajPoint;
    }
  }

  if(!(startCut.X() == endCut.X())){
    RecoCRTTrack recoCrtTrack = {id, tpc, startCut, endCut, time, complete};
    recoCrtTracks.push_back(recoCrtTrack);
  }

  return recoCrtTracks;

} // CreateRecoCRTTrack()


// Function to calculate if a CRTTrack crosses the TPC volume
bool CrossesTPC(sbnd::crt::CRTTrack track){

  geo::GeometryCore const* fGeometryService = lar::providerFrom<geo::Geometry>();
  // Check if particle enters the TPC
  bool enters = false;
  double xmin = -2.0 * fGeometryService->DetHalfWidth();
  double xmax = 2.0 * fGeometryService->DetHalfWidth();
  double ymin = -fGeometryService->DetHalfHeight();
  double ymax = fGeometryService->DetHalfHeight();
  double zmin = 0.;
  double zmax = fGeometryService->DetLength();
  if(track.complete){
    // Get track info
    TVector3 start(track.x1_pos, track.y1_pos, track.z1_pos);
    TVector3 end(track.x2_pos, track.y2_pos, track.z2_pos);
    TVector3 diff = end - start;
    // Loop over trajectory points
    int npts = 100;
    for (int traj_i = 0; traj_i < npts; traj_i++){
      TVector3 trajPoint = start + ((traj_i+1)/100.)*diff;
      // Check if point is within reconstructable volume
      if (trajPoint[0] >= xmin-5 && trajPoint[0] <= xmax+5 && trajPoint[1] >= ymin-5 && trajPoint[1] <= ymax+5 && trajPoint[2] >= zmin-5 && trajPoint[2] <= zmax+5){
        enters = true;
      }
    }
  }
  // If track just between top two planes
  else{
    //
    TVector3 start(track.x1_pos, track.y1_pos, track.z1_pos);
    TVector3 end(track.x2_pos, track.y2_pos, track.z2_pos);
    if(start.Y() < end.Y()) std::swap(start, end);
    TVector3 diff = (end - start).Unit();
    int npts = 100;
    for (int traj_i = 0; traj_i < npts; traj_i++){
      TVector3 trajPoint = start + (100.*(traj_i+1))*diff;
      // Check if point is within reconstructable volume
      if (trajPoint[0] >= xmin-5 && trajPoint[0] <= xmax+5 && trajPoint[1] >= ymin-5 && trajPoint[1] <= ymax+5 && trajPoint[2] >= zmin-5 && trajPoint[2] <= zmax+5){
        enters = true;
      }
    }
  }

  return enters;

} // CrossesTPC()


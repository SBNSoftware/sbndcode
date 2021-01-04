#include "CRTTrackRecoAlg.h"

namespace sbnd{

CRTTrackRecoAlg::CRTTrackRecoAlg(const Config& config)
  : hitAlg() {

  this->reconfigure(config);
  
}

CRTTrackRecoAlg::CRTTrackRecoAlg(double aveHitDist, double distLim)
  : hitAlg() {

  fAverageHitDistance = aveHitDist;
  fDistanceLimit = distLim;

}


CRTTrackRecoAlg::~CRTTrackRecoAlg(){

}


void CRTTrackRecoAlg::reconfigure(const Config& config){

  fTimeLimit          = config.TimeLimit();
  fAverageHitDistance = config.AverageHitDistance();
  fDistanceLimit      = config.DistanceLimit();

  return;

}


std::vector<std::vector<art::Ptr<sbn::crt::CRTHit>>> CRTTrackRecoAlg::CreateCRTTzeros(std::vector<art::Ptr<sbn::crt::CRTHit>> hits)
{

  std::vector<std::vector<art::Ptr<sbn::crt::CRTHit>>> crtTzeroVect;
  std::vector<int> iflag(hits.size(), 0);

  // Sort CRTHits by time
  std::sort(hits.begin(), hits.end(), [](auto& left, auto& right)->bool{
              return left->ts1_ns < right->ts1_ns;});

  // Loop over crt hits
  for(size_t i = 0; i<hits.size(); i++){
    if(iflag[i] == 0){
      std::vector<art::Ptr<sbn::crt::CRTHit>> crtTzero;
      double time_ns_A = hits[i]->ts1_ns;
      iflag[i]=1;
      crtTzero.push_back(hits[i]);

      // Sort into a Tzero collection
      // Loop over all the other CRT hits
      for(size_t j = i+1; j<hits.size(); j++){
        if(iflag[j] == 0){
          // If ts1_ns - ts1_ns < diff then put them in a vector
          double time_ns_B = hits[j]->ts1_ns; 
          double diff = std::abs(time_ns_B - time_ns_A) * 1e-3; // [us]
          if(diff < fTimeLimit){
            iflag[j] = 1;
            crtTzero.push_back(hits[j]);
          }
        }
      }

      crtTzeroVect.push_back(crtTzero);
    }
  }
  return crtTzeroVect;
}


// Function to make creating CRTTracks easier
sbn::crt::CRTTrack CRTTrackRecoAlg::FillCrtTrack(sbn::crt::CRTHit hit1, sbn::crt::CRTHit hit2, bool complete)
{

  sbn::crt::CRTTrack newtr;

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

} // CRTTrackRecoAlg::FillCrtTrack()


sbn::crt::CRTTrack CRTTrackRecoAlg::FillCrtTrack(sbn::crt::CRTHit hit1, sbn::crt::CRTHit hit2, size_t nhits)
{

  bool complete = true;
  if(nhits == 2){
    // Track is incomplete if just between 2 top planes
    if((hit1.tagger == "volTaggerTopHigh_0" && hit2.tagger == "volTaggerTopLow_0")
       || (hit2.tagger == "volTaggerTopHigh_0" && hit1.tagger == "volTaggerTopLow_0")) complete = false;
    return FillCrtTrack(hit1, hit2, complete);
  }
  else{
    // Project track on to the limits of the CRT volume TODO errors
    std::vector<double> crtLimits = fCrtGeo.CRTLimits();
    TVector3 min (crtLimits[0], crtLimits[1], crtLimits[2]);
    TVector3 max (crtLimits[3], crtLimits[4], crtLimits[5]);
    TVector3 start (hit1.x_pos, hit1.y_pos, hit1.z_pos);
    TVector3 end (hit2.x_pos, hit2.y_pos, hit2.z_pos);

    std::pair<TVector3, TVector3> intersection = CRTCommonUtils::CubeIntersection(min, max, start, end);
    if(intersection.first.X() == -99999) return FillCrtTrack(hit1, hit2, complete);

    hit1.x_pos = intersection.first.X();
    hit1.y_pos = intersection.first.Y();
    hit1.z_pos = intersection.first.Z();
    hit1.x_err = 0.;
    hit1.y_err = 0.;
    hit1.z_err = 0.;
    hit2.x_pos = intersection.second.X();
    hit2.y_pos = intersection.second.Y();
    hit2.z_pos = intersection.second.Z();
    hit1.x_err = 0.;
    hit1.y_err = 0.;
    hit1.z_err = 0.;
  }

  return FillCrtTrack(hit1, hit2, complete);
}


// Function to average hits within a certain distance of each other
std::vector<std::pair<sbn::crt::CRTHit, std::vector<int>>> CRTTrackRecoAlg::AverageHits(std::vector<art::Ptr<sbn::crt::CRTHit>> hits, std::map<art::Ptr<sbn::crt::CRTHit>, int> hitIds)
{

  std::vector<std::pair<sbn::crt::CRTHit, std::vector<int>>> returnHits;
  std::vector<art::Ptr<sbn::crt::CRTHit>> aveHits;
  std::vector<art::Ptr<sbn::crt::CRTHit>> spareHits;

  if (hits.size()>0){
    // loop over size of tx
    bool first = true;
    TVector3 middle(0., 0., 0.);
    for (size_t i = 0; i < hits.size(); i++){
      // Get the position of the hit
      TVector3 pos(hits[i]->x_pos, hits[i]->y_pos, hits[i]->z_pos);
      // If first then set average = hit pos
      if(first){
        middle = pos;
        first = false;
      }
      // If distance from average < limit then add to average
      if((pos-middle).Mag() < fAverageHitDistance){
        aveHits.push_back(hits[i]);
      }
      // Else add to another vector
      else{
        spareHits.push_back(hits[i]);
      }
    }

    sbn::crt::CRTHit aveHit = DoAverage(aveHits);
    std::vector<int> ids;
    for(size_t i = 0; i < aveHits.size(); i++){
      ids.push_back(hitIds[aveHits[i]]);
    }
    returnHits.push_back(std::make_pair(aveHit, ids));

    //Do this recursively
    std::vector<std::pair<sbn::crt::CRTHit, std::vector<int>>> moreHits = AverageHits(spareHits, hitIds);
    returnHits.insert(returnHits.end(), moreHits.begin(), moreHits.end());
    return returnHits;
  }
  else {
    return returnHits;
  }

} // CRTTrackRecoAlg::AverageHits()


std::vector<sbn::crt::CRTHit> CRTTrackRecoAlg::AverageHits(std::vector<art::Ptr<sbn::crt::CRTHit>> hits)
{

  std::map<art::Ptr<sbn::crt::CRTHit>, int> dummy;
  for(size_t i = 0; i < hits.size(); i++){
    dummy[hits[i]] = 0;
  }

  std::vector<std::pair<sbn::crt::CRTHit, std::vector<int>>> output = AverageHits(hits, dummy);

  std::vector<sbn::crt::CRTHit> returnHits;
  for(auto const& out : output){
    returnHits.push_back(out.first);
  }

  return returnHits;

} // CRTTrackRecoAlg::AverageHits()


  
// Take a list of hits and find average parameters
sbn::crt::CRTHit CRTTrackRecoAlg::DoAverage(std::vector<art::Ptr<sbn::crt::CRTHit>> hits)
{

  // Initialize values
  std::string tagger = hits[0]->tagger;
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
    xpos += hit->x_pos;
    ypos += hit->y_pos;
    zpos += hit->z_pos;
    ts1_ns += (double)(int)hit->ts1_ns;
    // For the errors get the maximum limits
    if(hit->x_pos + hit->x_err > xmax) xmax = hit->x_pos + hit->x_err;
    if(hit->x_pos - hit->x_err < xmin) xmin = hit->x_pos - hit->x_err;
    if(hit->y_pos + hit->y_err > ymax) ymax = hit->y_pos + hit->y_err;
    if(hit->y_pos - hit->y_err < ymin) ymin = hit->y_pos - hit->y_err;
    if(hit->z_pos + hit->z_err > zmax) zmax = hit->z_pos + hit->z_err;
    if(hit->z_pos - hit->z_err < zmin) zmin = hit->z_pos - hit->z_err;
    // Add all the unique IDs in the vector
    nhits++;
  }

  // Create a hit
  sbn::crt::CRTHit crtHit = hitAlg.FillCrtHit(hits[0]->feb_id, hits[0]->pesmap, hits[0]->peshit, 
                                  (ts1_ns/nhits)*1e-3, 0, xpos/nhits, (xmax-xmin)/2,
                                  ypos/nhits, (ymax-ymin)/2., zpos/nhits, (zmax-zmin)/2., tagger);

  return crtHit;

} // CRTTrackRecoAlg::DoAverage()


// Function to create tracks from tzero hit collections
std::vector<std::pair<sbn::crt::CRTTrack, std::vector<int>>> CRTTrackRecoAlg::CreateTracks(std::vector<std::pair<sbn::crt::CRTHit, std::vector<int>>> hits)
{

  std::vector<std::pair<sbn::crt::CRTTrack, std::vector<int>>> returnTracks;

  std::vector<std::vector<size_t>> trackCandidates;
  // Loop over all hits
  for(size_t i = 0; i < hits.size(); i++){

    // Loop over all unique pairs
    for(size_t j = i+1; j < hits.size(); j++){
      if(hits[i].first.tagger == hits[j].first.tagger) continue;

      // Draw a track between the two hits
      TVector3 start (hits[i].first.x_pos, hits[i].first.y_pos, hits[i].first.z_pos);
      TVector3 end (hits[j].first.x_pos, hits[j].first.y_pos, hits[j].first.z_pos);
      TVector3 diff = end - start;

      std::vector<size_t> candidate {i, j};

      // Loop over all other hits on different taggers and calculate DCA with variations
      for(size_t k = 0; k < hits.size(); k++){
        if(k == i || k == j || hits[k].first.tagger == hits[i].first.tagger 
           || hits[k].first.tagger == hits[j].first.tagger) continue;

        //  If hit within certain distance then add it to the track candidate
        if(CRTCommonUtils::DistToCrtHit(hits[k].first, start, end) < fDistanceLimit){
          candidate.push_back(k);
        }
      }
      trackCandidates.push_back(candidate);
    }
  }

  // Sort track candidates by number of hits
  std::sort(trackCandidates.begin(), trackCandidates.end(), [](auto& left, auto& right){
            return left.size() > right.size();});

  // Loop over track candidates
  std::vector<size_t> usedHits;
  for(auto const& candidate : trackCandidates){
    // Check if hits have been used
    bool used = false;
    for(size_t i = 0; i < candidate.size(); i++){
      //Check if any of the hits have been used
      if(std::find(usedHits.begin(), usedHits.end(), candidate[i]) != usedHits.end()) used = true;
    }
    if(used) continue;

    // Create track 
    if(candidate.size() < 2) continue;
    sbn::crt::CRTHit ihit = hits[candidate[0]].first;
    sbn::crt::CRTHit jhit = hits[candidate[1]].first;
    sbn::crt::CRTTrack crtTrack = FillCrtTrack(ihit, jhit, candidate.size());

    std::vector<int> ids;
    //TODO: Add charge matching for ambiguous cases
    // If nhits > 2 then record used hits
    for(size_t i = 0; i < candidate.size(); i++){
      ids.insert(ids.end(), hits[candidate[i]].second.begin(), hits[candidate[i]].second.end());
      if(candidate.size()>2) usedHits.push_back(candidate[i]);
    }

    returnTracks.push_back(std::make_pair(crtTrack, ids));

  }

  return returnTracks;

} // CRTTrackRecoAlg::CreateTracks()


std::vector<sbn::crt::CRTTrack> CRTTrackRecoAlg::CreateTracks(std::vector<sbn::crt::CRTHit> hits)
{

  std::vector<std::pair<sbn::crt::CRTHit, std::vector<int>>> input;
  for(auto const& hit : hits){
    std::vector<int> dummy;
    input.push_back(std::make_pair(hit, dummy));
  }

  std::vector<std::pair<sbn::crt::CRTTrack, std::vector<int>>> output = CreateTracks(input);

  std::vector<sbn::crt::CRTTrack> tracks;
  for(auto const& out : output){
    tracks.push_back(out.first);
  }

  return tracks;
  

} // CRTTrackRecoAlg::CreateTracks()



}

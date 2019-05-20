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

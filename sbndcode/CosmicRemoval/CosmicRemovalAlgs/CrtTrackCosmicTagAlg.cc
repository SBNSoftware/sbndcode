#include "CrtTrackCosmicTagAlg.h"

namespace sbnd{

CrtTrackCosmicTagAlg::CrtTrackCosmicTagAlg(const Config& config){

  this->reconfigure(config);

}


CrtTrackCosmicTagAlg::CrtTrackCosmicTagAlg(){

}


CrtTrackCosmicTagAlg::~CrtTrackCosmicTagAlg(){

}


void CrtTrackCosmicTagAlg::reconfigure(const Config& config){

  trackMatchAlg = config.TrackMatchAlg();
  fBeamTimeMin = config.BeamTimeLimits().BeamTimeMin();
  fBeamTimeMax = config.BeamTimeLimits().BeamTimeMax();

  return;
}


bool CrtTrackCosmicTagAlg::CrtTrackCosmicTag(recob::Track track, std::vector<crt::CRTTrack> crtTracks, int tpc){

  int crtID = trackMatchAlg.GetMatchedCRTTrackId(track, crtTracks, tpc);

  // If matching failed
  if(crtID == -99999) return false;

  // If track matched to a through going CRT track then it is a cosmic
  if(crtTracks.at(crtID).complete) return true;

  // If it matches a track through just the top planes make sure it is outside of the beam time
  double crtTime = ((double)(int)crtTracks.at(crtID).ts1_ns) * 1e-3; // [us]
  if(crtTime < fBeamTimeMin || crtTime > fBeamTimeMax) return true;

  return false;

}
 
}

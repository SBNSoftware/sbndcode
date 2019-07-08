#include "CrtHitCosmicIdAlg.h"

namespace sbnd{

CrtHitCosmicIdAlg::CrtHitCosmicIdAlg(const Config& config){

  this->reconfigure(config);

} //CrtHitCosmicIdAlg()


CrtHitCosmicIdAlg::CrtHitCosmicIdAlg(){

} //CrtHitCosmicIdAlg()


CrtHitCosmicIdAlg::~CrtHitCosmicIdAlg(){

} //~CrtHitCosmicIdAlg()


void CrtHitCosmicIdAlg::reconfigure(const Config& config){

  t0Alg = config.T0Alg();
  fBeamTimeMin = config.BeamTimeLimits().BeamTimeMin();
  fBeamTimeMax = config.BeamTimeLimits().BeamTimeMax();

  return;
} //reconfigure()


bool CrtHitCosmicIdAlg::CrtHitCosmicId(recob::Track track, std::vector<crt::CRTHit> crtHits, const art::Event& event){

  double crtHitTime = t0Alg.T0FromCRTHits(track, crtHits, event);

  if(crtHitTime != -99999 && (crtHitTime < fBeamTimeMin || crtHitTime > fBeamTimeMax)) return true;

  return false;

} //CrtHitCosmicId()
 
}

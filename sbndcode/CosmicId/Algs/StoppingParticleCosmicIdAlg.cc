#include "StoppingParticleCosmicIdAlg.h"

namespace sbnd{

StoppingParticleCosmicIdAlg::StoppingParticleCosmicIdAlg(const Config& config){

  this->reconfigure(config);

}


StoppingParticleCosmicIdAlg::StoppingParticleCosmicIdAlg(){

}


StoppingParticleCosmicIdAlg::~StoppingParticleCosmicIdAlg(){

}


void StoppingParticleCosmicIdAlg::reconfigure(const Config& config){

  fMinX = config.ContainmentCuts().MinX(); 
  fMinY = config.ContainmentCuts().MinY(); 
  fMinZ = config.ContainmentCuts().MinZ(); 
  fMaxX = config.ContainmentCuts().MaxX(); 
  fMaxY = config.ContainmentCuts().MaxY(); 
  fMaxZ = config.ContainmentCuts().MaxZ(); 
  fResRangeMin = config.ResRangeMin();
  fResRangeMax = config.ResRangeMax();
  fDEdxMax = config.DEdxMax();
  fStoppingChi2Limit = config.StoppingChi2Limit();

  return;
}

double StoppingParticleCosmicIdAlg::StoppingChiSq(geo::Point_t end, std::vector<art::Ptr<anab::Calorimetry>> calos){
  
  //Loop over residual range and dedx
  if(calos.size()==0) return false;
  size_t nhits = 0;
  art::Ptr<anab::Calorimetry> calo = calos[0];
  for( size_t i = calos.size(); i > 0; i--){
    if(calos[i-1]->dEdx().size() > nhits*1.5){
      nhits = calos[i-1]->dEdx().size();
      calo = calos[i-1];
    }
  }

  double distStart = (calo->XYZ()[0] - end).Mag2();
  double distEnd = (calo->XYZ()[nhits-1] - end).Mag2();

  double maxDedx = 0;
  double resrgStart = 0;
  std::vector<double> v_resrg;
  std::vector<double> v_dedx;

  for(size_t i = 0; i < nhits; i++){
    double dedx = calo->dEdx()[i];
    double resrg = calo->ResidualRange()[i];

    if(distStart < distEnd && calo->ResidualRange()[0] > calo->ResidualRange()[nhits-1]) resrg = calo->ResidualRange()[0] - calo->ResidualRange()[i];
    if(distStart > distEnd && calo->ResidualRange()[0] < calo->ResidualRange()[nhits-1]) resrg = calo->ResidualRange()[nhits-1] - calo->ResidualRange()[i];

    if(resrg < fResRangeMin && dedx > maxDedx && dedx < fDEdxMax){
      maxDedx = dedx;
      resrgStart = resrg;
    }

  }

  for(size_t i = 0; i < nhits; i++){
    double dedx = calo->dEdx()[i];
    double resrg = calo->ResidualRange()[i];

    if(distStart < distEnd && calo->ResidualRange()[0] > calo->ResidualRange()[nhits-1]) resrg = calo->ResidualRange()[0] - calo->ResidualRange()[i];
    if(distStart > distEnd && calo->ResidualRange()[0] < calo->ResidualRange()[nhits-1]) resrg = calo->ResidualRange()[nhits-1] - calo->ResidualRange()[i];

    if(resrg > resrgStart && resrg < resrgStart + fResRangeMax && dedx < fDEdxMax){
      v_resrg.push_back(resrg);
      v_dedx.push_back(dedx);
    }
  }

  if(v_dedx.size() < 10) return false;

  TGraph *gdedx = new TGraph(v_dedx.size(), &v_resrg[0], &v_dedx[0]);
  try{ gdedx->Fit("pol0", "Q"); } catch(...){ return false; }
  TF1* polfit = gdedx->GetFunction("pol0");
  double polchi2 = polfit->GetChisquare();

  try{ gdedx->Fit("expo", "Q"); } catch(...){ return false; }
  TF1* expfit = gdedx->GetFunction("expo");
  double expchi2 = expfit->GetChisquare();

  return polchi2/expchi2;

}

bool StoppingParticleCosmicIdAlg::StoppingEnd(geo::Point_t end, std::vector<art::Ptr<anab::Calorimetry>> calos){
  
  double chiSqRatio = StoppingChiSq(end, calos);

  if(chiSqRatio > fStoppingChi2Limit) return true;

  return false;
}

bool StoppingParticleCosmicIdAlg::StoppingParticleCosmicId(recob::Track track, std::vector<art::Ptr<anab::Calorimetry>> calos){

  bool startInFiducial = CosmicIdUtils::InFiducial(track.Vertex(), fMinX, fMinY, fMinZ, fMaxX, fMaxY, fMaxZ);
  bool endInFiducial = CosmicIdUtils::InFiducial(track.End(), fMinX, fMinY, fMinZ, fMaxX, fMaxY, fMaxZ);

  bool startStops = StoppingEnd(track.Vertex(), calos);
  bool endStops = StoppingEnd(track.End(), calos);

  if((startStops && !endInFiducial && track.End().Y() > 0) 
      || (endStops && !startInFiducial && track.Vertex().Y() > 0)){
    return true;
  }

  return false;

}

bool StoppingParticleCosmicIdAlg::StoppingParticleCosmicId(recob::Track track, recob::Track track2, std::vector<art::Ptr<anab::Calorimetry>> calos, std::vector<art::Ptr<anab::Calorimetry>> calos2){

  bool startInFiducial = CosmicIdUtils::InFiducial(track.End(), fMinX, fMinY, fMinZ, fMaxX, fMaxY, fMaxZ);
  bool endInFiducial = CosmicIdUtils::InFiducial(track.End(), fMinX, fMinY, fMinZ, fMaxX, fMaxY, fMaxZ);

  bool startStops = StoppingEnd(track.End(), calos);
  bool endStops = StoppingEnd(track2.End(), calos2);

  if((startStops && !endInFiducial && track2.End().Y() > 0) 
      || (endStops && !startInFiducial && track.End().Y() > 0)){
    return true;
  }

  return false;

}


}

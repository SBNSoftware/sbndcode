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

// Calculate the chi2 ratio of pol0 and exp fit to dE/dx vs residual range
double StoppingParticleCosmicIdAlg::StoppingChiSq(geo::Point_t end, std::vector<art::Ptr<anab::Calorimetry>> calos){

  // If calorimetry object is null then return 0
  if(calos.size()==0) return -99999;

  // Loop over planes (Y->V->U) and choose the next plane's calorimetry if there are 1.5x more points (collection plane more reliable)
  size_t nhits = 0;
  art::Ptr<anab::Calorimetry> calo = calos[0];
  for( size_t i = calos.size(); i > 0; i--){
    if(calos[i-1]->dEdx().size() > nhits*1.5){
      nhits = calos[i-1]->dEdx().size();
      calo = calos[i-1];
    }
  }

  // If there's something wrong with the calorimetry object return null
  if(calo->XYZ().size() != nhits || nhits < 1) return -99999;

  // Get the distance from the track point and the start/end of calo data
  double distStart = (calo->XYZ()[0] - end).Mag2();
  double distEnd = (calo->XYZ()[nhits-1] - end).Mag2();

  double maxDedx = 0;
  double resrgStart = 0;
  std::vector<double> v_resrg;
  std::vector<double> v_dedx;
  // Loop over plane's calorimetry data
  for(size_t i = 0; i < nhits; i++){
    double dedx = calo->dEdx()[i];
    double resrg = calo->ResidualRange()[i];

    // Flip the residual range if the track and calo objects don't match up
    if(distStart < distEnd && calo->ResidualRange()[0] > calo->ResidualRange()[nhits-1]) resrg = calo->ResidualRange()[0] - calo->ResidualRange()[i];
    if(distStart > distEnd && calo->ResidualRange()[0] < calo->ResidualRange()[nhits-1]) resrg = calo->ResidualRange()[nhits-1] - calo->ResidualRange()[i];

    // Find the maximum dE/dx within a limit and the corresponding res range
    if(resrg < fResRangeMin && dedx > maxDedx && dedx < fDEdxMax){
      maxDedx = dedx;
      resrgStart = resrg;
    }

  }

  // Loop over it again
  for(size_t i = 0; i < nhits; i++){
    double dedx = calo->dEdx()[i];
    double resrg = calo->ResidualRange()[i];

    // Flip the residual range if the track and calo objects don't match up
    if(distStart < distEnd && calo->ResidualRange()[0] > calo->ResidualRange()[nhits-1]) resrg = calo->ResidualRange()[0] - calo->ResidualRange()[i];
    if(distStart > distEnd && calo->ResidualRange()[0] < calo->ResidualRange()[nhits-1]) resrg = calo->ResidualRange()[nhits-1] - calo->ResidualRange()[i];

    // Record all dE/dx and residual ranges below limits
    if(resrg > resrgStart && resrg < resrgStart + fResRangeMax && dedx < fDEdxMax){
      v_resrg.push_back(resrg);
      v_dedx.push_back(dedx);
    }
  }

  // Return null value if not enough points to do fits
  if(v_dedx.size() < 10) return -99999;

  // Try to do a pol0 fit
  TGraph *gdedx = new TGraph(v_dedx.size(), &v_resrg[0], &v_dedx[0]);
  try{ gdedx->Fit("pol0", "Q"); } catch(...){ return -99999; }
  TF1* polfit = gdedx->GetFunction("pol0");
  double polchi2 = polfit->GetChisquare();

  // Try to do and exp fit
  try{ gdedx->Fit("expo", "Q"); } catch(...){ return -99999; }
  TF1* expfit = gdedx->GetFunction("expo");
  double expchi2 = expfit->GetChisquare();

  // Return the chi2 ratio
  return polchi2/expchi2;

}


// Determine if the track end looks like it stops
bool StoppingParticleCosmicIdAlg::StoppingEnd(geo::Point_t end, std::vector<art::Ptr<anab::Calorimetry>> calos){
  
  // Get the chi2 ratio
  double chiSqRatio = StoppingChiSq(end, calos);

  // If it is below a limit then tag it as stopping
  if(chiSqRatio > fStoppingChi2Limit) return true;

  return false;
}

// Determine if a track looks like a stopping cosmic
bool StoppingParticleCosmicIdAlg::StoppingParticleCosmicId(recob::Track track, std::vector<art::Ptr<anab::Calorimetry>> calos){

  // Check if start and end of track is inside the fiducial volume
  bool startInFiducial = fTpcGeo.InFiducial(track.Vertex(), fMinX, fMinY, fMinZ, fMaxX, fMaxY, fMaxZ);
  bool endInFiducial = fTpcGeo.InFiducial(track.End(), fMinX, fMinY, fMinZ, fMaxX, fMaxY, fMaxZ);

  // Check if the start and end of track stops
  bool startStops = StoppingEnd(track.Vertex(), calos);
  bool endStops = StoppingEnd(track.End(), calos);

  // If one end stops and the other end is outside of the FV then tag as stopping cosmic if the track is going downwards
  if((startStops && !endInFiducial && track.End().Y() > track.Vertex().Y()) 
      || (endStops && !startInFiducial && track.Vertex().Y() > track.End().Y())){
    return true;
  }

  return false;

}

// Determine if two tracks look like a stopping cosmic if they are merged
bool StoppingParticleCosmicIdAlg::StoppingParticleCosmicId(recob::Track track, recob::Track track2, std::vector<art::Ptr<anab::Calorimetry>> calos, std::vector<art::Ptr<anab::Calorimetry>> calos2){

  // Assume both tracks start from the same vertex so take end points as new start/end
  bool startInFiducial = fTpcGeo.InFiducial(track.End(), fMinX, fMinY, fMinZ, fMaxX, fMaxY, fMaxZ);
  bool endInFiducial = fTpcGeo.InFiducial(track.End(), fMinX, fMinY, fMinZ, fMaxX, fMaxY, fMaxZ);

  bool startStops = StoppingEnd(track.End(), calos);
  bool endStops = StoppingEnd(track2.End(), calos2);

  if((startStops && !endInFiducial && track2.End().Y() > 0) 
      || (endStops && !startInFiducial && track.End().Y() > 0)){
    return true;
  }

  return false;

}


// Calculate chi2 for particles which do not stop
std::vector<double> StoppingParticleCosmicIdAlg::FlatChi2(std::vector<art::Ptr<anab::Calorimetry>> calos){
  std::vector<double> chi2_vec;
  for(size_t i = 0; i < calos.size(); i++){
    std::vector<float> trkdedx = calos[i]->dEdx();
    std::vector<float> trkres = calos[i]->ResidualRange();
    double chi2 = 0;
    int npt = 0;
    for(size_t j = 0; j < trkdedx.size(); j++){
      if(j==0 || j==trkdedx.size()-1) continue;
      if(trkdedx[j] > 1000) continue;
      double errdedx = 0.04231+0.0001783*trkdedx[j]*trkdedx[j];
      errdedx *= trkdedx[j];
      chi2 += std::pow((trkdedx[i]-1.8)/std::sqrt(std::pow(0.5,2)+std::pow(errdedx,2)),2);
      npt++;
    }
    chi2_vec.push_back(chi2/npt);
  }
  return chi2_vec;
}


}

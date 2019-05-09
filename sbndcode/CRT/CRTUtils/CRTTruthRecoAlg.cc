#include "CRTTruthRecoAlg.h"

namespace sbnd{

CRTTruthRecoAlg::CRTTruthRecoAlg(const Config& config){

  this->reconfigure(config);
  
  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
  fAuxDetGeo = &(*fAuxDetGeoService);
  fAuxDetGeoCore = fAuxDetGeo->GetProviderPtr();

  nameToInd["volTaggerSideRight_0"] = 0;
  nameToInd["volTaggerSideLeft_0"] = 1;
  nameToInd["volTaggerBot_0"] = 2;
  nameToInd["volTaggerTopLow_0"] = 3;
  nameToInd["volTaggerTopHigh_0"] = 4;
  nameToInd["volTaggerFaceFront_0"] = 5;
  nameToInd["volTaggerFaceBack_0"] = 6;

}


CRTTruthRecoAlg::CRTTruthRecoAlg(){

  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); 
  fAuxDetGeo = &(*fAuxDetGeoService);
  fAuxDetGeoCore = fAuxDetGeo->GetProviderPtr();

  nameToInd["volTaggerSideRight_0"] = 0;
  nameToInd["volTaggerSideLeft_0"] = 1;
  nameToInd["volTaggerBot_0"] = 2;
  nameToInd["volTaggerTopLow_0"] = 3;
  nameToInd["volTaggerTopHigh_0"] = 4;
  nameToInd["volTaggerFaceFront_0"] = 5;
  nameToInd["volTaggerFaceBack_0"] = 6;

}


CRTTruthRecoAlg::~CRTTruthRecoAlg(){

}


void CRTTruthRecoAlg::reconfigure(const Config& config){

  return;
}
  
bool CRTTruthRecoAlg::CrossesTagger(const simb::MCParticle& particle, int tag_i){

  double tagCenter[3] = {0, 0, 208.25};
  tagCenter[fixCoord[tag_i*2]] = (crtPlanes[tag_i*2]+crtPlanes[tag_i*2+1])/2;
  double tagDim[3] = {0, 0, 0};
  if(tag_i==0 || tag_i==1){ tagDim[0] = 1.8; tagDim[1] = 360; tagDim[2] = 450; }
  if(tag_i==2){ tagDim[0] = 399.5; tagDim[1] = 1.8; tagDim[2] = 478; }
  if(tag_i==3 || tag_i==4){ tagDim[0] = 450; tagDim[1] = 1.8; tagDim[2] = 450; }
  if(tag_i==5 || tag_i==6){ tagDim[0] = 360; tagDim[1] = 360; tagDim[2] = 1.8; }
  bool crosses = false;
  // Get the trajectory of the true particle
  size_t npts = particle.NumberTrajectoryPoints();
  // Loop over particle trajectory
  for (size_t i = 0; i < npts; i++){
    double trajPoint[3] = {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
    // If the particle is inside the tagger volume then set to true.
    if(trajPoint[0]>tagCenter[0]-tagDim[0] && trajPoint[0]<tagCenter[0]+tagDim[0] &&
       trajPoint[1]>tagCenter[1]-tagDim[1] && trajPoint[1]<tagCenter[1]+tagDim[1] &&
       trajPoint[2]>tagCenter[2]-tagDim[2] && trajPoint[2]<tagCenter[2]+tagDim[2]) crosses = true;
  }
  return crosses;

}

bool CRTTruthRecoAlg::CrossesStrip(const simb::MCParticle& particle, int tag_i){

  double tagCenter[3] = {0, 0, 208.25};
  tagCenter[fixCoord[tag_i*2]] = (crtPlanes[tag_i*2]+crtPlanes[tag_i*2+1])/2;
  double tagDim[3] = {0, 0, 0};
  if(tag_i==0 || tag_i==1){ tagDim[0] = 1.8; tagDim[1] = 360; tagDim[2] = 450; }
  if(tag_i==2){ tagDim[0] = 399.5; tagDim[1] = 1.8; tagDim[2] = 478; }
  if(tag_i==3 || tag_i==4){ tagDim[0] = 450; tagDim[1] = 1.8; tagDim[2] = 450; }
  if(tag_i==5 || tag_i==6){ tagDim[0] = 360; tagDim[1] = 360; tagDim[2] = 1.8; }
  bool crosses = false;
  // Get the trajectory of the true particle
  size_t npts = particle.NumberTrajectoryPoints();
  // Loop over particle trajectory
  for (size_t i = 0; i < npts; i++){
    double trajPoint[3] = {particle.Vx(i), particle.Vy(i), particle.Vz(i)};
    // If the particle is inside the tagger volume then set to true.
    if(trajPoint[0]>tagCenter[0]-tagDim[0] && trajPoint[0]<tagCenter[0]+tagDim[0] &&
       trajPoint[1]>tagCenter[1]-tagDim[1] && trajPoint[1]<tagCenter[1]+tagDim[1] &&
       trajPoint[2]>tagCenter[2]-tagDim[2] && trajPoint[2]<tagCenter[2]+tagDim[2]){ 
        try{
          size_t adID, svID;
          fAuxDetGeoCore->PositionToAuxDetChannel(trajPoint,adID,svID);
          crosses = true;
        }catch(...){}
    }
  }
  return crosses;

}


TVector3 CRTTruthRecoAlg::TaggerCrossPoint(const simb::MCParticle& particle, int tag_i){

  double tagCenter[3] = {0, 0, 208.25};
  tagCenter[fixCoord[tag_i*2]] = (crtPlanes[tag_i*2]+crtPlanes[tag_i*2+1])/2;
  double tagDim[3] = {0, 0, 0};
  if(tag_i==0 || tag_i==1){ tagDim[0] = 1.8; tagDim[1] = 360; tagDim[2] = 450; }
  if(tag_i==2){ tagDim[0] = 399.5; tagDim[1] = 1.8; tagDim[2] = 478; }
  if(tag_i==3 || tag_i==4){ tagDim[0] = 450; tagDim[1] = 1.8; tagDim[2] = 450; }
  if(tag_i==5 || tag_i==6){ tagDim[0] = 360; tagDim[1] = 360; tagDim[2] = 1.8; }
  TVector3 start, end;
  bool first = true;
  // Get the trajectory of the true particle
  size_t npts = particle.NumberTrajectoryPoints();
  // Loop over particle trajectory
  for (size_t i = 0; i < npts; i++){
    TVector3 trajPoint(particle.Vx(i), particle.Vy(i), particle.Vz(i));
    // If the particle is inside the tagger volume then set to true.
    if(trajPoint[0]>tagCenter[0]-tagDim[0] && trajPoint[0]<tagCenter[0]+tagDim[0] &&
       trajPoint[1]>tagCenter[1]-tagDim[1] && trajPoint[1]<tagCenter[1]+tagDim[1] &&
       trajPoint[2]>tagCenter[2]-tagDim[2] && trajPoint[2]<tagCenter[2]+tagDim[2]){
      if(first) start = trajPoint;
      first = false;
      end = trajPoint;
    }
  }
  TVector3 crossPoint((start.X()+end.X())/2,(start.Y()+end.Y())/2,(start.Z()+end.Z())/2);
  return crossPoint;

}


bool CRTTruthRecoAlg::IsThroughGoing(const simb::MCParticle& particle){

  // Check if particle starts and ends outside the CRT planes
  bool startOutside = false;
  bool endOutside = false;
  TVector3 start(particle.Vx(), particle.Vy(), particle.Vz());
  TVector3 end(particle.EndX(), particle.EndY(), particle.EndZ());
  if(start[0]<crtPlanes[0] || start[0]>crtPlanes[3] ||
     start[1]<crtPlanes[4] || start[1]>crtPlanes[9] ||
     start[2]<crtPlanes[10] || start[2]>crtPlanes[13]) startOutside = true;
  if(end[0]<crtPlanes[0] || end[0]>crtPlanes[3] ||
     end[1]<crtPlanes[4] || end[1]>crtPlanes[9] ||
     end[2]<crtPlanes[10] || end[2]>crtPlanes[13]) endOutside = true;

  // Check if particle enters the TPC
  bool enters = false;

  // Loop over trajectory points
  int nTrajPoints = particle.NumberTrajectoryPoints();
  for (int traj_i = 0; traj_i < nTrajPoints; traj_i++){
    TVector3 trajPoint(particle.Vx(traj_i), particle.Vy(traj_i), particle.Vz(traj_i));
    // Check if point is within reconstructable volume
    if (trajPoint[0] >= crtPlanes[0] && trajPoint[0] <= crtPlanes[3] && trajPoint[1] >= crtPlanes[4] && trajPoint[1] <= -crtPlanes[4] && trajPoint[2] >= crtPlanes[10] && trajPoint[2] <= crtPlanes[13]){
      enters = true;
    }
  }
  // Only count as through going if particle starts and ends outside the CRT 
  // enclosed volume and enters the TPC
  if(startOutside && endOutside && enters) return true;
  return false;
}

void CRTTruthRecoAlg::DrawCube(TCanvas *c1, double *rmin, double *rmax, int col){

  c1->cd();
  TList *outline = new TList;
  TPolyLine3D *p1 = new TPolyLine3D(4);
  TPolyLine3D *p2 = new TPolyLine3D(4);
  TPolyLine3D *p3 = new TPolyLine3D(4);
  TPolyLine3D *p4 = new TPolyLine3D(4);
  p1->SetLineColor(col);
  p1->SetLineWidth(2);
  p1->Copy(*p2);
  p1->Copy(*p3);
  p1->Copy(*p4);
  outline->Add(p1);
  outline->Add(p2);
  outline->Add(p3);
  outline->Add(p4); 
  TPolyLine3D::DrawOutlineCube(outline, rmin, rmax);
  p1->Draw();
  p2->Draw();
  p3->Draw();
  p4->Draw();

}

// Function to calculate the CRT crossing points of a true particle
std::pair<TVector3, TVector3> CRTTruthRecoAlg::TpcCrossPoints(simb::MCParticle const& particle){

  double xmin = fTpcGeo.MinX();
  double xmax = fTpcGeo.MaxX();
  double ymin = fTpcGeo.MinY();
  double ymax = fTpcGeo.MaxY();
  double zmin = fTpcGeo.MinZ();
  double zmax = fTpcGeo.MaxZ();
  TVector3 start, end;

  bool first = true;
  // Get the trajectory of the true particle
  size_t npts = particle.NumberTrajectoryPoints();

  // Loop over particle trajectory
  for (size_t i = 0; i < npts; i++){
    TVector3 trajPoint(particle.Vx(i), particle.Vy(i), particle.Vz(i));
    // If the particle is inside the tagger volume then set to true.
    if(trajPoint[0]>xmin && trajPoint[0]<xmax &&
       trajPoint[1]>ymin && trajPoint[1]<ymax &&
       trajPoint[2]>zmin && trajPoint[2]<zmax){
      if(first) start = trajPoint;
      first = false;
      end = trajPoint;
    }
  }

  return std::make_pair(start, end);

} // CRTTrackMatchingAna::TpcCrossPoints()

// Function to calculate the CRT crossing points of a true particle
double CRTTruthRecoAlg::TpcLength(simb::MCParticle const& particle){

  double xmin = fTpcGeo.MinX();
  double xmax = fTpcGeo.MaxX();
  double ymin = fTpcGeo.MinY();
  double ymax = fTpcGeo.MaxY();
  double zmin = fTpcGeo.MinZ();
  double zmax = fTpcGeo.MaxZ();
  
  TVector3 start, end;
  TVector3 disp;
  double length = 0.;

  bool first = true;
  // Get the trajectory of the true particle
  size_t npts = particle.NumberTrajectoryPoints();

  // Loop over particle trajectory
  for (size_t i = 0; i < npts; i++){
    TVector3 trajPoint(particle.Vx(i), particle.Vy(i), particle.Vz(i));
    // If the particle is inside the tagger volume then set to true.
    if(trajPoint[0]>xmin && trajPoint[0]<xmax &&
       trajPoint[1]>ymin && trajPoint[1]<ymax &&
       trajPoint[2]>zmin && trajPoint[2]<zmax){
      if(first) start = trajPoint;
      else{
        disp -= trajPoint;
        length += disp.Mag();
      }
      first = false;
      disp = trajPoint;
      end = trajPoint;
    }
  }

  return length;

} // CRTTrackMatchingAna::TpcLength()

// Function to project a track position on to a tagger
TVector3 CRTTruthRecoAlg::T0ToXYZPosition(TVector3 position, TVector3 direction, std::string tagger, int tpc, double t0){

  //Here crt_i is index of tagger, so runs from 0 to 6
  TVector3 returnVal(-99999, -99999, -99999);
  int crt_i = nameToInd[tagger];

  // Convert the t0 into an x shift
  double shift = t0 * fDetectorProperties->DriftVelocity();
  // Apply the shift depending on which TPC the track is in
  if (tpc == 1) position[0] += shift;
  if (tpc == 0) position[0] -= shift;

  // Calculate the step to the CRT plane
  double step = (crtPlanes[crt_i*2]- position[fixCoord[crt_i*2]])/direction[fixCoord[crt_i*2]];

  // If the step is < 0 return a null position
  if (step < 0) return returnVal;

  // Calculate the CRT crossing point of the output coordinate
  returnVal[lenCoord[crt_i*2]] = position[lenCoord[crt_i*2]] + step*direction[lenCoord[crt_i*2]];
  returnVal[widthCoord[crt_i*2]] = position[widthCoord[crt_i*2]] + step*direction[widthCoord[crt_i*2]];
  returnVal[fixCoord[crt_i*2]] = crtPlanes[crt_i*2];

  return returnVal;

} // CRTTruthRecoAlg::T0ToXYZPosition()

}

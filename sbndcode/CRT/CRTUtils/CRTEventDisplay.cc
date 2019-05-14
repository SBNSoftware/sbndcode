#include "CRTEventDisplay.h"

namespace sbnd{

CRTEventDisplay::CRTEventDisplay(const Config& config){

  this->reconfigure(config);
  
}


CRTEventDisplay::CRTEventDisplay(){

}


CRTEventDisplay::~CRTEventDisplay(){

}


void CRTEventDisplay::reconfigure(const Config& config){

  fSimLabel = config.SimLabel();
  fCRTDataLabel = config.CRTDataLabel();
  fCRTHitLabel = config.CRTHitLabel();
  fCRTTrackLabel = config.CRTTrackLabel();
  fTPCTrackLabel = config.TPCTrackLabel();

  fDrawTaggers = config.DrawTaggers();
  fDrawModules = config.DrawModules();
  fDrawTpc = config.DrawTpc();
  fDrawCrtHits = config.DrawCrtHits();
  fDrawCrtTracks = config.DrawCrtTracks();
  fDrawIncompleteTracks = config.DrawIncompleteTracks();
  fDrawTpcTracks = config.DrawTpcTracks();
  fDrawTrueTracks = config.DrawTrueTracks();

  fTaggerColour = config.TaggerColour();
  fTpcColour = config.TpcColour();
  fCrtHitColour = config.CrtHitColour();
  fCrtTrackColour = config.CrtTrackColour();
  fTpcTrackColour = config.TpcTrackColour();
  fTrueTrackColour = config.TrueTrackColour();

  fUseTrueID = config.UseTrueID();
  fTrueID = config.TrueID();

  fPrint = config.Print();

  fLineWidth             = config.LineWidth();
  fIncompleteTrackLength = config.IncompleteTrackLength();
  fMinTime               = config.MinTime();
  fMaxTime               = config.MaxTime();

  return;

}
 
void CRTEventDisplay::SetDrawTaggers(bool tf){
  fDrawTaggers = tf;
}
void CRTEventDisplay::SetDrawTpc(bool tf){
  fDrawTpc = tf;
}
void CRTEventDisplay::SetDrawCrtHits(bool tf){
  fDrawCrtHits = tf;
}
void CRTEventDisplay::SetDrawCrtTracks(bool tf){
  fDrawCrtTracks = tf;
}
void CRTEventDisplay::SetDrawTpcTracks(bool tf){
  fDrawTpcTracks = tf;
}
void CRTEventDisplay::SetDrawTrueTracks(bool tf){
  fDrawTrueTracks = tf;
}

bool CRTEventDisplay::IsVisible(const simb::MCParticle& particle){
  int pdg = std::abs(particle.PdgCode());
  double momentum = particle.P();
  double momLimit = 0.05;
  if(momentum < momLimit) return false;
  if(pdg == 13) return true;
  if(pdg == 11) return true;
  if(pdg == 2212) return true;
  if(pdg == 211) return true;
  if(pdg == 321) return true;
  return false;
}

void CRTEventDisplay::DrawCube(TCanvas *c1, double *rmin, double *rmax, int colour){

  c1->cd();
  TList *outline = new TList;
  TPolyLine3D *p1 = new TPolyLine3D(4);
  TPolyLine3D *p2 = new TPolyLine3D(4);
  TPolyLine3D *p3 = new TPolyLine3D(4);
  TPolyLine3D *p4 = new TPolyLine3D(4);
  p1->SetLineColor(colour);
  p1->SetLineWidth(fLineWidth);
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

void CRTEventDisplay::Draw(const art::Event& event){
  // Create a canvas 
  TCanvas *c1 = new TCanvas("c1","",700,700);

  if(fDrawTaggers){
    for(size_t i = 0; i < fCrtGeo.NumTaggers(); i++){
      double rmin[3] = {fCrtGeo.GetTagger(i).minX, 
                        fCrtGeo.GetTagger(i).minY, 
                        fCrtGeo.GetTagger(i).minZ};
      double rmax[3] = {fCrtGeo.GetTagger(i).maxX, 
                        fCrtGeo.GetTagger(i).maxY, 
                        fCrtGeo.GetTagger(i).maxZ};
      DrawCube(c1, rmin, rmax, fTaggerColour);
    }
  }

  if(fDrawModules){
    for(size_t i = 0; i < fCrtGeo.NumModules(); i++){
      double rmin[3] = {fCrtGeo.GetModule(i).minX, 
                        fCrtGeo.GetModule(i).minY, 
                        fCrtGeo.GetModule(i).minZ};
      double rmax[3] = {fCrtGeo.GetModule(i).maxX, 
                        fCrtGeo.GetModule(i).maxY, 
                        fCrtGeo.GetModule(i).maxZ};
      DrawCube(c1, rmin, rmax, fTaggerColour);
    }
  }

  if(fDrawTpc){
    double rmin[3] = {fTpcGeo.MinX(), 
                      fTpcGeo.MinY(), 
                      fTpcGeo.MinZ()};
    double rmax[3] = {-fTpcGeo.CpaWidth(), 
                      fTpcGeo.MaxY(), 
                      fTpcGeo.MaxZ()};
    DrawCube(c1, rmin, rmax, fTpcColour);
    double rmin2[3] = {fTpcGeo.CpaWidth(), 
                      fTpcGeo.MinY(), 
                      fTpcGeo.MinZ()};
    double rmax2[3] = {fTpcGeo.MaxX(), 
                      fTpcGeo.MaxY(), 
                      fTpcGeo.MaxZ()};
    DrawCube(c1, rmin2, rmax2, fTpcColour);
  }

  if(fDrawCrtHits){
    auto crtHitHandle = event.getValidHandle<std::vector<crt::CRTHit>>(fCRTHitLabel);
    for(auto const& hit : (*crtHitHandle)){
      double time = (double)(int)hit.ts1_ns * 1e-3;
      if(fMinTime == fMaxTime || (time > fMinTime && time < fMaxTime)){
      double rmin[3] = {hit.x_pos - hit.x_err,
                        hit.y_pos - hit.y_err,
                        hit.z_pos - hit.z_err};
      double rmax[3] = {hit.x_pos + hit.x_err,
                        hit.y_pos + hit.y_err,
                        hit.z_pos + hit.z_err};
      DrawCube(c1, rmin, rmax, fCrtHitColour);
      }
    }
  }

  if(fDrawCrtTracks){
    auto crtTrackHandle = event.getValidHandle<std::vector<crt::CRTTrack>>(fCRTTrackLabel);
    for(auto const& track : (*crtTrackHandle)){
      double time = (double)(int)track.ts1_ns * 1e-3; 
      if(fMinTime == fMaxTime || (time > fMinTime && time < fMaxTime)){
      TPolyLine3D *line = new TPolyLine3D(2);
      line->SetPoint(0, track.x1_pos, track.y1_pos, track.z1_pos);
      line->SetPoint(1, track.x2_pos, track.y2_pos, track.z2_pos);
      line->SetLineColor(fCrtTrackColour);
      line->SetLineWidth(fLineWidth);
      if(track.complete) line->Draw();
      else if(fDrawIncompleteTracks){
        TVector3 start(track.x1_pos, track.y1_pos, track.z1_pos);
        TVector3 end(track.x2_pos, track.y2_pos, track.z2_pos);
        if(start.Y() < end.Y()) std::swap(start, end);
        TVector3 diff = (end - start).Unit();
        TVector3 newEnd = start + fIncompleteTrackLength * diff;
        line->SetPoint(0, start.X(), start.Y(), start.Z());
        line->SetPoint(1, newEnd.X(), newEnd.Y(), newEnd.Z());
        line->Draw();
      }
      }
    }
  }

  if(fDrawTpcTracks){
    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    for(auto const& track : (*tpcTrackHandle)){
      size_t npts = track.NumberTrajectoryPoints();
      TPolyLine3D *line = new TPolyLine3D(npts);
      int ipt = 0;
      for(size_t i = 0; i < npts; i++){
        auto& pos = track.LocationAtPoint(i);
        if(pos.X() == -999 || (pos.X() == 0 && pos.Y() == 0)) continue; 
        line->SetPoint(ipt, pos.X(), pos.Y(), pos.Z());
        ipt++;
      }
      line->SetLineColor(fTpcTrackColour);
      line->SetLineWidth(fLineWidth);
      line->Draw();
    }
  }

  if(fDrawTrueTracks){
    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimLabel);
    std::vector<double> crtLims = fCrtGeo.CRTLimits();
    for(auto const& part : (*particleHandle)){
      if(fUseTrueID && part.TrackId() != fTrueID) continue;
      double time = part.T() * 1e-3;
      if(fMinTime == fMaxTime || (time > fMinTime && time < fMaxTime)){
      if(!IsVisible(part)) continue;
      if(!fCrtGeo.EntersVolume(part)) continue;
      size_t npts = part.NumberTrajectoryPoints();
      TPolyLine3D *line = new TPolyLine3D(npts);
      int ipt = 0;
      for(size_t i = 0; i < npts; i++){
        double px = part.Vx(i);
        double py = part.Vy(i);
        double pz = part.Vz(i);
        if(px < crtLims[0] || px > crtLims[3] || py < crtLims[1] 
           || py > crtLims[4] || pz < crtLims[2] || pz > crtLims[5]) continue;
        line->SetPoint(ipt, px, py, pz);
        ipt++;
      }
      line->SetLineColor(fTrueTrackColour);
      line->SetLineWidth(fLineWidth);
      line->Draw();
      }
    }
  }

  c1->SaveAs("crtEventDisplay.root");

}

}

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
  fClockSpeedCRT = config.ClockSpeedCRT();

  fDrawTaggers = config.DrawTaggers();
  fDrawModules = config.DrawModules();
  fDrawTpc = config.DrawTpc();
  fDrawCrtData = config.DrawCrtData();
  fDrawCrtHits = config.DrawCrtHits();
  fDrawCrtTracks = config.DrawCrtTracks();
  fDrawIncompleteTracks = config.DrawIncompleteTracks();
  fDrawTpcTracks = config.DrawTpcTracks();
  fDrawTrueTracks = config.DrawTrueTracks();

  fTaggerColour = config.TaggerColour();
  fTpcColour = config.TpcColour();
  fCrtDataColour = config.CrtDataColour();
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
  
  fCrtBackTrack          = config.CrtBackTrack();

  return;

}
 
void CRTEventDisplay::SetDrawTaggers(bool tf){
  fDrawTaggers = tf;
}
void CRTEventDisplay::SetDrawTpc(bool tf){
  fDrawTpc = tf;
}
void CRTEventDisplay::SetDrawCrtData(bool tf){
  fDrawCrtData = tf;
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
void CRTEventDisplay::SetPrint(bool tf){
  fPrint = tf;
}

void CRTEventDisplay::SetTrueId(int id){
  fUseTrueID = true;
  fTrueID = id;
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

void CRTEventDisplay::Draw(detinfo::DetectorClocksData const& clockData,
                           const art::Event& event){
  // Create a canvas 
  TCanvas *c1 = new TCanvas("c1","",700,700);

  // Draw the CRT taggers
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

  // Draw individual CRT modules
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

  // Draw the TPC with central cathode
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

  // Draw the CRT data in the event
  if(fDrawCrtData){

    if(fPrint) std::cout<<"\nCRT data in event:\n";

    auto crtDataHandle = event.getValidHandle<std::vector<sbnd::crt::CRTData>>(fCRTDataLabel);
    //detinfo::DetectorClocks const* fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>();
    //detinfo::ElecClock fTrigClock = fDetectorClocks->TriggerClock();
    for(auto const& data : (*crtDataHandle)){

      // Skip if outside specified time window if time window used
      //fTrigClock.SetTime(data.T0());
      //double time = fTrigClock.Time();
      double time = (double)(int)data.T0()/fClockSpeedCRT; // [tick -> us]
      if(!(fMinTime == fMaxTime || (time > fMinTime && time < fMaxTime))) continue;

      // Skip if it doesn't match the true ID if true ID is used
      int trueId = fCrtBackTrack.TrueIdFromTotalEnergy(event, data);
      if(fUseTrueID && trueId != fTrueID) continue;

      std::string stripName = fCrtGeo.ChannelToStripName(data.Channel());
      double rmin[3] = {fCrtGeo.GetStrip(stripName).minX, 
                        fCrtGeo.GetStrip(stripName).minY, 
                        fCrtGeo.GetStrip(stripName).minZ};
      double rmax[3] = {fCrtGeo.GetStrip(stripName).maxX, 
                        fCrtGeo.GetStrip(stripName).maxY, 
                        fCrtGeo.GetStrip(stripName).maxZ};
      DrawCube(c1, rmin, rmax, fCrtDataColour);

      if(fPrint) std::cout<<"->True ID: "<<trueId<<", channel = "<<data.Channel()<<", tagger = "
                          <<fCrtGeo.GetModule(fCrtGeo.GetStrip(stripName).module).tagger<<", time = "<<time<<"\n";
    }
  }

  // Draw the CRT hits in the event
  if(fDrawCrtHits){

    if(fPrint) std::cout<<"\nCRT hits in event:\n";

    auto crtHitHandle = event.getValidHandle<std::vector<sbn::crt::CRTHit>>(fCRTHitLabel);
    for(auto const& hit : (*crtHitHandle)){

      // Skip if outside specified time window if time window used
      double time = (double)(int)hit.ts1_ns * 1e-3;
      if(!(fMinTime == fMaxTime || (time > fMinTime && time < fMaxTime))) continue;

      // Skip if it doesn't match the true ID if true ID is used
      int trueId = fCrtBackTrack.TrueIdFromTotalEnergy(event, hit);
      if(fUseTrueID && trueId != fTrueID) continue;

      double rmin[3] = {hit.x_pos - hit.x_err,
                        hit.y_pos - hit.y_err,
                        hit.z_pos - hit.z_err};
      double rmax[3] = {hit.x_pos + hit.x_err,
                        hit.y_pos + hit.y_err,
                        hit.z_pos + hit.z_err};
      DrawCube(c1, rmin, rmax, fCrtHitColour);

      if(fPrint) std::cout<<"->True ID: "<<trueId<<", position = ("<<hit.x_pos<<", "
                          <<hit.y_pos<<", "<<hit.z_pos<<"), time = "<<time<<"\n";
    }
  }

  // Draw CRT tracks in the event
  if(fDrawCrtTracks){

    if(fPrint) std::cout<<"\nCRT tracks in event:\n";

    auto crtTrackHandle = event.getValidHandle<std::vector<sbn::crt::CRTTrack>>(fCRTTrackLabel);
    for(auto const& track : (*crtTrackHandle)){

      // Skip if outside specified time window if time window used
      double time = (double)(int)track.ts1_ns * 1e-3; 
      if(!(fMinTime == fMaxTime || (time > fMinTime && time < fMaxTime))) continue;

      // Skip if it doesn't match the true ID if true ID is used
      int trueId = fCrtBackTrack.TrueIdFromTotalEnergy(event, track);
      if(fUseTrueID && trueId != fTrueID) continue;

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

      if(fPrint) std::cout<<"->True ID: "<<trueId<<", start = ("<<track.x1_pos<<", "
                          <<track.y1_pos<<", "<<track.z1_pos<<"), end = ("<<track.x2_pos
                          <<", "<<track.y2_pos<<", "<<track.z2_pos<<"), time = "<<time<<"\n";
    }
  }

  // Draw reconstructed TPC tracks in the event
  if(fDrawTpcTracks){

    if(fPrint) std::cout<<"\nTPC tracks in event:\n";

    auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTPCTrackLabel);
    art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTPCTrackLabel);
    for(auto const& track : (*tpcTrackHandle)){

      // Skip if it doesn't match the true ID if true ID is used
      std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(track.ID());
      int trueId = RecoUtils::TrueParticleIDFromTotalRecoHits(clockData, hits, false);
      if(fUseTrueID && trueId != fTrueID) continue;
      
      size_t npts = track.NumberTrajectoryPoints();
      TPolyLine3D *line = new TPolyLine3D(npts);
      int ipt = 0;
      bool first = true;
      geo::Point_t start {0,0,0};
      geo::Point_t end {0,0,0};
      for(size_t i = 0; i < npts; i++){
        auto& pos = track.LocationAtPoint(i);

        // Don't draw invalid points
        if(!track.HasValidPoint(i)) continue;
        //if(pos.X() == -999 || (pos.X() == 0 && pos.Y() == 0)) continue; 

        if(first){
          first = false;
          start = pos;
        }
        end = pos;

        line->SetPoint(ipt, pos.X(), pos.Y(), pos.Z());
        ipt++;
      }
      line->SetLineColor(fTpcTrackColour);
      line->SetLineWidth(fLineWidth);
      line->Draw();

      if(fPrint) std::cout<<"->True ID: "<<trueId<<", start = ("<<start.X()<<", "<<start.Y()<<", "
                          <<start.Z()<<"), end = ("<<end.X()<<", "<<end.Y()<<", "<<end.Z()<<")\n";
    }
  }

  // Draw true track trajectories for visible particles that cross the CRT
  if(fDrawTrueTracks){

    if(fPrint) std::cout<<"\nTrue tracks in event:\n";

    auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimLabel);
    std::vector<double> crtLims = fCrtGeo.CRTLimits();
    for(auto const& part : (*particleHandle)){

      // Skip if it doesn't match the true ID if true ID is used
      if(fUseTrueID && part.TrackId() != fTrueID) continue;

      // Skip if outside specified time window if time window used
      double time = part.T() * 1e-3;
      if(!(fMinTime == fMaxTime || (time > fMinTime && time < fMaxTime))) continue;

      // Skip if particle isn't visible
      if(!IsVisible(part)) continue;

      // Skip if particle doesn't cross the boundary enclosed by the CRTs
      if(!fCrtGeo.EntersVolume(part)) continue;

      size_t npts = part.NumberTrajectoryPoints();
      TPolyLine3D *line = new TPolyLine3D(npts);
      int ipt = 0;
      bool first = true;
      geo::Point_t start {0,0,0};
      geo::Point_t end {0,0,0};
      for(size_t i = 0; i < npts; i++){
        geo::Point_t pos {part.Vx(i), part.Vy(i), part.Vz(i)};

        // Don't draw trajectories outside of the CRT volume
        if(pos.X() < crtLims[0] || pos.X() > crtLims[3] || pos.Y() < crtLims[1] 
           || pos.Y() > crtLims[4] || pos.Z() < crtLims[2] || pos.Z() > crtLims[5]) continue;

        if(first){
          first = false;
          start = pos;
        }
        end = pos;

        line->SetPoint(ipt, pos.X(), pos.Y(), pos.Z());
        ipt++;
      }
      line->SetLineColor(fTrueTrackColour);
      line->SetLineWidth(fLineWidth);
      line->Draw();

      if(fPrint) std::cout<<"->True ID: "<<part.TrackId()<<", start = ("<<start.X()<<", "<<start.Y()<<", "
                          <<start.Z()<<"), end = ("<<end.X()<<", "<<end.Y()<<", "<<end.Z()<<")\n";
    }
  }

  c1->SaveAs("crtEventDisplay.root");

}

}

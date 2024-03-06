#include "CRTEventDisplayAlg.h"

namespace sbnd::crt {
  
  CRTEventDisplayAlg::CRTEventDisplayAlg(const Config& config)
    : fCRTGeoAlg(config.GeoAlgConfig())
    , fCRTBackTrackerAlg(config.BackTrackerAlgConfig())
  {
    this->reconfigure(config);
  }
  
  CRTEventDisplayAlg::CRTEventDisplayAlg(){}
  
  CRTEventDisplayAlg::~CRTEventDisplayAlg(){}

  void CRTEventDisplayAlg::reconfigure(const Config& config)
  {
    fSimLabel = config.SimLabel();
    fSimDepositLabel = config.SimDepositLabel();
    fStripHitLabel = config.StripHitLabel();
    fClusterLabel = config.ClusterLabel();
    fSpacePointLabel = config.SpacePointLabel();
    fTrackLabel = config.TrackLabel();

    fDrawTaggers = config.DrawTaggers();
    fDrawModules = config.DrawModules();
    fDrawFEBs = config.DrawFEBs();
    fDrawFEBEnds = config.DrawFEBEnds();
    fDrawStrips = config.DrawStrips();
    fDrawTpc = config.DrawTpc();
    fDrawTrueTracks = config.DrawTrueTracks();
    fDrawSimDeposits = config.DrawSimDeposits();
    fDrawStripHits = config.DrawStripHits();
    fDrawClusters = config.DrawClusters();
    fDrawSpacePoints = config.DrawSpacePoints();
    fDrawTracks = config.DrawTracks();

    fChoseTaggers = config.ChoseTaggers();
    fChosenTaggers = config.ChosenTaggers();

    fTaggerColour = config.TaggerColour();
    fFEBColour = config.FEBColour();
    fFEBEndColour = config.FEBEndColour();
    fTpcColour = config.TpcColour();
    fTrueTrackColour = config.TrueTrackColour();
    fSimDepositColour = config.SimDepositColour();
    fStripHitColour = config.StripHitColour();
    fClusterStartingColour = config.ClusterStartingColour();
    fClusterColourInterval = config.ClusterColourInterval();
    fSpacePointColour = config.SpacePointColour();
    fTrackColour = config.TrackColour();

    fMinTime = config.MinTime();
    fMaxTime = config.MaxTime();

    fPrint = config.Print();

    fLineWidth = config.LineWidth();

    return;
  }
 
  void CRTEventDisplayAlg::SetDrawTaggers(bool tf)
  {
    fDrawTaggers = tf;
  }

  void CRTEventDisplayAlg::SetDrawTpc(bool tf)
  {
    fDrawTpc = tf;
  }

  void CRTEventDisplayAlg::SetDrawTrueTracks(bool tf)
  {
    fDrawTrueTracks = tf;
  }

  void CRTEventDisplayAlg::SetDrawSimDeposits(bool tf)
  {
    fDrawSimDeposits = tf;
  }

  void CRTEventDisplayAlg::SetDrawStripHits(bool tf)
  {
    fDrawStripHits = tf;
  }

  void CRTEventDisplayAlg::SetDrawClusters(bool tf)
  {
    fDrawClusters = tf;
  }

  void CRTEventDisplayAlg::SetPrint(bool tf)
  {
    fPrint = tf;
  }

  void CRTEventDisplayAlg::DrawCube(TCanvas *c1, double *rmin, double *rmax, int colour)
  {
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

  void CRTEventDisplayAlg::Draw(detinfo::DetectorClocksData const& clockData,
                                const art::Event& event)
  {
    fCRTBackTrackerAlg.SetupMaps(event);

    double G4RefTime(clockData.G4ToElecTime(0) * 1e3);
    if(fPrint) std::cout << "G4RefTime: " << G4RefTime << std::endl;

    // Create a canvas 
    TCanvas *c1 = new TCanvas("c1","",700,700);
    
    std::vector<double> crtLims = fCRTGeoAlg.CRTLimits();
    crtLims[0] -= 100; crtLims[1] -= 100; crtLims[2] -= 100;
    crtLims[3] += 100; crtLims[4] += 100; crtLims[5] += 100;


    // Draw the CRT taggers
    if(fDrawTaggers)
      {
        for(auto const &[name, tagger] : fCRTGeoAlg.GetTaggers())
          {
            if(fChoseTaggers && std::find(fChosenTaggers.begin(), fChosenTaggers.end(), CRTCommonUtils::GetTaggerEnum(name)) == fChosenTaggers.end())
              continue;
               
            double rmin[3] = {tagger.minX, 
                              tagger.minY, 
                              tagger.minZ};
            double rmax[3] = {tagger.maxX, 
                              tagger.maxY, 
                              tagger.maxZ};
            DrawCube(c1, rmin, rmax, fTaggerColour);
          }
      }
    
    // Draw individual CRT modules
    if(fDrawModules)
      {
        for(auto const &[name, module] : fCRTGeoAlg.GetModules())
          {
            if(fChoseTaggers && std::find(fChosenTaggers.begin(), fChosenTaggers.end(), CRTCommonUtils::GetTaggerEnum(module.taggerName)) == fChosenTaggers.end())
              continue;

            double rmin[3] = {module.minX, 
                              module.minY, 
                              module.minZ};
            double rmax[3] = {module.maxX, 
                              module.maxY, 
                              module.maxZ};
            DrawCube(c1, rmin, rmax, fTaggerColour);

            if(fDrawFEBs)
              {
                const std::array<double, 6> febPos = fCRTGeoAlg.FEBWorldPos(module);
                
                double rmin[3] = {febPos[0],
                                  febPos[2],
                                  febPos[4]};

                double rmax[3] = {febPos[1],
                                  febPos[3],
                                  febPos[5]};

                DrawCube(c1, rmin, rmax, fFEBColour);

                if(fDrawFEBEnds)
                  {
                    const std::array<double, 6> febCh0Pos = fCRTGeoAlg.FEBChannel0WorldPos(module);

                    double rminCh0[3] = {febCh0Pos[0],
                                         febCh0Pos[2],
                                         febCh0Pos[4]};

                    double rmaxCh0[3] = {febCh0Pos[1],
                                         febCh0Pos[3],
                                         febCh0Pos[5]};

                    DrawCube(c1, rminCh0, rmaxCh0, fFEBEndColour);
                  }
              }
          }
      }
    
    // Draw individual CRT strips
    if(fDrawStrips)
      {
        for(auto const &[name, strip] : fCRTGeoAlg.GetStrips())
          {
            if(fChoseTaggers && std::find(fChosenTaggers.begin(), fChosenTaggers.end(), fCRTGeoAlg.ChannelToTaggerEnum(strip.channel0)) == fChosenTaggers.end())
              continue;

            double rmin[3] = {strip.minX, 
                              strip.minY, 
                              strip.minZ};
            double rmax[3] = {strip.maxX, 
                              strip.maxY, 
                              strip.maxZ};
            DrawCube(c1, rmin, rmax, fTaggerColour);
          }
      }
    
    // Draw the TPC with central cathode
    if(fDrawTpc)
      {
        double rmin[3] = {fTPCGeoAlg.MinX(), 
                          fTPCGeoAlg.MinY(), 
                          fTPCGeoAlg.MinZ()};
        double rmax[3] = {-fTPCGeoAlg.CpaWidth(), 
                          fTPCGeoAlg.MaxY(), 
                          fTPCGeoAlg.MaxZ()};
        DrawCube(c1, rmin, rmax, fTpcColour);
        double rmin2[3] = {fTPCGeoAlg.CpaWidth(), 
                           fTPCGeoAlg.MinY(), 
                           fTPCGeoAlg.MinZ()};
        double rmax2[3] = {fTPCGeoAlg.MaxX(), 
                           fTPCGeoAlg.MaxY(), 
                           fTPCGeoAlg.MaxZ()};
        DrawCube(c1, rmin2, rmax2, fTpcColour);
      }
    
    // Draw true track trajectories for visible particles that cross the CRT
    if(fDrawTrueTracks)
      { 
        auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimLabel);

        for(auto const& part : *particleHandle)
          {
            if(part.T() < fMinTime || part.T() > fMaxTime)
              continue;

            size_t npts = part.NumberTrajectoryPoints();
            TPolyLine3D *line = new TPolyLine3D(npts);
            int ipt = 0;
            bool first = true;
            geo::Point_t start {0,0,0};
            geo::Point_t end {0,0,0};
            for(size_t i = 0; i < npts; i++)
              {
                geo::Point_t pos {part.Vx(i), part.Vy(i), part.Vz(i)};
              
                // Don't draw trajectories outside of the CRT volume
                if(pos.X() < crtLims[0] || pos.X() > crtLims[3] || pos.Y() < crtLims[1] 
                   || pos.Y() > crtLims[4] || pos.Z() < crtLims[2] || pos.Z() > crtLims[5]) continue;
              
                if(first)
                  {
                    first = false;
                    start = pos;
                  }
                end = pos;
              
                line->SetPoint(ipt, pos.X(), pos.Y(), pos.Z());
                ipt++;
              }
            line->SetLineColor(fTrueTrackColour);
            line->SetLineWidth(fLineWidth+4);
            line->SetLineStyle(9);
            line->Draw();
          
            if(fPrint) std::cout<<"MCParticle, Track ID: " << part.TrackId() << " PDG: " << part.PdgCode() << ", traj points: "<<npts<<", start = ("<<start.X()<<", "<<start.Y()<<", "
                                <<start.Z()<<"), end = ("<<end.X()<<", "<<end.Y()<<", "<<end.Z()<<")\n";
          }
      }

    if(fDrawSimDeposits)
      {
        auto simDepositsHandle = event.getValidHandle<std::vector<sim::AuxDetSimChannel>>(fSimDepositLabel);

        for(auto const& simDep : *simDepositsHandle)
          {
            for(auto const& ide : simDep.AuxDetIDEs())
              {
                double x = (ide.entryX + ide.exitX) / 2.;
                double y = (ide.entryY + ide.exitY) / 2.;
                double z = (ide.entryZ + ide.exitZ) / 2.;
                double t = (ide.entryT + ide.exitT) / 2.;

                if(t < fMinTime || t > fMaxTime)
                  continue;

                double ex = std::abs(ide.entryX - ide.exitX) / 2.;
                double ey = std::abs(ide.entryY - ide.exitY) / 2.;
                double ez = std::abs(ide.entryZ - ide.exitZ) / 2.;

                ex = std::max(ex, 1.);
                ey = std::max(ey, 1.);
                ez = std::max(ez, 1.);

                double rmin[3] = { x - ex, y - ey, z - ez};
                double rmax[3] = { x + ex, y + ey, z + ez};

                if(fPrint)
                  std::cout << "Sim Energy Deposit: (" << x << ", " << y << ", " << z 
                            << ")  +/- (" << ex << ", " << ey << ", " << ez << ") by trackID: " 
                            << ide.trackID << " at t = " << t << std::endl;

                DrawCube(c1, rmin, rmax, fSimDepositColour);
              }
          }     
      }

    if(fDrawStripHits)
      {
        auto stripHitsHandle = event.getValidHandle<std::vector<CRTStripHit>>(fStripHitLabel);
        std::vector<art::Ptr<CRTStripHit>> stripHitsVec;
        art::fill_ptr_vector(stripHitsVec, stripHitsHandle);

        for(auto const stripHit : stripHitsVec)
          {
            if(stripHit->Ts1() - G4RefTime < fMinTime || stripHit->Ts1() - G4RefTime > fMaxTime)
              continue;

            CRTStripGeo strip = fCRTGeoAlg.GetStrip(stripHit->Channel());

            double rmin[3] = {strip.minX, strip.minY, strip.minZ};
            double rmax[3] = {strip.maxX, strip.maxY, strip.maxZ};

            CRTBackTrackerAlg::TruthMatchMetrics truthMatch = fCRTBackTrackerAlg.TruthMatching(event, stripHit);

            if(fPrint)
              std::cout << "Strip Hit: (" 
                        << rmin[0] << ", " << rmin[1] << ", " << rmin[2] << ") --> ("
                        << rmax[0] << ", " << rmax[1] << ", " << rmax[2] << ") at t1 = " << stripHit->Ts1()
                        << " (" << stripHit->Ts1() - G4RefTime << ")"
                        << "\t Matches to trackID: " << truthMatch.trackid 
                        << " with completeness: " << truthMatch.completeness 
                        << " and purity: " << truthMatch.purity << std::endl;

            DrawCube(c1, rmin, rmax, fStripHitColour);
          }
      }

    if(fDrawClusters || fDrawSpacePoints)
      {
        auto clustersHandle = event.getValidHandle<std::vector<CRTCluster>>(fClusterLabel);
        std::vector<art::Ptr<CRTCluster>> clustersVec;
        art::fill_ptr_vector(clustersVec, clustersHandle);
        art::FindManyP<CRTStripHit> clustersToStripHits(clustersHandle, event, fClusterLabel);
        art::FindManyP<CRTSpacePoint> clustersToSpacePoints(clustersHandle, event, fSpacePointLabel);

        int colour = fClusterStartingColour;

        for(auto const cluster : clustersVec)
          {
            if(cluster->Ts1() - G4RefTime < fMinTime || cluster->Ts1() - G4RefTime > fMaxTime)
              continue;

            auto stripHitVec   = clustersToStripHits.at(cluster.key());
            auto spacePointVec = clustersToSpacePoints.at(cluster.key());

            if(fDrawClusters)
              {
                for(auto stripHit : stripHitVec)
                  {
                    CRTStripGeo strip = fCRTGeoAlg.GetStrip(stripHit->Channel());
                    
                    double rmin[3] = {strip.minX, strip.minY, strip.minZ};
                    double rmax[3] = {strip.maxX, strip.maxY, strip.maxZ};
                    
                    DrawCube(c1, rmin, rmax, colour);
                  }

                CRTBackTrackerAlg::TruthMatchMetrics truthMatch = fCRTBackTrackerAlg.TruthMatching(event, cluster);

                if(fPrint)
                  std::cout << "Cluster of " << cluster->NHits() << " hits at t1 = " << cluster->Ts1()
                            << " (" << cluster->Ts1() - G4RefTime << ")"
                            << "\t Matches to trackID: " << truthMatch.trackid
                            << " with completeness: " << truthMatch.completeness
                            << " and purity: " << truthMatch.purity << std::endl;

                colour += fClusterColourInterval;
              }
            
            if(fDrawSpacePoints)
              {
                if(spacePointVec.size() == 1)
                  {
                    const art::Ptr<CRTSpacePoint> spacepoint = spacePointVec[0];
                    const geo::Point_t pos = spacepoint->Pos();
                    const geo::Point_t err = spacepoint->Err();

                    double rmin[3] = {pos.X() - err.X(), pos.Y() - err.Y(), pos.Z() - err.Z()};
                    double rmax[3] = {pos.X() + err.X(), pos.Y() + err.Y(), pos.Z() + err.Z()};

                    DrawCube(c1, rmin, rmax, fSpacePointColour);

                    if(fPrint)
                      std::cout << "Space Point: (" 
                                << rmin[0] << ", " << rmin[1] << ", " << rmin[2] << ") --> ("
                                << rmax[0] << ", " << rmax[1] << ", " << rmax[2] << ") at t1 = "
                                << spacepoint->Time() << " (" << spacepoint->Time() - G4RefTime << ")"
                                << " with PE " << spacepoint->PE()
                                << std::endl;
                  }
                else if(spacePointVec.size() != 0)
                  std::cout << "What an earth is going on here then..." << std::endl;
              }
          }
      }

    if(fDrawTracks)
      {
        auto tracksHandle = event.getValidHandle<std::vector<CRTTrack>>(fTrackLabel);
        std::vector<art::Ptr<CRTTrack>> tracksVec;
        art::fill_ptr_vector(tracksVec, tracksHandle);

        for(auto track : tracksVec)
          {
            if(track->Time() - G4RefTime < fMinTime || track->Time() - G4RefTime > fMaxTime)
              continue;

            const geo::Point_t start = track->Start();
            const geo::Vector_t dir  = track->Direction();

            TPolyLine3D *line = new TPolyLine3D(2);
            geo::Point_t a {0,0,0};
            geo::Point_t b {0,0,0};

            int i = 0;
            do
              {
                a = start + i * dir;
                ++i;
              }
            while(IsPointInsideBox(crtLims, a));

            i = 0;
            do
              {
                b = start + i * dir;
                --i;
              }
            while(IsPointInsideBox(crtLims, b));

            line->SetPoint(0, a.X(), a.Y(), a.Z());
            line->SetPoint(1, b.X(), b.Y(), b.Z());

            line->SetLineColor(fTrackColour);
            line->SetLineWidth(fLineWidth);
            line->Draw();

            if(fPrint)
              std::cout << "Track at (" << start.X() << ", " << start.Y() << ", " << start.Z() << ")\n"
                        << "\twith direction (" << dir.X() << ", " << dir.Y() << ", " << dir.Z() << ")\n"
                        << "\tdrawn between (" << a.X() << ", " << a.Y() << ", " << a.Z() << ")\n"
                        << "\tand (" << b.X() << ", " << b.Y() << ", " << b.Z() << ")\n"
                        << "\tat time " << track->Time() << " (" << track->Time() - G4RefTime << ")\n"
                        << "\tfrom three hits? " << track->Triple() << std::endl;

          }
      }

    c1->SaveAs(Form("crtEventDisplayEvent%d.root", event.event()));
    delete c1;
  }

  bool CRTEventDisplayAlg::IsPointInsideBox(const std::vector<double> &lims, const geo::Point_t &p)
  {
    return (p.X() > lims[0] && p.X() < lims[3])
      && (p.Y() > lims[1] && p.Y() < lims[4])
      && (p.Z() > lims[2] && p.Z() < lims[5]);
  }
}

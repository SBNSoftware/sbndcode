#include "CRTEventDisplayAlg.h"

namespace sbnd{
  
  CRTEventDisplayAlg::CRTEventDisplayAlg(const Config& config)
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

    fDrawTaggers = config.DrawTaggers();
    fDrawModules = config.DrawModules();
    fDrawStrips = config.DrawStrips();
    fDrawTpc = config.DrawTpc();
    fDrawTrueTracks = config.DrawTrueTracks();
    fDrawSimDeposits = config.DrawSimDeposits();
    fDrawStripHits = config.DrawStripHits();

    fTaggerColour = config.TaggerColour();
    fTpcColour = config.TpcColour();
    fTrueTrackColour = config.TrueTrackColour();
    fSimDepositColour = config.SimDepositColour();
    fStripHitColour = config.StripHitColour();

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
    // Create a canvas 
    TCanvas *c1 = new TCanvas("c1","",700,700);
    
    // Draw the CRT taggers
    if(fDrawTaggers)
      {
        for(auto const &[name, tagger] : fCrtGeo.GetTaggers())
          {
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
        for(auto const &[name, module] : fCrtGeo.GetModules())
          {
            double rmin[3] = {module.minX, 
                              module.minY, 
                              module.minZ};
            double rmax[3] = {module.maxX, 
                              module.maxY, 
                              module.maxZ};
            DrawCube(c1, rmin, rmax, fTaggerColour);
          }
      }
    
    // Draw individual CRT strips
    if(fDrawStrips)
      {
        for(auto const &[name, strip] : fCrtGeo.GetStrips())
          {
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
    
    // Draw true track trajectories for visible particles that cross the CRT
    if(fDrawTrueTracks)
      { 
        if(fPrint) std::cout<<"\nTrue tracks in event:\n";
        
        auto particleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fSimLabel);
        std::vector<double> crtLims = fCrtGeo.CRTLimits();
	crtLims[0] -= 100; crtLims[2] -= 100; crtLims[4] -= 100;
	crtLims[1] += 100; crtLims[3] -= 100; crtLims[5] -= 100;

        for(auto const& part : *particleHandle)
          {
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
            line->SetLineWidth(fLineWidth);
            line->Draw();
          
            if(fPrint) std::cout<<"->True ID: "<<part.TrackId()<<", start = ("<<start.X()<<", "<<start.Y()<<", "
                                <<start.Z()<<"), end = ("<<end.X()<<", "<<end.Y()<<", "<<end.Z()<<")\n";
          }
      }

    if(fDrawSimDeposits)
      {
	auto simDepositsHandle = event.getValidHandle<std::vector<sim::AuxDetSimChannel>>(fSimDepositLabel);

	for(auto const simDep : *simDepositsHandle)
	  {
	    for(auto const ide : simDep.AuxDetIDEs())
	      {
 		double x = (ide.entryX + ide.exitX) / 2.;
		double y = (ide.entryY + ide.exitY) / 2.;
		double z = (ide.entryZ + ide.exitZ) / 2.;

 		double ex = std::abs(ide.entryX - ide.exitX) / 2.;
		double ey = std::abs(ide.entryY - ide.exitY) / 2.;
		double ez = std::abs(ide.entryZ - ide.exitZ) / 2.;

		double rmin[3] = { x - ex, y - ey, z - ez};
		double rmax[3] = { x + ex, y + ey, z + ez};

		if(fPrint)
		  std::cout << "Sim Energy Deposit: (" 
			    << x << ", " << y << ", " << z << ")" << std::endl;
		DrawCube(c1, rmin, rmax, fSimDepositColour);
	      }
	  }	
      }

    if(fDrawStripHits)
      {
	auto stripHitsHandle = event.getValidHandle<std::vector<sbnd::crt::CRTStripHit>>(fStripHitLabel);

	for(auto const stripHit : *stripHitsHandle)
	  {
	    TVector3 xyz  = stripHit.XYZ();
	    TVector3 exyz = stripHit.XYZ_Error();

	    double rmin[3] = {xyz.X() - exyz.X(), xyz.Y() - exyz.Y(), xyz.Z() - exyz.Z()};
	    double rmax[3] = {xyz.X() + exyz.X(), xyz.Y() + exyz.Y(), xyz.Z() + exyz.Z()};

	    if(fPrint)
	      std::cout << "Strip Hit: (" 
			<< rmin[0] << ", " << rmin[1] << ", " << rmin[2] << ") --> ("
			<< rmax[0] << ", " << rmax[1] << ", " << rmax[2] << ")" << std::endl;
    
	    DrawCube(c1, rmin, rmax, fStripHitColour);
	  }
      }

    c1->SaveAs(Form("crtEventDisplayEvent%d.root", event.event()));
    delete c1;
  }
}

#include "TPCEventDisplayAlg.h"

namespace sbnd::crt {
  
  TPCEventDisplayAlg::TPCEventDisplayAlg(const Config& config)
  {
    this->reconfigure(config);
  }
  
  TPCEventDisplayAlg::TPCEventDisplayAlg(){}
  
  TPCEventDisplayAlg::~TPCEventDisplayAlg(){}

  void TPCEventDisplayAlg::reconfigure(const Config& config)
  {
    fSliceLabel = config.SliceLabel();

    fDrawTPC = config.DrawTPC();
    fDrawSlices = config.DrawSlices();

    fTPCColour = config.TPCColour();
    fSliceStartingColour = config.SliceStartingColour();
    fSliceColourInterval = config.SliceColourInterval();

    fLineWidth = config.LineWidth();

    return;
  }
 
  void TPCEventDisplayAlg::DrawCube(TCanvas *c1, double *rmin, double *rmax, int colour, int lineWidth)
  {
    c1->cd();
    TList *outline = new TList;
    TPolyLine3D *p1 = new TPolyLine3D(4);
    TPolyLine3D *p2 = new TPolyLine3D(4);
    TPolyLine3D *p3 = new TPolyLine3D(4);
    TPolyLine3D *p4 = new TPolyLine3D(4);
    p1->SetLineColor(colour);

    if(lineWidth == -1)
      p1->SetLineWidth(fLineWidth);
    else
      p1->SetLineWidth(lineWidth);

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

  void TPCEventDisplayAlg::Draw(detinfo::DetectorClocksData const& clockData,
                                const art::Event& event, const TString& saveName)
  {
    // Create a canvas 
    TCanvas *c1 = new TCanvas(Form("c%s",saveName.Data()),"",700,700);
    
    // Draw the TPC with central cathode
    if(fDrawTPC)
      {
        double rmin[3] = {fTPCGeoAlg.MinX(), 
                          fTPCGeoAlg.MinY(), 
                          fTPCGeoAlg.MinZ()};
        double rmax[3] = {-fTPCGeoAlg.CpaWidth(), 
                          fTPCGeoAlg.MaxY(), 
                          fTPCGeoAlg.MaxZ()};
        DrawCube(c1, rmin, rmax, fTPCColour);
        double rmin2[3] = {fTPCGeoAlg.CpaWidth(), 
                           fTPCGeoAlg.MinY(), 
                           fTPCGeoAlg.MinZ()};
        double rmax2[3] = {fTPCGeoAlg.MaxX(), 
                           fTPCGeoAlg.MaxY(), 
                           fTPCGeoAlg.MaxZ()};

        DrawCube(c1, rmin2, rmax2, fTPCColour);
      }

    if(fDrawSlices)
      {
        art::ValidHandle<std::vector<recob::Slice>> sliceHandle = event.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
        art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle = event.getValidHandle<std::vector<recob::PFParticle>>(fSliceLabel);

        std::vector<art::Ptr<recob::Slice>> sliceVec;
        art::fill_ptr_vector(sliceVec, sliceHandle);

        art::FindManyP<recob::PFParticle> slicesToPFParticles(sliceHandle, event, fSliceLabel);
        art::FindManyP<recob::SpacePoint> pfpsToSpacePoints(pfpHandle, event, fSliceLabel);

        int colour = fSliceStartingColour;

        for(auto const& slc : sliceVec)
          {
            std::vector<art::Ptr<recob::PFParticle>> pfpVec = slicesToPFParticles.at(slc.key());

            for(auto const& pfp : pfpVec)
              {
                std::vector<art::Ptr<recob::SpacePoint>> spacePointVec = pfpsToSpacePoints.at(pfp.key());

                for(auto const& sp : spacePointVec)
                  {
                    const double* pos = sp->XYZ();

                    double rmin[3] = { pos[0] - 0.15,
                                       pos[1] - 0.15,
                                       pos[2] - 0.15
                    };

                    double rmax[3] = { pos[0] + 0.15,
                                       pos[1] + 0.15,
                                       pos[2] + 0.15
                    };

                    DrawCube(c1, rmin, rmax, colour);
                  }
              }

            colour += fSliceColourInterval;
          }
      }

    c1->SaveAs(Form("%s.root", saveName.Data()));


    TView3D *view = (TView3D*) TView::CreateView(1);

    double c[3] = { 0, 0, 250 };
    double s[3] = { 350, -350, 350 };

    view->SetRange(-300, -300, -50, 300, 300, 550);

    view->DefineViewDirection(s, c,
                              .2, .8,
                              .5, .5,
                              1, 0,
                              view->GetTnorm(),
                              view->GetTback());

    c1->SaveAs(Form("%s.pdf", saveName.Data()));
    c1->SaveAs(Form("%s.png", saveName.Data()));

    view->DefineViewDirection(s, c,
                              0, 1,
                              1, 0,
                              1, 0,
                              view->GetTnorm(),
                              view->GetTback());

    c1->SaveAs(Form("%s_front.png", saveName.Data()));
    c1->SaveAs(Form("%s_front.pdf", saveName.Data()));

    view->DefineViewDirection(s, c,
                              0, 1,
                              0, 1,
                              1, 0,
                              view->GetTnorm(),
                              view->GetTback());

    c1->SaveAs(Form("%s_top.png", saveName.Data()));
    c1->SaveAs(Form("%s_top.pdf", saveName.Data()));

    view->DefineViewDirection(s, c,
                              1, 0,
                              0, 1,
                              0, 1,
                              view->GetTnorm(),
                              view->GetTback());

    c1->SaveAs(Form("%s_side.png", saveName.Data()));
    c1->SaveAs(Form("%s_side.pdf", saveName.Data()));
  }
}

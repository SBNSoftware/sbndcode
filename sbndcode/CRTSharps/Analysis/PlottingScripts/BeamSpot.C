void BeamSpot()
{
  const TString saveDir = "/sbnd/data/users/hlay/crt/sharps/fancy_plots";
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  const TString run_name = "run4500";
  gSystem->Exec("mkdir -p " + saveDir + "/" + run_name);

  TChain *tree = new TChain("crtana/tree");
  tree->Add("/pnfs/sbnd/scratch/users/hlay/crt_sharps_data/" + run_name + "/crtana_sbnd.root");

  TCanvas *cBeamSpot = new TCanvas("cBeamSpot", "cBeamSpot");
  cBeamSpot->cd();

  gStyle->SetNdivisions(505, "x");
  gStyle->SetPalette(kBlueRedYellow);

  TH2F* sighist  = new TH2F("sighist", ";Hit -x position (cm);Hit y position (cm);CRT 3D Hits",
			    15, -50, 200, 15, -120, 130);
  TH2F* backhist  = new TH2F("backhist", ";Hit -x position (cm);Hit y position (cm);CRT 3D Hits",
			     15, -50, 200, 15, -120, 130);

  cBeamSpot->SetRightMargin(0.25);

  tree->Draw("chit_y:-chit_x>>sighist", "chit_z < 0 && chit_t1 > 332900 && chit_t1 < 334400", "colz");
  tree->Draw("chit_y:-chit_x>>backhist", "chit_z < 0 && chit_t1 > 330000 && chit_t1 < 331500", "colz");

  sighist->GetYaxis()->SetTitleOffset(1.25);
  sighist->GetZaxis()->SetTitleOffset(1.3);

  TH2F* subhist = (TH2F*) sighist->Clone();
  subhist->Add(backhist, -1.);

  subhist->Draw("colz");

  if(save)
    {
      cBeamSpot->SaveAs(saveDir + "/crt_sharps_" + run_name + "_beam_spot.png");
      cBeamSpot->SaveAs(saveDir + "/crt_sharps_" + run_name + "_beam_spot.pdf");
    }
}

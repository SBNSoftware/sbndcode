void TopHat()
{
  const TString saveDir = "/sbnd/data/users/hlay/crt/plots_for_michelle";
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  const TString run_name = "run2100";

  TChain *tree = new TChain("crtana/tree");
  tree->Add("/pnfs/sbnd/scratch/users/hlay/crt_sharps_data/" + run_name + "/crtana_sbnd.root");

  TCanvas *cT1BeamAcceptance = new TCanvas("cT1BeamAcceptance", "cT1BeamAcceptance");
  cT1BeamAcceptance->cd();
  
  TH1F* hT1BeamAcceptance  = new TH1F("hT1BeamAcceptance", ";CRT - BES (#kern[.1]{#mus});CRT 3D Hits", 100, 319, 344);
  hT1BeamAcceptance->SetLineColor(kMagenta+2);
  hT1BeamAcceptance->SetLineWidth(3);
  hT1BeamAcceptance->GetYaxis()->SetTitleOffset(1);
  hT1BeamAcceptance->SetMinimum(0);

  tree->Draw("chit_t1 * 1e-3 >> hT1BeamAcceptance", "", "hist");

  TPaveText *pt = new TPaveText(.28,.85,.48,.9, "NDC");
  pt->AddText("SBND CRT## Run2100 Data");
  pt->SetFillColor(kWhite);
  pt->SetTextColor(kGray+2);
  pt->SetTextSize(0.035);
  pt->Draw("same");

  if(save)
    {
      cT1BeamAcceptance->SaveAs(saveDir + "/crt_sharps_" + run_name + "_beam_peak.png");
      cT1BeamAcceptance->SaveAs(saveDir + "/crt_sharps_" + run_name + "_beam_peak.pdf");
      cT1BeamAcceptance->SaveAs(saveDir + "/crt_sharps_" + run_name + "_beam_peak.C");
    }
}

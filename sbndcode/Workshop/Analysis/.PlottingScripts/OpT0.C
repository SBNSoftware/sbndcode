void OpT0(const std::string fileName)
{
  gStyle->SetOptStat(0);

  TFile *file = new TFile(fileName.c_str(), "READ");

  TTree* tree = (TTree*)file->Get("ana/tree");

  TH1D *hOpT0 = new TH1D("hOpT0", ";T0 (#mu s); Slices", 50, 1.602, 1.608);

  TCanvas *c = new TCanvas("c", "c", 1400, 1000);
  c->cd();

  tree->Draw("opT0>>hOpT0");

  hOpT0->SetLineColor(kGreen+3);
  hOpT0->SetLineWidth(2);

  hOpT0->Draw("hist");
}

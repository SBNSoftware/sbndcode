void TrackLength(const std::string fileName)
{
  gStyle->SetOptStat(0);

  TFile *file = new TFile(fileName.c_str(), "READ");

  TTree* tree = (TTree*)file->Get("ana/tree");

  TH1D *hTLLong = new TH1D("hTLLong", ";Track Length (cm); Tracks", 70, 0, 350);
  TH1D *hTLNotLong = new TH1D("hTLNotLong", ";Track Length (cm); Tracks", 70, 0, 350);

  TCanvas *c = new TCanvas("c", "c", 1400, 1000);
  c->cd();

  tree->Draw("childTrackLengths>>hTLLong", "childTrackIsLongest");
  tree->Draw("childTrackLengths>>hTLNotLong", "!childTrackIsLongest");

  hTLLong->SetLineColor(kMagenta+2);
  hTLNotLong->SetLineColor(kOrange+9);

  hTLLong->SetLineWidth(2);
  hTLNotLong->SetLineWidth(2);

  hTLNotLong->Draw("hist");
  hTLLong->Draw("histsame");

  TLegend *leg = new TLegend(.54, .63, .83, .78);
  leg->AddEntry("hTLLong", "Longest Track", "l");
  leg->AddEntry("hTLNotLong", "Other Tracks", "l");
  leg->SetLineWidth(0);
  leg->Draw();
}

void dEdxResRange(const std::string fileName)
{
  gStyle->SetOptStat(0);

  TFile *file = new TFile(fileName.c_str(), "READ");

  TTree* tree = (TTree*)file->Get("ana/tree");

  TH2D *hdEdxRRLong = new TH2D("hdEdxRRLong", ";Residual Range (cm); dE/dx (MeV/cm)", 100, 0, 50, 100, 0, 30);
  TH2D *hdEdxRRNotLong = new TH2D("hdEdxRRNotLong", ";Residual Range (cm); dE/dx (MeV/cm)", 100, 0, 50, 100, 0, 30);

  TCanvas *c = new TCanvas("c", "c", 1400, 1000);
  c->cd();

  tree->Draw("childTrackdEdx:childTrackResRange>>hdEdxRRLong", "childTrackIsLongest");
  tree->Draw("childTrackdEdx:childTrackResRange>>hdEdxRRNotLong", "!childTrackIsLongest");

  hdEdxRRLong->SetMarkerColor(kMagenta+2);
  hdEdxRRNotLong->SetMarkerColor(kOrange+9);

  hdEdxRRLong->SetMarkerStyle(8);
  hdEdxRRNotLong->SetMarkerStyle(8);

  hdEdxRRLong->SetMarkerSize(.5);
  hdEdxRRNotLong->SetMarkerSize(.5);

  hdEdxRRNotLong->Draw("hist");
  hdEdxRRLong->Draw("histsame");

  TH2D *hdEdxRRLongLeg = (TH2D*)hdEdxRRLong->Clone("hdEdxRRLongLeg");
  TH2D *hdEdxRRNotLongLeg = (TH2D*)hdEdxRRNotLong->Clone("hdEdxRRNotLongLeg");
  hdEdxRRLongLeg->SetMarkerSize(1.5);
  hdEdxRRNotLongLeg->SetMarkerSize(1.5);

  TLegend *leg = new TLegend(.54, .63, .83, .78);
  leg->AddEntry(hdEdxRRLongLeg, "Longest Track", "p");
  leg->AddEntry(hdEdxRRNotLongLeg, "Other Tracks", "p");
  leg->SetLineWidth(0);
  leg->Draw();
}

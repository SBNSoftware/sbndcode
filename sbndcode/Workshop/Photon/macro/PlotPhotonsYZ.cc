void PlotPhotonsYZ(TString file_name, int event)
{
  TFile * file = new TFile(file_name,"READ");
  TTree * tree = (TTree *)file->Get("opanalyzer/PhotonsPerOpDet");

  // Cuts to select only PMTs and particular event
  TCut c1 = "isPMT";
  TCut c2 = Form("EventID == %d", event);
  
  //2D Histogram with SBND PDS dimensions
  TH2F * hPhotocatode = new TH2F("hPhotocathode", ";Z [cm];Y [cm]; #PE", 50, 0, 500, 40, -200, 200);
  hPhotocatode->SetStats(0);
  tree->Draw("OpDetY:OpDetZ >> hPhotocathode", (c1&&c2)*"CountAll", "COLZ");
}

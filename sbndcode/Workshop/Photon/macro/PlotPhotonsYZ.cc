void PlotPhotonsYZ(TString file_name, int event)
{
  TFile * file = new TFile(file_name,"READ");
  TTree * tree = (TTree *)file->Get("opanalyzer/PhotonsPerOpDet");

  // Cut to select event and only PMTs
  TCut c1 = Form("EventID == %d && isPMT == 1", event);

  //2D Histogram with SBND PDS dimensions
  TH2F * hPhotocatode = new TH2F("hPhotocathode", ";Z [cm];Y [cm]; #PE", 50, 0, 500, 40, -200, 200);
  hPhotocatode->SetStats(0);
  tree->Draw("OpDetY:OpDetZ >> hPhotocathode", (c1)*"CountAll", "COLZ");
}

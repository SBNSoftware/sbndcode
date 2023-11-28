void InterPlaneDriftBugPlots()
{
  const TString saveDir = "/exp/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv10/interplanedriftbugplots";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFileOld = "/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv8/NCPiZeroAv8_rockbox.root";
  const TString rockboxFileNew = "/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv9/NCPiZeroAv9_rockbox.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxEventsOld = new TChain("ncpizeroana/events");
  rockboxEventsOld->Add(rockboxFileOld);

  TChain *rockboxEventsNew = new TChain("ncpizeroana/events");
  rockboxEventsNew->Add(rockboxFileNew);

  TCanvas *cUnclusteredHitsDiff = new TCanvas("cUnclusteredHitsDiff", "cUnclusteredHitsDiff");
  cUnclusteredHitsDiff->cd();
  cUnclusteredHitsDiff->SetLogy();

  TH1F *hUnclusteredHitsOld = new TH1F("hUnclusteredHitsOld", ";Unclusted Hits;Slices (normalised)", 50, 0, 4000);
  rockboxEventsOld->Draw("slc_n_hits-slc_n_used_hits>>hUnclusteredHitsOld");

  TH1F *hUnclusteredHitsNew = new TH1F("hUnclusteredHitsNew", ";Unclusted Hits;Slices (normalised)", 50, 0, 4000);
  rockboxEventsNew->Draw("slc_n_hits-slc_n_used_hits>>hUnclusteredHitsNew");

  hUnclusteredHitsOld->SetLineColor(kRed+2);
  hUnclusteredHitsOld->GetXaxis()->SetNdivisions(505);
  hUnclusteredHitsNew->SetLineColor(kBlue+2);

  hUnclusteredHitsOld->DrawNormalized("hist");
  hUnclusteredHitsNew->DrawNormalized("histsame");

  cUnclusteredHitsDiff->SaveAs(saveDir + "/unclustered_hits.png");
  cUnclusteredHitsDiff->SaveAs(saveDir + "/unclustered_hits.pdf");
}

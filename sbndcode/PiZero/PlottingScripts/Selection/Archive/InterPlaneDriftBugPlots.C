void InterPlaneDriftBugPlots()
{
  const TString saveDir = "/exp/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv10/interplanedriftbugplots";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFileOld = "/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv8/NCPiZeroAv8_rockbox.root";
  const TString rockboxFileNew = "/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv10/NCPiZeroAv10_rockbox.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxEventsOld = new TChain("ncpizeroana/events");
  rockboxEventsOld->Add(rockboxFileOld);

  TChain *rockboxEventsNew = new TChain("ncpizeroana/events");
  rockboxEventsNew->Add(rockboxFileNew);

  TLegend *legend = new TLegend(.7, .7, .9, .8);
  legend->SetBorderSize(0);

  TCanvas *cUnclusteredHits = new TCanvas("cUnclusteredHits", "cUnclusteredHits");
  cUnclusteredHits->cd();
  cUnclusteredHits->SetLogy();

  TH1F *hUnclusteredHitsOld = new TH1F("hUnclusteredHitsOld", ";Unclustered Hits;Slices (normalised)", 50, 0, 4000);
  rockboxEventsOld->Draw("slc_n_hits-slc_n_used_hits>>hUnclusteredHitsOld");

  TH1F *hUnclusteredHitsNew = new TH1F("hUnclusteredHitsNew", ";Unclustered Hits;Slices (normalised)", 50, 0, 4000);
  rockboxEventsNew->Draw("slc_n_hits-slc_n_used_hits>>hUnclusteredHitsNew");

  hUnclusteredHitsOld->SetLineColor(kRed+2);
  hUnclusteredHitsOld->GetXaxis()->SetNdivisions(505);
  hUnclusteredHitsNew->SetLineColor(kBlue+2);

  legend->AddEntry(hUnclusteredHitsOld, "Unfixed", "l");
  legend->AddEntry(hUnclusteredHitsNew, "Fixed", "l");

  hUnclusteredHitsOld->DrawNormalized("hist");
  hUnclusteredHitsNew->DrawNormalized("histsame");

  legend->Draw();

  cUnclusteredHits->SaveAs(saveDir + "/unclustered_hits.png");
  cUnclusteredHits->SaveAs(saveDir + "/unclustered_hits.pdf");

  TCanvas *cPFPComp = new TCanvas("cPFPComp", "cPFPComp");
  cPFPComp->cd();

  TH1F *hPFPCompOld = new TH1F("hPFPCompOld", ";Completeness;PFPs (normalised)", 40, 0, 1);
  rockboxEventsOld->Draw("slc_pfp_comp>>hPFPCompOld", "!slc_is_clear_cosmic && slc_pfp_comp>0 && slc_pfp_primary_child");

  TH1F *hPFPCompNew = new TH1F("hPFPCompNew", ";Completeness;PFPs (normalised)", 40, 0, 1);
  rockboxEventsNew->Draw("slc_pfp_comp>>hPFPCompNew", "!slc_is_clear_cosmic && slc_pfp_comp>0 && slc_pfp_primary_child");

  hPFPCompOld->SetLineColor(kRed+2);
  hPFPCompOld->GetXaxis()->SetNdivisions(505);
  hPFPCompNew->SetLineColor(kBlue+2);

  hPFPCompOld->SetMaximum(1.2 * hPFPCompOld->GetMaximum());

  hPFPCompOld->DrawNormalized("hist][");
  hPFPCompNew->DrawNormalized("hist][same");

  legend->Draw();

  cPFPComp->SaveAs(saveDir + "/pfp_comp.png");
  cPFPComp->SaveAs(saveDir + "/pfp_comp.pdf");

  TCanvas *cPFPPur = new TCanvas("cPFPPur", "cPFPPur");
  cPFPPur->cd();

  TH1F *hPFPPurOld = new TH1F("hPFPPurOld", ";Purity;PFPs (normalised)", 40, 0, 1);
  rockboxEventsOld->Draw("slc_pfp_pur>>hPFPPurOld", "!slc_is_clear_cosmic && slc_pfp_comp>0 && slc_pfp_primary_child");

  TH1F *hPFPPurNew = new TH1F("hPFPPurNew", ";Purity;PFPs (normalised)", 40, 0, 1);
  rockboxEventsNew->Draw("slc_pfp_pur>>hPFPPurNew", "!slc_is_clear_cosmic && slc_pfp_comp>0 && slc_pfp_primary_child");

  hPFPPurOld->SetLineColor(kRed+2);
  hPFPPurOld->GetXaxis()->SetNdivisions(505);
  hPFPPurNew->SetLineColor(kBlue+2);

  hPFPPurOld->DrawNormalized("hist][");
  hPFPPurNew->DrawNormalized("hist][same");

  legend->Draw();

  cPFPPur->SaveAs(saveDir + "/pfp_pur.png");
  cPFPPur->SaveAs(saveDir + "/pfp_pur.pdf");
}

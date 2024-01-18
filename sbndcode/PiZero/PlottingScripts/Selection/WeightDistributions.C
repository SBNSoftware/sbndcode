#include "Selections.h"
#include "Common.C"
#include "WeightNames.h"
#include "Enumerate.h"

void WeightDistributions(const TString productionVersion, const TString saveDirExt, const std::vector<Cut> &signals,
                         std::vector<std::string> weight_names)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/weight_distributions/" + saveDirExt;
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_rockbox.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxEvents = new TChain("ncpizeroana/events");
  rockboxEvents->Add(rockboxFile);

  for(auto&& [weight_i, weight] : enumerate(weight_names))
    {
      TCanvas *canvas = new TCanvas(Form("canvas%lu", weight_i), Form("canvas%lu", weight_i));
      canvas->cd();
      canvas->SetTopMargin(0.12);

      rockboxEvents->Draw(Form("nu_weight_%s", weight.c_str()));
      TH1D *hist = new TH1D(*((TH1D*) gPad->GetPrimitive("htemp")));
      hist->SetName(Form("hist%lu", weight_i));
      hist->SetTitle(Form("%s;Weight;#nus x Univs (A.U.)", weight.c_str()));
      hist->GetYaxis()->SetTitleOffset(1.2);
      hist->SetLineColor(kBlue+2);
      hist->DrawNormalized("hist");

      canvas->SaveAs(saveDir + "/" + weight.c_str() + ".png");
      canvas->SaveAs(saveDir + "/" + weight.c_str() + ".pdf");

      TCanvas *signalCanvas = new TCanvas(Form("signalCanvas%lu", weight_i), Form("signalCanvas%lu", weight_i));
      signalCanvas->cd();
      signalCanvas->SetTopMargin(0.12);

      const int nbins    = hist->GetNbinsX();
      const double xlow  = hist->GetBinLowEdge(1);
      const double xhigh = hist->GetBinLowEdge(nbins) + hist->GetBinWidth(nbins);

      TH1D *histIncl = new TH1D(Form("histIncl%lu", weight_i), Form("%s;Weight;#nus x Univs (A.U.)", weight.c_str()), nbins, xlow, xhigh);
      rockboxEvents->Draw(Form("nu_weight_%s>>histIncl%lu", weight.c_str(), weight_i), "nu_event_type_incl==0");
      histIncl->SetLineColor(kMagenta+2);
 
      TH1D *hist0p0pi = new TH1D(Form("hist0p0pi%lu", weight_i), Form("%s;Weight;#nus x Univs (A.U.)", weight.c_str()), nbins, xlow, xhigh);
      rockboxEvents->Draw(Form("nu_weight_%s>>hist0p0pi%lu", weight.c_str(), weight_i), "nu_event_type_0p0pi==0");
      hist0p0pi->SetLineColor(kRed+2);

      TH1D *histNp0pi = new TH1D(Form("histNp0pi%lu", weight_i), Form("%s;Weight;#nus x Univs (A.U.)", weight.c_str()), nbins, xlow, xhigh);
      rockboxEvents->Draw(Form("nu_weight_%s>>histNp0pi%lu", weight.c_str(), weight_i), "nu_event_type_Np0pi==0");
      histNp0pi->SetLineColor(kGreen+2);

      const double max = std::max({hist->GetMaximum() / hist->GetEntries(), histIncl->GetMaximum() / histIncl->GetEntries(),
            hist0p0pi->GetMaximum() / hist0p0pi->GetEntries(), histNp0pi->GetMaximum() / histNp0pi->GetEntries() });

      hist->Scale(max / hist->GetMaximum());
      hist->SetMaximum(1.2 * max);
      histIncl->Scale(max/(histIncl->GetMaximum()));
      hist0p0pi->Scale(max/(hist0p0pi->GetMaximum()));
      histNp0pi->Scale(max/(histNp0pi->GetMaximum()));
      
      hist->Draw("hist");
      histIncl->Draw("histsame");
      hist0p0pi->Draw("histsame");
      histNp0pi->Draw("histsame");

      TLegend *leg = new TLegend(0.25, 0.85, 0.94, 0.8);
      leg->SetNColumns(4);
      leg->AddEntry(hist, "All", "l");
      leg->AddEntry(histIncl, "NC#pi^{0}", "l");
      leg->AddEntry(hist0p0pi, "NC#pi^{0}0p0#pi^{#pm}", "l");
      leg->AddEntry(histNp0pi, "NC#pi^{0}1p0#pi^{#pm}", "l");

      leg->Draw();

      signalCanvas->SaveAs(saveDir + "/" + weight.c_str() + "_signal.png");
      signalCanvas->SaveAs(saveDir + "/" + weight.c_str() + "_signal.pdf");
    }
}

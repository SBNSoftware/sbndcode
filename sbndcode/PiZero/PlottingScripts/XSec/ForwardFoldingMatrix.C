#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "XSecPlot.h"
#include "Selection.h"
#include "Common.C"

void ForwardFoldingMatrix(const TString productionVersion)
{
  gROOT->SetStyle("henrySBND");
  gStyle->SetPaintTextFormat("1.2g");
  gROOT->ForceStyle();

  const TString saveDir = baseSaveDir + "/" + productionVersion + "/forwardfoldingmatrices";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString ncpizeroFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_ncpizero.root";

  TChain *ncpizeroSlices = new TChain("ncpizeroxsectrees/slices");
  ncpizeroSlices->Add(ncpizeroFile);

  std::deque<bool> sel(selections.size(), false);
  std::vector<int> event_type(selections.size(), -1);
  double pizero_mom, cos_theta_pizero, reco_pizero_mom, reco_cos_theta_pizero;

  ncpizeroSlices->SetBranchStatus("*", 0);

  ncpizeroSlices->SetBranchAddress("pizero_mom", &pizero_mom);
  ncpizeroSlices->SetBranchAddress("cos_theta_pizero", &cos_theta_pizero);
  ncpizeroSlices->SetBranchAddress("reco_pizero_mom", &reco_pizero_mom);
  ncpizeroSlices->SetBranchAddress("reco_cos_theta_pizero", &reco_cos_theta_pizero);

  for(auto&& [ selection_i, selection ] : enumerate(selections))
    {
      ncpizeroSlices->SetBranchAddress(selection.signal, &event_type[selection_i]);
      ncpizeroSlices->SetBranchAddress(selection.cut, &sel[selection_i]);
    }

  TFile *outfile = new TFile(saveDir + "/forwardfoldingmatrices.root", "RECREATE");

  const int N = ncpizeroSlices->GetEntries();

  const double piZeroMomBins[9]       = { 0., 60., 120., 180., 240., 300., 400., 600., 1000. };
  const double cosThetaPiZeroBins[10] = { -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. };

  TH1D *hPiZeroMom      = new TH1D("hPiZeroMom", "", 8, piZeroMomBins);
  TH1D *hCosThetaPiZero = new TH1D("hCosThetaPiZero", "", 9, cosThetaPiZeroBins);

  for(auto&& [ selection_i, selection ] : enumerate(selections))
    {
      TH2D *hForwardFoldPiZeroMom      = new TH2D(Form("hForwardFoldPiZeroMom%lu", selection_i), ";True Bin;Reco Bin", 10, -0.5, 9.5, 10, -0.5, 9.5);
      TH2D *hForwardFoldCosThetaPiZero = new TH2D(Form("hForwardFoldCosThetaPiZero%lu", selection_i), ";True Bin;Reco Bin", 11, -0.5, 10.5, 11, -0.5, 10.5);

      for(int slc = 0; slc < N; ++slc)
        {
          ncpizeroSlices->GetEntry(slc);

          if(sel[selection_i] && event_type[selection_i] == 0)
            {
              const int trueMomBin = hPiZeroMom->FindBin(pizero_mom);
              const int recoMomBin = hPiZeroMom->FindBin(reco_pizero_mom);

              const int trueCosThetaBin = hCosThetaPiZero->FindBin(cos_theta_pizero);
              const int recoCosThetaBin = hCosThetaPiZero->FindBin(reco_cos_theta_pizero);

              hForwardFoldPiZeroMom->Fill(trueMomBin, recoMomBin);
              hForwardFoldCosThetaPiZero->Fill(trueCosThetaBin, recoCosThetaBin);
            }
        }

      TCanvas *piZeroMomCanvas = new TCanvas(Form("piZeroMomCanvas%lu", selection_i), Form("piZeroMomCanvas%lu", selection_i));
      piZeroMomCanvas->cd();

      gPad->SetTopMargin(0.12);
      gPad->SetRightMargin(0.2);

      NormaliseEntriesByXTotal(hForwardFoldPiZeroMom);

      hForwardFoldPiZeroMom->GetXaxis()->SetBinLabel(1, "U/F");
      hForwardFoldPiZeroMom->GetXaxis()->SetBinLabel(10, "O/F");
      hForwardFoldPiZeroMom->GetYaxis()->SetBinLabel(1, "U/F");
      hForwardFoldPiZeroMom->GetYaxis()->SetBinLabel(10, "O/F");

      for(int i = 2; i < 10; ++i)
        {
          hForwardFoldPiZeroMom->GetXaxis()->SetBinLabel(i, std::to_string(i-1).c_str());
          hForwardFoldPiZeroMom->GetYaxis()->SetBinLabel(i, std::to_string(i-1).c_str());
        }

      hForwardFoldPiZeroMom->Draw("colztext");
      hForwardFoldPiZeroMom->Write("hForwardFoldPiZeroMom_" + selection.name);

      piZeroMomCanvas->SaveAs(saveDir + "/" + selection.name + "_pizero_mom.png");
      piZeroMomCanvas->SaveAs(saveDir + "/" + selection.name + "_pizero_mom.pdf");

      TCanvas *cosThetaPiZeroCanvas = new TCanvas(Form("cosThetaPiZeroCanvas%lu", selection_i), Form("cosThetaPiZeroCanvas%lu", selection_i));
      cosThetaPiZeroCanvas->cd();

      gPad->SetTopMargin(0.12);
      gPad->SetRightMargin(0.2);

      NormaliseEntriesByXTotal(hForwardFoldCosThetaPiZero);

      hForwardFoldCosThetaPiZero->GetXaxis()->SetBinLabel(1, "U/F");
      hForwardFoldCosThetaPiZero->GetXaxis()->SetBinLabel(11, "O/F");
      hForwardFoldCosThetaPiZero->GetYaxis()->SetBinLabel(1, "U/F");
      hForwardFoldCosThetaPiZero->GetYaxis()->SetBinLabel(11, "O/F");

      for(int i = 2; i < 11; ++i)
        {
          hForwardFoldCosThetaPiZero->GetXaxis()->SetBinLabel(i, std::to_string(i-1).c_str());
          hForwardFoldCosThetaPiZero->GetYaxis()->SetBinLabel(i, std::to_string(i-1).c_str());
        }

      hForwardFoldCosThetaPiZero->Draw("colztext");
      hForwardFoldCosThetaPiZero->Write("hForwardFoldCosThetaPiZero_" + selection.name);

      cosThetaPiZeroCanvas->SaveAs(saveDir + "/" + selection.name + "_cos_theta_pizero.png");
      cosThetaPiZeroCanvas->SaveAs(saveDir + "/" + selection.name + "_cos_theta_pizero.pdf");
    }
}

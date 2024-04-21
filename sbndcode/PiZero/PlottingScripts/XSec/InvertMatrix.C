#include "Common.C"

void InvertMatrix(const TString productionVersion)
{
  gROOT->SetStyle("henrySBND");
  gStyle->SetPaintTextFormat("1.2g");
  gStyle->SetPadTopMargin(0.12);
  gStyle->SetPadRightMargin(0.2);
  gROOT->ForceStyle();

  const TString saveDir = baseSaveDir + "/" + productionVersion + "/unfoldingmatrices";
  gSystem->Exec("mkdir -p " + saveDir);

  TFile* forwardFoldFile = new TFile("/exp/sbnd/data/users/hlay/ncpizero/plots/" + productionVersion + "/forwardfoldingmatrices/forwardfoldingmatrices.root", "READ");
  TFile *outfile = new TFile(saveDir + "/unfoldingmatrices.root", "RECREATE");

  for (auto&& keyAsObj : *forwardFoldFile->GetListOfKeys())
    {
      auto key = (TKey*) keyAsObj;
      TString name = key->GetName();

      TCanvas *canvas = new TCanvas("c" + name, "c" + name);
      canvas->cd();

      TH2D* hist = (TH2D*) forwardFoldFile->Get(name);
      name.ReplaceAll("ForwardFold", "Unfold");

      const int size   = strstr(key->GetName(), "cos_theta") ? hist->GetNbinsX() - 2 : hist->GetNbinsX() - 1;
      const int lim    = strstr(key->GetName(), "cos_theta") ? hist->GetNbinsX() : hist->GetNbinsX() + 1;
      const int start  = strstr(key->GetName(), "cos_theta") ? 2 : 2;

      TMatrixD forwardMatrix = THilbertMatrixD(size, size);

      for(int i = start; i < lim; ++i)
        {
          for(int j = start; j < lim; ++j)
            {
              forwardMatrix(i - 2, j - 2) = hist->GetBinContent(i, j);
            }
        }

      forwardMatrix.Invert();

      const int rawSize = hist->GetNbinsX();

      TH2D *invHist = new TH2D(name, "", rawSize, -0.5, rawSize - 0.5, rawSize, -0.5, rawSize - 0.5);

      for(int i = start; i < lim; ++i)
        {
          for(int j = start; j < lim; ++j)
            {
              invHist->SetBinContent(i, j, forwardMatrix(i - 2, j - 2));
              invHist->SetBinError(i, j, 0);
            }
        }

      invHist->Draw("colzetext");
      invHist->Write(name);

      canvas->SaveAs(saveDir + "/" + name + ".png");
      canvas->SaveAs(saveDir + "/" + name + ".pdf");
    }
}

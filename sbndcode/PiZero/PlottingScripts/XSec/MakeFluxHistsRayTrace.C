#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Common.C"

void MakeFluxHistsRayTrace(const TString productionVersion)
{
  TString saveDir = baseSaveDir + "/" + productionVersion + "/flux_hists_ray_trace";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString file = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_flux_configL_*.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();
  gStyle->SetPadTopMargin(.1);

  TChain *nus = new TChain("fluxana/tree");
  nus->Add(file);

  const int N = nus->GetEntries();

  int nu_pdg;
  float nu_e;
  double ray_fv_length; 

  nus->SetBranchStatus("*", 0);
  nus->SetBranchStatus("nu_pdg", 1);
  nus->SetBranchStatus("nu_e", 1);
  nus->SetBranchStatus("ray_fv_length", 1);

  nus->SetBranchAddress("nu_pdg", &nu_pdg);
  nus->SetBranchAddress("nu_e", &nu_e);
  nus->SetBranchAddress("ray_fv_length", &ray_fv_length);

  TH1F *hRayFVLength = new TH1F("hRayFVLength", ";#nu Ray Length in FV (cm);#nu", 100, -20, 460);

  TH1F *hNuEnergyNuMu  = new TH1F("hNuEnergyNuMu", ";True E_{#nu} (GeV);#nu_{#mu}", 200, 0, 10);
  TH1F *hNuEnergyANuMu = new TH1F("hNuEnergyANuMu", ";True E_{#nu} (GeV);#bar{#nu_{#mu}}", 200, 0, 10);
  TH1F *hNuEnergyNuE   = new TH1F("hNuEnergyNuE", ";True E_{#nu} (GeV);#nu_{e}", 200, 0, 10);
  TH1F *hNuEnergyANuE  = new TH1F("hNuEnergyANuE", ";True E_{#nu} (GeV);#bar{#nu_{e}}", 200, 0, 10);

  double mean_ray_length = 0.;
  int N_considered = 0;

  for(int i = 0; i < N; ++i)
    {
      if(!(i%1000000))
        std::cout << i << " / " << N << " (" << (100. * i) / N << "%)" << std::endl;

      nus->GetEntry(i);

      if(ray_fv_length > 0)
        {
          hRayFVLength->Fill(ray_fv_length);
          mean_ray_length += ray_fv_length;
          ++N_considered;
        }
    }

  mean_ray_length /= N_considered;

  for(int i = 0; i < N; ++i)
    {
      if(!(i%1000000))
        std::cout << i << " / " << N << " (" << (100. * i) / N << "%)" << std::endl;

      nus->GetEntry(i);

      if(ray_fv_length > 0)
        {
          const double w = ray_fv_length / 440;//mean_ray_length;

          switch(nu_pdg)
            {
            case 14:
              hNuEnergyNuMu->Fill(nu_e, w);
              break;
            case -14:
              hNuEnergyANuMu->Fill(nu_e, w);
              break;
            case 12:
              hNuEnergyNuE->Fill(nu_e, w);
              break;
            case -12:
              hNuEnergyANuE->Fill(nu_e, w);
              break;
            }
        }
    }

  TFile *outfile = new TFile(saveDir + "/sbnd_flux.root", "RECREATE");

  hNuEnergyNuMu->Write("flux_sbnd_numu");
  hNuEnergyANuMu->Write("flux_sbnd_anumu");
  hNuEnergyNuE->Write("flux_sbnd_nue");
  hNuEnergyANuE->Write("flux_sbnd_anue");

  outfile->Close();

  TCanvas *cRayLength = new TCanvas("cRayLength", "cRayLength");
  cRayLength->cd();
  cRayLength->SetLogy();

  hRayFVLength->SetLineColor(kViolet-5);
  hRayFVLength->Draw("hist");

  cRayLength->SaveAs(saveDir + "/ray_fv_length.png");
  cRayLength->SaveAs(saveDir + "/ray_fv_length.pdf");

  TCanvas *cRayLengthCumul = new TCanvas("cRayLengthCumul", "cRayLengthCumul");
  cRayLengthCumul->cd();

  TH1* hRayFVLengthCumul = hRayFVLength->GetCumulative();
  hRayFVLengthCumul->GetYaxis()->SetTitle("Cumulative #nu");
  hRayFVLengthCumul->SetLineColor(kPink-5);
  hRayFVLengthCumul->Draw("hist");

  cRayLengthCumul->SaveAs(saveDir + "/ray_fv_length_cumul.png");
  cRayLengthCumul->SaveAs(saveDir + "/ray_fv_length_cumul.pdf");

  TCanvas *cNuMu = new TCanvas("cNuMu", "cNuMu");
  cNuMu->cd();
  hNuEnergyNuMu->Draw("hist");
  cNuMu->SaveAs(saveDir + "/numu_ray_traced_flux_energy.png");
  cNuMu->SaveAs(saveDir + "/numu_ray_traced_flux_energy.pdf");

  TCanvas *cANuMu = new TCanvas("cANuMu", "cANuMu");
  cANuMu->cd();
  hNuEnergyANuMu->Draw("hist");
  cANuMu->SaveAs(saveDir + "/anumu_ray_traced_flux_energy.png");
  cANuMu->SaveAs(saveDir + "/anumu_ray_traced_flux_energy.pdf");

  TCanvas *cNuE = new TCanvas("cNuE", "cNuE");
  cNuE->cd();
  hNuEnergyNuE->Draw("hist");
  cNuE->SaveAs(saveDir + "/nue_ray_traced_flux_energy.png");
  cNuE->SaveAs(saveDir + "/nue_ray_traced_flux_energy.pdf");

  TCanvas *cANuE = new TCanvas("cANuE", "cANuE");
  cANuE->cd();
  hNuEnergyANuE->Draw("hist");
  cANuE->SaveAs(saveDir + "/anue_ray_traced_flux_energy.png");
  cANuE->SaveAs(saveDir + "/anue_ray_traced_flux_energy.pdf");
}

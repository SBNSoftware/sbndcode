#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Common.C"

void MakeFluxHists(const TString productionVersion)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/flux_hists";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString file = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_flux_configI.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *nus = new TChain("fluxana/tree");
  nus->Add(file);

  const double fv_face_area = 175 * 360 * 2.;

  const int N = nus->GetEntries();

  int nu_pdg;
  float nu_x, nu_y, nu_e;

  nus->SetBranchAddress("nu_pdg", &nu_pdg);
  nus->SetBranchAddress("nu_x", &nu_x);
  nus->SetBranchAddress("nu_y", &nu_y);
  nus->SetBranchAddress("nu_e", &nu_e);

  TH1F *hNuEnergyNuMu  = new TH1F("hNuEnergyNuMu", ";True E_{#nu} (GeV);#nu_{#mu}", 200, 0, 10);
  TH1F *hNuEnergyANuMu = new TH1F("hNuEnergyANuMu", ";True E_{#nu} (GeV);#bar{#nu_{#mu}}", 100, 0, 10);
  TH1F *hNuEnergyNuE   = new TH1F("hNuEnergyNuE", ";True E_{#nu} (GeV);#nu_{e}", 50, 0, 10);
  TH1F *hNuEnergyANuE  = new TH1F("hNuEnergyANuE", ";True E_{#nu} (GeV);#bar{#nu_{e}}", 25, 0, 10);

  for(int i = 0; i < N; ++i)
    {
      nus->GetEntry(i);

      if(abs(nu_x) > 180 || abs(nu_x) < 5 || abs(nu_y) > 180)
        continue;

      switch(nu_pdg)
        {
        case 14:
          hNuEnergyNuMu->Fill(nu_e);
          break;
        case -14:
          hNuEnergyANuMu->Fill(nu_e);
          break;
        case 12:
          hNuEnergyNuE->Fill(nu_e);
          break;
        case -12:
          hNuEnergyANuE->Fill(nu_e);
          break;
        }
    }

  TFile *outfile = new TFile(saveDir + "/sbnd_flux.root", "RECREATE");

  hNuEnergyNuMu->Write("flux_sbnd_numu");
  hNuEnergyANuMu->Write("flux_sbnd_anumu");
  hNuEnergyNuE->Write("flux_sbnd_nue");
  hNuEnergyANuE->Write("flux_sbnd_anue");
}

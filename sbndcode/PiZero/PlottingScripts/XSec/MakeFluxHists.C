#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Common.C"

void Extrapolate(float &nu_x, float &nu_y, const float nu_z, const float &nu_other_x,
                 const float &nu_other_y, const float nu_other_z, const float extrap_z);

void MakeFluxHists(const TString productionVersion, const bool back, const bool fv_centre, const bool eff_z);

void MakeFluxHists(const TString productionVersion)
{
  MakeFluxHists(productionVersion, false, false, false);
  MakeFluxHists(productionVersion, true, false, false);
  MakeFluxHists(productionVersion, false, true, false);
  MakeFluxHists(productionVersion, false, false, true);
}

void MakeFluxHists(const TString productionVersion, const bool back, const bool fv_centre, const bool eff_z)
{
  if(back + fv_centre + eff_z > 1)
    throw std::runtime_error("Poorly configured, cannot consider back face & FV centre simultaneously");

  TString saveDir = baseSaveDir + "/" + productionVersion + "/flux_hists";

  if(back)
    saveDir += "_back_face";

  if(fv_centre)
    saveDir += "_fv_centre";

  if(eff_z)
    saveDir += "_eff_z";

  gSystem->Exec("mkdir -p " + saveDir);

  const TString file = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_flux_configL_*.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *nus = new TChain("fluxana/tree");
  nus->Add(file);

  const int N = nus->GetEntries();

  int nu_pdg;
  float nu_x, nu_y, nu_other_x, nu_other_y, nu_e;

  nus->SetBranchStatus("*", 0);
  nus->SetBranchStatus("nu_x", 1);
  nus->SetBranchStatus("nu_y", 1);
  nus->SetBranchStatus("nu_other_x", 1);
  nus->SetBranchStatus("nu_other_y", 1);
  nus->SetBranchStatus("nu_pdg", 1);
  nus->SetBranchStatus("nu_e", 1);

  nus->SetBranchAddress("nu_x", &nu_x);
  nus->SetBranchAddress("nu_y", &nu_y);
  nus->SetBranchAddress("nu_other_x", &nu_other_x);
  nus->SetBranchAddress("nu_other_y", &nu_other_y);
  nus->SetBranchAddress("nu_pdg", &nu_pdg);
  nus->SetBranchAddress("nu_e", &nu_e);

  TH1F *hNuEnergyNuMu  = new TH1F("hNuEnergyNuMu", ";True E_{#nu} (GeV);#nu_{#mu}", 200, 0, 10);
  TH1F *hNuEnergyANuMu = new TH1F("hNuEnergyANuMu", ";True E_{#nu} (GeV);#bar{#nu_{#mu}}", 200, 0, 10);
  TH1F *hNuEnergyNuE   = new TH1F("hNuEnergyNuE", ";True E_{#nu} (GeV);#nu_{e}", 200, 0, 10);
  TH1F *hNuEnergyANuE  = new TH1F("hNuEnergyANuE", ";True E_{#nu} (GeV);#bar{#nu_{e}}", 200, 0, 10);

  for(int i = 0; i < N; ++i)
    {
      if(!(i%1000000))
        std::cout << i << " / " << N << " (" << (100. * i) / N << "%)" << std::endl;

      nus->GetEntry(i);

      if(back)
        {
          nu_x = nu_other_x;
          nu_y = nu_other_y;
        }

      if(fv_centre)
        Extrapolate(nu_x, nu_y, 0, nu_other_x, nu_other_y, 500, 230);

      if(eff_z)
        Extrapolate(nu_x, nu_y, 0, nu_other_x, nu_other_y, 500, effbaseline - 11000);

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

void Extrapolate(float &nu_x, float &nu_y, const float nu_z, const float &nu_other_x,
                 const float &nu_other_y, const float nu_other_z, const float extrap_z)
{
  const TVector3 start(nu_x, nu_y, nu_z);
  const TVector3 end(nu_other_x, nu_other_y, nu_other_z);

  const float k = extrap_z / (nu_other_z - nu_z);

  const TVector3 extrap = start + k * (end - start);

  nu_x = extrap.X();
  nu_y = extrap.Y();
}

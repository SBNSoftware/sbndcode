#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

constexpr float BF    = 0.98823; // PDG branching fraction for pi0 -> gamma gamma
constexpr float units = 1e38;    // TO have cross section in 10^-38 cm^-2

void NuisanceXSecExtract(const TString productionVersion, const TString flavour, const TString cfg);

void NuisanceXSecExtract(const TString productionVersion, const TString cfg)
{
  NuisanceXSecExtract(productionVersion, "numu", cfg);
  NuisanceXSecExtract(productionVersion, "anumu", cfg);
  NuisanceXSecExtract(productionVersion, "nue", cfg);
  NuisanceXSecExtract(productionVersion, "anue", cfg);
}

void NuisanceXSecExtract(const TString productionVersion, const TString flavour, const TString cfg)
{
  const TString baseDir = "/exp/sbnd/data/users/hlay/ncpizero/generators";

  TChain* events = new TChain("FlatTree_VARS");
  events->Add(baseDir + "/genie/" + productionVersion + "/" + cfg + "/genie_" + flavour + "_prep.flat.root");

  Int_t Mode;

  Char_t cc;
  Int_t  PDGnu, nfsp;

  Int_t    pdg[40];
  Float_t  px[40], py[4], pz[40];

  Float_t  weight;
  Double_t scale_factor;

  events->SetBranchStatus("*", 0);
  events->SetBranchAddress("Mode", &Mode);
  events->SetBranchAddress("cc", &cc);
  events->SetBranchAddress("PDGnu", &PDGnu);
  events->SetBranchAddress("nfsp", &nfsp);
  events->SetBranchAddress("pdg", &pdg);
  events->SetBranchAddress("px", &px);
  events->SetBranchAddress("py", &py);
  events->SetBranchAddress("pz", &pz);
  events->SetBranchAddress("Weight", &weight);
  events->SetBranchAddress("fScaleFactor", &scale_factor);

  gSystem->Exec("mkdir -p " + baseDir + "/plots");
  TFile* outfile = new TFile(baseDir + "/genie_xsec_" + flavour + ".root", "RECREATE");
  if(!outfile->IsOpen())
    std::runtime_error("Couldn't open save file");

  const int N = events->GetEntries();

  float piZeroMomBins[9]       = { 0. - std::numeric_limits<float>::epsilon(), 60., 120., 180., 240., 300., 400., 600., 1000. + std::numeric_limits<float>::epsilon() };
  float cosThetaPiZeroBins[10] = { -1. - std::numeric_limits<float>::epsilon(), -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. + std::numeric_limits<float>::epsilon() };

  TH1F *hTotalIncl           = new TH1F("hTotalIncl", "", 1, -10, 1e10);
  TH1F *hPiZeroMomIncl       = new TH1F("hPiZeroMomIncl", "", 8, piZeroMomBins);
  TH1F *hCosThetaPiZeroIncl  = new TH1F("hCosThetaPiZeroIncl", "", 9, cosThetaPiZeroBins);
  TH1F *hTotal0p0pi          = new TH1F("hTotal0p0pi", "", 1, -10, 1e10);
  TH1F *hPiZeroMom0p0pi      = new TH1F("hPiZeroMom0p0pi", "", 8, piZeroMomBins);
  TH1F *hCosThetaPiZero0p0pi = new TH1F("hCosThetaPiZero0p0pi", "", 9, cosThetaPiZeroBins);
  TH1F *hTotalNp0pi          = new TH1F("hTotalNp0pi", "", 1, -10, 1e10);
  TH1F *hPiZeroMomNp0pi      = new TH1F("hPiZeroMomNp0pi", "", 8, piZeroMomBins);
  TH1F *hCosThetaPiZeroNp0pi = new TH1F("hCosThetaPiZeroNp0pi", "", 9, cosThetaPiZeroBins);

  double maxMomIncl = hPiZeroMomIncl->GetBinLowEdge(9),
    maxMom0p0pi = hPiZeroMom0p0pi->GetBinLowEdge(9), maxMomNp0pi = hPiZeroMomNp0pi->GetBinLowEdge(9);

  for(int evt = 0; evt < N; ++evt)
    {
      events->GetEntry(evt);

      int npizeros = 0, nchargedpi = 0, nprotons = 0;
      float mom, cosTheta;

      for(int fsp = 0; fsp < nfsp; ++fsp)
        {
          TVector3 momVec(px[fsp], py[fsp], pz[fsp]);

          if(pdg[fsp] == 111)
            {
              ++npizeros;
              mom      = momVec.Mag();
              cosTheta = momVec.Z()/ mom;

              mom *= 1000;
            }
          else if(abs(pdg[fsp]) == 211 && momVec.Mag() > .15)
            ++nchargedpi;
          else if(pdg[fsp] == 2212 && momVec.Mag() > .4)
            ++nprotons;
        }

      const double bfeffects = BF * TMath::Power(1 - BF, npizeros - 1);

      const double w = weight * scale_factor * units * bfeffects;

      if(cc == 0 && npizeros > 0)
        {
          if(abs(Mode) == 36)
            continue;

          hTotalIncl->Fill(0., w);
          hPiZeroMomIncl->Fill(mom, w);
          hCosThetaPiZeroIncl->Fill(cosTheta, w);

          if(mom > maxMomIncl)
            maxMomIncl = mom;

          if(nchargedpi == 0)
            {
              if(nprotons == 0)
                {
                  hTotal0p0pi->Fill(0., w);
                  hPiZeroMom0p0pi->Fill(mom, w);
                  hCosThetaPiZero0p0pi->Fill(cosTheta, w);

                  if(mom > maxMom0p0pi)
                    maxMom0p0pi = mom;
                }
              else
                {
                  hTotalNp0pi->Fill(0., w);
                  hPiZeroMomNp0pi->Fill(mom, w);
                  hCosThetaPiZeroNp0pi->Fill(cosTheta, w);

                  if(mom > maxMomNp0pi)
                    maxMomNp0pi = mom;
                }
            }
        }
    }

  std::cout << "\nTotal XSec\n"
            << "\t Incl:  " << hTotalIncl->Integral() << " " << hPiZeroMomIncl->Integral() << " " << hCosThetaPiZeroIncl->Integral() << '\n'
            << "\t 0p0pi: " << hTotal0p0pi->Integral() << " " << hPiZeroMom0p0pi->Integral() << " " << hCosThetaPiZero0p0pi->Integral() << '\n'
            << "\t Np0pi: " << hTotalNp0pi->Integral() << " " << hPiZeroMomNp0pi->Integral() << " " << hCosThetaPiZeroNp0pi->Integral() << '\n'
            << std::endl;

  TCanvas *cTotalIncl = new TCanvas("cTotalIncl", "cTotalIncl");
  cTotalIncl->cd();

  hTotalIncl->Draw("histe");
  hTotalIncl->Write("ncpizero_incl_" + flavour);

  cTotalIncl->SaveAs(baseDir + "/plots/total_incl_" + flavour + ".png");
  cTotalIncl->SaveAs(baseDir + "/plots/total_incl_" + flavour + ".pdf");

  TCanvas *cPiZeroMomIncl = new TCanvas("cPiZeroMomIncl", "cPiZeroMomIncl");
  cPiZeroMomIncl->cd();

  NormaliseEntriesByBinWidth(hPiZeroMomIncl);
  hPiZeroMomIncl->SetBinContent(9, hPiZeroMomIncl->GetBinContent(9) / (maxMomIncl - hPiZeroMomIncl->GetBinLowEdge(9)));
  hPiZeroMomIncl->Draw("histe");
  hPiZeroMomIncl->Write("pizero_mom_ncpizero_incl_" + flavour);

  cPiZeroMomIncl->SaveAs(baseDir + "/plots/pizero_mom_incl_" + flavour + ".png");
  cPiZeroMomIncl->SaveAs(baseDir + "/plots/pizero_mom_incl_" + flavour + ".pdf");

  TCanvas *cCosThetaPiZeroIncl = new TCanvas("cCosThetaPiZeroIncl", "cCosThetaPiZeroIncl");
  cCosThetaPiZeroIncl->cd();

  NormaliseEntriesByBinWidth(hCosThetaPiZeroIncl);
  hCosThetaPiZeroIncl->Draw("histe");
  hCosThetaPiZeroIncl->Write("cos_theta_pizero_ncpizero_incl_" + flavour);

  cCosThetaPiZeroIncl->SaveAs(baseDir + "/plots/cos_theta_pizero_incl_" + flavour + ".png");
  cCosThetaPiZeroIncl->SaveAs(baseDir + "/plots/cos_theta_pizero_incl_" + flavour + ".pdf");

  delete hTotalIncl;
  delete hPiZeroMomIncl;
  delete hCosThetaPiZeroIncl;

  delete cTotalIncl;
  delete cPiZeroMomIncl;
  delete cCosThetaPiZeroIncl;

  TCanvas *cTotal0p0pi = new TCanvas("cTotal0p0pi", "cTotal0p0pi");
  cTotal0p0pi->cd();

  hTotal0p0pi->Draw("histe");
  hTotal0p0pi->Write("ncpizero_0p0pi_" + flavour);

  cTotal0p0pi->SaveAs(baseDir + "/plots/total_0p0pi_" + flavour + ".png");
  cTotal0p0pi->SaveAs(baseDir + "/plots/total_0p0pi_" + flavour + ".pdf");

  TCanvas *cPiZeroMom0p0pi = new TCanvas("cPiZeroMom0p0pi", "cPiZeroMom0p0pi");
  cPiZeroMom0p0pi->cd();

  NormaliseEntriesByBinWidth(hPiZeroMom0p0pi);
  hPiZeroMom0p0pi->SetBinContent(9, hPiZeroMom0p0pi->GetBinContent(9) / (maxMom0p0pi - hPiZeroMom0p0pi->GetBinLowEdge(9)));
  hPiZeroMom0p0pi->Draw("histe");
  hPiZeroMom0p0pi->Write("pizero_mom_ncpizero_0p0pi_" + flavour);

  cPiZeroMom0p0pi->SaveAs(baseDir + "/plots/pizero_mom_0p0pi_" + flavour + ".png");
  cPiZeroMom0p0pi->SaveAs(baseDir + "/plots/pizero_mom_0p0pi_" + flavour + ".pdf");

  TCanvas *cCosThetaPiZero0p0pi = new TCanvas("cCosThetaPiZero0p0pi", "cCosThetaPiZero0p0pi");
  cCosThetaPiZero0p0pi->cd();

  NormaliseEntriesByBinWidth(hCosThetaPiZero0p0pi);
  hCosThetaPiZero0p0pi->Draw("histe");
  hCosThetaPiZero0p0pi->Write("cos_theta_pizero_ncpizero_0p0pi_" + flavour);

  cCosThetaPiZero0p0pi->SaveAs(baseDir + "/plots/cos_theta_pizero_0p0pi_" + flavour + ".png");
  cCosThetaPiZero0p0pi->SaveAs(baseDir + "/plots/cos_theta_pizero_0p0pi_" + flavour + ".pdf");

  delete hTotal0p0pi;
  delete hPiZeroMom0p0pi;
  delete hCosThetaPiZero0p0pi;

  delete cTotal0p0pi;
  delete cPiZeroMom0p0pi;
  delete cCosThetaPiZero0p0pi;

  TCanvas *cTotalNp0pi = new TCanvas("cTotalNp0pi", "cTotalNp0pi");
  cTotalNp0pi->cd();

  hTotalNp0pi->Draw("histe");
  hTotalNp0pi->Write("ncpizero_Np0pi_" + flavour);

  cTotalNp0pi->SaveAs(baseDir + "/plots/total_Np0pi_" + flavour + ".png");
  cTotalNp0pi->SaveAs(baseDir + "/plots/total_Np0pi_" + flavour + ".pdf");

  TCanvas *cPiZeroMomNp0pi = new TCanvas("cPiZeroMomNp0pi", "cPiZeroMomNp0pi");
  cPiZeroMomNp0pi->cd();

  NormaliseEntriesByBinWidth(hPiZeroMomNp0pi);
  hPiZeroMomNp0pi->SetBinContent(9, hPiZeroMomNp0pi->GetBinContent(9) / (maxMomNp0pi - hPiZeroMomNp0pi->GetBinLowEdge(9)));
  hPiZeroMomNp0pi->Draw("histe");
  hPiZeroMomNp0pi->Write("pizero_mom_ncpizero_Np0pi_" + flavour);

  cPiZeroMomNp0pi->SaveAs(baseDir + "/plots/pizero_mom_Np0pi_" + flavour + ".png");
  cPiZeroMomNp0pi->SaveAs(baseDir + "/plots/pizero_mom_Np0pi_" + flavour + ".pdf");

  TCanvas *cCosThetaPiZeroNp0pi = new TCanvas("cCosThetaPiZeroNp0pi", "cCosThetaPiZeroNp0pi");
  cCosThetaPiZeroNp0pi->cd();

  NormaliseEntriesByBinWidth(hCosThetaPiZeroNp0pi);
  hCosThetaPiZeroNp0pi->Draw("histe");
  hCosThetaPiZeroNp0pi->Write("cos_theta_pizero_ncpizero_Np0pi_" + flavour);

  cCosThetaPiZeroNp0pi->SaveAs(baseDir + "/plots/cos_theta_pizero_Np0pi_" + flavour + ".png");
  cCosThetaPiZeroNp0pi->SaveAs(baseDir + "/plots/cos_theta_pizero_Np0pi_" + flavour + ".pdf");

  delete hTotalNp0pi;
  delete hPiZeroMomNp0pi;
  delete hCosThetaPiZeroNp0pi;

  delete cTotalNp0pi;
  delete cPiZeroMomNp0pi;
  delete cCosThetaPiZeroNp0pi;
}

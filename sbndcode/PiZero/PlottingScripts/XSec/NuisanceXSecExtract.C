#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

constexpr float BF    = 0.98823; // PDG branching fraction for pi0 -> gamma gamma
constexpr float units = 1e38;    // TO have cross section in 10^-38 cm^-2

void NuisanceXSecExtract(const TString flavour);

void NuisanceXSecExtract()
{
  NuisanceXSecExtract("numu");
  NuisanceXSecExtract("anumu");
  NuisanceXSecExtract("nue");
  NuisanceXSecExtract("anue");
}

void NuisanceXSecExtract(const TString flavour)
{
  const TString baseDir = "/exp/sbnd/data/users/hlay/ncpizero/generators";

  TChain* events = new TChain("FlatTree_VARS");
  events->Add(baseDir + "/genie/genie_" + flavour + "_prep.flat.root");

  Char_t cc;
  Int_t  PDGnu, nfsp;

  Int_t    pdg[40];
  Float_t  px[40], py[4], pz[40];

  Float_t  weight;
  Double_t scale_factor;

  events->SetBranchStatus("*", 0);
  events->SetBranchAddress("cc", &cc);
  events->SetBranchAddress("PDGnu", &PDGnu);
  events->SetBranchAddress("nfsp", &nfsp);
  events->SetBranchAddress("pdg", &pdg);
  events->SetBranchAddress("px", &px);
  events->SetBranchAddress("py", &py);
  events->SetBranchAddress("pz", &pz);
  events->SetBranchAddress("Weight", &weight);
  events->SetBranchAddress("fScaleFactor", &scale_factor);

  TFile* outfile = new TFile(baseDir + "/genie_xsec_" + flavour + ".root", "RECREATE");
  if(!outfile->IsOpen())
    std::runtime_error("Couldn't open save file");

  const int N = events->GetEntries();

  float piZeroMomBins[9]       = { 0., 60., 120., 180., 240., 300., 400., 600., 1000. };
  float cosThetaPiZeroBins[10] = { -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. };

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
            << "\t Incl:  " << hPiZeroMomIncl->Integral() << " " << hCosThetaPiZeroIncl->Integral() << '\n'
            << "\t 0p0pi: " << hPiZeroMom0p0pi->Integral() << " " << hCosThetaPiZero0p0pi->Integral() << '\n'
            << "\t Np0pi: " << hPiZeroMomNp0pi->Integral() << " " << hCosThetaPiZeroNp0pi->Integral() << '\n'
            << std::endl;

  hTotalIncl->Write("ncpizero_incl_" + flavour);

  NormaliseEntriesByBinWidth(hPiZeroMomIncl);
  hPiZeroMomIncl->SetBinContent(9, hPiZeroMomIncl->GetBinContent(9) / (maxMomIncl - hPiZeroMomIncl->GetBinLowEdge(9)));
  hPiZeroMomIncl->Draw("histe");
  hPiZeroMomIncl->Write("pizero_mom_ncpizero_incl_" + flavour);

  NormaliseEntriesByBinWidth(hCosThetaPiZeroIncl);
  hCosThetaPiZeroIncl->Draw("hist");
  hCosThetaPiZeroIncl->Write("cos_theta_pizero_ncpizero_incl_" + flavour);

  delete hTotalIncl;
  delete hPiZeroMomIncl;
  delete hCosThetaPiZeroIncl;

  hTotal0p0pi->Write("ncpizero_0p0pi_" + flavour);

  NormaliseEntriesByBinWidth(hPiZeroMom0p0pi);
  hPiZeroMom0p0pi->SetBinContent(9, hPiZeroMom0p0pi->GetBinContent(9) / (maxMom0p0pi - hPiZeroMom0p0pi->GetBinLowEdge(9)));
  hPiZeroMom0p0pi->Draw("histe");
  hPiZeroMom0p0pi->Write("pizero_mom_ncpizero_0p0pi_" + flavour);

  NormaliseEntriesByBinWidth(hCosThetaPiZero0p0pi);
  hCosThetaPiZero0p0pi->Draw("hist");
  hCosThetaPiZero0p0pi->Write("cos_theta_pizero_ncpizero_0p0pi_" + flavour);

  delete hTotal0p0pi;
  delete hPiZeroMom0p0pi;
  delete hCosThetaPiZero0p0pi;

  hTotalNp0pi->Write("ncpizero_Np0pi_" + flavour);

  NormaliseEntriesByBinWidth(hPiZeroMomNp0pi);
  hPiZeroMomNp0pi->SetBinContent(9, hPiZeroMomNp0pi->GetBinContent(9) / (maxMomNp0pi - hPiZeroMomNp0pi->GetBinLowEdge(9)));
  hPiZeroMomNp0pi->Draw("histe");
  hPiZeroMomNp0pi->Write("pizero_mom_ncpizero_Np0pi_" + flavour);

  NormaliseEntriesByBinWidth(hCosThetaPiZeroNp0pi);
  hCosThetaPiZeroNp0pi->Draw("hist");
  hCosThetaPiZeroNp0pi->Write("cos_theta_pizero_ncpizero_Np0pi_" + flavour);

  delete hTotalNp0pi;
  delete hPiZeroMomNp0pi;
  delete hCosThetaPiZeroNp0pi;
}

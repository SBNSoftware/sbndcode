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

  TFile* outfile = new TFile(baseDir + "/genie_xsec_" + flavour + ".root", "RECREATE");

  TChain* events = new TChain("FlatTree_VARS");
  events->Add(baseDir + "/genie/tmp/genie_" + flavour + "_prep.flat.root");

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

  const int N = events->GetEntries();

  float piZeroMomBins[9]       = { 0., 60., 120., 180., 240., 300., 400., 600., 1000. };
  float cosThetaPiZeroBins[10] = { -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. };

  TH1F *hPiZeroMom      = new TH1F("hPiZeroMom", "", 8, piZeroMomBins);
  TH1F *hCosThetaPiZero = new TH1F("hCosThetaPiZero", "", 9, cosThetaPiZeroBins);

  for(int evt = 0; evt < N; ++evt)
    {
      events->GetEntry(evt);

      int npizeros = 0;
      float mom, cosTheta;

      for(int fsp = 0; fsp < nfsp; ++fsp)
        {
          if(pdg[fsp] == 111)
            {
              ++npizeros;
              TVector3 momVec(px[fsp], py[fsp], pz[fsp]);
              mom      = momVec.Mag();
              cosTheta = momVec.Z()/ mom;

              mom *= 1000;
            }
        }

      const double w = weight * scale_factor * units * BF;

      if(cc == 0 && npizeros == 1)
        {
          hPiZeroMom->Fill(mom, w);
          hCosThetaPiZero->Fill(cosTheta, w);
        }
    }

  NormaliseEntriesByBinWidth(hPiZeroMom);
  hPiZeroMom->Write("pizero_mom_" + flavour);

  NormaliseEntriesByBinWidth(hCosThetaPiZero);
  hCosThetaPiZero->Write("cos_theta_pizero_" + flavour);

  delete hPiZeroMom;
  delete hCosThetaPiZero;

  outfile->Close();
}

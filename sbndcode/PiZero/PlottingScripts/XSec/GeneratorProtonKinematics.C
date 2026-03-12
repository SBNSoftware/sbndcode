#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

constexpr float BF    = 0.98823; // PDG branching fraction for pi0 -> gamma gamma

void FillFromTree(TTree *tree, TH1D *mult, TH1D *multThresh, TH1D *protonMom);
TH1D* Combine(const std::map<TString, double> &flavours, std::map<TString, TH1D*> &hists);

void GeneratorProtonKinematics()
{
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  const TString baseDir = "/exp/sbnd/data/users/hlay/ncpizero/generators";
  gSystem->Exec("mkdir -p " + baseDir + "/comparison_plots");

  const std::map<TString, double> flavours = { { "numu", 0.92536980 },
					       { "anumu", 0.068619804 },
					       { "nue", 0.0054254241 },
					       { "anue", 0.00058497585 },
  };

  std::map<TString, TH1D*> hProtonMultiplicityGENIE, hProtonMultiplicityThreshGENIE, hProtonMomentumGENIE,
    hProtonMultiplicityNuWro, hProtonMultiplicityThreshNuWro, hProtonMomentumNuWro;

  for(auto const& [ flavour, frac ] : flavours)
    {
      hProtonMultiplicityGENIE[flavour]       = new TH1D(Form("hProtonMultiplicityGENIE_%s", flavour.Data()), "", 8, -.5, 7.5);
      hProtonMultiplicityThreshGENIE[flavour] = new TH1D(Form("hProtonMultiplicityThreshGENIE_%s", flavour.Data()), "", 8, -.5, 7.5);
      hProtonMomentumGENIE[flavour]           = new TH1D(Form("hProtonMomentumGENIE_%s", flavour.Data()), "", 25, 0, 1000);

      hProtonMultiplicityNuWro[flavour]       = new TH1D(Form("hProtonMultiplicityNuWro_%s", flavour.Data()), "", 8, -.5, 7.5);
      hProtonMultiplicityThreshNuWro[flavour] = new TH1D(Form("hProtonMultiplicityThreshNuWro_%s", flavour.Data()), "", 8, -.5, 7.5);
      hProtonMomentumNuWro[flavour]           = new TH1D(Form("hProtonMomentumNuWro_%s", flavour.Data()), "", 25, 0, 1000);

      TChain* genieEvents = new TChain("FlatTree_VARS");
      genieEvents->Add(baseDir + "/genie/NCPiZeroBv2/v3_04_02/" + flavour + "/genie_" + flavour + "_prep.flat.root");

      TChain* nuwroEvents = new TChain("FlatTree_VARS");
      nuwroEvents->Add(baseDir + "/nuwro/NCPiZeroBv2/NC/" + flavour + "/nuwro_" + flavour + "_prep.flat.root");

      FillFromTree(genieEvents, hProtonMultiplicityGENIE[flavour], hProtonMultiplicityThreshGENIE[flavour], hProtonMomentumGENIE[flavour]);
      FillFromTree(nuwroEvents, hProtonMultiplicityNuWro[flavour], hProtonMultiplicityThreshNuWro[flavour], hProtonMomentumNuWro[flavour]);
    }

  TCanvas *cProtonMultiplicity = new TCanvas("cProtonMultiplicity", "cProtonMultiplicity");
  cProtonMultiplicity->cd();

  TH1D *combinedProtonMultiplicityGENIE = Combine(flavours, hProtonMultiplicityGENIE);
  TH1D *combinedProtonMultiplicityNuWro = Combine(flavours, hProtonMultiplicityNuWro);

  combinedProtonMultiplicityGENIE->SetLineColor(kOrange+2);
  combinedProtonMultiplicityNuWro->SetLineColor(kGreen+2);

  combinedProtonMultiplicityGENIE->SetTitle(";N Protons;Interactions (A.U.)");
  combinedProtonMultiplicityGENIE->GetYaxis()->SetTitleOffset(1.1);
  combinedProtonMultiplicityGENIE->GetXaxis()->SetTitleOffset(1.1);

  combinedProtonMultiplicityGENIE->SetMarkerStyle(0);
  combinedProtonMultiplicityNuWro->SetMarkerStyle(0);

  combinedProtonMultiplicityGENIE->DrawNormalized("hist");
  combinedProtonMultiplicityNuWro->DrawNormalized("histsame");

  TLegend *leg = new TLegend(.5, .6, .75, .8);
  leg->AddEntry(combinedProtonMultiplicityGENIE, "GENIEv3 AR23_20i_00_000", "l");
  leg->AddEntry(combinedProtonMultiplicityNuWro, "NuWro v21.09.2", "l");
  leg->Draw();

  cProtonMultiplicity->SaveAs(baseDir + "/comparison_plots/proton_multiplicity.png");
  cProtonMultiplicity->SaveAs(baseDir + "/comparison_plots/proton_multiplicity.pdf");

  TCanvas *cProtonMultiplicityThresh = new TCanvas("cProtonMultiplicityThresh", "cProtonMultiplicityThresh");
  cProtonMultiplicityThresh->cd();

  TH1D *combinedProtonMultiplicityThreshGENIE = Combine(flavours, hProtonMultiplicityThreshGENIE);
  TH1D *combinedProtonMultiplicityThreshNuWro = Combine(flavours, hProtonMultiplicityThreshNuWro);

  combinedProtonMultiplicityThreshGENIE->SetLineColor(kOrange+2);
  combinedProtonMultiplicityThreshNuWro->SetLineColor(kGreen+2);

  combinedProtonMultiplicityThreshNuWro->SetTitle(";N Protons, p_{p}>400 MeV;Interactions (A.U.)");

  combinedProtonMultiplicityThreshGENIE->GetYaxis()->SetTitleOffset(1.1);
  combinedProtonMultiplicityThreshGENIE->GetXaxis()->SetTitleOffset(1.2);

  combinedProtonMultiplicityThreshGENIE->SetMarkerStyle(0);
  combinedProtonMultiplicityThreshNuWro->SetMarkerStyle(0);

  combinedProtonMultiplicityThreshNuWro->DrawNormalized("hist");
  combinedProtonMultiplicityThreshGENIE->DrawNormalized("histsame");

  leg->Draw();

  cProtonMultiplicityThresh->SaveAs(baseDir + "/comparison_plots/proton_multiplicity_thresh.png");
  cProtonMultiplicityThresh->SaveAs(baseDir + "/comparison_plots/proton_multiplicity_thresh.pdf");

  TCanvas *cProtonMomentum = new TCanvas("cProtonMomentum", "cProtonMomentum");
  cProtonMomentum->cd();

  TH1D *combinedProtonMomentumGENIE = Combine(flavours, hProtonMomentumGENIE);
  TH1D *combinedProtonMomentumNuWro = Combine(flavours, hProtonMomentumNuWro);

  combinedProtonMomentumGENIE->SetLineColor(kOrange+2);
  combinedProtonMomentumNuWro->SetLineColor(kGreen+2);

  combinedProtonMomentumGENIE->SetTitle(";p_{p} (MeV);Protons (A.U.)");

  combinedProtonMomentumGENIE->SetMarkerStyle(0);
  combinedProtonMomentumNuWro->SetMarkerStyle(0);

  combinedProtonMomentumGENIE->GetYaxis()->SetTitleOffset(1.1);
  combinedProtonMomentumGENIE->GetXaxis()->SetTitleOffset(1.1);

  combinedProtonMomentumGENIE->DrawNormalized("histe");
  combinedProtonMomentumNuWro->DrawNormalized("histesame");

  leg->Draw();

  cProtonMomentum->SaveAs(baseDir + "/comparison_plots/proton_momentum.png");
  cProtonMomentum->SaveAs(baseDir + "/comparison_plots/proton_momentum.pdf");
}

void FillFromTree(TTree *tree, TH1D *mult, TH1D *multThresh, TH1D *protonMom)
{
  Char_t cc;
  Int_t  PDGnu, nfsp;

  Int_t    pdg[40];
  Float_t  px[40], py[40], pz[40];
      
  tree->SetBranchStatus("*", 0);
  tree->SetBranchAddress("cc", &cc);
  tree->SetBranchAddress("PDGnu", &PDGnu);
  tree->SetBranchAddress("nfsp", &nfsp);
  tree->SetBranchAddress("pdg", &pdg);
  tree->SetBranchAddress("px", &px);
  tree->SetBranchAddress("py", &py);
  tree->SetBranchAddress("pz", &pz);

  const int N_genie = tree->GetEntries();

  for(int evt = 0; evt < N_genie; ++evt)
    {
      tree->GetEntry(evt);

      int npizeros = 0, nprotons = 0, nprotonsthresh = 0;

      for(int fsp = 0; fsp < nfsp; ++fsp)
	{
	  TVector3 momVec(px[fsp], py[fsp], pz[fsp]);
	      
	  if(pdg[fsp] == 111)
	    ++npizeros;
	      
	  if(pdg[fsp] == 2212 && momVec.Mag() > .4)
	    ++nprotonsthresh;
	      
	  if(pdg[fsp] == 2212)
	    ++nprotons;
	}

      const double bfeffects = BF * TMath::Power(1 - BF, npizeros - 1);

      if(cc == 0 && npizeros > 0)
	{
	  mult->Fill(nprotons, bfeffects);
	  multThresh->Fill(nprotonsthresh, bfeffects);

	  for(int fsp = 0; fsp < nfsp; ++fsp)
	    {
	      TVector3 momVec(px[fsp], py[fsp], pz[fsp]);
	      
	      if(pdg[fsp] == 2212)
		protonMom->Fill(momVec.Mag() * 1000, bfeffects);
	    }
	}
    }   
}

TH1D* Combine(const std::map<TString, double> &flavours, std::map<TString, TH1D*> &hists)
{
  TH1D* h;
  int counter = 0;

  for(auto const& [ flavour, frac ] : flavours)
    {
      if(counter == 0)
	{
	  h = (TH1D*) hists[flavour]->Clone(Form("%s_%s", hists[flavour]->GetName(), "_COMBINED"));
	  h->Scale(frac);
	}
      else
	h->Add(hists[flavour], frac);

      ++counter;
    }

  return h;
}

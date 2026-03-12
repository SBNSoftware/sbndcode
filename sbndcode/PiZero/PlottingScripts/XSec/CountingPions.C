#include "XSecCommon.C"

void CountingPions(const TString &productionVersion)
{
  gROOT->SetStyle("henrySBND");
  gStyle->SetPaintTextFormat("1.2g");
  gStyle->SetPadTopMargin(.1);
  gStyle->SetTitleOffset(1.35, "y");
  gROOT->ForceStyle();

  TString saveDir = baseSaveDir + "/" + productionVersion + "/counting_pions";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile  = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_rockbox.root";

  TChain *rockboxEvents = new TChain("ncpizeroana/events");
  rockboxEvents->Add(rockboxFile);
  TChain *rockboxSubruns = new TChain("ncpizeroana/subruns");
  rockboxSubruns->Add(rockboxFile);

  const double rockboxPOT     = GetPOT(rockboxSubruns);
  const double rockboxScaling = goalPOT / rockboxPOT;

  const int N = rockboxEvents->GetEntries();

  std::vector<int> *nu_pdg = 0, *nu_ccnc = 0, *nu_n_neutral_pions = 0, *nu_n_dalitz_neutral_pions = 0;
  std::vector<bool> *nu_av = 0, *nu_fv = 0;

  rockboxEvents->SetBranchStatus("*", 0);
  rockboxEvents->SetBranchStatus("nu_pdg", 1);
  rockboxEvents->SetBranchStatus("nu_ccnc", 1);
  rockboxEvents->SetBranchStatus("nu_av", 1);
  rockboxEvents->SetBranchStatus("nu_fv", 1);
  rockboxEvents->SetBranchStatus("nu_n_neutral_pions", 1);
  rockboxEvents->SetBranchStatus("nu_n_dalitz_neutral_pions", 1);

  rockboxEvents->SetBranchAddress("nu_pdg", &nu_pdg);
  rockboxEvents->SetBranchAddress("nu_ccnc", &nu_ccnc);
  rockboxEvents->SetBranchAddress("nu_av", &nu_av);
  rockboxEvents->SetBranchAddress("nu_fv", &nu_fv);
  rockboxEvents->SetBranchAddress("nu_n_neutral_pions", &nu_n_neutral_pions);
  rockboxEvents->SetBranchAddress("nu_n_dalitz_neutral_pions", &nu_n_dalitz_neutral_pions);

  int n_pizeros = 0, n_pizeros_ccnumu = 0, n_pizeros_nc = 0, n_pizeros_ccnue = 0;
  int n_pizero_evts = 0, n_pizero_evts_ccnumu = 0, n_pizero_evts_nc = 0, n_pizero_evts_ccnue = 0;

  int n_numu = 0, n_numucc = 0, n_numunc = 0, n_nue = 0, n_nuecc = 0, n_nuenc = 0;

  TH1D *hPiZeroMultiplicity       = new TH1D("hPiZeroMultiplicity", ";#pi^{0} multiplicity;AV Interactions", 5, -0.5, 4.5);
  TH1D *hPiZeroMultiplicityCCNuMu = new TH1D("hPiZeroMultiplicityCCNuMu", ";#pi^{0} multiplicity;AV CC#nu_{#mu} Interactions", 5, -0.5, 4.5);
  TH1D *hPiZeroMultiplicityCCNuE  = new TH1D("hPiZeroMultiplicityCCNuE", ";#pi^{0} multiplicity;AV CC#nu_{e} Interactions", 5, -0.5, 4.5);
  TH1D *hPiZeroMultiplicityNC     = new TH1D("hPiZeroMultiplicityNC", ";#pi^{0} multiplicity;AV NC Interactions", 5, -0.5, 4.5);

  TH1D *hPiZeroMultiplicityNoZero       = new TH1D("hPiZeroMultiplicityNoZero", ";#pi^{0} multiplicity;AV Interactions", 5, 0.5, 5.5);
  TH1D *hPiZeroMultiplicityNoZeroCCNuMu = new TH1D("hPiZeroMultiplicityNoZeroCCNuMu", ";#pi^{0} multiplicity;AV CC#nu_{#mu} Interactions", 5, 0.5, 5.5);
  TH1D *hPiZeroMultiplicityNoZeroCCNuE  = new TH1D("hPiZeroMultiplicityNoZeroCCNuE", ";#pi^{0} multiplicity;AV CC#nu_{e} Interactions", 5, 0.5, 5.5);
  TH1D *hPiZeroMultiplicityNoZeroNC     = new TH1D("hPiZeroMultiplicityNoZeroNC", ";#pi^{0} multiplicity;AV NC Interactions", 5, 0.5, 5.5);

  for(int evt = 0; evt < N; ++evt)
    {
      rockboxEvents->GetEntry(evt);

      for(int nu = 0; nu < nu_pdg->size(); ++nu)
	{
	  if(!nu_av->at(nu))
	    continue;

	  if(abs(nu_pdg->at(nu)) == 14)
	    {
	      ++n_numu;

	      if(nu_ccnc->at(nu) == 0)
		++n_numucc;
	      else
		++n_numunc;
	    }
	  if(abs(nu_pdg->at(nu)) == 12)
	    {
	      ++n_nue;

	      if(nu_ccnc->at(nu) == 0)
		++n_nuecc;
	      else
		++n_nuenc;
	    }

	  const int total_event_pizero = nu_n_neutral_pions->at(nu) + nu_n_dalitz_neutral_pions->at(nu);

	  n_pizeros += total_event_pizero;
	  hPiZeroMultiplicity->Fill(total_event_pizero);
	  hPiZeroMultiplicityNoZero->Fill(total_event_pizero);

	  if(nu_ccnc->at(nu) == 0)
	    {
	      if(abs(nu_pdg->at(nu)) == 14)
		{
		  n_pizeros_ccnumu += total_event_pizero;
		  hPiZeroMultiplicityCCNuMu->Fill(total_event_pizero);
		  hPiZeroMultiplicityNoZeroCCNuMu->Fill(total_event_pizero);
		}
	      if(abs(nu_pdg->at(nu)) == 12)
		{
		  n_pizeros_ccnue  += total_event_pizero;
		  hPiZeroMultiplicityCCNuE->Fill(total_event_pizero);
		  hPiZeroMultiplicityNoZeroCCNuE->Fill(total_event_pizero);
		}
	    }
	  else
	    {
	      n_pizeros_nc += total_event_pizero;
	      hPiZeroMultiplicityNC->Fill(total_event_pizero);
	      hPiZeroMultiplicityNoZeroNC->Fill(total_event_pizero);
	    }

	  if(total_event_pizero)
	    {
	      n_pizero_evts += 1;
	      if(nu_ccnc->at(nu) == 0)
		{
		  if(abs(nu_pdg->at(nu)) == 14)
		    n_pizero_evts_ccnumu += 1;
		  if(abs(nu_pdg->at(nu)) == 12)
		    n_pizero_evts_ccnue  += 1;
		}
	      else
		n_pizero_evts_nc += 1;
	    }
	}
    }

  std::cout << "\nRAW\n" << std::endl;

  std::cout << "Total Events with PiZeros: " << n_pizero_evts << '\n'
	    << "of which CCNuMu:           " << n_pizero_evts_ccnumu << '\n'
    	    << "of which CCNuE:            " << n_pizero_evts_ccnue << '\n'
    	    << "of which NC:               " << n_pizero_evts_nc << '\n'
	    << std::endl;

  std::cout << "Total PiZeros:   " << n_pizeros << '\n'
	    << "of which CCNuMu: " << n_pizeros_ccnumu << '\n'
    	    << "of which CCNuE:  " << n_pizeros_ccnue << '\n'
    	    << "of which NC:     " << n_pizeros_nc << '\n'
	    << std::endl;

  std::cout << "\nScaled to 1e21\n" << std::endl;

  std::cout << "Total Events with PiZeros: " << rockboxScaling * n_pizero_evts << '\n'
	    << "of which CCNuMu:           " << rockboxScaling * n_pizero_evts_ccnumu << '\n'
    	    << "of which CCNuE:            " << rockboxScaling * n_pizero_evts_ccnue << '\n'
    	    << "of which NC:               " << rockboxScaling * n_pizero_evts_nc << '\n'
	    << std::endl;

  std::cout << "Total PiZeros:   " << rockboxScaling * n_pizeros << '\n'
	    << "of which CCNuMu: " << rockboxScaling * n_pizeros_ccnumu << '\n'
    	    << "of which CCNuE:  " << rockboxScaling * n_pizeros_ccnue << '\n'
    	    << "of which NC:     " << rockboxScaling * n_pizeros_nc << '\n'
	    << std::endl;

  std::cout << "\nTotal Numbers\n" << std::endl;

  std::cout << "Total NuMu:  " << (rockboxScaling / 3.) * n_numu << '\n'
	    << "of which CC: " << (rockboxScaling / 3.) * n_numucc << '\n'
	    << "of which NC: " << (rockboxScaling / 3.) * n_numunc << '\n'
	    << "Total NuE:   " << (rockboxScaling / 3.) * n_nue << '\n'
	    << "of which CC: " << (rockboxScaling / 3.) * n_nuecc << '\n'
	    << "of which NC: " << (rockboxScaling / 3.) * n_nuenc << '\n'
	    << std::endl;
  
  TCanvas *cPiZeroMultiplicity = new TCanvas("cPiZeroMultiplicity", "cPiZeroMultiplicity");
  cPiZeroMultiplicity->cd();

  hPiZeroMultiplicity->Scale(rockboxScaling);
  hPiZeroMultiplicity->SetLineColor(kGray+2);
  hPiZeroMultiplicity->Draw("histe");

  cPiZeroMultiplicity->SaveAs(saveDir + "/pizero_av_multiplicity.png");
  cPiZeroMultiplicity->SaveAs(saveDir + "/pizero_av_multiplicity.pdf");

  TCanvas *cPiZeroMultiplicityCCNuMu = new TCanvas("cPiZeroMultiplicityCCNuMu", "cPiZeroMultiplicityCCNuMu");
  cPiZeroMultiplicityCCNuMu->cd();

  hPiZeroMultiplicityCCNuMu->Scale(rockboxScaling);
  hPiZeroMultiplicityCCNuMu->SetLineColor(kGreen+2);
  hPiZeroMultiplicityCCNuMu->Draw("histe");

  cPiZeroMultiplicityCCNuMu->SaveAs(saveDir + "/pizero_av_multiplicity_ccnumu.png");
  cPiZeroMultiplicityCCNuMu->SaveAs(saveDir + "/pizero_av_multiplicity_ccnumu.pdf");

  TCanvas *cPiZeroMultiplicityCCNuE = new TCanvas("cPiZeroMultiplicityCCNuE", "cPiZeroMultiplicityCCNuE");
  cPiZeroMultiplicityCCNuE->cd();

  hPiZeroMultiplicityCCNuE->Scale(rockboxScaling);
  hPiZeroMultiplicityCCNuE->SetLineColor(kCyan+2);
  hPiZeroMultiplicityCCNuE->Draw("histe");

  cPiZeroMultiplicityCCNuE->SaveAs(saveDir + "/pizero_av_multiplicity_ccnue.png");
  cPiZeroMultiplicityCCNuE->SaveAs(saveDir + "/pizero_av_multiplicity_ccnue.pdf");

  TCanvas *cPiZeroMultiplicityNC = new TCanvas("cPiZeroMultiplicityNC", "cPiZeroMultiplicityNC");
  cPiZeroMultiplicityNC->cd();

  hPiZeroMultiplicityNC->Scale(rockboxScaling);
  hPiZeroMultiplicityNC->SetLineColor(kMagenta+2);
  hPiZeroMultiplicityNC->Draw("histe");

  cPiZeroMultiplicityNC->SaveAs(saveDir + "/pizero_av_multiplicity_nc.png");
  cPiZeroMultiplicityNC->SaveAs(saveDir + "/pizero_av_multiplicity_nc.pdf");

  TCanvas *cPiZeroMultiplicityNoZero = new TCanvas("cPiZeroMultiplicityNoZero", "cPiZeroMultiplicityNoZero");
  cPiZeroMultiplicityNoZero->cd();

  hPiZeroMultiplicityNoZero->Scale(rockboxScaling);
  hPiZeroMultiplicityNoZero->SetLineColor(kGray+2);
  hPiZeroMultiplicityNoZero->Draw("histe");

  cPiZeroMultiplicityNoZero->SaveAs(saveDir + "/pizero_av_multiplicity_nozero.png");
  cPiZeroMultiplicityNoZero->SaveAs(saveDir + "/pizero_av_multiplicity_nozero.pdf");

  TCanvas *cPiZeroMultiplicityNoZeroCCNuMu = new TCanvas("cPiZeroMultiplicityNoZeroCCNuMu", "cPiZeroMultiplicityNoZeroCCNuMu");
  cPiZeroMultiplicityNoZeroCCNuMu->cd();

  hPiZeroMultiplicityNoZeroCCNuMu->Scale(rockboxScaling);
  hPiZeroMultiplicityNoZeroCCNuMu->SetLineColor(kGreen+2);
  hPiZeroMultiplicityNoZeroCCNuMu->Draw("histe");

  cPiZeroMultiplicityNoZeroCCNuMu->SaveAs(saveDir + "/pizero_av_multiplicity_nozero_ccnumu.png");
  cPiZeroMultiplicityNoZeroCCNuMu->SaveAs(saveDir + "/pizero_av_multiplicity_nozero_ccnumu.pdf");

  TCanvas *cPiZeroMultiplicityNoZeroCCNuE = new TCanvas("cPiZeroMultiplicityNoZeroCCNuE", "cPiZeroMultiplicityNoZeroCCNuE");
  cPiZeroMultiplicityNoZeroCCNuE->cd();

  hPiZeroMultiplicityNoZeroCCNuE->Scale(rockboxScaling);
  hPiZeroMultiplicityNoZeroCCNuE->SetLineColor(kCyan+2);
  hPiZeroMultiplicityNoZeroCCNuE->Draw("histe");

  cPiZeroMultiplicityNoZeroCCNuE->SaveAs(saveDir + "/pizero_av_multiplicity_nozero_ccnue.png");
  cPiZeroMultiplicityNoZeroCCNuE->SaveAs(saveDir + "/pizero_av_multiplicity_nozero_ccnue.pdf");

  TCanvas *cPiZeroMultiplicityNoZeroNC = new TCanvas("cPiZeroMultiplicityNoZeroNC", "cPiZeroMultiplicityNoZeroNC");
  cPiZeroMultiplicityNoZeroNC->cd();

  hPiZeroMultiplicityNoZeroNC->Scale(rockboxScaling);
  hPiZeroMultiplicityNoZeroNC->SetLineColor(kMagenta+2);
  hPiZeroMultiplicityNoZeroNC->Draw("histe");

  cPiZeroMultiplicityNoZeroNC->SaveAs(saveDir + "/pizero_av_multiplicity_nozero_nc.png");
  cPiZeroMultiplicityNoZeroNC->SaveAs(saveDir + "/pizero_av_multiplicity_nozero_nc.pdf");
}

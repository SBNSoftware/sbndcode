#include "Common.C"

std::vector<int> *slc_true_event_type_incl = 0;
std::vector<float> *slc_comp = 0;
std::vector<bool> *slc_sel_incl = 0;
std::vector<double> *slc_best_pzc_pizero_mom = 0, *slc_best_pzc_invariant_mass = 0, *slc_best_pzc_cos_theta_pizero = 0;

std::vector<size_t> *slc_best_pzc_photon_0_id = 0, *slc_best_pzc_photon_1_id = 0;

std::vector<std::vector<double>> *slc_true_pz_pizero_mom = 0, *slc_true_pz_cos_theta_pizero = 0, *slc_pfp_shower_dir_x = 0,
  *slc_pfp_shower_dir_y = 0, *slc_pfp_shower_dir_z = 0, *slc_pfp_track_dir_x = 0, *slc_pfp_track_dir_y = 0, *slc_pfp_track_dir_z = 0,
  *slc_pfp_shower_energy = 0, *slc_pfp_true_energy = 0;

TFile* file = TFile::Open("/exp/sbnd/app/users/hlay/ncpizero/srcs/sbndcode/sbndcode/PiZero/ShowerEnergyCorrection/shower_energy_correction_hist_NCPiZeroAv12.root");
TProfile *fShowerEnergyCorrectionHist = (TProfile*) file->Get("hShowerEnergy2DRecoFractionalResolution_pfx");

void InitialiseTree(TChain *tree);

double CorrectEnergy(const double &energy);

class Chi2Func
{
private:
  TMatrixD _covInv;
  TMatrixD _alpha_0;

public:
  Chi2Func(TMatrixD &covInv, TMatrixD &alpha_0)
  {
    _covInv.ResizeTo(3, 3);
    _covInv = covInv;

    _alpha_0.ResizeTo(1, 3);
    _alpha_0 = alpha_0;
  }

  double Eval(const double* pars)
  {
    std::cout << "Call of Eval with: \n" 
	      << "\t Lambda: " << pars[0] << '\n'
	      << "\t En0:    " << pars[1] << '\n'
	      << "\t En1:    " << pars[2] << '\n'
	      << "\t Theta:  " << pars[3] << '\n'
	      << std::endl;

    /*
    ROOT::Math::Minimizer* min = new ROOT::Minuit2::Minuit2Minimizer("MIGRAD");
    min->SetMaxFunctionCalls(100000);
    min->SetTolerance(0.01);
    min->SetPrintLevel(0);

    ROOT::Math::Functor f(this, &Chi2Func::Eval2, 4);
    min->SetFunction(f);

    double step[4]  = { 0., 0.01 * pars[1], 0.01 * pars[2], 0.01 * pars[3] };

    min->SetFixedVariable(0, "lambda", pars[0]);
    min->SetLimitedVariable(1, "en0", pars[1], step[1], 0, 134.9769);
    min->SetLimitedVariable(2, "en1", pars[2], step[2], 0, 134.9769);
    min->SetLimitedVariable(3, "theta", pars[3], step[3], 0, TMath::Pi());

    min->Minimize();
    */

    TMatrixD alpha = TMatrixD(1, 3);
    alpha(0, 0) = pars[1];
    alpha(0, 1) = pars[2];
    alpha(0, 2) = pars[3];

    double H = abs(2 * pars[1] * pars[2] * (1 - cos(pars[3])) - TMath::Power(134.9769, 2));

    /*
    if(min->Status() == 0)
      {
	alpha(0, 0) = min->X()[1];
	alpha(0, 1) = min->X()[2];
	alpha(0, 2) = min->X()[3];

	H = abs(2 * min->X()[1] * min->X()[2] * (1 - cos(min->X()[3])) - TMath::Power(134.9769, 2));
      }
    */

    TMatrixD a = ((alpha - _alpha_0) * _covInv) * (alpha - _alpha_0).T();

    return abs(pars[0] * H + a(0, 0));
  }

  double Eval2(const double* pars)
  {
    std::cout << "Call of Eval2 with: \n" 
	      << "\t Lambda: " << pars[0] << '\n'
	      << "\t En0:    " << pars[1] << '\n'
	      << "\t En1:    " << pars[2] << '\n'
	      << "\t Theta:  " << pars[3] << '\n'
	      << std::endl;

    TMatrixD alpha = TMatrixD(1, 3);
    alpha(0, 0) = pars[1];
    alpha(0, 1) = pars[2];
    alpha(0, 2) = pars[3];

    TMatrixD a = ((alpha - _alpha_0) * _covInv) * (alpha - _alpha_0).T();
    double H   = abs(2 * pars[1] * pars[2] * (1 - cos(pars[3])) - TMath::Power(134.9769, 2));

    return abs(pars[0] * H + a(0, 0));
  }

  double Eval3(const double* pars)
  {
    std::cout << "Call of Eval3 with: \n" 
	      << "\t En0:     " << pars[0] << '\n'
	      << "\t En1:     " << pars[1] << '\n'
	      << "\t Theta:   " << pars[2] << '\n'
	      << "\t InvMass: " << TMath::Sqrt(2 * pars[0] * pars[1] * (1 - cos(pars[2]))) << '\n'
	      << "\t Diff:    " << abs(2 * pars[0] * pars[1] * (1 - cos(pars[2])) - TMath::Power(134.9769, 2)) << '\n'
	      << std::endl;

    return abs(2 * pars[0] * pars[1] * (1 - cos(pars[2])) - TMath::Power(134.9769, 2));
  }
};

void FittedObservables(const TString productionVersion)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/observables_resolution_2";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString ncpizeroFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_ncpizero.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *ncpizeroEvents = new TChain("ncpizeroana/events");
  ncpizeroEvents->Add(ncpizeroFile);

  InitialiseTree(ncpizeroEvents);

  const double pizeroMomBins[9] = { 0., 60., 120., 180., 240., 300., 400., 600., 1000. };
  const double cosThetaBins[10] = { -1., -0.5, 0., 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1. };

  const int N = ncpizeroEvents->GetEntries();

  TH2F *covMatrix = new TH2F("covMatrix", "", 3, 0, 3, 3, 0, 3);

  struct PiZero
  {
    double en0;
    double en1;
    double theta;
    double trueEn0;
    double trueEn1;
    double trueTheta;
  };

  std::vector<PiZero> pizeros;

  for(int ev_i = 0; ev_i < N; ++ev_i)
    {
      ncpizeroEvents->GetEntry(ev_i);

      for(int slc_i = 0; slc_i < slc_true_event_type_incl->size(); ++slc_i)
	{
          if(slc_true_event_type_incl->at(slc_i) == 0 && slc_comp->at(slc_i) > .5 && slc_sel_incl->at(slc_i))
            {
              const double shwEn0 = slc_pfp_shower_energy->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i));
              const double shwEn1 = slc_pfp_shower_energy->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i));

              const double corrEn0 = CorrectEnergy(shwEn0);
              const double corrEn1 = CorrectEnergy(shwEn1);

              const TVector3 trkDir0 = TVector3(slc_pfp_track_dir_x->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)),
                                                slc_pfp_track_dir_y->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)),
                                                slc_pfp_track_dir_z->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)));

              const TVector3 trkDir1 = TVector3(slc_pfp_track_dir_x->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)),
                                                slc_pfp_track_dir_y->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)),
                                                slc_pfp_track_dir_z->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)));

              const double cosineThetaGammaGammaTrk = trkDir0.Dot(trkDir1) / (trkDir0.Mag() * trkDir1.Mag());
              const double thetaGammaGammaTrk       = acos(trkDir0.Dot(trkDir1) / (trkDir0.Mag() * trkDir1.Mag()));

              const double trueEn0   = slc_pfp_true_energy->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i));
              const double trueEn1   = slc_pfp_true_energy->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i));
	      const double trueTheta = acos(slc_true_pz_cos_theta_pizero->at(slc_i).at(0));

	      //	      std::cout << corrEn0 << " " << corrEn1 << "\n" << trueEn0 << " " << trueEn1 << '\n' << std::endl;

	      pizeros.push_back({ corrEn0, corrEn1, thetaGammaGammaTrk, trueEn0 * 1e3, trueEn1 * 1e3, trueTheta });
	    }
	}
    }

  const int npizeros = pizeros.size();
  double aveEn0 = 0., aveEn1 = 0., aveTheta = 0.;
  double varEn0En0 = 0., varEn0En1 = 0., varEn0Theta = 0.,
    varEn1En1 = 0., varEn1Theta = 0., varThetaTheta = 0.;

  for(const PiZero pizero : pizeros)
    {
      aveEn0   += pizero.en0;
      aveEn1   += pizero.en1;
      aveTheta += pizero.theta;

      varEn0En0     += (pizero.en0 - pizero.trueEn0) * (pizero.en0 - pizero.trueEn0);
      varEn0En1     += (pizero.en0 - pizero.trueEn0) * (pizero.en1 - pizero.trueEn1);
      varEn0Theta   += (pizero.en0 - pizero.trueEn0) * (pizero.theta - pizero.trueTheta);
      varEn1En1     += (pizero.en1 - pizero.trueEn1) * (pizero.en1 - pizero.trueEn1);
      varEn1Theta   += (pizero.en1 - pizero.trueEn1) * (pizero.theta - pizero.trueTheta);
      varThetaTheta += (pizero.theta - pizero.trueTheta) * (pizero.theta - pizero.trueTheta);
    }

  aveEn0   /= npizeros;
  aveEn1   /= npizeros;
  aveTheta /= npizeros;

  varEn0En0     /= npizeros;
  varEn0En1     /= npizeros;
  varEn0Theta   /= npizeros;
  varEn1En1     /= npizeros;
  varEn1Theta   /= npizeros;
  varThetaTheta /= npizeros;

  covMatrix->SetBinContent(1, 1, varEn0En0);
  covMatrix->SetBinContent(1, 2, varEn0En1);
  covMatrix->SetBinContent(1, 3, varEn0Theta);
  covMatrix->SetBinContent(2, 1, varEn0En1);
  covMatrix->SetBinContent(2, 2, varEn1En1);
  covMatrix->SetBinContent(2, 3, varEn1Theta);
  covMatrix->SetBinContent(3, 1, varEn0Theta);
  covMatrix->SetBinContent(3, 2, varEn1Theta);
  covMatrix->SetBinContent(3, 3, varThetaTheta);

  TCanvas *cCovMatrix = new TCanvas("cCovMatrix", "cCovMatrix");
  cCovMatrix->cd();

  cCovMatrix->SetRightMargin(.15);
  covMatrix->Draw("colztext");

  TCanvas *cCorrMatrix = new TCanvas("cCorrMatrix", "cCorrMatrix");
  cCorrMatrix->cd();

  TH2F *corrMatrix = (TH2F*) covMatrix->Clone("corrMatrix");
  
  for(int i = 1; i < 4; ++i)
    {
      for(int j = 1; j < 4; ++j)
	{
	  const double sigma_i = TMath::Sqrt(covMatrix->GetBinContent(i, i));
	  const double sigma_j = TMath::Sqrt(covMatrix->GetBinContent(j, j));

	  corrMatrix->SetBinContent(i, j, covMatrix->GetBinContent(i, j)/(sigma_i * sigma_j));
	}
    }

  cCorrMatrix->SetRightMargin(.15);
  corrMatrix->Draw("colztext");

  TMatrixD alpha_0 = TMatrixD(1, 3);
  alpha_0(0, 0) = aveEn0;
  alpha_0(0, 1) = aveEn1;
  alpha_0(0, 2) = aveTheta;

  TMatrixD cov(3, 3);
  cov(0, 0) = varEn0En0;
  cov(0, 1) = varEn0En1;
  cov(0, 2) = varEn0Theta;
  cov(1, 0) = varEn0En1;
  cov(1, 1) = varEn1En1;
  cov(1, 2) = varEn1Theta;
  cov(2, 0) = varEn0Theta;
  cov(2, 1) = varEn1Theta;
  cov(2, 2) = varThetaTheta;

  TMatrixD covInv = cov.Invert();

  Chi2Func chi2Func(covInv, alpha_0);

  ROOT::Math::Minimizer* min = new ROOT::Minuit2::Minuit2Minimizer("MIGRAD");

  min->SetMaxFunctionCalls(100000);
  min->SetTolerance(0.01);
  min->SetPrintLevel(1);

  //  ROOT::Math::Functor f(&chi2Func, &Chi2Func::Eval3, 3);
  ROOT::Math::Functor f(&chi2Func, &Chi2Func::Eval, 4);
  min->SetFunction(f);

  for(const PiZero pizero : pizeros)
    {
      double variable[4] = { 0.1, pizero.en0, pizero.en1, pizero.theta };
      double step[4]     = { 0.001, 0.1, 0.1, 0.1 };

      min->SetFixedVariable(0, "lambda", variable[0]);//, step[0]);
      min->SetLimitedVariable(1, "en0", variable[1], step[1], 0, 134.9769);
      min->SetLimitedVariable(2, "en1", variable[2], step[2], 0, 134.9769);
      min->SetLimitedVariable(3, "theta", variable[3], step[3], 0, TMath::Pi());

      /*
      double variable[3] = { pizero.en0, pizero.en1, pizero.theta };

      min->SetLimitedVariable(0, "en0", variable[0], 0.1 * variable[0], 0, 134.9769);
      min->SetLimitedVariable(1, "en1", variable[1], 0.1 * variable[1], 0, 134.9769);
      min->SetFixedVariable(2, "theta", variable[2]);
      */

      min->Minimize();

      break;
    }
}

void InitialiseTree(TChain *tree)
{
  tree->SetBranchStatus("*", 0);

  tree->SetBranchAddress("slc_true_event_type_incl", &slc_true_event_type_incl);
  tree->SetBranchAddress("slc_comp", &slc_comp);
  tree->SetBranchAddress("slc_sel_incl", &slc_sel_incl);
  tree->SetBranchAddress("slc_best_pzc_pizero_mom", &slc_best_pzc_pizero_mom);
  tree->SetBranchAddress("slc_best_pzc_invariant_mass", &slc_best_pzc_invariant_mass);
  tree->SetBranchAddress("slc_best_pzc_cos_theta_pizero", &slc_best_pzc_cos_theta_pizero);
  tree->SetBranchAddress("slc_best_pzc_photon_0_id", &slc_best_pzc_photon_0_id);
  tree->SetBranchAddress("slc_best_pzc_photon_1_id", &slc_best_pzc_photon_1_id);
  tree->SetBranchAddress("slc_true_pz_pizero_mom", &slc_true_pz_pizero_mom);
  tree->SetBranchAddress("slc_true_pz_cos_theta_pizero", &slc_true_pz_cos_theta_pizero);
  tree->SetBranchAddress("slc_pfp_shower_dir_x", &slc_pfp_shower_dir_x);
  tree->SetBranchAddress("slc_pfp_shower_dir_y", &slc_pfp_shower_dir_y);
  tree->SetBranchAddress("slc_pfp_shower_dir_z", &slc_pfp_shower_dir_z);
  tree->SetBranchAddress("slc_pfp_track_dir_x", &slc_pfp_track_dir_x);
  tree->SetBranchAddress("slc_pfp_track_dir_y", &slc_pfp_track_dir_y);
  tree->SetBranchAddress("slc_pfp_track_dir_z", &slc_pfp_track_dir_z);
  tree->SetBranchAddress("slc_pfp_shower_energy", &slc_pfp_shower_energy);
  tree->SetBranchAddress("slc_pfp_true_energy", &slc_pfp_true_energy);
}

double CorrectEnergy(const double &energy)
{
  const int bin = fShowerEnergyCorrectionHist->FindBin(energy);

  return energy * (1 - fShowerEnergyCorrectionHist->GetBinContent(bin));
}

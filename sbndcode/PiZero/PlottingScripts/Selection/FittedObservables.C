#include "Common.C"

constexpr double kPiZeroMass = 134.9769;

std::vector<int> *slc_true_event_type_incl = 0;
std::vector<float> *slc_comp = 0;
std::vector<bool> *slc_sel_incl = 0;
std::vector<double> *slc_best_pzc_pizero_mom = 0, *slc_best_pzc_invariant_mass = 0, *slc_best_pzc_cos_theta_pizero = 0;

std::vector<size_t> *slc_best_pzc_photon_0_id = 0, *slc_best_pzc_photon_1_id = 0;
std::vector<int> *slc_best_pzc_photon_0_true_trackid = 0, *slc_best_pzc_photon_1_true_trackid = 0;
std::vector<float> *slc_best_pzc_photon_0_comp = 0, *slc_best_pzc_photon_1_comp = 0, *slc_best_pzc_photon_0_pur = 0, *slc_best_pzc_photon_1_pur = 0;

std::vector<std::vector<double>> *slc_true_pz_pizero_mom = 0, *slc_true_pz_cos_theta_pizero = 0, *slc_pfp_shower_dir_x = 0,
  *slc_pfp_shower_dir_y = 0, *slc_pfp_shower_dir_z = 0, *slc_pfp_track_dir_x = 0, *slc_pfp_track_dir_y = 0, *slc_pfp_track_dir_z = 0,
  *slc_pfp_shower_energy = 0, *slc_pfp_true_energy = 0, *slc_true_pz_open_angle = 0, *slc_true_pz_gamma0_energy = 0, *slc_true_pz_gamma1_energy = 0;
std::vector<std::vector<int>> *slc_true_pz_gamma0_trackid = 0, *slc_true_pz_gamma1_trackid = 0;

TFile* file = TFile::Open("/exp/sbnd/app/users/hlay/ncpizero/srcs/sbndcode/sbndcode/PiZero/ShowerEnergyCorrection/shower_energy_correction_hist_NCPiZeroAv12.root");
TProfile *fShowerEnergyCorrectionHist = (TProfile*) file->Get("hShowerEnergy2DRecoFractionalResolution_pfx");

void InitialiseTree(TChain *tree);

double CorrectEnergy(const double &energy);

class KinFit{

public:
  static void SetPar(const double par[], const int npar)
  {
    if(npar != fNpar - 1)
      {
        printf("npar mismatch! %d %d\n", npar, fNpar);
        exit(1);
      }

    fLambda   = par[0];
    fOptEn0   = par[1];
    fOptEn1   = par[2];
    fOptTheta = par[3];
  }

  static double Constraint()
  {
    return abs(2 * fOptEn0 * fOptEn1 * (1 - cos(fOptTheta)) - TMath::Power(kPiZeroMass, 2));
  }

  static double FullLikelihood()
  {
    return CoreLikelihood() + fLambda * Constraint();
  }

  static void IniCoreMIN(TMinuit * mnt, const double inputl)
  {
    mnt->DefineParameter(0, "lambda", inputl, 1e-2, -1e6, 1e6);
    mnt->FixParameter(0);

    mnt->DefineParameter(1, "en0", fOptEn0, 50, 0.8 * fOptEn0, 1500);
    mnt->DefineParameter(2, "en1", fOptEn1, 50, 0.8 * fOptEn1, 1200);

    double upper = 1.2 * fOptTheta;
    if(upper > TMath::Pi())
      upper = TMath::Pi();

    mnt->DefineParameter(3, "theta", fOptTheta, 1e-1, 0.8 * fOptTheta, upper);
  }

  static bool IsConstraintGood()
  {
    const double eps = 10;
    return (TMath::Abs(Constraint())<eps);
  }

  static double GetOptEn0() { return fOptEn0; }

  static double GetOptEn1() { return fOptEn1; }

  static double GetOptTheta() { return fOptTheta; }

  static double GetEn0() { return fEn0; }

  static double GetEn1() { return fEn1; }

  static double GetTheta() { return fTheta; }

  static double GetLambda(){ return fLambda; }

  static double GetNpar() { return fNpar; }

  static void SetVars(const double iniVarEn0, const double iniVarEn1, const double iniVarTheta)
  {
    fEn0   = iniVarEn0;
    fEn1   = iniVarEn1;
    fTheta = iniVarTheta;
  }

  static void SetOptVars(const double varEn0, const double varEn1, const double varTheta)
  {
    fOptEn0   = varEn0;
    fOptEn1   = varEn1;
    fOptTheta = varTheta;
  }

  static void SetCVM(const vector<double> V)
  {
    int size = V.size();
    int dim  = sqrt(size);

    for(int i = 0; i < size; i++)
      {
        int x = floor(i / dim);
        int y = i % dim;

        fCovMatrix[x][y] = V[i];
      }
  }

private:
  static double fLambda;
  static double fOptEn0;
  static double fOptEn1;
  static double fOptTheta;
  static const int fNpar;

  static double fEn0;
  static double fEn1;
  static double fTheta;

  static TMatrixD fCovMatrix;

  static double CoreLikelihood()
  {
    TMatrixD CovMatrixtmp     = fCovMatrix;
    TMatrixD CovMatrixInverse = CovMatrixtmp.Invert();

    int nparameters = fNpar-1;

    TMatrixD Diff(nparameters,1);
    Diff[0][0] = fOptEn0 - fEn0;
    Diff[1][0] = fOptEn1 - fEn1;
    Diff[2][0] = fOptTheta - fTheta;

    TMatrixD DiffT(1,nparameters);
    DiffT.Transpose(Diff);

    TMatrixD Chi2 = (DiffT*CovMatrixInverse*Diff);

    if(Chi2[0][0]<0)
      throw std::runtime_error("Get out");

    return Chi2[0][0];
  }
};

double KinFit::fLambda   = -999;
double KinFit::fOptEn0   = -999;
double KinFit::fOptEn1   = -999;
double KinFit::fOptTheta = -999;
const int KinFit::fNpar  = 4;

double KinFit::fEn0   = -999;
double KinFit::fEn1   = -999;
double KinFit::fTheta = -999;

const int npars = KinFit::GetNpar() - 1;

TMatrixD KinFit::fCovMatrix(npars, npars);

void CoreFCN(int &npars, double *grad, double &value, double *par, int flag)
{
  KinFit::SetPar(par, npars);

  value = KinFit::FullLikelihood();
}

void LambdaFCN(int &npars, double *grad, double &value, double *par, int flag)
{
  TMinuit * CoreMIN = new TMinuit(3);
  CoreMIN->SetPrintLevel(-1);

  CoreMIN->SetFCN(CoreFCN);

  KinFit::IniCoreMIN(CoreMIN, par[0]);

  int flagL = CoreMIN->Command("MIGRAD");

  int irun = 1;
  const int maxnrun = 2;

  while(flagL!=0)
    {
      ++irun;
      flagL = CoreMIN->Command("MIGRAD");

      if(irun>=maxnrun)
        break;
    }

  if(flagL==0)
    value = TMath::Abs(KinFit::Constraint());
  else
    value = 1E50;

  delete CoreMIN;
}

void SetIniValues(const vector<double> iniVar, const vector<double> CVM)
{
  int VarsSize = iniVar.size();
  int CVMDim   = sqrt(CVM.size());

  if(VarsSize != CVMDim)
    {
      cout <<  "Variable size is not the same as CVM dimension! " << VarsSize << " " << CVMDim << endl;
      exit(1);
    }

  KinFit::SetVars(iniVar[0], iniVar[1], iniVar[2]);
  KinFit::SetCVM(CVM);
}

bool DoubleMin(const double iniLambda, const double lmin, const double lmax, vector<double> iniVar, vector<double> CVM)
{
  SetIniValues(iniVar, CVM);

  TMinuit * LambdaMIN = new TMinuit(1);
  LambdaMIN->SetPrintLevel(-1);

  LambdaMIN->SetFCN(LambdaFCN);

  LambdaMIN->DefineParameter(0, "lambda", iniLambda, 0.1 * iniLambda, lmin, lmax);

  int flag = LambdaMIN->Command("MIGRAD");

  int irun = 1;
  const int maxnrun=20;
  std::srand(std::time(nullptr));

  while(flag != 0 || ! KinFit::IsConstraintGood())
    {
      irun++;
      KinFit::SetOptVars(KinFit::GetEn0(), KinFit::GetEn1(), KinFit::GetTheta());

      flag = LambdaMIN->Command("MIGRAD");

      if(irun>=maxnrun)
        break;
    }

  delete LambdaMIN;

  if(flag==0 && KinFit::IsConstraintGood())
    return true;
  else
    return false;
}

vector<double> DoKF(const double &LdShowerEnergyRaw, const double &SlShowerEnergyRaw, const double &OpenAngle, vector<double> CVM, bool &GoodFit)
{
  vector<double> IniVars;
  vector<double> FittedVars;
  IniVars.push_back(LdShowerEnergyRaw);
  IniVars.push_back(SlShowerEnergyRaw);
  IniVars.push_back(OpenAngle);
  KinFit::SetOptVars(LdShowerEnergyRaw,SlShowerEnergyRaw,OpenAngle);

  GoodFit = DoubleMin(1e-5, -1e-4, 1e-4, IniVars, CVM);

  if(GoodFit)
    {
      FittedVars.push_back(KinFit::GetOptEn0());
      FittedVars.push_back(KinFit::GetOptEn1());
      FittedVars.push_back(KinFit::GetOptTheta());
    }
  else
    {
      FittedVars.push_back(KinFit::GetEn0());
      FittedVars.push_back(KinFit::GetEn1());
      FittedVars.push_back(KinFit::GetTheta());
    }

  IniVars.clear();

  return FittedVars;
}

std::vector<double> GetCovMatrix(const TString productionVersion, const bool diag = false)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/observables_resolution_2";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString ncpizeroFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_ncpizero.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *ncpizeroEvents = new TChain("ncpizeroana/events");
  ncpizeroEvents->Add(ncpizeroFile);

  InitialiseTree(ncpizeroEvents);

  const int N = ncpizeroEvents->GetEntries();

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
          bool signal = slc_true_event_type_incl->at(slc_i) == 0 && slc_comp->at(slc_i) > .5;

          if(signal && slc_sel_incl->at(slc_i))
            {
              const bool goodReco = ((slc_best_pzc_photon_0_true_trackid->at(slc_i) == slc_true_pz_gamma0_trackid->at(slc_i).at(0) &&
                                      slc_best_pzc_photon_1_true_trackid->at(slc_i) == slc_true_pz_gamma1_trackid->at(slc_i).at(0)) ||
                                     (slc_best_pzc_photon_0_true_trackid->at(slc_i) == slc_true_pz_gamma1_trackid->at(slc_i).at(0) &&
                                      slc_best_pzc_photon_1_true_trackid->at(slc_i) == slc_true_pz_gamma0_trackid->at(slc_i).at(0)))
                && slc_best_pzc_photon_0_comp->at(slc_i) > .8 && slc_best_pzc_photon_1_comp->at(slc_i) > .8
                && slc_best_pzc_photon_0_pur->at(slc_i) > .8 && slc_best_pzc_photon_1_pur->at(slc_i) > .8;

              /*
                if(!goodReco)
                continue;
              */

              const double shwEn0 = slc_pfp_shower_energy->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i));
              const double shwEn1 = slc_pfp_shower_energy->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i));

              double corrEn0 = CorrectEnergy(shwEn0);
              double corrEn1 = CorrectEnergy(shwEn1);

              const TVector3 trkDir0 = TVector3(slc_pfp_track_dir_x->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)),
                                                slc_pfp_track_dir_y->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)),
                                                slc_pfp_track_dir_z->at(slc_i).at(slc_best_pzc_photon_0_id->at(slc_i)));

              const TVector3 trkDir1 = TVector3(slc_pfp_track_dir_x->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)),
                                                slc_pfp_track_dir_y->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)),
                                                slc_pfp_track_dir_z->at(slc_i).at(slc_best_pzc_photon_1_id->at(slc_i)));

              const double cosineThetaGammaGammaTrk = trkDir0.Dot(trkDir1) / (trkDir0.Mag() * trkDir1.Mag());
              const double thetaGammaGammaTrk       = acos(trkDir0.Dot(trkDir1) / (trkDir0.Mag() * trkDir1.Mag()));

              double trueEn0   = slc_true_pz_gamma0_energy->at(slc_i).at(0);
              double trueEn1   = slc_true_pz_gamma1_energy->at(slc_i).at(0);
              const double trueTheta = TMath::DegToRad() * slc_true_pz_open_angle->at(slc_i).at(0);

              if(corrEn0 < corrEn1)
                std::swap(corrEn0, corrEn1);
              if(corrEn0 < corrEn1)
                throw std::runtime_error("Reco energy order");
              if(trueEn0 < trueEn1)
                std::swap(trueEn0, trueEn1);
              if(trueEn0 < trueEn1)
                throw std::runtime_error("True energy order");

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

  varEn0En0     /= npizeros;
  varEn0En1     /= npizeros;
  varEn0Theta   /= npizeros;
  varEn1En1     /= npizeros;
  varEn1Theta   /= npizeros;
  varThetaTheta /= npizeros;

  std::vector<double> covVec = { varEn0En0,   varEn0En1,   varEn0Theta,
                                 varEn0En1,   varEn1En1,   varEn1Theta,
                                 varEn0Theta, varEn1Theta, varThetaTheta };

  std::vector<double> diagCovVec = { varEn0En0, 0,         0,
                                     0,         varEn1En1, 0,
                                     0,         0,         varThetaTheta };

  return diag ? diagCovVec : covVec;
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
  tree->SetBranchAddress("slc_best_pzc_photon_0_true_trackid", &slc_best_pzc_photon_0_true_trackid);
  tree->SetBranchAddress("slc_best_pzc_photon_1_true_trackid", &slc_best_pzc_photon_1_true_trackid);
  tree->SetBranchAddress("slc_best_pzc_photon_0_comp", &slc_best_pzc_photon_0_comp);
  tree->SetBranchAddress("slc_best_pzc_photon_1_comp", &slc_best_pzc_photon_1_comp);
  tree->SetBranchAddress("slc_best_pzc_photon_0_pur", &slc_best_pzc_photon_0_pur);
  tree->SetBranchAddress("slc_best_pzc_photon_1_pur", &slc_best_pzc_photon_1_pur);
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
  tree->SetBranchAddress("slc_true_pz_open_angle", &slc_true_pz_open_angle);
  tree->SetBranchAddress("slc_true_pz_gamma0_energy", &slc_true_pz_gamma0_energy);
  tree->SetBranchAddress("slc_true_pz_gamma1_energy", &slc_true_pz_gamma1_energy);
  tree->SetBranchAddress("slc_true_pz_gamma0_trackid", &slc_true_pz_gamma0_trackid);
  tree->SetBranchAddress("slc_true_pz_gamma1_trackid", &slc_true_pz_gamma1_trackid);
}

double CorrectEnergy(const double &energy)
{
  const int bin = fShowerEnergyCorrectionHist->FindBin(energy);

  return energy * (1 - fShowerEnergyCorrectionHist->GetBinContent(bin));
}

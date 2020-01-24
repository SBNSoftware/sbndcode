#include "NueRecoOptimiser.hh"

//C++ Includes 
#include <iostream>
#include <typeinfo>

//Root Includes
#include "TFile.h"

optimiser::NueRecoOptimiser::NueRecoOptimiser(TTree* signaltree, TTree* backgroundtree, float bkpot, float sigpot){

  SignalTree.AddBranch<std::vector<float>>("number_of_showers_per_neutrino",signaltree);
  SignalTree.AddBranch<std::vector<float>>("vertex_recoK",signaltree);
  SignalTree.AddBranch<std::vector<float>>("vertex_trueK",signaltree);
  SignalTree.AddBranch<std::vector<float>>("vertex_reco",signaltree);
  SignalTree.AddBranch<std::vector<float>>("nu_reco_energy",signaltree);
  SignalTree.AddBranch<std::vector<float>>("nu_truth_energy",signaltree);
  SignalTree.AddBranch<std::vector<float>>("nu_interaction_type",signaltree);
  SignalTree.AddBranch<std::vector<float>>("nu_mode",signaltree);
  SignalTree.AddBranch<std::vector<float>>("nu_E",signaltree);
  SignalTree.AddBranch<std::vector<float>>("nu_E_numtrue",signaltree);
  SignalTree.AddBranch<std::vector<float>>("nu_distance",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("truepionE",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("trueprotonE",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("truekaonE",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("truetrackE",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_energy",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("truth_pid",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("true_energy",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_coversion_gap",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_residual_dist",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_length",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_length_perp",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_perp",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_3D",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_perp",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_ratio",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_ratio",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_sq",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_perp_sq",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_ratio_sq",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_dEdx",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("track_lengths",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("track_E",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("track_trueE",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("track_pdg",signaltree); 
  SignalTree.AddBranch<std::vector<std::vector<float> >>("track_resE",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_new",signaltree);
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_open_angle",signaltree);
  
  SignalTree.AddBranch<std::vector<float> >("in_FV",signaltree);
  SignalTree.AddBranch<std::vector<float> >("nu_osc_prob",signaltree);
  SignalTree.AddBranch<std::vector<float> >("nu_pdg",signaltree);

  SignalTree.AddBranch<int>("trueShower_num",signaltree);
  SignalTree.AddBranch<int>("numtrueVtx_branch",signaltree);
  SignalTree.AddBranch<float>("POT",signaltree);

  BackgroundTree.AddBranch<std::vector<float>>("number_of_showers_per_neutrino",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("vertex_recoK",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("vertex_trueK",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("vertex_reco",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("nu_reco_energy",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("nu_truth_energy",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("nu_interaction_type",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("nu_mode",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("nu_E",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("nu_E_numtrue",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float>>("nu_distance",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("truepionE",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("trueprotonE",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("truekaonE",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("truetrackE",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_energy",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("truth_pid",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("true_energy",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_coversion_gap",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_residual_dist",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_length",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_length_perp",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_perp",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_3D",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_perp",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_ratio",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_ratio",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_perp_sq",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_sq",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_ratio_sq",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_dEdx",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("track_lengths",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("track_E",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("track_trueE",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("track_pdg",backgroundtree); 
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("track_resE",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_density_grad_new",backgroundtree);
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_open_angle",backgroundtree);

  BackgroundTree.AddBranch<std::vector<float> >("in_FV",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float> >("nu_osc_prob",backgroundtree);
  BackgroundTree.AddBranch<std::vector<float> >("nu_pdg",backgroundtree);


  BackgroundTree.AddBranch<int>("trueShower_num",backgroundtree);
  BackgroundTree.AddBranch<int>("numtrueVtx_branch",backgroundtree);
  BackgroundTree.AddBranch<float>("POT",backgroundtree);

  MVAFile = new TFile("MVAFile.root","RECREATE");
  MVAFile->cd();
  BackgroundTreeMVA = new TTree("BackgroundTreeMVA","BackgroundTreeMVA");
  SignalTreeMVA     = new TTree("SignalTreeMVA","SignalTreeMVA");

  BackgroundTreeMVA->Branch("weight",BackgrndWeight);
  SignalTreeMVA->Branch("weight",SignalWeight);


  TotalSigPOT = sigpot;
  TotalBKPOT  = bkpot;

}

void optimiser::NueRecoOptimiser::InitialiseMetrics(){
  
  //MVA Analsysis
  ApplyOscWht = true;
  ApplyPOTWht = true;
  InitialiseMetric("dEdxMVA","shower_dEdx",100,0,10,"StandardShowerCutFinder",true,true,true,-1);
  InitialiseMetric("CoversionGapMVA","shower_coversion_gap",100,0,10,"StandardShowerCutFinder",true,true,true,-1);
  InitialiseMetric("MaxTrackLenghMVA","track_lengths",50,0,200,"StandardShowerCutFinder",true,true,true,0); 
  InitialiseMetric("ShowerEnergyMVA","shower_energy",150,0,3000,"StandardShowerCutFinder",true,true,true,-100);
  InitialiseMetric("ShowerLengthOverEMVA","shower_length",50,0,0.5,"StandardOverEDaughterCutFinder",true,true,true,-0.2);
  InitialiseMetric("ShowerDensityOverEMVA","shower_density",50,0,0.1,"StandardOverEDaughterCutFinder",true,true,true,-0.1);
  InitialiseMetric("ShowerLengthPerpOverEMVA","shower_length_perp",50,0,0.2,"StandardOverEDaughterCutFinder",true,true, true,-0.1);
  InitialiseMetric("ShowerDensityGradPerpSqOverEMVA","shower_density_grad_perp_sq",50,0,0.1,"StandardOverEDaughterCutFinder",true,true,true,-0.1);
  InitialiseMetric("ShowerDensityGradSqOverEMVA","shower_density_grad_sq",50,0,0.1,"StandardOverEDaughterCutFinder",true,true,true,-0.1);
  InitialiseMetric("ShowerDensityRatioOverEMVA","shower_density_ratio",50,0,0.2,"StandardOverEDaughterCutFinder",true,true,true,-0.1);
  InitialiseMetric("NumberOfNeutrinosMVA","nu_reco_energy",10,0,10,"OneShowerSizeNeutrinoAnalysisCutFinder",true,true,true,-1);
  InitialiseMetric("ShowerResidualNumShowersMVA","shower_residual_dist",10,0,10,"ShowerResidualNumShowersCutFinder",true,true,true,-1);
  

  InitialiseMetric("dEdxOneShowerCut","shower_dEdx",100,0,10,"OneShower1DCutFinder");
  InitialiseMetric("CoversionGapOneShowerCut","shower_coversion_gap",50,0,10,"OneShower1DCutFinder");
  InitialiseMetric("MaxTrackLenghOneShowerCut","track_lengths",25,0,200,"OneShower1DCutFinder"); 
  
  //Sandard Cut Analysis
  InitialiseMetric("dEdx","shower_dEdx",100,0,10);
  InitialiseMetric("dEdxPOT","shower_dEdx",100,0,10,"StandardShowerCutFinder",true,true);

  InitialiseMetric("CoversionGap","shower_coversion_gap",100,0,10);
  InitialiseMetric("MaxTrackLengh","track_lengths",25,0,200); 
  InitialiseMetric("ShowerEnergy","shower_energy",150,0,3000);
  InitialiseMetric("zBDT","nu_reco_energy",200,-1,1,"MVAAnalysisCutFinder",true,true,false,-999,"BDT");
  InitialiseMetric("zBDTG","nu_reco_energy",200,-1,1,"MVAAnalysisCutFinder",true,true,false,-999,"BDTG");

  InitialiseMetric("ShowerTrueEnergy","true_energy",150,0,3000);
  InitialiseMetric("ShowerResidual","shower_residual_dist",150,0,100);
  InitialiseMetric("ShowerLength","shower_length",50,0,200);
  InitialiseMetric("ShowerDensity","shower_density",50,0,50);
  InitialiseMetric("ShowerLengthPerp","shower_length_perp",50,0,200);
  InitialiseMetric("ShowerDensityPerp","shower_density_perp",50,0,50);
  InitialiseMetric("ShowerDensity3D","shower_density_3D",50,0,2);
  InitialiseMetric("ShowerDensityGradPerp","shower_density_grad_perp",50,-20,20);
  InitialiseMetric("ShowerDensityGrad","shower_density_grad",50,-10,10);
  InitialiseMetric("ShowerDensityRatio","shower_density_ratio",50,0,10);
  InitialiseMetric("ShowerDensityGradRatio","shower_density_grad_ratio",100,-100,100);
  InitialiseMetric("ShowerDensityGradPerpSq","shower_density_grad_perp_sq",100,0,1);
  InitialiseMetric("ShowerDensityGradSq","shower_density_grad_sq",100,0,1);
  InitialiseMetric("ShowerDensityGradRatioSq","shower_density_grad_ratio_sq",100,0,1);

  InitialiseMetric("ShowerEnergy2D","shower_energy",150,0,3000,"true_energy",150,0,3000);
  InitialiseMetric("ShowerTrueEnergy2D","true_energy",150,0,3000,"true_energy",150,0,3000);
  InitialiseMetric("ShowerResidual2D","shower_residual_dist",150,0,100);
  InitialiseMetric("ShowerLength2D","shower_length",50,0,200,"true_energy",150,0,3000);
  InitialiseMetric("ShowerDensity2D","shower_density",50,0,50,"true_energy",150,0,3000);
  InitialiseMetric("ShowerLengthPerp2D","shower_length_perp",50,0,200,"true_energy",150,0,3000);
  InitialiseMetric("ShowerDensityPerp2D","shower_density_perp",50,0,50,"true_energy",150,0,3000);
  InitialiseMetric("ShowerDensity3D2D","shower_density_3D",50,0,2,"true_energy",150,0,3000);
  InitialiseMetric("ShowerDensityGradPerp2D","shower_density_grad_perp",50,-20,20,"true_energy",150,0,3000);
  InitialiseMetric("ShowerDensityGrad2D","shower_density_grad",50,-10,10,"true_energy",150,0,3000);
  InitialiseMetric("ShowerDensityRatio2D","shower_density_ratio",50,0,10,"true_energy",150,0,3000);
  InitialiseMetric("ShowerDensityGradRatio2D","shower_density_grad_ratio",100,-100,100,"true_energy",150,0,3000);
  InitialiseMetric("ShowerDensityGradPerpSq2D","shower_density_grad_perp_sq",100,0,1,"true_energy",150,0,3000);
  InitialiseMetric("ShowerDensityGradSq2D","shower_density_grad_sq",100,0,1,"true_energy",150,0,3000);
  InitialiseMetric("ShowerDensityGradRatioSq2D","shower_density_grad_ratio_sq",100,0,1,"true_energy",150,0,3000);


  InitialiseMetric("ShowerLengthOverE","shower_length",50,0,2,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensityOverE","shower_density",50,0,0.5,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerLengthPerpOverE","shower_length_perp",50,0,0.2,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensityPerpOverE","shower_density_perp",50,0,0.5,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensity3DOverE","shower_density_3D",50,0,0.01,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensityGradPerpOverE","shower_density_grad_perp",50,-0.05,0.5,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensityGradOverE","shower_density_grad",200,-0.05,0.05,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensityRatioOverE","shower_density_ratio",50,0,1,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensityGradRatioOverE","shower_density_grad_ratio",100,-0.05,0.05,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensityGradPerpSqOverE","shower_density_grad_perp_sq",100,0,0.001,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensityGradSqOverE","shower_density_grad_sq",200,0,0.001,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensityGradRatioSqOverE","shower_density_grad_ratio_sq",100,0,0.001,"StandardOverEDaughterCutFinder");

  InitialiseMetric("ShowerLengthOverE2D","shower_length",50,0,2,"shower_energy",150,0,3000,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensityOverE2D","shower_density",50,0,0.5,"shower_energy",150,0,3000,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerLengthPerpOverE2D","shower_length_perp",50,0,0.2,"shower_energy",150,0,3000,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensityPerpOverE2D","shower_density_perp",50,0,0.5,"shower_energy",150,0,3000,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensity3DOverE2D","shower_density_3D",50,0,0.01,"shower_energy",150,0,3000,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensityGradPerpOverE2D","shower_density_grad_perp",50,-0.05,0.5,"shower_energy",150,0,3000,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensityGradOverE2D","shower_density_grad",50,-0.05,0.05,"shower_energy",150,0,3000,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensityRatioOverE2D","shower_density_ratio",50,0,1,"shower_energy",150,0,3000,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensityGradRatioOverE2D","shower_density_grad_ratio",100,-0.05,0.05,"shower_energy",150,0,3000,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensityGradPerpSqOverE2D","shower_density_grad_perp_sq",100,0,0.001,"shower_energy",150,0,3000,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensityGradSqOverE2D","shower_density_grad_sq",100,0,0.001,"shower_energy",150,0,3000,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerDensityGradRatioSqOverE2D","shower_density_grad_ratio_sq",100,0,0.001,"shower_energy",150,0,3000,"StandardOverEDaughterCutFinder");

  InitialiseMetric("ShowerDensityGradNew","shower_density_grad_new",200,-1000,1000);
  InitialiseMetric("ShowerDensityGradOverNew","shower_density_grad_new",200,0,10,"StandardOverEDaughterCutFinder");
  InitialiseMetric("ShowerOpenAngle","shower_open_angle",50,0,2);



  InitialiseMetric("dEdxWithResCut","shower_dEdx",100,0,10,"ShowerResidualCutFinder"); 
 
  InitialiseMetric("dEdxVsOneShowerCut","shower_energy",100,0,10,"shower_energy",50,0,3000,"OneShowerDaughterAnalysisCutFinder");

  //Neutrino Cut analysis
  InitialiseMetric("NumberOfShowers","number_of_showers_per_neutrino",10,0,10,"StandardNeutrinoCutFinder");

  InitialiseMetric("NumberOfNeutrinos","nu_reco_energy",10,0,10,"OneShowerSizeNeutrinoAnalysisCutFinder");


  //2D analysis
  InitialiseMetric("NumberOfShowersEnergy","shower_energy",150,0,3000,"number_of_showers_per_neutrino",10,0,10,"StandardNeutrinoDaughterCutFinder");
  InitialiseMetric("CoversionGapdEdx","shower_coversion_gap",100,0,10,"shower_dEdx",100,0,10);
  InitialiseMetric("CoversionGapVertexEnergy","shower_coversion_gap",100,0,10,"vertex_recoK",50,0,500,"StandardNeutrinoDaughterCutFinder");
  InitialiseMetric("dEdxVsEnergy","shower_dEdx",100,0,10,"shower_energy",150,0,3000,"DaughterComparision");

  //Special analyses that don't fit my stupid way.
  InitialiseMetric("NumberOfShowersIntegral","shower_energy",50,0,1000,"number_of_showers_per_neutrino",10,0,10,"StandardNumShowerAnalysis");
  //  InitialiseMetric("OneShowerdEdxAnalysis","shower_energy",50,0,3000,"shower_dEdx",100,0,10,"OneShower2DEnergyAnalysis");
  InitialiseMetric("ShowerResidualAnalysis","shower_residual_dist",100,0,50,"number_of_showers_per_neutrino",10,0,10,"ShowerResidualAnalysis");

  SetLogic("CoversionGapVertexEnergy",true,true,false);
  SetLogic("ShowerEnergy",false,false,false);
  SetLogic("ShowerTrueEnergy",false,false,false);
  SetLogic("NumberOfShowersEnergy",false,true,true);
  SetLogic("NumberOfShowers",false,false,false);

  SetLogic("ShowerLength",false,false,false);
  SetLogic("ShowerDensity",false,false,false);
  SetLogic("ShowerLengthPerp",false,false,false);
  SetLogic("ShowerDensityPerp",false,false,false);
  SetLogic("ShowerDensity3D",false,false,false);
  SetLogic("ShowerDensityGradPerp",false,false,false);
  SetLogic("ShowerDensityGrad",false,false,false);
  SetLogic("ShowerDensityRatio",false,false,false);
  SetLogic("ShowerDensityGradRatio",false,false,false);
  SetLogic("ShowerDensityGradNew",false,false,false);


  SetLogic("ShowerLength2D",false,true,false);
  SetLogic("ShowerDensity2D",false,true,false);
  SetLogic("ShowerLengthPerp2D",false,true,false);
  SetLogic("ShowerDensityPerp2D",false,true,false);
  SetLogic("ShowerDensity3D2D",false,true,false);
  SetLogic("ShowerDensityGradPerp2D",false,true,false);
  SetLogic("ShowerDensityGrad2D",false,true,false);
  SetLogic("ShowerDensityRatio2D",false,true,false);
  SetLogic("ShowerDensityGradRatio2D",false,true,false);


  //InitialiseMetric("StandardRecoEfficiency","track_lengths",50,0,200); 

  //Big Boi Full selection 
  InitialiseMetric("ShowerAllsCut","ShowerAllsCut");
  InitialiseMetric("zBDTAllsCut","MVAAllsCut",true,true,false,-999,"BDT");
  InitialiseMetric("zBDTGAllsCut","MVAAllsCut",true,true,false,-999,"BDTG");
  SetLogic("zBDT",false,false,false);  
  SetLogic("zBDTG",false,false,false);  

}

void optimiser::NueRecoOptimiser::InitialiseMetric(std::string Name, std::string branchname,int numbins, float xmin, float xmax, std::string AnalysisName, bool ApplyPOT, 
						   bool ApplyOscProb, bool fFillMVA, float MVAErrorVal, TString MVAName){
  MetricMap[Name] = optimiser::MetricHolder(Name,branchname,numbins,xmin,xmax,TotalSigPOT,TotalBKPOT,ApplyPOT,ApplyOscProb,AnalysisName,fFillMVA, MVAName,MVAErrorVal);
  if(fFillMVA){
    std::string OverF = Name + "/F";
    BackgroundTreeMVA->Branch(Name.c_str(),&MetricMap[Name].GetMVAValBackground(),OverF.c_str());
    SignalTreeMVA->Branch(Name.c_str(),&MetricMap[Name].GetMVAValSignal(),OverF.c_str());
  }
  return;
}

void optimiser::NueRecoOptimiser::InitialiseMetric(std::string Name, std::string AnalysisName, bool ApplyPOT,bool ApplyOscProb, bool fFillMVA, float MVAErrorVal, TString MVAName){
  MetricMap[Name] = optimiser::MetricHolder(Name,TotalSigPOT,TotalBKPOT,ApplyPOT,ApplyOscProb,AnalysisName,fFillMVA,MVAName,MVAErrorVal);
  if(fFillMVA){
    std::string OverF = Name + "/F";
    BackgroundTreeMVA->Branch(Name.c_str(),&MetricMap[Name].GetMVAValBackground(),OverF.c_str());
    SignalTreeMVA->Branch(Name.c_str(),&MetricMap[Name].GetMVAValSignal(),OverF.c_str());
  }
  return;
}


void optimiser::NueRecoOptimiser::InitialiseMetric(std::string Name, std::string branchname,int numbins, float xmin, float xmax, std::string metricpartner,int numbins_partner, float xmin_partner, float xmax_partner, std::string AnalysisName, bool ApplyPOT, bool ApplyOscProb, bool fFillMVA, float MVAErrorVal, TString MVAName){
  MetricMap[Name] = optimiser::MetricHolder(Name,branchname,numbins,xmin,xmax,metricpartner,numbins_partner,xmin_partner,xmax_partner,TotalSigPOT,TotalBKPOT,ApplyPOT,ApplyOscProb,AnalysisName,fFillMVA,MVAName,MVAErrorVal);

  if(fFillMVA){
    std::string OverF = Name + "/F";
    BackgroundTreeMVA->Branch(Name.c_str(),&MetricMap[Name].GetMVAValBackground(),OverF.c_str());
    SignalTreeMVA->Branch(Name.c_str(),&MetricMap[Name].GetMVAValSignal(),OverF.c_str());
  }
  return;
}

void optimiser::NueRecoOptimiser::SetLogic(std::string Name, bool ls, bool pls, bool andor){
  MetricMap[Name].SetLogic(ls,pls,andor);
  return;
}

float optimiser::NueRecoOptimiser::GetOscProb(int& iter, std::string sample, bool applyoscprob){

  if(!applyoscprob){return 1;}

  std::string branchname = "nu_osc_prob";
  std::vector<float>* oscprob;  
  int err = GetBranchPointer(branchname,sample,oscprob);
  if(err){
    std::cerr << "Osclillation probability broken" << std::endl;
    return 0;
  }

  if(oscprob->size() <= iter){
    //    std::cerr << "Osclillation probability broken" << std::endl;
    return 0;
  }
  return oscprob->at(iter);
}

void optimiser::NueRecoOptimiser::FillData(TTree* tree_signal, TTree* tree_background){

  //Fill The Background
  for (Long64_t evt = 0; evt < tree_background->GetEntries(); evt++) {
    tree_background->GetEntry(evt);
    for(auto& Metric: MetricMap){
      
      std::vector<float> vals;
      std::vector<std::pair<float,float> > TwoDvals;
      TString MVAName = Metric.second.GetMVAName();
      int err = PerformCuts(Metric.second,"background",vals,TwoDvals,MVAName);
      if(err){return;}
      for(int i=0; i<vals.size(); ++i){
	if(vals[i] == Metric.second.GetErrorVal()){continue;}
	double oscprob = GetOscProb(i,"background", Metric.second.ApplyOscProb());
	double potweight = 1;
        if(Metric.second.ApplyPOT()){potweight = Metric.second.ReturnPOT("background");}
	else{potweight = potweight/tree_background->GetEntries();}

	double oscpot = oscprob*potweight;  
	Metric.second.FillBackground(vals[i],oscpot);
      }
      for(int i=0; i<TwoDvals.size(); ++i){
	if(TwoDvals[i].first == Metric.second.GetErrorVal() || TwoDvals[i].second == Metric.second.GetErrorVal()){continue;}
	double oscprob = GetOscProb(i,"background", Metric.second.ApplyOscProb());
	double potweight = 1;
        if(Metric.second.ApplyPOT()){potweight = Metric.second.ReturnPOT("background");}
	else{potweight = potweight/tree_background->GetEntries();}

	double oscpot = oscprob*potweight;
	Metric.second.FillBackground(TwoDvals[i].first,TwoDvals[i].second,oscpot);
      }

      if(Metric.second.FillForMVA()){
	Metric.second.MVAValBackgroundVec(vals);
      }
    }
    
    int MVAVectorSize = -999;
    for(auto& Metric: MetricMap){
      
      if(!Metric.second.FillForMVA()){continue;}
      if(MVAVectorSize < Metric.second.GetMVAVectorBackgroundSize()){
	MVAVectorSize =  Metric.second.GetMVAVectorBackgroundSize();
      }
    }

    for(int neutrino=0; neutrino<MVAVectorSize; ++neutrino){
      for(auto& Metric: MetricMap){

	if(!Metric.second.FillForMVA()){continue;}
	float val = Metric.second.GetMVAValBackground(neutrino);
	Metric.second.FillMVAValBackground(val);
      }
      if(ApplyOscWht){BackgrndWeight =  GetOscProb(neutrino,"background", true);}
      if(ApplyPOTWht){BackgrndWeight *= TotalBKPOT;}
      else{BackgrndWeight /= tree_background->GetEntries();}
      BackgroundTreeMVA->Fill();
    }

  }

  //Fill The Signal
  for (Long64_t evt = 0; evt < tree_signal->GetEntries(); evt++) {
    tree_signal->GetEntry(evt);
    for(auto& Metric: MetricMap){

      std::vector<float> vals;
      std::vector<std::pair<float,float> > TwoDvals;
      TString MVAName = Metric.second.GetMVAName();
      int err = PerformCuts(Metric.second,"signal",vals,TwoDvals,MVAName);
      if(err){return;}
      for(int i=0; i<vals.size(); ++i){
        if(vals[i] == Metric.second.GetErrorVal()){continue;}
	double oscprob = GetOscProb(i,"signal",Metric.second.ApplyOscProb());
	double potweight = 1;
	if(Metric.second.ApplyPOT()){potweight = Metric.second.ReturnPOT("signal");}
	else{potweight = potweight/tree_signal->GetEntries();}

	double oscpot = oscprob*potweight;  
	Metric.second.FillSignal(vals[i],oscpot);
      }
      for(int i=0; i<TwoDvals.size(); ++i){
	if(TwoDvals[i].first == Metric.second.GetErrorVal() || TwoDvals[i].second == Metric.second.GetErrorVal()){continue;}
	double oscprob = GetOscProb(i,"signal",Metric.second.ApplyOscProb());
	double potweight = 1;  
	if(Metric.second.ApplyPOT()){potweight = Metric.second.ReturnPOT("signal");}
	else{potweight = potweight/tree_signal->GetEntries();}
	double oscpot = oscprob*potweight;
	Metric.second.FillSignal(TwoDvals[i].first,TwoDvals[i].second,oscpot);
      }
      
      if(Metric.second.FillForMVA()){
	Metric.second.FillMVAValSignalVec(vals);
      }
    }

    int MVAVectorSize = -999;
    for(auto& Metric: MetricMap){
      if(!Metric.second.FillForMVA()){continue;}
      if(MVAVectorSize < Metric.second.GetMVAVectorSignalSize()){
	MVAVectorSize = Metric.second.GetMVAVectorSignalSize();
      }
    }
    
    for(int neutrino=0; neutrino<MVAVectorSize; ++neutrino){
      for(auto& Metric: MetricMap){

	if(!Metric.second.FillForMVA()){continue;}
	float val = Metric.second.GetMVAValSignal(neutrino);
	Metric.second.FillMVAValSignal(val);
      }
      if(ApplyOscWht){SignalWeight =  GetOscProb(neutrino,"signal", true);}
      if(ApplyPOTWht){SignalWeight *= TotalSigPOT;}
      else{SignalWeight /= tree_signal->GetEntries();}

      SignalTreeMVA->Fill();
    }
  }
  return;
} 

void optimiser::NueRecoOptimiser::AnalyseData(){

  TFile* file = new TFile("CutFile.root", "RECREATE");

  //Make new directory for plots 
  for(auto const& Metric: MetricMap){
    gDirectory->mkdir((Metric.first).c_str());
  }

  //Loop over metrics and make the graphs.
  for(auto& Metric: MetricMap){

    file->cd((Metric.first).c_str());

    std::cout << "##################" << std::endl;
    std::cout << "On Metric: " << Metric.first << std::endl;
    std::cout << "##################" << std::endl;

    if(Metric.second.GetAnalysisName().find("CutFinder") != std::string::npos){
      //Perform the 1D analysis
      if(Metric.second.TotalEntries() != 0){Metric.second.Perform1DCutFinder();}

      //Perform the 2D analysis if any.
      if(Metric.second.GetPartnerMetricName() != "" 
	 && Metric.second.Total2DEntries() != 0){
	Metric.second.Perform2DCutFinder();
      }
    }
    
    //Now for the special analyses
    if(Metric.second.GetAnalysisName() == "StandardNumShowerAnalysis"){
      float totalsig = Metric.second.GetSignalHistPartner()->GetEntries()/(Metric.second.GetSignalHistPartner()->GetNbinsX());
      float totalbk  = Metric.second.GetBackgroundHistPartner()->GetEntries()/(Metric.second.GetSignalHistPartner()->GetNbinsX());
      Metric.second.MakeStandardNumShowerAnalysisGraphs(totalsig,totalbk);
    }
    if(Metric.second.GetAnalysisName() == "OneShower2DEnergyAnalysis"){

      float totalsig = Metric.second.GetSignalHistPartner()->GetEntries()/(Metric.second.GetSignalHistPartner()->GetNbinsX()*Metric.second.GetSignalHistPartner()->GetNbinsY());
      float totalbk  = Metric.second.GetBackgroundHistPartner()->GetEntries()/(Metric.second.GetSignalHistPartner()->GetNbinsX()*Metric.second.GetSignalHistPartner()->GetNbinsY());
      Metric.second.MakeStandardNumShowerAnalysisGraphs(totalsig,totalbk);
    }
    if(Metric.second.GetAnalysisName() == "ShowerResidualAnalysis"){
      float totalsig = Metric.second.GetSignalHistPartner()->GetEntries()/(Metric.second.GetSignalHistPartner()->GetNbinsX());
      float totalbk  = Metric.second.GetBackgroundHistPartner()->GetEntries()/(Metric.second.GetSignalHistPartner()->GetNbinsX());
      Metric.second.MakeStandardNumShowerAnalysisGraphs(totalsig,totalbk);
    }
    if(Metric.second.GetAnalysisName().find("Comparision") != std::string::npos){
      
      if(Metric.second.TotalEntries() != 0){
	Metric.second.Plot1D();
      }

      if(Metric.second.GetPartnerMetricName() != "" 
	 && Metric.second.Total2DEntries() != 0){
	Metric.second.Plot2D();
      }
    }

    if(Metric.second.GetAnalysisName() == "ShowerAllsCut"){
      Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram("NeutrinoRecoEAfter","signal"),Metric.second.GetHistogram("NeutrinoRecoEBefore","signal"),Metric.second.GetHistogram("NeutrinoRecoEAfterExtra","signal"));

      Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram("NeutrinoRecoEAfter","background"),Metric.second.GetHistogram("NeutrinoRecoEBefore","background"),Metric.second.GetHistogram("NeutrinoRecoEAfterExtra","background"));
      Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram("NeutrinoTrueEAfter","signal"),Metric.second.GetHistogram("NeutrinoTrueEBefore","signal"),Metric.second.GetHistogram("NeutrinoTrueEAfterExtra","signal"));
      Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram("NeutrinoTrueEAfter","background"),Metric.second.GetHistogram("NeutrinoTrueEBefore","background"),Metric.second.GetHistogram("NeutrinoTrueEAfterExtra","background"));
    }

    if(Metric.second.GetAnalysisName() == "MVAAllsCut"){

      TString MVAMethod = Metric.second.GetMVAName();

      std::string NeutrinoTrueEBeforeName = "NeutrinoTrueEBefore" + (std::string) MVAMethod; 
      std::string NeutrinoRecoEBeforeName = "NeutrinoRecoEBeforeName" +(std::string) MVAMethod;
      std::string NeutrinoRecoEBeforeExtraName = "NeutrinoRecoEBeforeExtraName" +(std::string) MVAMethod;
      std::string NeutrinoRecoEAfterName = "NeutrinoRecoEAfterName" +(std::string) MVAMethod;
      std::string NeutrinoTrueEAfterName = "NeutrinoTrueEAfterName" +(std::string) MVAMethod;
      std::string NeutrinoRecoEAfterExtraName = "NeutrinoRecoEAfterExtraName" +(std::string) MVAMethod;
      std::string NeutrinoTrueEAfterExtraName = "NeutrinoTrueEAfterExtraName" +(std::string) MVAMethod;
      

      Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram(NeutrinoRecoEAfterName,"signal"),Metric.second.GetHistogram(NeutrinoRecoEBeforeName,"signal"),Metric.second.GetHistogram(NeutrinoRecoEAfterExtraName,"signal"));
      
      Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram(NeutrinoRecoEAfterName,"background"),Metric.second.GetHistogram(NeutrinoRecoEBeforeName,"background"),Metric.second.GetHistogram(NeutrinoRecoEAfterExtraName,"background"));
      Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram(NeutrinoTrueEAfterName,"signal"),Metric.second.GetHistogram(NeutrinoTrueEBeforeName,"signal"),Metric.second.GetHistogram(NeutrinoTrueEAfterExtraName,"signal"));
      Metric.second.MakeEfficiencyPlots(Metric.second.GetHistogram(NeutrinoTrueEAfterName,"background"),Metric.second.GetHistogram(NeutrinoTrueEBeforeName,"background"),Metric.second.GetHistogram(NeutrinoTrueEAfterExtraName,"background"));
    }


  }
  
  file->Close();
  return;
}

int optimiser::NueRecoOptimiser::PerformCuts(optimiser::MetricHolder& Metric, std::string sample,std::vector<float>& vals, std::vector<std::pair<float,float> >& TwoDvals, TString& MVAName){

  int err=0;

  //Make Efficiency Hists
  //Make Resolutions Hits

  if(Metric.GetAnalysisName() == "StandardShowerCutFinder"){err = StandardDaughterAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "StandardNeutrinoCutFinder"){err = StandardNeutrinoAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "StandardNeutrinoDaughterCutFinder"){err = StandardNeutrinoDaughterAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "StandardNumShowerAnalysis"){err = StandardNumShowerAnalysis(Metric,sample,TwoDvals);}
  else if(Metric.GetAnalysisName() == "OneShower2DEnergyAnalysis"){err = OneShower2DEnergyAnalysis(Metric,sample,TwoDvals);}  
  else if(Metric.GetAnalysisName() == "ShowerResidualAnalysis"){ err = ShowerResidualAnalysis(Metric,sample,TwoDvals);}
  else if(Metric.GetAnalysisName() == "DaughterComparision"){err = StandardDaughterAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "NeutrinoComparision"){err = StandardNeutrinoAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "NeutrinoDaughterComparision"){err = StandardNeutrinoDaughterAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "ShowerResidualCutFinder"){err = ShowerResidualCut(Metric,sample,vals);}
  else if(Metric.GetAnalysisName() == "ShowerAllsCut"){err = ShowerAllsCut(Metric,sample);}
  else if(Metric.GetAnalysisName() == "OneShowerSizeNeutrinoAnalysisCutFinder"){err =StandardSizeNeutrinoAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "StandardOverEDaughterCutFinder"){err = StandardOverEDaughterAnalysis(Metric,sample,vals,TwoDvals);}
  else if(Metric.GetAnalysisName() == "ShowerResidualNumShowersCutFinder"){err = ShowerResidualNumShowers(Metric,sample,vals);}
  else if(Metric.GetAnalysisName() == "MVAAnalysisCutFinder"){err = MVAAnalysis(Metric,sample,vals,MVAName);}
  else if(Metric.GetAnalysisName() == "MVAAllsCut"){err = MVAAnalysisAllsCut(Metric,sample,MVAName);}
  else if(Metric.GetAnalysisName() == "OneShower1DCutFinder"){err = OneShower1DCutAnalysis(Metric,sample,vals);}
  return err;
}

template <class T>
int optimiser::NueRecoOptimiser::GetBranchPointer(std::string& branchname, std::string& sample, T* &vals){

  optimiser::BranchTypeBase* branchbase;

  //Decide if it is a signal or background branch
  if(sample == "signal"){
    if(SignalTree.Tree.find(branchname) == SignalTree.Tree.end()){
      std::cerr << "Signal Tree branch doesn't exist. You probably mistypes the name in the metric initialser" << std::endl;
      return 1;
    }
    branchbase = SignalTree.Tree[branchname];
  }
  else{
    if(BackgroundTree.Tree.find(branchname) == BackgroundTree.Tree.end()){
      std::cerr << "Background Tree branch doesn't exist. You probably mistypes the name in the metric initialser" << std::endl;
      return 1;
    }
    branchbase = BackgroundTree.Tree[branchname];
  }

  //Get the branch.
  optimiser::BranchType<T>* branch  =  dynamic_cast<optimiser::BranchType<T>*>(branchbase);
  if(branch == NULL){
    std::cerr << "could not find branch: " << branchname << " with type given " << typeid(*vals).name() << std::endl;
    return 1;
  }
  //Actual values.
  vals = branch->GetPointer();
  return 0;
}

//#######################################
//######### Standard Daughter ###########
//#######################################


int optimiser::NueRecoOptimiser::StandardDaughterAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals,  std::vector<std::pair<float,float> >& TwoDvals){

  int err = StandardDaughter1DAnalysis(Metric,sample,vals);
  if(err){
    std::cerr << "Error in 1D daughter analysis and bailing" << std::endl;
    return 1;
  }

  if(Metric.GetPartnerMetricName() == ""){return 0;}

  err =  StandardDaughter2DAnalysis(Metric,sample,TwoDvals);
  if(err){
    std::cerr << "Error in 2D daughter analysis and bailing" << std::endl;
    return 1;
  }
  return 0;
}


int optimiser::NueRecoOptimiser::StandardDaughter1DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;  
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Loop over the neutrinos
  for(auto const& neutrino: *branch_vals){
    //Add the largest shower in the event. The metric are ordered this way
    if(neutrino.size() > 0){
      float maxval = -999;
      for(auto const& val: neutrino){
	if(maxval < val){
	  maxval = val;
	}
      }
      vals.push_back(maxval);
    }
    else{
      vals.push_back(Metric.GetErrorVal());
    }
  }
  return 0;
}

int optimiser::NueRecoOptimiser::StandardDaughter2DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername = Metric.GetPartnerMetricName();
  std::vector<std::vector<float> >* partner_vals;
  err = GetBranchPointer(partnername,sample,partner_vals);
  if(err){return 1;}

  if(partner_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for daughter and daughter info. Please Check" << std::endl;
  }

  //Loop over the neutrinos
  for(int iter=0; iter<branch_vals->size(); ++iter){
    //Add the largest shower in the event. The metric are ordered this way
    if(branch_vals->at(iter).size() > 0 && partner_vals->at(iter).size() > 0){
      std::pair<float,float> vals = std::make_pair(branch_vals->at(iter).at(0),partner_vals->at(iter).at(0));
      TwoDvals.push_back(vals);
    }
    else{
      std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
      TwoDvals.push_back(vals);
    }

  }
  return 0;
}

//#############################################
//######### Standard OverE Daughter ###########
//#############################################


int optimiser::NueRecoOptimiser::StandardOverEDaughterAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals, std::vector<std::pair<float,float> >& TwoDvals){

  int err = StandardOverEDaughter1DAnalysis(Metric,sample,vals);
  if(err){
    std::cerr << "Error in 1D daughter analysis and bailing" << std::endl;
    return 1;
  }

  if(Metric.GetPartnerMetricName() == ""){return 0;}

  err =  StandardOverEDaughter2DAnalysis(Metric,sample,TwoDvals);
  if(err){
    std::cerr << "Error in 2D daughter analysis and bailing" << std::endl;
    return 1;
  }
  return 0;
}


int optimiser::NueRecoOptimiser::StandardOverEDaughter1DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;  
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the Energy
  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* shower_energy_vals;
  err = GetBranchPointer(shower_energy,sample,shower_energy_vals);
  if(err){return 1;}

  //Loop over the neutrinos
  for(int iter=0; iter<branch_vals->size(); ++iter){
    //Add the largest shower in the event. The metric are ordered this way
    if(branch_vals->at(iter).size() > 0){
      vals.push_back(branch_vals->at(iter).at(0)/shower_energy_vals->at(iter).at(0));
    }
    else{
      vals.push_back(Metric.GetErrorVal());
    }
  }
  return 0;
}

int optimiser::NueRecoOptimiser::StandardOverEDaughter2DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername = Metric.GetPartnerMetricName();
  std::vector<std::vector<float> >* partner_vals;
  err = GetBranchPointer(partnername,sample,partner_vals);
  if(err){return 1;}

  //Get the Energy
  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* shower_energy_vals;
  err = GetBranchPointer(shower_energy,sample,shower_energy_vals);
  if(err){return 1;}


  if(partner_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for daughter and daughter info. Please Check" << std::endl;
  }

  //Loop over the neutrinos
  for(int iter=0; iter<branch_vals->size(); ++iter){
    //Add the largest shower in the event. The metric are ordered this way
    if(branch_vals->at(iter).size() > 0 && partner_vals->at(iter).size() > 0){
      std::pair<float,float> vals = std::make_pair(branch_vals->at(iter).at(0)/shower_energy_vals->at(iter).at(0),partner_vals->at(iter).at(0));
      TwoDvals.push_back(vals);
    }
    else{
      std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());

      TwoDvals.push_back(vals);
    }
  }
  return 0;
}


//#######################################
//######### Standard Neutrino ###########
//#######################################

int optimiser::NueRecoOptimiser::StandardNeutrinoAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals, std::vector<std::pair<float,float> >& TwoDvals){

  int err = StandardNeutrino1DAnalysis(Metric,sample,vals);
  if(err){
    std::cerr << "Error or in 1D neutrino analysis and bailing" << std::endl;
    return 1;
  }

  if(Metric.GetPartnerMetricName() == ""){return 0;}

  err =  StandardNeutrino2DAnalysis(Metric,sample,TwoDvals);
  if(err){
    std::cerr << "Error or in 2D neutrino analysis and bailing" << std::endl;
    return 1;
  }
  return 0;
}


int optimiser::NueRecoOptimiser::StandardNeutrino1DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  //Get the values
  std::string branchname = Metric.GetBranchName();
  std::vector<float>* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}
  
  //Loop over the neutrinos
  for(auto const& neutrino: *branch_vals){
    vals.push_back(neutrino);
  }
  return 0;
}

int optimiser::NueRecoOptimiser::StandardNeutrino2DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<float>* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername = Metric.GetPartnerMetricName();
  std::vector<float>* partner_vals;
  err = GetBranchPointer(partnername,sample,partner_vals);
  if(err){return 1;}

  if(partner_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for neutrino and neutrino info. Please Check" << std::endl;
  }


  //Loop over the neutrinos
  for(int iter=0; iter<branch_vals->size(); ++iter){
    std::pair<float,float> vals = std::make_pair(branch_vals->at(iter),partner_vals->at(iter));
    TwoDvals.push_back(vals);
  }
  return 0;
}

//############################################
//######### Standard Size Neutrino ###########
//############################################

int optimiser::NueRecoOptimiser::StandardSizeNeutrinoAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals, std::vector<std::pair<float,float> >& TwoDvals){

  int err = StandardSizeNeutrino1DAnalysis(Metric,sample,vals);
  if(err){
    std::cerr << "Error or in 1D neutrino analysis and bailing" << std::endl;
    return 1;
  }

  if(Metric.GetPartnerMetricName() == ""){return 0;}

  err =  StandardSizeNeutrino2DAnalysis(Metric,sample,TwoDvals);
  if(err){
    std::cerr << "Error or in 2D neutrino analysis and bailing" << std::endl;
    return 1;
  }
  return 0;
}


int optimiser::NueRecoOptimiser::StandardSizeNeutrino1DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  //Get the values
  std::string branchname = Metric.GetBranchName();
  std::vector<float>* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Loop over the neutrinos
  vals.push_back(branch_vals->size());

  return 0;
}

int optimiser::NueRecoOptimiser::StandardSizeNeutrino2DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<float>* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername = Metric.GetPartnerMetricName();
  std::vector<float>* partner_vals;
  err = GetBranchPointer(partnername,sample,partner_vals);
  if(err){return 1;}

  //Loop over the neutrinos
  std::pair<float,float> vals = std::make_pair(branch_vals->size(),partner_vals->size());
  TwoDvals.push_back(vals);

  return 0;
}


//##################################################
//######### Standard Neutrino & Daughter ###########
//##################################################

int optimiser::NueRecoOptimiser::StandardNeutrinoDaughterAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals, std::vector<std::pair<float,float> >& TwoDvals){

  if(Metric.GetPartnerMetricName() == ""){return 1;}

  int err =  StandardNeutrinoDaughter2DAnalysis(Metric,sample,TwoDvals);
  if(err){
    std::cerr << "Error or in 2D daughter+neutrino analysis and bailing" << std::endl;
    return 1;
  }
  return 0;
}


int optimiser::NueRecoOptimiser::StandardNeutrinoDaughter2DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername = Metric.GetPartnerMetricName();
  std::vector<float>* partner_vals;
  err = GetBranchPointer(partnername,sample,partner_vals);
  if(err){return 1;}

  if(partner_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for neutrino and neutrino daughter info. Please Check" << std::endl;
    return 1;
  }

  //Loop over the neutrinos
  for(int iter=0; iter<branch_vals->size(); ++iter){
    //Add the largest shower in the event. The metric are ordered this way
    if(branch_vals->at(iter).size() > 0){
      std::pair<float,float> vals = std::make_pair(branch_vals->at(iter).at(0),partner_vals->at(iter));
      TwoDvals.push_back(vals);
    }
    else{
      std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());

      TwoDvals.push_back(vals);
    }
  }
  return 0;
}

//######################
//### Num Shower Ana ###
//######################

int optimiser::NueRecoOptimiser::StandardNumShowerAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  int numshowerbins = Metric.GetSignalHist()->GetNbinsX();
  
  //Loop over the neutrinos
  for(int bin=0; bin<numshowerbins; ++bin){
  
    float binenergy = ((TAxis*)Metric.GetSignalHist()->GetXaxis())->GetBinCenter(bin);
  
    for(auto const& neutrino: *branch_vals){
      int  numshowers = 0;
      for(auto const& shower_energy: neutrino){
	if(shower_energy > binenergy){
	  ++numshowers;
	}
	// if(shower_energy < binenergy){
	//   break;
	// }
	// ++numshowers;
      }
      std::pair<float,float> vals = std::make_pair(binenergy,numshowers);
      TwoDvals.push_back(vals);
    }
  }
  return 0;
}

//######################
//### Num Shower Ana ###
//######################

int optimiser::NueRecoOptimiser::OneShower1DCutAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  float fEnergyCut = 210;

  //Get the shower energy 
  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* shower_energy_vals;
  int err = GetBranchPointer(shower_energy,sample,shower_energy_vals);
  if(err){return 1;}
  
  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  for(int iter=0; iter<branch_vals->size(); ++iter){

    std::vector<float> shower_energy_vals_neutrino = shower_energy_vals->at(iter);
    int  numshowers = 0;
    for(auto const& shower_energy: shower_energy_vals_neutrino){
      if(shower_energy > fEnergyCut){
        ++numshowers;
      }
    }
    
    if(numshowers != 1){
      vals.push_back(Metric.GetErrorVal());
      continue;
    }
  
    //Add the largest shower in the event. The metric are ordered this way
    float maxmetric = -999;
    for(auto const& metric: branch_vals->at(iter)){
      if(metric > maxmetric){
	maxmetric = metric;
      }
    }
     
    if(maxmetric == -999){
      vals.push_back(Metric.GetErrorVal());
      continue;
    }
    vals.push_back(maxmetric);
    
    //    if(branch_vals->at(iter).size() > 0){
    //      vals.push_back(branch_vals->at(iter).at(0));
    //    }
  }
  return 0;



}

int optimiser::NueRecoOptimiser::OneShower2DEnergyAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Not optimised use sparingly

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  //Get the partner metric branch
  std::string partnername = Metric.GetPartnerMetricName();
  std::vector<std::vector<float> >* partner_vals;
  err = GetBranchPointer(partnername,sample,partner_vals);
  if(err){return 1;}

  if(partner_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for neutrino and neutrino daughter info. Please Check" << std::endl;
  }

  int numshowerbins = Metric.GetSignalHist()->GetNbinsX();
  int numshowerbins_partner = Metric.GetSignalHistPartner()->GetNbinsY();
  
  //Loop over the neutrinos
  for(int bin_p=0; bin_p<numshowerbins_partner; ++bin_p){
    float fCut = ((TAxis*)Metric.GetSignalHistPartner()->GetYaxis())->GetBinCenter(bin_p);
    for(int bin=0; bin<numshowerbins; ++bin){
  
      float binenergy = ((TAxis*)Metric.GetSignalHist()->GetXaxis())->GetBinCenter(bin);
 
      for(int neutrino_iter=0; neutrino_iter<branch_vals->size(); ++neutrino_iter){
      
	//No shower no cut continue;
	if(branch_vals->at(neutrino_iter).size() == 0){
	  std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	  TwoDvals.push_back(vals);
	  continue;
	}

	//Count the number of showers.
	int numshowers = 0;
	int shower_iter_to_ana = -999;
	for(int shower_iter=0; shower_iter<branch_vals->at(neutrino_iter).size(); ++shower_iter){
	  float shower_energy = branch_vals->at(neutrino_iter).at(shower_iter);
	  if(shower_energy > binenergy){
	    ++numshowers;
	    shower_iter_to_ana = shower_iter;
	  }
	}

	//Only Keep one shower events
	if(numshowers > 1){
	  std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	  TwoDvals.push_back(vals);
	  continue;
	}

	//Only analyse match showers
	if(shower_iter_to_ana == -999 || partner_vals->at(neutrino_iter).size() < shower_iter_to_ana){
	  std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	  TwoDvals.push_back(vals);
	  continue;
	}

	//Cut is unfortately done here 
	if(partner_vals->at(neutrino_iter).at(shower_iter_to_ana) > fCut){
	  std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	  TwoDvals.push_back(vals);
	  continue;
	}
      
	//Add the highest energy shower and its counter part 
	std::pair<float,float> vals = std::make_pair(binenergy,fCut);
	TwoDvals.push_back(vals);
      }
    }
  }
  return 0;
}


//######################
//### Shower Res Ana ###
//######################

int optimiser::NueRecoOptimiser::ShowerResidualAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  int numshowerbins = Metric.GetSignalHist()->GetNbinsX();

  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* energy_vals;
  err = GetBranchPointer(shower_energy,sample,energy_vals);
  if(err){return 1;}

  if(energy_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for neutrino and neutrino daughter info. Please Check" << std::endl;
  }

  
  //Loop over the neutrinos
  for(int bin=0; bin<numshowerbins; ++bin){
  
    float binresidual = ((TAxis*)Metric.GetSignalHist()->GetXaxis())->GetBinCenter(bin);
  
    for(int n_iter=0; n_iter<branch_vals->size(); ++n_iter){
      int  numshowers = 0;
      
      if(energy_vals->at(n_iter).size() != branch_vals->at(n_iter).size()){
	std::cerr << "energy_vals  size: " << energy_vals->at(n_iter).size() << " branch_vals->at(n_iter).size(): " << branch_vals->at(n_iter).size() << std::endl;
	std::pair<float,float> vals = std::make_pair(Metric.GetErrorVal(),Metric.GetErrorVal());
	TwoDvals.push_back(vals);
	continue;
      }

      
      for(int s_iter=0; s_iter<branch_vals->at(n_iter).size(); ++s_iter){
	//    for(auto const& neutrino: *branch_vals){
	//      for(auto const& shower_residual: neutrino){
	//if(shower_residual > binresidual){ //|| shower_residual<1e-5){
	if(branch_vals->at(n_iter).at(s_iter) > binresidual && energy_vals->at(n_iter).at(s_iter) > 50){
	  ++numshowers;
	}
      }
      
      //Only Keep one shower events
      // if(numshowers > 1){
      // 	std::pair<float,float> vals = std::make_pair(-999,-999);
      // 	TwoDvals.push_back(vals);
      // 	continue;
      // }
      
      std::pair<float,float> vals = std::make_pair(binresidual,numshowers);
      TwoDvals.push_back(vals);
    }
  }
  return 0;
}

int optimiser::NueRecoOptimiser::ShowerResidualCut(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  int numshowerbins = Metric.GetSignalHist()->GetNbinsX();

  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* energy_vals;
  err = GetBranchPointer(shower_energy,sample,energy_vals);
  if(err){return 1;}

  std::string shower_residual_dist = "shower_residual_dist";
  std::vector<std::vector<float> >* shower_residual_vals;
  err = GetBranchPointer(shower_residual_dist,sample,shower_residual_vals);
  if(err){return 1;}


  if(energy_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for neutrino and neutrino daughter info. Please Check" << std::endl;
  }

  
  //Loop over the neutrinos
  for(int n_iter=0; n_iter<branch_vals->size(); ++n_iter){
    int  numshowers = 0;
    
    if(energy_vals->at(n_iter).size() != branch_vals->at(n_iter).size()){
      std::cerr << "energy_vals  size: " << energy_vals->at(n_iter).size() << " branch_vals->at(n_iter).size(): " << branch_vals->at(n_iter).size() << std::endl;
      vals.push_back(Metric.GetErrorVal());
      continue;
    }
    
    if(branch_vals->at(n_iter).size() == 0){continue;}

    
    for(int s_iter=0; s_iter<branch_vals->at(n_iter).size(); ++s_iter){
      if(shower_residual_vals->at(n_iter).at(s_iter) > 1.25 && energy_vals->at(n_iter).at(s_iter) > 50){
	++numshowers;
      }
    }
    
    //Only Keep one shower events
    if(numshowers > 1){
      vals.push_back(Metric.GetErrorVal());
      continue;
    }
    
    vals.push_back(branch_vals->at(n_iter).at(0));
  }

return 0;
}

int optimiser::NueRecoOptimiser::ShowerResidualNumShowers(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  std::vector<std::vector<float> >* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}

  int numshowerbins = Metric.GetSignalHist()->GetNbinsX();

  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* energy_vals;
  err = GetBranchPointer(shower_energy,sample,energy_vals);
  if(err){return 1;}

  std::string shower_residual_dist = "shower_residual_dist";
  std::vector<std::vector<float> >* shower_residual_vals;
  err = GetBranchPointer(shower_residual_dist,sample,shower_residual_vals);
  if(err){return 1;}


  if(energy_vals->size() != branch_vals->size()){
    std::cerr << "expecting same vector length for neutrino and neutrino daughter info. Please Check" << std::endl;
  }

  
  //Loop over the neutrinos
  for(int n_iter=0; n_iter<branch_vals->size(); ++n_iter){
    int  numshowers = 0;
    
    if(energy_vals->at(n_iter).size() != branch_vals->at(n_iter).size()){
      std::cerr << "energy_vals  size: " << energy_vals->at(n_iter).size() << " branch_vals->at(n_iter).size(): " << branch_vals->at(n_iter).size() << std::endl;
      vals.push_back(Metric.GetErrorVal());
      continue;
    }
    
    if(branch_vals->at(n_iter).size() == 0){
      vals.push_back(Metric.GetErrorVal());
      continue;
    }

    
    for(int s_iter=0; s_iter<branch_vals->at(n_iter).size(); ++s_iter){
      if(shower_residual_vals->at(n_iter).at(s_iter) > 1.25 && energy_vals->at(n_iter).at(s_iter) > 50){
	++numshowers;
      }
    }
    
    vals.push_back(numshowers);
  }
return 0;
}



int optimiser::NueRecoOptimiser::ShowerAllsCut(optimiser::MetricHolder& Metric, std::string& sample){

  //Cuts 
  const float fResidualCutVal          = 1.25;
  const float fdEdxCutVal              = 3.1;
  const float fShowerTrackLengthCutVal = 40; 
  const float fConversionGapVal        = 2.8;
  const float fSingleShowerEnergyVal   = 210;
  const float fMinEnergyCutVal         = 70;

  //Bools to swithc cut on 
  const bool fResidualCut          = true;
  const bool fSingleShowerCut      = false;
  const bool fdEdxCut              = true;
  const bool fShowerTrackLengthCut = true;
  const bool fConversionGap        = true;
  const bool fMinEnergyCut         = true; 

  if(!Metric.CheckHistogram("NeutrinoTrueEBefore",sample)){Metric.SetHistogram("NeutrinoTrueEBefore",30,0,3000,sample);}
  if(!Metric.CheckHistogram("NeutrinoRecoEBefore",sample)){Metric.SetHistogram("NeutrinoRecoEBefore",30,0,3000,sample);}
  if(!Metric.CheckHistogram("NeutrinoRecoEBeforeExtra",sample)){Metric.SetHistogram("NeutrinoRecoEBeforeExtra",30,0,3000,sample);}
  if(!Metric.CheckHistogram("NeutrinoRecoEAfter",sample)){Metric.SetHistogram("NeutrinoRecoEAfter",30,0,3000,sample);}
  if(!Metric.CheckHistogram("NeutrinoTrueEAfter",sample)){Metric.SetHistogram("NeutrinoTrueEAfter",30,0,3000,sample);}
  if(!Metric.CheckHistogram("NeutrinoRecoEAfterExtra",sample)){Metric.SetHistogram("NeutrinoRecoEAfterExtra",30,0,3000,sample);}
  if(!Metric.CheckHistogram("NeutrinoTrueEAfterExtra",sample)){Metric.SetHistogram("NeutrinoTrueEAfterExtra",30,0,3000,sample);}

  //Get the shower energy 
  std::string shower_energy = "shower_energy";
  std::vector<std::vector<float> >* shower_energy_vals;
  int err = GetBranchPointer(shower_energy,sample,shower_energy_vals);
  if(err){return 1;}

  //Get the residual distance 
  std::string shower_residual_dist = "shower_residual_dist";
  std::vector<std::vector<float> >* shower_residual_vals;
  err = GetBranchPointer(shower_residual_dist,sample,shower_residual_vals);
  if(err){return 1;}

  //Get the dEdx
  std::string shower_dEdx = "shower_dEdx";
  std::vector<std::vector<float> >* shower_dEdx_vals;
  err = GetBranchPointer(shower_dEdx,sample,shower_dEdx_vals);
  if(err){return 1;}

  //Get the conversion distance
  std::string shower_coversion_gap = "shower_coversion_gap";
  std::vector<std::vector<float> >* shower_coversion_gap_vals;
  err = GetBranchPointer(shower_coversion_gap,sample,shower_coversion_gap_vals);
  if(err){return 1;}

  //Get the max track length 
  std::string track_lengths = "track_lengths";
  std::vector<std::vector<float> >* shower_track_length_vals;
  err = GetBranchPointer(track_lengths,sample,shower_track_length_vals);
  if(err){return 1;}

  //Get the nu reco energy  
  std::string nu_reco_energy = "nu_reco_energy";
  std::vector<float>* nu_reco_energy_vals;
  err = GetBranchPointer(nu_reco_energy,sample,nu_reco_energy_vals);
  if(err){return 1;}

  //Get the nu truth energy  
  std::string nu_truth_energy = "nu_truth_energy";
  std::vector<float> * nu_truth_energy_vals;
  err = GetBranchPointer(nu_truth_energy,sample,nu_truth_energy_vals);
  if(err){return 1;}

  //Get the nu truth energy  
  std::string nu_E = "nu_E";
  std::vector<float> * nu_E_vals;
  err = GetBranchPointer(nu_E,sample,nu_E_vals);
  if(err){return 1;}


  //Loop over the neutrinos and get their true energy
  for(int n_iter=0; n_iter<nu_E_vals->size(); ++n_iter){
    Metric.GetHistogram("NeutrinoTrueEBefore",sample)->Fill(nu_truth_energy_vals->at(n_iter)*1000);
  }

  //Loop over the reco neutrinos and idnetify the reco neutrino with the most energy
  std::map<float,std::vector<float> > NeutrinoE;
  std::map<float,float > MaxNeutrinoE;
  std::map<float, int >  MaxNeutrinoE_iter;
  for(int n_iter=0; n_iter<nu_E_vals->size(); ++n_iter){
    MaxNeutrinoE[nu_E_vals->at(n_iter)]      = -99999;
    MaxNeutrinoE_iter[nu_E_vals->at(n_iter)] = -99999;   
  }

  for(int n_iter=0; n_iter<nu_reco_energy_vals->size(); ++n_iter){
    //check if the true neutrino is already be identified 
    if(MaxNeutrinoE[nu_truth_energy_vals->at(n_iter)] < nu_reco_energy_vals->at(n_iter)){
      MaxNeutrinoE[nu_truth_energy_vals->at(n_iter)] = nu_reco_energy_vals->at(n_iter);
      MaxNeutrinoE_iter[nu_truth_energy_vals->at(n_iter)] = n_iter;
    }
  }
  //Add the max reco values to the efficiency and the others to the additional graph
  for(int n_iter=0; n_iter<nu_reco_energy_vals->size(); ++n_iter){
    if(MaxNeutrinoE[nu_truth_energy_vals->at(n_iter)] == nu_reco_energy_vals->at(n_iter)){
      Metric.GetHistogram("NeutrinoRecoEBefore",sample)->Fill(nu_reco_energy_vals->at(n_iter));
    }
    else{
      Metric.GetHistogram("NeutrinoRecoEBeforeExtra",sample)->Fill(nu_reco_energy_vals->at(n_iter)); 
    }
  }

  //Re initialise the maps
  for(int n_iter=0; n_iter<nu_E_vals->size(); ++n_iter){
    MaxNeutrinoE[nu_E_vals->at(n_iter)]      = -99999;
    MaxNeutrinoE_iter[nu_E_vals->at(n_iter)] = -99999;   
  }


  //Loop over the neutrinos
  for(int n_iter=0; n_iter<nu_reco_energy_vals->size(); ++n_iter){
    
    if(shower_energy_vals->at(n_iter).size() == 0){
      continue;
    }
    
    //Residual Cut
    if(fResidualCut){
      int  numshowers = 0;
      for(int s_iter=0; s_iter<shower_residual_vals->at(n_iter).size(); ++s_iter){
	if(shower_residual_vals->at(n_iter).at(s_iter) > fResidualCut && shower_energy_vals->at(n_iter).at(s_iter) > 25){
	  ++numshowers;
	}
      }
      //Only Keep one shower events
      if(numshowers > 1 && fResidualCut){
	continue;
      }
    }
    else if(fSingleShowerCut){
      int  numshowers = 0;
      for(int s_iter=0; s_iter<shower_residual_vals->at(n_iter).size(); ++s_iter){
	if(shower_energy_vals->at(n_iter).at(s_iter) > fSingleShowerEnergyVal){
	  ++numshowers;
	}
      }
      //Only Keep one shower events
      if(numshowers != 1 && fSingleShowerCut){
	continue;
      }
    }
  	  
    //Get the biggest shower 
    float MaxEnergy = -99999;
    int MaxEnergy_iter = -999;
    for(int s_iter=0; s_iter<shower_energy_vals->at(n_iter).size(); ++s_iter){
      if(shower_energy_vals->at(n_iter).at(s_iter) > MaxEnergy){
	MaxEnergy      = shower_energy_vals->at(n_iter).at(s_iter);
	MaxEnergy_iter = s_iter;
      }
    }
    //If there is no shower continue;
    if(MaxEnergy_iter == -999){continue;}


    //Energy Cut
    if(shower_energy_vals->at(n_iter).at(MaxEnergy_iter) < fMinEnergyCutVal && fMinEnergyCut){
      continue;
    } 

    //dEdx Cut
    if(shower_dEdx_vals->at(n_iter).at(MaxEnergy_iter) > fdEdxCutVal && fdEdxCut){
      continue;
    }

    //Conversion Gap Cut
    if(shower_coversion_gap_vals->at(n_iter).at(MaxEnergy_iter) > fConversionGapVal && fConversionGap){
      continue;
    }

    //Track length cut
    float max_track_length = -99999;
    for(auto const& shower_track_length_val: shower_track_length_vals->at(n_iter)){
      if(max_track_length < shower_track_length_val){
	max_track_length = shower_track_length_val;
      }
    }
    if(max_track_length > fShowerTrackLengthCutVal && fShowerTrackLengthCut){continue;}
    
    //Save the Neutrino E
    NeutrinoE[nu_truth_energy_vals->at(n_iter)].push_back(nu_reco_energy_vals->at(n_iter));
  }  


  //Identify if this neutrino is the largest
  for(auto const& trueNeutrino: NeutrinoE){
    float MaxNeutrinoEval      = -99999;
    int   MaxNeutrinoEval_iter = -99999; 
    //Identify the remianing reco neutrino with the max energy
    for(int n_iter=0; n_iter<(trueNeutrino.second.size()); ++n_iter){
      if(MaxNeutrinoEval < trueNeutrino.second.at(n_iter)){
	MaxNeutrinoEval = trueNeutrino.second.at(n_iter);
	MaxNeutrinoEval_iter = n_iter;
      }
    }
    for(int n_iter=0; n_iter<(trueNeutrino.second.size()); ++n_iter){ 
      //Add the biggest to account for the efficiency.
      if(n_iter == MaxNeutrinoEval_iter){
	Metric.GetHistogram("NeutrinoRecoEAfter",sample)->Fill(trueNeutrino.second.at(n_iter));
	Metric.GetHistogram("NeutrinoTrueEAfter",sample)->Fill(trueNeutrino.first*1000);
      }
      //Add the extra reco
      else{
	Metric.GetHistogram("NeutrinoRecoEAfterExtra",sample)->Fill(trueNeutrino.second.at(n_iter));
	Metric.GetHistogram("NeutrinoTrueEAfterExtra",sample)->Fill(trueNeutrino.first*1000);
      }
    }
  }
return 0;
}

int optimiser::NueRecoOptimiser::TrainMVA(){

  MVAFile->cd();
  SignalTreeMVA->Write();
  BackgroundTreeMVA->Write();

  std::cout << "sig entries: " << SignalTreeMVA->GetEntries() << " bk entries: " << BackgroundTreeMVA->GetEntries() << std::endl;

  //Load the libaries 
  TMVA::Tools::Instance();

  //Create the output label
  TString outfileName( "TMVA.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  //Create the factory 
  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
					      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

  TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");

  //Define the variables.
  for(auto& Metric: MetricMap){
    if(!Metric.second.FillForMVA()){continue;}
    dataloader->AddVariable(Metric.first.c_str(),Metric.first.c_str(),"units",'F');
  }

  double signalWeight = 1;
  double backgroundWeight = 1;
  
  // You can add an arbitrary number of signal or background trees
  dataloader->AddSignalTree    (SignalTreeMVA,     signalWeight);
  dataloader->AddBackgroundTree(BackgroundTreeMVA, backgroundWeight);

  //event weights 
  dataloader->SetSignalWeightExpression("weight");
  dataloader->SetBackgroundWeightExpression("weight");


  TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
  TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

  dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
					  "nTrain_Signal=5000:nTrain_Background=5000:SplitMode=Random:NormMode=NumEvents:!V" );

  //State the training methods. 

  //Standard Cut Based
  factory->BookMethod( dataloader, TMVA::Types::kCuts, "Cuts",
		       "!H:!V:FitMethod=MC:EffSel:SampleSize=10000:VarProp=FSmart" );
  //Some whack deep learning
  //  factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

  //Mr BDT
  factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
		       "!H:!V:NTrees=500:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
  factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTSimp",
		       "!H:!V:NTrees=500:MaxDepth=2" );



  factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
		       "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2");
  
  // Now you can tell the factory to train, test, and evaluate the MVAs
  //
  // Train MVAs using the set of training events
  factory->TrainAllMethods();
  // Evaluate all MVAs using the set of test events
  factory->TestAllMethods();
  // Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();
  // --------------------------------------------------------------
  // Save the output
  outputFile->Close();
  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;
  delete factory;
  delete dataloader;
  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );
  return 0;
}

int optimiser::NueRecoOptimiser::SetupMVA(){

  reader = new TMVA::Reader( "!Color:!Silent" );

  //Define the variables.
  for(auto& Metric: MetricMap){
    if(!Metric.second.FillForMVA()){continue;}
    std::cout << "adding variable: " << Metric.first << std::endl;
    reader->AddVariable(Metric.first.c_str(),&Metric.second.GetMVAVal());
  }
  
  std::vector<TString> methods = {"BDT","BDTG"};
  
  for(auto const method: methods){
    TString dir    = "dataset/weights/";
    TString prefix = "TMVAClassification";
    TString weightfile = dir + prefix + TString("_") + method + TString(".weights.xml");
    reader->BookMVA( method, weightfile);
  }

  return 0; 
}

int optimiser::NueRecoOptimiser::MVAAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals, TString& MVAMethod){

  //Get the values
  std::string branchname = Metric.GetBranchName();
  std::vector<float>* branch_vals;
  int err = GetBranchPointer(branchname,sample,branch_vals);
  if(err){return 1;}
  
  //Loop over the neutrinos
  vals.push_back(reader->EvaluateMVA(MVAMethod));

  return 0;
}


int optimiser::NueRecoOptimiser::MVAAnalysisAllsCut(optimiser::MetricHolder& Metric, std::string& sample, TString& MVAMethod){

  //Get the nu reco energy  
  std::string nu_reco_energy = "nu_reco_energy";
  std::vector<float>* nu_reco_energy_vals;
  int err = GetBranchPointer(nu_reco_energy,sample,nu_reco_energy_vals);
  if(err){return 1;}

  //Get the nu truth energy  
  std::string nu_truth_energy = "nu_truth_energy";
  std::vector<float> * nu_truth_energy_vals;
  err = GetBranchPointer(nu_truth_energy,sample,nu_truth_energy_vals);
  if(err){return 1;}

  //Get the nu truth energy  
  std::string nu_E = "nu_E";
  std::vector<float> * nu_E_vals;
  err = GetBranchPointer(nu_E,sample,nu_E_vals);
  if(err){return 1;}

  std::map<float,std::vector<float> > NeutrinoE;
  std::map<float,float > MaxNeutrinoE;
  std::map<float, int >  MaxNeutrinoE_iter;
  
  
  //Re initialise the maps
  for(int n_iter=0; n_iter<nu_E_vals->size(); ++n_iter){
    MaxNeutrinoE[nu_E_vals->at(n_iter)]      = -99999;
    MaxNeutrinoE_iter[nu_E_vals->at(n_iter)] = -99999;   
  }

  std::string NeutrinoTrueEBeforeName = "NeutrinoTrueEBefore" + (std::string) MVAMethod; 
  std::string NeutrinoRecoEBeforeName = "NeutrinoRecoEBeforeName" +(std::string) MVAMethod;
  std::string NeutrinoRecoEBeforeExtraName = "NeutrinoRecoEBeforeExtraName" +(std::string) MVAMethod;
  std::string NeutrinoRecoEAfterName = "NeutrinoRecoEAfterName" +(std::string) MVAMethod;
  std::string NeutrinoTrueEAfterName = "NeutrinoTrueEAfterName" +(std::string) MVAMethod;
  std::string NeutrinoRecoEAfterExtraName = "NeutrinoRecoEAfterExtraName" +(std::string) MVAMethod;
  std::string NeutrinoTrueEAfterExtraName = "NeutrinoTrueEAfterExtraName" +(std::string) MVAMethod;

  if(!Metric.CheckHistogram(NeutrinoTrueEBeforeName,sample)){Metric.SetHistogram(NeutrinoTrueEBeforeName,30,0,3000,sample);}
  if(!Metric.CheckHistogram(NeutrinoRecoEBeforeName,sample)){Metric.SetHistogram(NeutrinoRecoEBeforeName,30,0,3000,sample);}
  if(!Metric.CheckHistogram(NeutrinoRecoEBeforeExtraName,sample)){Metric.SetHistogram(NeutrinoRecoEBeforeExtraName,30,0,3000,sample);}
  if(!Metric.CheckHistogram(NeutrinoRecoEAfterName,sample)){Metric.SetHistogram(NeutrinoRecoEAfterName,30,0,3000,sample);}
  if(!Metric.CheckHistogram(NeutrinoTrueEAfterName,sample)){Metric.SetHistogram(NeutrinoTrueEAfterName,30,0,3000,sample);}
  if(!Metric.CheckHistogram(NeutrinoRecoEAfterExtraName,sample)){Metric.SetHistogram(NeutrinoRecoEAfterExtraName,30,0,3000,sample);}
  if(!Metric.CheckHistogram(NeutrinoTrueEAfterExtraName,sample)){Metric.SetHistogram(NeutrinoTrueEAfterExtraName,30,0,3000,sample);}
  
  
  //Loop over the neutrinos and get their true energy
  for(int n_iter=0; n_iter<nu_E_vals->size(); ++n_iter){
    Metric.GetHistogram(NeutrinoTrueEBeforeName,sample)->Fill(nu_truth_energy_vals->at(n_iter)*1000);
  }

  for(int n_iter=0; n_iter<nu_reco_energy_vals->size(); ++n_iter){
    //check if the true neutrino is already be identified 
    if(MaxNeutrinoE[nu_truth_energy_vals->at(n_iter)] < nu_reco_energy_vals->at(n_iter)){
      MaxNeutrinoE[nu_truth_energy_vals->at(n_iter)] = nu_reco_energy_vals->at(n_iter);
      MaxNeutrinoE_iter[nu_truth_energy_vals->at(n_iter)] = n_iter;
    }
  }
  //Add the max reco values to the efficiency and the others to the additional graph
  for(int n_iter=0; n_iter<nu_reco_energy_vals->size(); ++n_iter){
    if(MaxNeutrinoE[nu_truth_energy_vals->at(n_iter)] == nu_reco_energy_vals->at(n_iter)){
      Metric.GetHistogram(NeutrinoRecoEBeforeName,sample)->Fill(nu_reco_energy_vals->at(n_iter));
    }
    else{
      Metric.GetHistogram(NeutrinoRecoEBeforeExtraName,sample)->Fill(nu_reco_energy_vals->at(n_iter)); 
    }
  }

  //Perform the BDT Cut
  
  //Loop over the neutrinos
  for(int n_iter=0; n_iter<nu_reco_energy_vals->size(); ++n_iter){

    //Perform the BDT Cut
    if(reader->EvaluateMVA(MVAMethod) < 0){continue;}

    //Save the Neutrino E
    NeutrinoE[nu_truth_energy_vals->at(n_iter)].push_back(nu_reco_energy_vals->at(n_iter));
  }

  //Identify if this neutrino is the largest
  for(auto const& trueNeutrino: NeutrinoE){
    float MaxNeutrinoEval      = -99999;
    int   MaxNeutrinoEval_iter = -99999; 
    //Identify the remianing reco neutrino with the max energy
    for(int n_iter=0; n_iter<(trueNeutrino.second.size()); ++n_iter){
      if(MaxNeutrinoEval < trueNeutrino.second.at(n_iter)){
	MaxNeutrinoEval = trueNeutrino.second.at(n_iter);
	MaxNeutrinoEval_iter = n_iter;
      }
    }
    for(int n_iter=0; n_iter<(trueNeutrino.second.size()); ++n_iter){ 
      //Add the biggest to account for the efficiency.
      if(n_iter == MaxNeutrinoEval_iter){
	Metric.GetHistogram(NeutrinoRecoEAfterName,sample)->Fill(trueNeutrino.second.at(n_iter));
	Metric.GetHistogram(NeutrinoTrueEAfterName,sample)->Fill(trueNeutrino.first*1000);
      }
      //Add the extra reco
      else{
	Metric.GetHistogram(NeutrinoRecoEAfterExtraName,sample)->Fill(trueNeutrino.second.at(n_iter));
	Metric.GetHistogram(NeutrinoTrueEAfterExtraName,sample)->Fill(trueNeutrino.first*1000);
      }
    }
  }
return 0;
}

  

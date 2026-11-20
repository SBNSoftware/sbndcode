#include "FRANSAnaCommandLineParser.h"

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

int MacroFRANSTraining(std::string fInputFileName="", double nTrainFrac = -1, int view = 2, std::string fBGLabel = "", std::string fTreeDirName = "FRANSCheatedVx/", std::string fTreeName = "FRANSTree")
{

  //--------- Configuration Parameters
  bool fUseMultiplicityVars = false; //fUseMultiplicityVars = true;
  bool fApplyFiducialCut = true; //fApplyFiducialCut = false;
  std::string fVw = std::to_string(view);

  //--------- Input file
  TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
  fFile->ls();
  TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());

  //--------- Output file
  gSystem->Exec( "rm -rf TMVAResults" );
  gSystem->Exec( "mkdir TMVAResults" );
  std::string fOutputTMVAROOtFileName = "TMVAResults/TMVAResults.root";
  TFile* outputFile = TFile::Open( fOutputTMVAROOtFileName.c_str(), "CREATE" );

  //---------  MVA methods to be trained
  std::map<std::string,int> Use;
  // Rectangular cut optimisation
  Use["Cuts"]            = 0;
  // 1-dimensional likelihood ("naive Bayes estimator")
  Use["Likelihood"]      = 0;
  // Linear Discriminant Analysis
  Use["Fisher"]          = 0;
  // Boosted Decision Trees
  Use["BDT"]             = 1;
  // Neural Networks (all are feed-forward Multilayer Perceptrons)
  Use["MLP"]             = 0;

  //--------- Variables to use
  std::map<std::string,int> UseVars;
  UseVars["Eta"] = 1;
  UseVars["Delta"] = 1;
  UseVars["Alpha"] = 0;
  UseVars["Iota"] = 1;
  UseVars["FitScore"] = 1;
  // get number of vafriables to use
  int nVars = 0;
  for(auto it = UseVars.begin(); it != UseVars.end(); ++it){
    if(it->second) nVars++;
  }

  //---------  Load TMVA
  TMVA::Tools::Instance();
  std::string PreTransformations = nVars>1 ? "I;D;P;G,D" : "I;D;G,D";
  TMVA::Factory *factory = new TMVA::Factory( "FRANSSelectionTMVA", outputFile,
                                              "!V:!Silent:Color:DrawProgressBar:Transformations="+PreTransformations+":AnalysisType=Classification" );
  TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");


  //---------  Add variables to the dataloader
  if(UseVars["Eta"]) dataloader->AddVariable( "FRANSObj"+fVw+".fEta", "#eta", "", 'D' );
  if(UseVars["Delta"]) dataloader->AddVariable( "FRANSObj"+fVw+".fDelta", "#Delta", "", 'D' );
  if(UseVars["Alpha"]) dataloader->AddVariable( "FRANSObj"+fVw+".fAlpha", "#alpha", "", 'D' );
  if(UseVars["Iota"]) dataloader->AddVariable( "FRANSObj"+fVw+".fIota", "#iota", "", 'D' );
  if(UseVars["FitScore"]) dataloader->AddVariable( "FRANSObj"+fVw+".fFitScore", "r", "", 'D' );
  if(fUseMultiplicityVars){
    dataloader->AddVariable( "FRANSObj"+fVw+".fNOrigins", "N", "", 'I' );
    dataloader->AddVariable( "FRANSObj"+fVw+".fNOriginsM1", "N^{1}", "", 'I' );
    dataloader->AddVariable( "FRANSObj"+fVw+".fNOriginsM2", "N^{2}", "", 'I' );
    //dataloader->AddVariable( "FRANSObj2.NOriginsM3_C", "N^{>3}_C", "", 'I' );
    dataloader->AddVariable( "FRANSObj"+fVw+".fHitDensity", "d", "", 'D' );
  }

  //--------- Add spectators
  dataloader->AddSpectator( "Gap", "Gap", "", 'D' );
  dataloader->AddSpectator( "ProtonKE", "ProtonKE", "", 'D' );
  dataloader->AddSpectator( "PionKE", "PionKE", "", 'D' );

  //--------- Add the signal and background trees
  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;
  dataloader->AddSignalTree    ( fTree,     signalWeight );
  dataloader->AddBackgroundTree( fTree, backgroundWeight );

  //-------- Signal and BG definitions
  TCut signalCut  = "IsSignal==1";
  TCut bgCut   = "IsSignal!=1";
  TCut bgCutQE = "IsSignal!=1 && IntMode==0";
  TCut bgCutRES = "IsSignal!=1  && IntMode==1";
  TCut bgCutDIS = "IsSignal!=1  && IntMode==2";
  if(fBGLabel=="QE") bgCut = bgCutQE;
  else if(fBGLabel=="RES") bgCut = bgCutRES;
  else if(fBGLabel=="DIS") bgCut = bgCutDIS;

  //--------- Fiducial cut
  TCut fiducialCut = "abs(TrueVx)<190 && abs(TrueVy)<190 && TrueVz>10 && TrueVz<490";
  if(fApplyFiducialCut){
    signalCut = signalCut && fiducialCut;
    bgCut = bgCut && fiducialCut;
  }

  //--------- Set the training and testing fractions
  int nMaxSignal = fTree->Draw(">>selectedEntries", signalCut, "entrylist")-1;
  int nMaxBg = fTree->Draw(">>selectedEntries", bgCut, "entrylist")-1;
  std::cout<<"nMaxSignal: "<<nMaxSignal<<std::endl;
  std::cout<<"nMaxBg: "<<nMaxBg<<std::endl;

  //--------- Prepare training and test trees
  std::string tmva_options = "";
  if(nTrainFrac<0) tmva_options = "nTrain_Signal="+std::to_string(nMaxSignal)+":nTrain_Background="+to_string(nMaxBg);
  else tmva_options = "nTrain_Signal="+std::to_string((int)nTrainFrac*nMaxSignal)+":nTrain_Background="+to_string((int)nTrainFrac*nMaxBg);
  tmva_options+=":SplitMode=Random:NormMode=NumEvents:!V";
  dataloader->PrepareTrainingAndTestTree( signalCut, bgCut, tmva_options );

  //--------- Book MVA methods
  // Cut optimisation
  if (Use["Cuts"])
    factory->BookMethod( dataloader, TMVA::Types::kCuts, "Cuts",
                         "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );
  // Likelihood ("naive Bayes estimator")
  if (Use["Likelihood"])
    factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "Likelihood",
                         "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );
  // Fisher discriminant (same as LD)
  if (Use["Fisher"])
    factory->BookMethod( dataloader, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );
  // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
  if (Use["MLP"])
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
  // Boosted Decision Tree
  if (Use["BDT"])  // Adaptive Boost
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
                         "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

  //-------- Train MVAs using the set of training events
  factory->TrainAllMethods();

  //--------  Evaluate all MVAs using the set of test events
  factory->TestAllMethods();

  //--------  Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();


  //--------  Save the output
  outputFile->Close();

  delete factory;
  delete dataloader;

  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVA::TMVAGui( fOutputTMVAROOtFileName.c_str() );

  return 0;
}

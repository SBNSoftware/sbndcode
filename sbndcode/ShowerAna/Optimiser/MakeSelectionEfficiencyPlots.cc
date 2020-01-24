//C++ Inlcudes
#include <string> 
#include <getopt.h>
#include <iostream> 

//Root Includes 
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

//Framework inludes
#include "BranchType.hh"
#include "MetricHolder.hh"
#include "NueRecoOptimiser.hh"

int main(int argc, char* argv[]){
  
  //Read The files.
  std::string signalfile;
  std::string backgroundfile;
  bool RunMVATrainng = false;

  static struct option long_options[] = {
    {"signalfile"     , required_argument,  0,  's' },
    {"backgroundfiles", required_argument,  0,  'b' },
    {"MVA", no_argument,  0,  'M' },
    {"help", no_argument,  0,  'h' },
    {0, 0}
  };

  if(argc == 1){
    std::cout << "Please give a signal file with a -s myfile_sig.root and a background file -b myfile_bk.root" << std::endl; 
    return 1;
  }
  
  int opt = 0;
  int c = 0;
  while ((opt = getopt_long_only(argc, argv,"", long_options, &c )) != -1) {
    switch (opt) {
    case 's' : signalfile = optarg;
      break;
    case 'b' : backgroundfile = optarg;
      break;
    case 'M': RunMVATrainng = true;
      break;
    case 'h': 
      std::cout << "Please give a signal file with a -s myfile_sig.root and a background file -b myfile_bk.root" << std::endl;
    default:
      std::cout << "Please give a signal file with a -s myfile_sig.root and a background file -b myfile_bk.root" << std::endl;
      return 1;
    }
  }

  //Root Dictionaries cos I'm lazy.
  gROOT->ProcessLine("#include <vector>");
  gInterpreter->GenerateDictionary("vector<float>", "vector");
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gStyle->SetPalette(1);

  //Get the histogram for the signal
  TFile *f_eff = new TFile(signalfile.c_str());
  TDirectory* dir_eff = f_eff->GetDirectory("efftwo");
  TTree *tree_sig = (TTree*)dir_eff->Get("RecoEffMetricTree");  
  TTree *POT_sigtree = (TTree*)dir_eff->Get("RecoEffPOTTree");

  //Get the histogram for the background
  TFile *f_bk = new TFile(backgroundfile.c_str());
  TDirectory* dir_bk = f_bk->GetDirectory("efftwo");
  TTree *tree_bk = (TTree*)dir_bk->Get("RecoEffMetricTree");  
  TTree *POT_bktree = (TTree*)dir_bk->Get("RecoEffPOTTree");

  //Get The POT
  float TotalSigPOT = 0;
  float TotalBKPOT  = 0;
  float SignalPOT;
  float BackgrndPOT;
  
  POT_sigtree->SetBranchAddress("POT",&SignalPOT);
  POT_bktree->SetBranchAddress("POT",&BackgrndPOT);

  for (Long64_t evt = 0; evt < POT_sigtree->GetEntries(); evt++) {
    POT_sigtree->GetEntry(evt);
    TotalSigPOT += SignalPOT;
  }

  for (Long64_t evt = 0; evt < POT_bktree->GetEntries(); evt++) {
    POT_bktree->GetEntry(evt);
    TotalBKPOT += BackgrndPOT;
  }

  //Make New Optimiser Object. Setup up the branches
  optimiser::NueRecoOptimiser NueOptimiser(tree_sig,tree_bk,TotalSigPOT,TotalBKPOT);

  //Set Up the metrics.
  NueOptimiser.InitialiseMetrics();

  NueOptimiser.SetupMVA();

  //Read the data
  NueOptimiser.FillData(tree_sig,tree_bk);

  //Run the MVA Training
  if(RunMVATrainng){
    NueOptimiser.TrainMVA();
  }

  TFile* file = new TFile("CutFile.root", "RECREATE");

  //Make the pretty graphs 
  NueOptimiser.AnalyseData();

  return 0;
}

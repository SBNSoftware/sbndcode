//Framework Includes
#include "BranchType.hh"
#include "MetricHolder.hh"

//Root Inlcues
#include "TTree.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/Reader.h"

//C++ Includes 
#include <map>
#include <string> 
#include <vector>

namespace optimiser { 
  class NueRecoOptimiser;
}

class optimiser::NueRecoOptimiser {
  
public:

  NueRecoOptimiser(TTree* signaltree, TTree* backgroundtree, float bkpot, float sigpot);

  int TrainMVA();
  int SetupMVA();
 
  //Branches from Reco Efficiency Finder.
  struct EfficiencyTree {
    std::map<std::string, optimiser::BranchTypeBase*> Tree;
    
    template <class T> 
    void AddBranch(std::string Name,TTree* tree){
      Tree[Name] = new optimiser::BranchType<T>();
      BranchType<T>* Branch = dynamic_cast<BranchType<T>*>(Tree[Name]);
      tree->SetBranchAddress(Name.c_str(),&(Branch->GetPointer()));
      return;
    }
  };

  void InitialiseMetrics();
  
  void InitialiseMetric(std::string Name, 
			std::string branchname,
			int numbins, 
			float xmin, 
			float xmax,
			std::string AnalysisName="StandardShowerCutFinder",
			bool ApplyPOT = false,
		        bool ApplyOscProb = false,
			bool fFillMVA=false,
			float MVAErrorVal= -999,
			TString MVAName=""
			);

  void InitialiseMetric(std::string Name, 
			std::string AnalysisName="StandardShowerCutFinder",
			bool ApplyPOT = false,
		        bool ApplyOscProb = false,
			bool fFillMVA=false,
			float MVAErrorVal= -999,
			TString MVAName="");


  void InitialiseMetric(std::string Name, 
			std::string branchname,
			int numbins, 
			float xmin, 
			float xmax, 
			std::string metricpartner,
			int numbins_partner, 
			float xmin_partner, 
			float xmax_partner,
			std::string AnalysisName="StandardShowerCutFinder",
			bool ApplyPOT = false,
		        bool ApplyOscProb = false,
			bool fFillMVA=false,
			float MVAErrorVal= -999,
			TString MVAName="");

  void SetLogic(std::string Name, bool ls, bool pls, bool andor);
  
  void FillData(TTree* tree_signal, TTree* tree_background);
  
  int PerformCuts(optimiser::MetricHolder& Metric, 
		  std::string sample,
		  std::vector<float>& vals, 
		  std::vector<std::pair<float,float> >& TwoDvals,
		  TString& MVAName);
 
  template <class T>
  int GetBranchPointer(std::string& branchname, std::string& sample, T* &vals); 
  
  void AnalyseData();
  
  int StandardDaughterAnalysis(optimiser::MetricHolder& Metric, 
			       std::string& sample, 
			       std::vector<float>& vals, 
			       std::vector<std::pair<float,float> >& TwoDvals
			       );
  
  int StandardDaughter1DAnalysis(optimiser::MetricHolder& Metric,
				 std::string& sample, 
				 std::vector<float>& vals
				 );
  
  int StandardDaughter2DAnalysis(optimiser::MetricHolder& Metric, 
				 std::string& sample, 
				 std::vector<std::pair<float,float> >& TwoDvals
				 );
  
  int StandardOverEDaughterAnalysis(optimiser::MetricHolder& Metric, 
				    std::string& sample, 
				    std::vector<float>& vals, 
				    std::vector<std::pair<float,float> >& TwoDvals
				    );
  
  int StandardOverEDaughter1DAnalysis(optimiser::MetricHolder& Metric,
				      std::string& sample, 
				      std::vector<float>& vals
				      );
  
  int StandardOverEDaughter2DAnalysis(optimiser::MetricHolder& Metric, 
				      std::string& sample, 
				      std::vector<std::pair<float,float> >& TwoDvals
				      );


  int StandardNeutrinoAnalysis(optimiser::MetricHolder& Metric, 
			       std::string& sample, 
			       std::vector<float>& vals, 
			       std::vector<std::pair<float,float> >& TwoDvals
			       );
  
  int StandardNeutrino1DAnalysis(optimiser::MetricHolder& Metric, 
				 std::string& sample, 
				 std::vector<float>& vals
				 );

  int StandardNeutrino2DAnalysis(optimiser::MetricHolder& Metric, 
				 std::string& sample,
				 std::vector<std::pair<float,float> >& TwoDvals
				 );
  
  int StandardSizeNeutrinoAnalysis(optimiser::MetricHolder& Metric, 
				   std::string& sample, 
				   std::vector<float>& vals, 
				   std::vector<std::pair<float,float> >& TwoDvals
				   );
  
  int StandardSizeNeutrino1DAnalysis(optimiser::MetricHolder& Metric, 
				     std::string& sample, 
				     std::vector<float>& vals
				     );

  int StandardSizeNeutrino2DAnalysis(optimiser::MetricHolder& Metric, 
				     std::string& sample,
				     std::vector<std::pair<float,float> >& TwoDvals
				     );


  int StandardNeutrinoDaughterAnalysis(optimiser::MetricHolder& Metric, 
				       std::string& sample, 
				       std::vector<float>& vals, 
				       std::vector<std::pair<float,float> >& TwoDvals
				       );
  
  int StandardNeutrinoDaughter2DAnalysis(optimiser::MetricHolder& Metric, 
					 std::string& sample, 
					 std::vector<std::pair<float,float> >& TwoDvals
					 );

  int StandardNumShowerAnalysis(optimiser::MetricHolder& Metric, 
				std::string& sample, 
				std::vector<std::pair<float,float> >& TwoDvals
				);

  int OneShower2DEnergyAnalysis(optimiser::MetricHolder& Metric, 
				std::string& sample, 
				std::vector<std::pair<float,float> >& TwoDvals
				); 
 
  int OneShower1DCutAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals);

 
  int ShowerResidualAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals);

  int ShowerResidualCut(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals);
  
  int ShowerResidualNumShowers(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals);

  int ShowerAllsCut(optimiser::MetricHolder& Metric, std::string& sample);

  int MVAAnalysis(optimiser::MetricHolder& Metric, 
		  std::string& sample, 
		  std::vector<float>& vals,
		  TString& MVAMethod);

  int MVAAnalysisAllsCut(optimiser::MetricHolder& Metric, 
			 std::string& sample, 
			 TString& MVAMethod);

  float GetOscProb(int& iter, std::string sample, bool applyoscprob);


  EfficiencyTree SignalTree;
  EfficiencyTree BackgroundTree;

  
  TFile* MVAFile;
  TTree* BackgroundTreeMVA;
  TTree* SignalTreeMVA;
  TMVA::Reader *reader;

  double BackgrndWeight = 1;
  double SignalWeight  = 1;

  double TotalSigPOT = 1;
  double TotalBKPOT  = 1;

  bool ApplyOscWht;
  bool ApplyPOTWht;

  std::map<std::string,optimiser::MetricHolder> MetricMap;

  //Analysis cuts 
  float fEnergyCut = 200;

};

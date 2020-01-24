//Class to hold the histograms of the metrics analysed by NueOptimiser.
//This performs the optimisation process.

//C++ Includes 
#include <string>
#include <map> 
#include <iostream>

//Root Includes 
#include "TString.h"
#include "TH1.h"
#include "TH2.h"

#ifndef MetricHolder_HH
#define MetricHolder_HH

namespace optimiser { 
  class MetricHolder;
}


class optimiser::MetricHolder{

public:

  MetricHolder(){
  };

  MetricHolder(std::string name, std::string branchname, int numbins, float xmin, float xmax, float sigpot, float bkpot, bool ApplyPOT, bool ApplyOscProb, std::string analysisname, bool FillMVA,  TString mvaname, float mvaerror){
    MetricName = name;
    BranchName = branchname;
    AnalysisName  = analysisname;
    MetricPartner = "";
    std::string SignalName     = name + "_signal"; 
    std::string BackgroundName = name + "_background";
    SignalHist = new TH1D(SignalName.c_str(),SignalName.c_str(),numbins,xmin,xmax);
    BackgroundHist = new TH1D(BackgroundName.c_str(),BackgroundName.c_str(),numbins,xmin,xmax);

    fLessThan = true;
    fPartnerLessThan = true;
    fAND = true;
    fFillMVA = FillMVA;
    MVAName = mvaname;
    MVAErrorVal = mvaerror;
    MVAMinCap = xmin;
    MVAMaxCap = xmax;
    
    TotalSigPOT = sigpot;
    TotalBKPOT  = bkpot;

    applyoscprob = ApplyOscProb;
    applypot     = ApplyPOT; 
  };

  MetricHolder(std::string name, float sigpot, float bkpot, bool ApplyPOT, bool ApplyOscProb, std::string analysisname, bool FillMVA, TString mvaname, float mvaerror){
    MetricName = name;
    BranchName = "";
    AnalysisName  = analysisname;
    MetricPartner = "";
    fLessThan = true;
    fPartnerLessThan = true;
    fAND = true;
    fFillMVA = FillMVA;
    MVAName =  mvaname;
    MVAErrorVal = mvaerror;

    TotalSigPOT = sigpot;
    TotalBKPOT  = bkpot;

    applyoscprob =  ApplyOscProb;
    applypot   = ApplyPOT;
  };

  


  MetricHolder(std::string name, std::string branchname, int numbins, float xmin, float xmax, std::string metricpartner,int numbins_partner, float xmin_partner, float xmax_partner,float sigpot, float bkpot, bool ApplyPOT, bool ApplyOscProb, std::string analysisname, bool FillMVA, TString mvaname, float mvaerror){
    MetricName = name;
    BranchName = branchname;
    MetricPartner = metricpartner;
    AnalysisName  = analysisname;
    std::string SignalName     = name + "_signal"; 
    std::string BackgroundName = name + "_background";
    SignalHist = new TH1D(SignalName.c_str(),SignalName.c_str(),numbins,xmin,xmax);
    BackgroundHist = new TH1D(BackgroundName.c_str(),BackgroundName.c_str(),numbins,xmin,xmax);

    std::string SignalName2D     = name + "_signal_2D"; 
    std::string BackgroundName2D = name + "_background_2D";
    SignalHistPartner = new TH2D(SignalName2D.c_str(),SignalName2D.c_str(),numbins,xmin,xmax,numbins_partner,xmin_partner,xmax_partner);
    BackgroundHistPartner = new TH2D(BackgroundName2D.c_str(),BackgroundName2D.c_str(),numbins,xmin,xmax,numbins_partner,xmin_partner,xmax_partner);

    fLessThan = true;
    fPartnerLessThan = true;
    fAND = true;
    fFillMVA = FillMVA;
    MVAName = mvaname;
    MVAErrorVal = mvaerror;
    MVAMinCap = xmin;
    MVAMaxCap = xmax;

    TotalSigPOT = sigpot;
    TotalBKPOT  = bkpot;


    applyoscprob =  ApplyOscProb;
    applypot   = ApplyPOT;
  };

  void SetLogic(bool lt, bool plt, bool andor){
    fLessThan = lt;
    fPartnerLessThan = plt;
    fAND = andor;
  }

  void FillSignal(float& val,double& oscprob){
    SignalHist->Fill(val,oscprob);
  }
  void FillBackground(float& val, double& oscprob){
    BackgroundHist->Fill(val,oscprob);
  }

  void FillSignal(float& x,float& y, double& oscprob){
    SignalHistPartner->Fill(x,y,oscprob);
  }
  void FillBackground(float& x, float& y, double& oscprob){
    BackgroundHistPartner->Fill(x,y,oscprob);
  }

  int TotalEntries(){
    return SignalHist->GetEntries() + BackgroundHist->GetEntries();
  }

  int Total2DEntries(){
    return SignalHistPartner->GetEntries() + BackgroundHistPartner->GetEntries();
  }

  TH1D* GetSignalHist(){
    return SignalHist;
  }
  TH1D* GetBackgroundHist(){
    return BackgroundHist;
  }
  TH2D* GetSignalHistPartner(){
    return SignalHistPartner;
  }
  TH2D* GetBackgroundHistPartner(){
    return BackgroundHistPartner;
  }

  void Perform1DCutFinder();
  void Perform1DCutFinder(TH1D* signalhist, TH1D* backgroundhist);
  void Perform2DCutFinder();  

  void Plot1D(TH1D* signalhist, TH1D* backgroundhist);
  void Plot1D();
  void Plot2D(TH2D* signalhist, TH2D* backgroundhist);
  void Plot2D();

  void MakeEfficiencyPlots(TH1D* postcut, TH1D* precut, TH1D* extra);

  void MakeStandardNumShowerAnalysisGraphs(float totalsig, float totalbk);

  std::string GetBranchName(){return BranchName;}
  std::string GetMetricName(){return MetricName;}
  std::string GetPartnerMetricName(){return MetricPartner;}
  std::string GetAnalysisName(){return AnalysisName;}

  TH1D* GetHistogram(std::string Name, std::string signal){
    if(signal == "signal"){
      Name += "_signal";
      return HistSignalmap[Name];
    }
    else{
      Name += "_background";
      return HistBackgroundmap[Name];
    }
  }

  void SetHistogram(std::string Name,int bins,float xmin, float xmax, std::string signal){
    if(signal == "signal"){
      Name += "_signal";
      HistSignalmap[Name] = new TH1D(Name.c_str(),Name.c_str(),bins,xmin,xmax);
    }
    else{
      Name += "_background";
      HistBackgroundmap[Name] = new TH1D(Name.c_str(),Name.c_str(),bins,xmin,xmax);
    }
  };

  bool CheckHistogram(std::string Name, std::string signal){
    if(signal == "signal"){
      Name += "_signal";
      if(HistSignalmap.find(Name) != HistSignalmap.end()){return true;}
      else return false;
    }
    else{
      Name += "_background";
      if(HistBackgroundmap.find(Name) != HistBackgroundmap.end()){return true;}
      else return false;
    }
  };

  int GetMVAVectorSignalSize(){
    return MVAVectorSignal.size();
  }

  int GetMVAVectorBackgroundSize(){
    return MVAVectorBackground.size();
  }

  void FillMVAValSignal(float& val){
    MVAValSignal = val;
    MVAVal = val;
  }

  void FillMVAValBackground(float& val){
    MVAValBackground = val;
    MVAVal = val;
  }

  void FillMVAValSignalVec(std::vector<float>& val){
    MVAVectorSignal = val;
  }
  
  void MVAValBackgroundVec(std::vector<float>& val){
    MVAVectorBackground = val;
  }

  bool FillForMVA(){
    return fFillMVA;
  }

  float GetMVAValSignal(int ind){
    if(MVAVectorSignal.size() <= ind){
      return MVAErrorVal;
    }
    if(MVAVectorSignal.at(ind) == -999){
      return MVAErrorVal;
    }
    if(MVAVectorSignal.at(ind) < MVAMinCap){
      return MVAMinCap;
    }
    if(MVAVectorSignal.at(ind) > MVAMaxCap){
      return MVAMaxCap;
    }
    return MVAVectorSignal.at(ind);
  }

  float GetMVAValBackground(int ind){
    if(MVAVectorBackground.size() <= ind){
      return MVAErrorVal;
    }
    if(MVAVectorBackground.at(ind) == -999){
      return MVAErrorVal;
    }
    if(MVAVectorBackground.at(ind) < MVAMinCap){
      return MVAMinCap;
    }
    if(MVAVectorBackground.at(ind) > MVAMaxCap){
      return MVAMaxCap;
    }
    return MVAVectorBackground.at(ind);
  }

  float& GetMVAValSignal(){
    return MVAValSignal;
  }

  float& GetMVAValBackground(){
    return MVAValBackground;
  }

  float& GetMVAVal(){
    return MVAVal;
  }

  TString GetMVAName(){
    return MVAName;
  }

  float GetErrorVal(){
    return MVAErrorVal;
  }

  bool ApplyOscProb(){
    return applyoscprob;
  }

  bool ApplyPOT(){
    return applypot;
  }

  double ReturnPOT(std::string signal){
    if(signal == "signal"){return 6.6e20/TotalSigPOT;}
    else if(signal =="background"){return 6.6e20/TotalBKPOT;}

    std::cerr << "POT not return. BUG!" << std::endl;
    return -999;
  }
    
private:

  bool applyoscprob;
  bool applypot;

  //MVA Training 
  std::vector<float> MVAVectorSignal;
  float MVAValSignal;
  std::vector<float> MVAVectorBackground;
  float MVAValBackground;
  bool fFillMVA;
  float MVAErrorVal;
  float MVAMinCap;
  float MVAMaxCap;

  //MVA Usage 
  float MVAVal;
  TString MVAName;

  std::string MetricName;
  std::string BranchName;
  std::string MetricPartner;
  std::string AnalysisName;
  TString units;
  TString axisTitle;
  TH1D* SignalHist;
  TH1D* BackgroundHist;

  TH2D* SignalHistPartner;
  TH2D* BackgroundHistPartner;

  std::map<std::string,TH1D*> HistSignalmap;
  std::map<std::string,TH1D*> HistBackgroundmap;
  bool fLessThan;
  bool fPartnerLessThan;
  bool fAND;

  float TotalSigPOT;
  float TotalBKPOT;
  
};

#endif

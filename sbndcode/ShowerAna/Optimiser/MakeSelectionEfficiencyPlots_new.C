#include "TH1.h"
#include "TTree.h"
#include <vector>
#include "TGraph.h"
#include "TH2.h"
#include "TGraph2D.h"
#include "TMath.h"

#include <string>
#include <vector>
#include <iostream> 

namespace optimiser { 
  class NueRecoOptimiser;
  class MetricHolder;
  class BranchType<class T>;
}


template <class T>
class optimiser::BranchType{

public:

  BranchType(){
    Branch = 0;
  }

  T*& GetPointer(){return Branch;}

private:

  T* Branch;

}

class optimiser::MetricHolder{

public:

  MetricHolder(std::string name, std::string branchname, int numbins, float xmin, float xmax){
    MetricName = name;
    BranchName = branchname;
    MetricPartner = "";
    std::string SignalName     = name + "_signal"; 
    std::string BackgroundName = name + "_background";
    SignalHist = new TH1D(SignalName.c_str(),SignalName.c_str(),numbins,xmin,xmax);
    BackgroundHist = new TH1D(BackgroundName.c_str(),BackgroundName.c_str(),numbins,xmin,xmax);
  };

  MetricHolder(std::string name, std::string branchname, int numbins, float xmin, float xmax, std::string metricpartner,int numbins_partner, float xmin_partner, float xmax_partner){
    MetricName = name;
    BranchName = branchname;
    MetricPartner = metricpartner;
    std::string SignalName     = name + "_signal"; 
    std::string BackgroundName = name + "_background";
    SignalHist = new TH1D(SignalName.c_str(),SignalName.c_str(),numbins,xmin,xmax);
    BackgroundHist = new TH1D(BackgroundName.c_str(),BackgroundName.c_str(),numbins,xmin,xmax);

    std::string SignalName2D     = name + "_signal_2D"; 
    std::string BackgroundName2D = name + "_background_2D";
    SignalHistPartner = new TH2D(SignalName2D.c_str(),SignalName2D.c_str(),numbins,xmin,xmax,numbins_partner,xmin_partner,xmax_partner);
    BackgroundHistPartner = new TH2D(BackgroundName2D.c_str(),BackgroundName2D.c_str(),numbins,xmin,xmax,numbins_partner,xmin_partner,xmax_partner);

  };


  void FillSignal(float& val){
    SignalHist->Fill(val);
  }
  void FillBackground(float& val){
    BackgroundHist->Fill(val);
  }

  void FillSignal(float& x,float& y){
    SignalHistPartner->Fill(x,y);
  }
  void FillBackground(float& x, float& y){
    BackgroundHistPartner->Fill(x,y);
  }

  void Perform1DCutFinder();
  void Perform2DCutFinder();  

  std::string GetBranchName(){return BranchName;}
  std::string GetMetricName(){return MetricName;}
  std::string GetPartnerMetricName(){return MetricPartner;}

private:

  std::string MetricName;
  std::string BranchName;
  std::string MetricPartner;
  TString units;
  TString axisTitle;
  TH1D* SignalHist;
  TH1D* BackgroundHist;

  TH2D* SignalHistPartner;
  TH2D* BackgroundHistPartner;

};


void optimiser::MetricHolder::Perform2DCutFinder(){
  return;
}

void optimiser::MetricHolder::Perform1DCutFinder(){

  // TH1::AddDirectory(kFALSE);
  SignalHist->Scale(1 / SignalHist->GetEntries());
  BackgroundHist->Scale(1 /BackgroundHist->GetEntries());
  
  //Find the min and max start positions
  float minsig = SignalHist->GetBinCenter(0);
  float maxsig = SignalHist->GetBinCenter(SignalHist->GetNbinsX());
  float minbk = BackgroundHist->GetBinCenter(0);
  float maxbk = BackgroundHist->GetBinCenter(BackgroundHist->GetNbinsX());

  if (minsig != minbk) {
    std::cout << "min point does not match" << std::endl;
    return;
  }
  if (maxbk != maxbk) {
    std::cout << "max point does not match" << std::endl;
    return;
  }
  if (SignalHist->GetNbinsX() != BackgroundHist->GetNbinsX()) {
    std::cout << "bin numbers sig " << SignalHist->GetNbinsX() << " and background " << BackgroundHist->GetNbinsX() << "do not match" << std::endl;
    return;
  }

  int max_ind = -999;
  int maxsigbk_ind = -999;
  int n = SignalHist->GetNbinsX();
  double bktot_err = 0;
  double bktot = BackgroundHist->IntegralAndError(1, n+1, bktot_err);
  double sigtot_err = 0;
  double sigtot = SignalHist->IntegralAndError(1, n+1, sigtot_err);
  double maxeffpur = -999;
  double maxsigbk = -999;
  double efficiency[n];
  double efficencyErr[n];
  double purity[n];
  double purityErr[n];
  double bkRej[n];
  double bkRejErr[n];
  double effpur[n];
  double effpurErr[n];
  double sigbk[n];
  double sigbkErr[n];
  double xval[n];
  double xvalErr[n];

  //Get the values.
  for (int i = 1; i < n+1; ++i) {

    double presig_err = 0;
    double prebk_err = 0;
    double presig = SignalHist->IntegralAndError(1, i, presig_err);
    double prebk = BackgroundHist->IntegralAndError(1, i, prebk_err);

    //    std::cout << "presig: " << presig << " prebk: " << prebk << std::endl;

    efficiency[i] = presig / sigtot;
    efficencyErr[i] = efficiency[i] * TMath::Sqrt(((presig_err / presig) * (presig_err / presig)) + (sigtot_err / sigtot) * (sigtot_err / sigtot));

    bkRej[i] = 1 - prebk / bktot;
    bkRejErr[i] = bkRej[i] * TMath::Sqrt((prebk_err / prebk) * (prebk_err / prebk) + (bktot_err / bktot) * (bktot_err / bktot));

    if (presig + prebk != 0) {
      purity[i] = presig / (presig + prebk);
      double purityErr_dem = TMath::Sqrt((presig_err * presig_err) + (prebk_err * prebk_err));
      purityErr[i] = purity[i] * TMath::Sqrt(((presig_err / presig) * (presig_err / presig)) + (purityErr_dem / (presig + prebk)) * (purityErr_dem / (presig + prebk)));
    } else {
      purity[i] = 0;
      purityErr[i] = 0;
    }

    effpur[i] = efficiency[i] * purity[i];
    if (efficiency[i] == 0 || purity[i] == 0) {
      effpurErr[i] = 0;
    } else {
      effpurErr[i] = effpur[i] * TMath::Sqrt(((efficencyErr[i] / efficiency[i]) * (efficencyErr[i] / efficiency[i])) + (purityErr[i] / purity[i]) * (purityErr[i] / purity[i]));
    }

    sigbk[i] = efficiency[i] * bkRej[i];
    if (efficiency[i] == 0 || bkRej[i] == 0) {
      sigbkErr[i] = 0;
    } else {
      sigbkErr[i] = sigbk[i] * TMath::Sqrt(((efficencyErr[i] / efficiency[i]) * (efficencyErr[i] / efficiency[i])) + (bkRejErr[i] / bkRej[i]) * (bkRejErr[i] / bkRej[i]));
    }

    xval[i] = SignalHist->GetBinCenter(i);
    xvalErr[i] = SignalHist->GetBinWidth(i);

    //    std::cout << "xval " << xval[i] << " efficiency: " << efficiency[i] << " purity: " << purity[i] << " eff * pur "
    //          << effpur[i] << "background Rejection: " << bkRej[i] << " and sig*bkrej: " << sigbk[i] << std::endl;

    //Keep the largest value.
    if (effpur[i] > maxeffpur) {
      max_ind = i;
      maxeffpur = effpur[i];
    }
    if (sigbk[i] > maxsigbk) {
      maxsigbk = sigbk[i];
      maxsigbk_ind = i;
    }
  }

  //Draw the curves
  auto effpur_canvas = new TCanvas("effpur_canvas", "effpur_canvas", 600, 400);
  auto sigbk_canvas = new TCanvas("sigbk_canvas", "sigbk_canvas", 600, 400);

  TGraphErrors* efficiencyPlot = new TGraphErrors(n, xval, efficiency, xvalErr, efficencyErr);
  efficiencyPlot->GetYaxis()->SetTitle("Efficency");
  efficiencyPlot->SetTitle("Efficency");
  efficiencyPlot->SetLineColor(kBlue);
  efficiencyPlot->Draw("a4");
  efficiencyPlot->Write();

  TGraphErrors* bkRejPlot = new TGraphErrors(n, xval, bkRej, xvalErr, bkRejErr);
  bkRejPlot->GetYaxis()->SetTitle("Background Rejection");
  bkRejPlot->SetTitle("Background Rejection");
  bkRejPlot->SetLineColor(kRed);
  bkRejPlot->Draw("a4");
  bkRejPlot->Write();

  TGraphErrors* purityPlot = new TGraphErrors(n, xval, purity, xvalErr, purityErr);
  purityPlot->GetYaxis()->SetTitle("Purity");
  purityPlot->SetTitle("Purity");
  purityPlot->SetLineColor(kRed);
  purityPlot->Draw("a4");
  purityPlot->Write();

  TGraphErrors* effxpurPlot = new TGraphErrors(n, xval, effpur, xvalErr, effpurErr);
  effxpurPlot->GetYaxis()->SetTitle("Effxpur");
  effxpurPlot->SetTitle("Eff*pur");
  effxpurPlot->SetLineColor(kGreen);
  effxpurPlot->Draw("a4");
  effxpurPlot->Write();

  TGraphErrors* sigbkPlot = new TGraphErrors(n, xval, sigbk, xvalErr, sigbkErr);
  sigbkPlot->GetYaxis()->SetTitle("Signal*Background Rejection");
  sigbkPlot->SetTitle("Signal*Background Rejection");
  sigbkPlot->SetLineColor(kGreen);
  sigbkPlot->Draw("a4");
  sigbkPlot->Write();

  //Make the graph to present
  effpur_canvas->cd();
  auto mg1 = new TMultiGraph("mg1", "mg1");
  mg1->Add(efficiencyPlot);
  mg1->Add(purityPlot);
  mg1->Add(effxpurPlot);
  mg1->GetYaxis()->SetRangeUser(0, 1);
  mg1->Write();
  mg1->Draw("AP");
  effpur_canvas->BuildLegend();
  effpur_canvas->Write();

  //Make the graph to present
  sigbk_canvas->cd();
  auto mg2 = new TMultiGraph("mg2", "mg2");
  mg2->Add(efficiencyPlot);
  mg2->Add(bkRejPlot);
  mg2->Add(sigbkPlot);
  mg2->GetYaxis()->SetRangeUser(0, 1);
  mg2->Write();

  mg2->Draw("AP");
  sigbk_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
  sigbk_canvas->Write();

  double maxScale = TMath::Max(BackgroundHist->GetBinContent(BackgroundHist->GetMaximumBin()),
			       SignalHist->GetBinContent(SignalHist->GetMaximumBin()));

  SignalHist->SetLineColor(kBlue);
  SignalHist->SetFillStyle(3003);
  SignalHist->SetFillColor(6);
  SignalHist->Scale(1. / maxScale);
  SignalHist->SetTitle("Signal");
  SignalHist->Write();
  BackgroundHist->SetLineColor(kRed);
  BackgroundHist->SetFillStyle(3003);
  BackgroundHist->SetFillColor(42);
  BackgroundHist->Scale(1. / maxScale);
  BackgroundHist->SetTitle("Background");
  BackgroundHist->Write();
  
  auto sigbk_cut_canvas = new TCanvas("sigbk_cut_canvas", "sigbk_cut_canvas", 600, 400);
  sigbk_cut_canvas->cd();
  mg2->GetHistogram()->SetTitle(Form("Cut at %0.1f " + units + " with efficiency %0.2f and background rejection %0.2f", xval[maxsigbk_ind], efficiency[maxsigbk_ind], bkRej[maxsigbk_ind]));
  mg2->GetXaxis()->SetTitle(axisTitle);
  mg2->Draw("AP");
  SignalHist->Draw("HIST SAME");
  BackgroundHist->Draw("HIST SAME");
  sigbk_cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
  TLine* l1 = new TLine(xval[maxsigbk_ind], 0, xval[maxsigbk_ind], 1);
  l1->SetLineColor(kBlack);
  l1->Draw("same");

  sigbk_cut_canvas->Write();

  //Draw the overlayed samples with the cut.
  auto cut_canvas = new TCanvas("cut_canvas", "cut_canvas", 600, 400);
  cut_canvas->cd();
  mg1->GetHistogram()->SetTitle(Form("Cut at %0.1f " + units + " with efficiency %0.2f and purity %0.2f", xval[max_ind], efficiency[max_ind], purity[max_ind]));
  mg1->GetXaxis()->SetTitle(axisTitle);
  mg1->Draw("AP");
  SignalHist->Draw("HIST SAME");
  BackgroundHist->Draw("HIST SAME");
  cut_canvas->Update();
  cut_canvas->BuildLegend();

  TLine* l = new TLine(xval[max_ind], 0, xval[max_ind], 1);
  l->SetLineColor(kBlack);
  l->Draw("same");

  cut_canvas->Write();

  std::cout << "Best Cut is at: " << xval[max_ind] << std::endl;
  std::cout << "Efficiency: " << efficiency[max_ind] << " +- " << efficencyErr[max_ind] << std::endl;
  std::cout << "Purity    : " << purity[max_ind] << " +- " << purityErr[max_ind] << std::endl;
  std::cout << "Eff*Purity: " << effpur[max_ind] << " +- " << effpurErr[max_ind] << std::endl;

  std::cout << "Best sig*bkrej Cut is at: " << xval[maxsigbk_ind] << std::endl;
  std::cout << "Efficiency:     " << efficiency[maxsigbk_ind] << " +- " << efficencyErr[maxsigbk_ind] << std::endl;
  std::cout << "Background Rej: " << bkRej[maxsigbk_ind] << " +- " << bkRejErr[maxsigbk_ind] << std::endl;
  std::cout << "Eff*Purity:     " << sigbk[maxsigbk_ind] << " +- " << sigbkErr[maxsigbk_ind] << std::endl;

  delete efficiencyPlot;
  delete bkRejPlot;
  delete purityPlot;
  delete effxpurPlot;
  delete sigbkPlot;

  delete mg1;
  delete mg2;

  delete effpur_canvas;
  delete sigbk_canvas;
  delete sigbk_cut_canvas;
  delete cut_canvas;
}


class optimiser::NueRecoOptimiser {
  
public:

  //Branches from Reco Efficiency Finder.
  struct EfficiencyTree {
    std::map<std::string, optimiser::BranchType*> Tree;
    
    template <class T> 
    void AddBranch(std::string Name){
      Tree[Name] = new optimiser::BranchType<T>();
      return;
    }

    void SetBranchAddresses(TTree* tree){
      for(auto const& Branch: Tree){
	tree_signal->SetBranchAddress(Branch.first,&(Branch.second->GetPointer()));
      }
      return;
    }
  };


  EfficiencyTree SignalTree;
  EfficiencyTree BackgroundTree;

  std::map<std::string,optimiser::MetricHolder> MetricMap;
  

};

optimiser::NueRecoOptimiser::NueRecoOptimiser(){

  SignalTree.AddBranch<std::vector<float>>("number_of_showers_per_neutrino");
  SignalTree.AddBranch<std::vector<float>>("vertex_recoKE");
  SignalTree.AddBranch<std::vector<float>>("vertex_trueKE");
  SignalTree.AddBranch<std::vector<float>>("vertex_reco");
  SignalTree.AddBranch<std::vector<float>>("nu_reco_energy");
  SignalTree.AddBranch<std::vector<float>>("nu_truth_energy");
  SignalTree.AddBranch<std::vector<float>>("nu_interaction_type");
  SignalTree.AddBranch<std::vector<float>>("nu_mode");
  SignalTree.AddBranch<std::vector<float>>("nu_E");
  SignalTree.AddBranch<std::vector<float>>("nu_E_numtrue");
  SignalTree.AddBranch<std::vector<float>>("nu_distance");
  SignalTree.AddBranch<std::vector<std::vector<float> >>("truepionE");
  SignalTree.AddBranch<std::vector<std::vector<float> >>("trueprotonE");
  SignalTree.AddBranch<std::vector<std::vector<float> >>("truekaonE");
  SignalTree.AddBranch<std::vector<std::vector<float> >>("truetrackE");
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_energy");
  SignalTree.AddBranch<std::vector<std::vector<float> >>("truth_pid");
  SignalTree.AddBranch<std::vector<std::vector<float> >>("true_energy");
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_coversion_gap");
  SignalTree.AddBranch<std::vector<std::vector<float> >>("shower_dEdx");
  SignalTree.AddBranch<std::vector<std::vector<float> >>("track_lengths");
  SignalTree.AddBranch<std::vector<std::vector<float> >>("track_E");
  SignalTree.AddBranch<std::vector<std::vector<float> >>("track_trueE");
  SignalTree.AddBranch<std::vector<std::vector<float> >>("track_pdg"); 
  SignalTree.AddBranch<std::vector<std::vector<float> >>("track_resE");
  SignalTree.AddBranch<int>("trueShower_num");
  SignalTree.AddBranch<int>("numtrueVtx_branch");
  SignalTree.AddBranch<float>("POT");

  BackgroundTree.AddBranch<std::vector<float>>("number_of_showers_per_neutrino");
  BackgroundTree.AddBranch<std::vector<float>>("vertex_recoKE");
  BackgroundTree.AddBranch<std::vector<float>>("vertex_trueKE");
  BackgroundTree.AddBranch<std::vector<float>>("vertex_reco");
  BackgroundTree.AddBranch<std::vector<float>>("nu_reco_energy");
  BackgroundTree.AddBranch<std::vector<float>>("nu_truth_energy");
  BackgroundTree.AddBranch<std::vector<float>>("nu_interaction_type");
  BackgroundTree.AddBranch<std::vector<float>>("nu_mode");
  BackgroundTree.AddBranch<std::vector<float>>("nu_E");
  BackgroundTree.AddBranch<std::vector<float>>("nu_E_numtrue");
  BackgroundTree.AddBranch<std::vector<float>>("nu_distance");
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("truepionE");
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("trueprotonE");
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("truekaonE");
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("truetrackE");
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_energy");
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("truth_pid");
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("true_energy");
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_coversion_gap");
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("shower_dEdx");
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("track_lengths");
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("track_E");
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("track_trueE");
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("track_pdg"); 
  BackgroundTree.AddBranch<std::vector<std::vector<float> >>("track_resE");
  BackgroundTree.AddBranch<int>("trueShower_num");
  BackgroundTree.AddBranch<int>("numtrueVtx_branch");
  BackgroundTree.AddBranch<float>("POT");

}

void optimiser::NueRecoOptimiser::InitialiseMetric(std::string Name; std::string branchname,int numbins, float xmin, float xmax){
  MetricMap[Name] = optimiser::MetricHolder(Name,branchname,numbins,xmin,xmax);
  return;
}

void optimiser::NueRecoOptimiser::SetBranchAddresses(TTree* tree_signal, TTree* tree_background){
  SignalTree.SetBranchAddresses(tree_signal);
  BackgroundTree.SetBranchAddresses(tree_background);
}

void optimiser::NueRecoOptimiser::FillData(TTree* tree_signal, TTree* tree_background){

  //Fill The Background
  for (Long64_t evt = 0; evt < tree_background->GetEntries(); evt++) {
    tree_background->GetEntry(evt);
    for(auto const& Metric: MetricMap){
      std::vector<float> vals;
      std::vector<std::pair<float,float> > TwoDvals;
      int err = PerformCuts(Metric.first,"background",vals,TwoDvals);
      if(err){return;}
      for(auto const& val: vals){
	Metric.second.FillBackground(val);
      }
      for(auto const& TwoDval: TwoDvals){
	Metric.second.FillBackground(TwoDval.first,TwoDval.second);
      }
    }
  }
  
  //Fill The Signal
  for (Long64_t evt = 0; evt < tree_signal->GetEntries(); evt++) {
    tree_signal->GetEntry(evt);
    for(auto const& Metric: MetricMap){
      std::vector<float> vals;
      std::vector<std::pair<float,float> > TwoDvals;
      int err = PerformCuts(Metric,"signal",vals,TwoDvals);
      if(err){return;}
      for(auto const& val: vals){
	Metric.second.FillSignal(val);
      }
      for(auto const& TwoDval: TwoDvals){
	Metric.second.FillSignal(TwoDval.first,TwoDval.second);
      }
    }
  }
  return;
} 

void optimiser::NueRecoOptimiser::AnalyseData(){

  TFile* file = new TFile("CutFile.root", "RECREATE");

  //Make new directory for plots 
  for(auto const& Metric: MetricMap){
    gDirectory->mkdir((Metric.GetMetricName()).c_str());
  }

  //Loop over metrics and make the graphs.
  for(auto const& Metric: MetricMap){
    file->cd((Metric.GetMetricName()).c_str());
    //Perform the 1D analysis
    Metric.Perform1DCutFinder();
    //Perform the 2D analysis if any.
    if(Metric.GetPartnerMetricName() == ""){return;}
    Metric.Perform2DCutFinder();
  }

  file->Close();
  return;
}

int optimiser::NueRecoOptimiser::PerformCuts(optimiser::MetricHolder& Metric, std::string& sample,std::vector<float>& vals, std::vector<std::pair<float,float> >& TwoDvals){

  int err=0;

  //Make Efficiency Hists
  //Make Resolutions Hits

  if(Metric.GetMetricName() == "VisibleEnergy"){err = StandardDaughterAnalysis(Metric,sample,vals);}
  else{err = StandardAnalysis(Metric,vals);}

  return err;
}

optimiser::BranchType* optimiser::NueRecoOptimiser::GetBranch(std::string& branchname, std::string& sample){

  optimiser::BranchType* branch_vals;

  //Decide if it is a signal or background branch
  if(sample == "signal"){
    if(SignalTree.find(branchname) == SignalTree.end()){
      std::cerr << "Signal Tree branch doesn't exist. You probably mistypes the name in the metric initialser" << std::endl;
      return 1;
    }
    branch_vals = SignalTree[branchname];
  }
  else{
    if(BackgroundTree.find(branchname) == BackgroundTree.end()){
      std::cerr << "Background Tree branch doesn't exist. You probably mistypes the name in the metric initialser" << std::endl;
      return 1;
    }
    branch_vals = BackgroundTree[branchname];
  }
  return branch_vals;
}


int optimiser::NueRecoOptimiser::StandardDaughterAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals, std::vector<std::pair<float,float> >& TwoDvals){

  int err = StandardDaughter1DAnalysis(Metric,sample,vals);
  if(err){
    std::cerr << "Error or in 1D daughter analysis and bailing" << std::endl;
    return 1;
  }

  if(Metric.GetPartnerMetricName() == ""){return 0;}

  err =  StandardDaughter2DAnalysis(Metric,sample,TwoDvals);
  if(err){
    std::cerr << "Error or in 2D daughter analysis and bailing" << std::endl;
    return 1;
  }
  return 0;
}


int optimiser::NueRecoOptimiser::StandardDaughter1DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<float>& vals){

  std::string branchname = Metric.GetBranchName();
  
  optimiser::BranchType* branch_vals = GetBranch(branchname,sample);

  //Loop over the neutrinos
  for(auto const& neutrino: metric){
    //Add the largest shower in the event. The metric are ordered this way
    if(neutrino.size() > 0){
      vals.push_back(neutrino[0]);
    }
  }
  
  return 0;
}

int optimiser::NueRecoOptimiser::StandardDaughter2DAnalysis(optimiser::MetricHolder& Metric, std::string& sample, std::vector<std::pair<float,float> >& TwoDvals){

  //Get the metric branch 
  std::string branchname = Metric.GetBranchName();
  optimiser::BranchType* branch_vals = GetBranch(branchname,sample);

  //Standard analysis assumes  
  std::vector<std::vector<float> > * metric  = dynamic_cast<std::vector<std::vector<float> >*> (*branch_vals->GetPointer());
  
  //Check the metric is of the type expected 
  if(metric == NULL){
    std::cerr << "Metric: " << Metric.GetMetricName() << " has gone under the standard analysis but is not the correct type" << std::endl;
    return 1;
  }

  //Get the partner metric branch
  std::string partnername = Metric.GetPartnerMetricName();
  optimiser::BranchType* partner_branch_vals = GetBranch(partnername,sample);

  //Standard analysis assumes  
  std::vector<std::vector<float> > * partner_metric  = dynamic_cast<std::vector<std::vector<float> >*> (*partner_branch_vals->GetPointer());
  
  //Check the metric is of the type expected 
  if(partner_metric == NULL){
    std::cerr << "Partner Metric: " << partnername << " has gone under the standard analysis but is not the correct type" << std::endl;
    return 1;
  }

  //Loop over the neutrinos
  for(int iter=0; iter<metric.size(); ++iter){
    //Add the largest shower in the event. The metric are ordered this way
    if(metric[iter].size() > 0 && partner_metric[iter].size() > 0){
      std::pair<float,float> vals = std::make_pair(metric[iter][0],partner_metric[iter][0]);
      TwoDvals.push_back(vals);
    }
  }
  return 0;
}


int main(int argc, char* argv[]){
  
  //Read The files.
 std:string signalfile;
  std::string backgroundfile;

  static struct option long_options[] = {
    {"signalfile"     , required_argument,  0,  's' },
    {"backgroundfiles", required_argument,  0,  'b' },
    {0, 0}
  };
  
  int opt = 0;
  int c = 0;
  while ((opt = getopt_long_only(argc, argv,"", long_options, &c )) != -1) {
    switch (opt) {
    case 's' : signalfile = optarg;
      break;
    case 'b' : backgroundfile = optarg;
      break;
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

  //FIX LATER
  TDirectory* dir_eff = f_eff->GetDirectory("eff");
  TTree *tree_sig = (TTree*)dir_eff->Get("RecoEffMetricTree");  

  //Get the histogram for the background
  TFile *f_bk = new TFile(backgroundfile.c_str());
  TDirectory* dir_bk = f_bk->GetDirectory("eff");
  TTree *tree_bk = (TTree*)dir_bk->Get("RecoEffMetricTree");  


  //Make New Optimiser Object
  optimiser::NueRecoOptimiser NueOptimiser();

  //Set Up the metrics.
  NueOptimiser.InitialiseMetric("VisibleEnergy","shower_energy_branch",100,0,3000);

  NueOptimiser.SetBranchAddress(tree_sig,tree_bk);

  //Read the data
  NueOptimiser.FillData(tree_sig,tree_bk);

  TFile* file = new TFile("CutFile.root", "RECREATE");

  //Make the pretty graphs 
  NueOptimiser.AnalyseData();

  return 0;
}



void MakeEfficiencyPlots(TH1F* postcut, TH1F* precut, TH1F* extra, TString postcutname, TString precutname, TString extraname, TString title, TString eff_name){

  float maxSignal = precut->GetBinContent(precut->GetMaximumBin());
  
  if(TEfficiency::CheckConsistency(*postcut,*precut)){
    //Make the eff plot.
    TEfficiency* reco_eff = new TEfficiency(*postcut,*precut);
    reco_eff->SetTitle(title);
    reco_eff->SetName(eff_name);
    reco_eff->Write();

    postcut->Scale(1./maxSignal);
    precut->Scale(1./maxSignal);
      
    postcut->SetLineColor(kBlue);
    postcut->SetFillStyle(3003);
    postcut->SetFillColor(6);
    postcut->SetTitle(postcutname);
    postcut->Write();
    precut->SetLineColor(kRed);
    precut->SetFillStyle(3003);
    precut->SetFillColor(42);
    precut->SetTitle(precutname);
    precut->Write();

    TString canavs_name = eff_name + "_canvas";
    auto cut_canvas = new TCanvas(canavs_name, canavs_name, 600, 400);
    precut->Draw("HIST SAME");
    cut_canvas->cd();
    reco_eff->Draw("same");
    postcut->Draw("HIST SAME");
    cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
    cut_canvas->Write();
    
      
    postcut->Scale(maxSignal);
    precut->Scale(maxSignal);
  }

  return;

  if(TEfficiency::CheckConsistency(*extra,*precut)){
    //Make the eff plot.
    TEfficiency* reco_eff = new TEfficiency(*extra,*precut);
    reco_eff->SetTitle(title);
    reco_eff->SetName(eff_name);
    reco_eff->Write();

    extra->Scale(1./maxSignal);
    precut->Scale(1./maxSignal);
      
    extra->SetLineColor(kBlue);
    extra->SetFillStyle(3003);
    extra->SetFillColor(6);
    extra->SetTitle(extraname);
    extra->Write();
    precut->SetLineColor(kRed);
    precut->SetFillStyle(3003);
    precut->SetFillColor(42);
    precut->SetTitle(precutname);
    precut->Write();

    TString canavs_name = eff_name + "_canvas_extra";
    auto cut_canvas = new TCanvas(canavs_name, canavs_name, 600, 400);
    precut->Draw("HIST SAME");
    cut_canvas->cd();
    reco_eff->Draw("same");
    extra->Draw("HIST SAME");
    cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
    cut_canvas->Write();
    
      
    extra->Scale(maxSignal);
    precut->Scale(maxSignal);
  }



}

void MakeResolutionGraphs(std::vector<float>& vals, std::vector<float>& energyval, float min, float max, float num_steps, float width,TString yTitle, TString Metric, TString name, float nbins){


  float step_size = (max-min)/(num_steps);

  float x[(int) num_steps];
  float dx[(int)num_steps];
  float y[(int)num_steps];
  float dy[(int)num_steps];
  float ddy[(int)num_steps];
  float dddy[(int)num_steps];

  TString TrueEnergy = "track_trueE";

  //Loop over the energy step.
  for(int i=0; i<num_steps; ++i){
    
    //Get the widths to access 
    float energylow = min + i*step_size - width/2;
    float energymax = min + i*step_size + width/2;

    TString tempHistname = Metric + "tempHist " + std::to_string(energylow);  
    TH1F* tempHist = new TH1F(tempHistname, tempHistname, nbins, -2, 2);
    for(int val=0; val<vals.size(); ++val){
      //      std::cout << "energyval[val]: " << energyval[val] << " energylow: " << energylow << " energymax: " << energymax  << std::endl;
      if(energyval[val] > energylow && energyval[val] < energymax){
	tempHist->Fill(vals[val]);
      }
    }
    x[i] = (energylow + energymax) / 2.0;
    dx[i] = (energymax - energylow) / 2.0;
    y[i] = tempHist->GetMean();
    dy[i] = tempHist->GetRMS();
    ddy[i] = tempHist->GetRMS();
    dddy[i] = dy[i]/TMath::Sqrt(2*(tempHist->Integral()-1));

    // if(tempHist->GetEntries() > 20){
    //   TF1* fit = new TF1("fit","gaus");
    //   tempHist->Fit("fit");
    //   y[i]  = fit->GetParameter(1);
    //   dy[i] = fit->GetParError(1);
    //   ddy[i]  = fit->GetParameter(2);
    //   dddy[i] = fit->GetParError(2);
      
    //   tempHist->Clear();
    //   for(int val=0; val<vals.size(); ++val){
    // 	//      std::cout << "energyval[val]: " << energyval[val] << " energylow: " << energylow << " energymax: " << energymax  << std::endl;
    // 	if(energyval[val] > energylow && energyval[val] < energymax){
    // 	  if(TMath::Abs(vals[val]-y[i]) <  ddy[i]){ 
    // 	    tempHist->Fill(vals[val]);
    // 	  }
    // 	}
    //   }
    //   delete fit;
    // }

    // if(tempHist->GetEntries() > 20){
    //   TF1* fit1 = new TF1("fit1","gaus");
    //   tempHist->Fit("fit1");
    //   y[i]  = fit1->GetParameter(1);
    //   dy[i] = fit1->GetParError(1);
    //   ddy[i]  = fit1->GetParameter(2);
    //   dddy[i] = fit1->GetParError(2);
    //   tempHist->Write();
    //   delete fit1;
    // }
    
    //    y[i] = tempHist->GetMean();
    //dy[i] = tempHist->GetRMS();
    //ddy[i] = dy[i]/TMath::Sqrt(2*(tempHist->Integral()-1));
    delete tempHist;

    std::cout << "x: " << x[i] << " y: " << y[i]  << " dy: " << dy[i] << " ddy: " << ddy[i] << " dddy[i]: " << dddy[i] << std::endl;
  }


  TString yTitle_mean = yTitle + " mean";
  TString mean_name = yTitle_mean + " " + name;
  TGraphErrors* tempGraph = new TGraphErrors(num_steps, x, y, dx, dy);
  //  tempGraphs->SetLineColor();
  //  tempGraphs->SetMarkerColor(pass + 1);
  tempGraph->GetXaxis()->SetTitle("True Track Energy (MeV)");
  tempGraph->GetYaxis()->SetTitle(yTitle_mean);
  tempGraph->SetTitle(mean_name);
  tempGraph->Write();
 
  TString yTitle_rms = yTitle + " rms"; 
  TString rms_name = yTitle_rms + " " + name;
  TGraphErrors* tempRMSGraph = new TGraphErrors(num_steps, x, ddy, dx, dddy);
  //  tempRMSGraphs->SetLineColor(pass + 1);
  //tempRMSGraphs->SetMarkerColor(pass + 1);
  tempRMSGraph->GetXaxis()->SetTitle("True True Energy (MeV)");
  tempRMSGraph->GetYaxis()->SetTitle(yTitle_rms);
  tempRMSGraph->SetTitle(rms_name);
  

  TString yTitle_mean_canvas = yTitle + " mean_canvas" + name;    
  TCanvas* cMeanGraph = new TCanvas(yTitle_mean_canvas);
  tempGraph->Draw();
  cMeanGraph->Write();

  TString yTitle_rms_canvas = yTitle + " rms_canvas" + name;    
  TCanvas* cRMSGraph = new TCanvas(yTitle_rms_canvas);
  tempRMSGraph->Draw();
  cRMSGraph->Write();
}




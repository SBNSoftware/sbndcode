//Quick Tool to finder a where a cut should be made

#include "TH1.h"
#include "TFile.h"
#include "TMath.h" 

void DoCutFinding(TH1D* signalhist, TH1D* backgroundhist){

  signalhist->Scale(1/signalhist->GetEntries());
  backgroundhist->Scale(1/backgroundhist->GetEntries());

  //Find the min and max start positions 
  float minsig  = signalhist->GetBinCenter(0); 
  float maxsig  = signalhist->GetBinCenter(signalhist->GetNbinsX());
  float minbk   = backgroundhist->GetBinCenter(0);
  float maxbk   = backgroundhist->GetBinCenter(backgroundhist->GetNbinsX());

  if(minsig != minbk){
    std::cout << "min point does not match" << std::endl;
    return;
  }
  if(maxbk != maxbk){
    std::cout << "max point does not match" << std::endl;
    return;
  }
  if(signalhist->GetNbinsX() != backgroundhist->GetNbinsX()){
    std::cout << "bin numbers sig " << signalhist->GetNbinsX() << " and background " << backgroundhist->GetNbinsX() << "do not match" << std::endl;
    return;
  }

  int max_ind = -999;
  int n = signalhist->GetNbinsX();
  double sigtot_err = 0; 
  double sigtot = signalhist->IntegralAndError(1,n,sigtot_err);
  double maxeffpur = -999;
  double efficiency[n];
  double efficencyErr[n];
  double purity[n];
  double purityErr[n];
  double effpur[n];
  double effpurErr[n];
  double xval[n];
  double xvalErr[n];

  //Get the values.
  for(int i=1; i<n; ++i){
    
    double presig_err = 0;
    double prebk_err  = 0;
    double presig = signalhist->IntegralAndError(1,i,presig_err);
    double prebk  = backgroundhist->IntegralAndError(1,i,prebk_err);

    //std::cout << "presig: " << presig << " prebk: " << prebk << std::endl;

    efficiency[i]   = presig/sigtot;
    efficencyErr[i] = efficiency[i]*TMath::Sqrt(((presig_err/presig)*(presig_err/presig)) + (sigtot_err/sigtot)*(sigtot_err/sigtot));

    if(presig+prebk != 0){
    purity[i]       = presig/(presig+prebk);
    double purityErr_dem   = TMath::Sqrt((presig_err*presig_err) + (prebk_err*prebk_err));
    purityErr[i]    = purity[i]*TMath::Sqrt(((presig_err/presig)*(presig_err/presig)) + (purityErr_dem/(presig+prebk))*(purityErr_dem/(presig+prebk)));
    }
    else{
      purity[i]  = 0;
      purityErr[i] = 0;
    }

    effpur[i]     = efficiency[i]*purity[i];
    if(efficiency[i] == 0 || purity[i] == 0){
      effpurErr[i] = 0;
    }
    else{
      effpurErr[i] = effpur[i]*TMath::Sqrt(((efficencyErr[i]/efficiency[i])*(efficencyErr[i]/efficiency[i])) + (purityErr[i]/purity[i])*(purityErr[i]/purity[i]));
    }
    xval[i]       = signalhist->GetBinCenter(i);
    xvalErr[i]    = signalhist->GetBinWidth(i);

    std::cout << "xval " << xval[i] << " efficiency: " << efficiency[i] << " purity: " << purity[i] << " eff * pur " << effpur[i] << std::endl;


    //Keep the largest value.
    if(effpur[i] > maxeffpur){
      max_ind = i;
      maxeffpur= effpur[i];
    }
  }
  
  TFile *file  = new TFile("CutFile.root", "RECREATE");
  
  //Draw the curves 
  auto roc_canvas = new TCanvas("roc_canvas","roc_canvas",600,400);

  TGraphErrors* efficiencyPlot = new TGraphErrors(n, xval, efficiency, xvalErr, efficencyErr);
  efficiencyPlot->GetYaxis()->SetTitle("Efficency");
  efficiencyPlot->SetTitle("Efficency");
  efficiencyPlot->SetLineColor(kBlue);
  efficiencyPlot->Draw("a4");
  efficiencyPlot->Write();

  TGraphErrors* purityPlot = new TGraphErrors(n, xval, purity,  xvalErr, purityErr);
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

  //Make the graph to present
  auto mg = new TMultiGraph();
  mg->Add(efficiencyPlot);
  mg->Add(purityPlot);
  mg->Add(effxpurPlot);
  mg->Write();

  mg->Draw("AP");
  roc_canvas->BuildLegend();
  roc_canvas->Write();


  //Draw the overlayed samples with the cut.
  auto cut_canvas = new TCanvas("cut_canvas","cut_canvas",600,400);
  cut_canvas->cd();
  
  //  signalhist->SetLineColor(kBlue);
  signalhist->SetFillStyle(3003);
  signalhist->SetFillColor(6);
  signalhist->GetYaxis()->SetRangeUser(0,1);
  signalhist->Draw("histsame");
  backgroundhist->SetLineColor(kRed);
  backgroundhist->SetFillStyle(3003);
  backgroundhist->SetFillColor(42);
  backgroundhist->GetYaxis()->SetRangeUser(0,1);
  backgroundhist->Draw("histsame");
  cut_canvas->Update();

  double ymax = signalhist->GetBinContent(signalhist->GetMaximumBin());
  if(backgroundhist->GetBinContent(backgroundhist->GetMaximumBin()) >
     signalhist->GetBinContent(signalhist->GetMaximumBin()))
    ymax = backgroundhist->GetBinContent(backgroundhist->GetMaximumBin());

  TLine *l=new TLine(xval[max_ind],0,xval[max_ind],ymax);
  l->SetLineColor(kBlack);
  l->Draw("same");
  
  efficiencyPlot->Draw("same");
  purityPlot->Draw("same");
  effxpurPlot->Draw("same");
  

  cut_canvas->BuildLegend();

  cut_canvas->Write();

  std::cout << "Best Cut is at: " << xval[max_ind] << std::endl;
  std::cout << "Efficiency: " << efficiency[max_ind] << " +- " << efficencyErr[max_ind] << std::endl;
  std::cout << "Purity    : " << purity[max_ind] << " +- " << purityErr[max_ind] << std::endl;
  std::cout << "Eff*Purity: " << effpur[max_ind] << " +- " << effpurErr[max_ind] << std::endl;

  file->Close();

}


void CutFinder(const char* signalfile_name, const char* backgroundfile_name, const char* signalhist_name, const char* backgroundhist_name){

  //Get the histogram for the signal
  auto signalfile = new TFile(signalfile_name);
  TH1D *signalhist  = (TH1D*)signalfile->Get(signalhist_name);


  //Get the histogram for the background
  auto backgroundfile = new TFile(backgroundfile_name);
  TH1D *backgroundhist  = (TH1D*)backgroundfile->Get(backgroundhist_name);

  if(signalhist == NULL || backgroundhist == NULL){
    return;
  }


  //Doing the Cut finding.
  DoCutFinding(signalhist,backgroundhist);

  delete signalhist;
  delete backgroundhist;
  signalfile->Close();
  backgroundfile->Close();

}


void CutFinder(const char* file_name, const char* signalhist_name, const char* backgroundhist_name){

  //Get the histogram for the signal
  auto  file  = new TFile(file_name);
  TH1D  *signalhist      = (TH1D*)file->Get(signalhist_name);
  TH1D  *backgroundhist  = (TH1D*)file->Get(backgroundhist_name);
  
  //Doing the Cut finding.
  DoCutFinding(signalhist,backgroundhist);

  //  delete signalhist;
  //delete backgroundhist;


}



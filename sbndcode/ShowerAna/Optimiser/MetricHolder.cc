//Framework includes
#include "MetricHolder.hh"

//Root includes
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TGraph2D.h"
#include "TProfile.h" 
#include "TEfficiency.h"

//C++ includes
#include <iostream>

void optimiser::MetricHolder::Perform2DCutFinder(){

  std::cout << "Perform2DCutFinder" << std::endl;

  //  SignalHistPartner->Scale(1/TotalSigPOT);
  //  BackgroundHistPartner->Scale(1/TotalBKPOT);
  
  int max_ind = -999;
  int maxsigbk_ind = -999;
  int n_x = SignalHistPartner->GetNbinsX();
  int n_y = SignalHistPartner->GetNbinsY();

  double bktot_err = 0;
  double sigtot_err = 0;
  double bktot = BackgroundHistPartner->IntegralAndError(1, n_x+1,1,n_y+1, bktot_err);
  double sigtot = SignalHistPartner->IntegralAndError(1, n_x+1,1,n_y+1, sigtot_err);
  double maxeffpur = -999;
  double maxsigbk = -999;
  double efficiency[n_x*n_y];
  double purity[n_x*n_y];
  double bkRej[n_x*n_y];
  double effpur[n_x*n_y];
  double effpurErr[n_x*n_y];
  double sigbk[n_x*n_y];
  double sigbkErr[n_x*n_y];
  double xval[n_x*n_y];
  double yval[n_x*n_y];

  int maxeffpur_i;
  int maxeffpur_j;
  int maxsigbk_i;
  int maxsigbk_j;

  double maxeffpur_err = 0;
  double maxsigbkErr = 0;
  double maxefficiency_effpur_err = 0;
  double maxpurity_err = 0;
  double maxefficiency_effbk_err = 0;
  double maxbackgnd_err = 0;

  //Get the values.
  for(int j=1; j<n_y+1; ++j){
    for(int i=1; i<n_x+1; ++i){
      
      int ipj = i-1 + n_x*(j-1);

      double presig_err = 0;
      double prebk_err = 0;
      double presig = 0;
      double prebk =  0;
      
      if(fLessThan && fPartnerLessThan && fAND){
	presig = SignalHistPartner->IntegralAndError(1, i, 1, j, presig_err);
	prebk  = BackgroundHistPartner->IntegralAndError(1, i, 1, j, prebk_err);
      }
      else if(!fLessThan && fPartnerLessThan && fAND){
	presig = SignalHistPartner->IntegralAndError(i,n_x+1, 1, j, presig_err);
	prebk  = BackgroundHistPartner->IntegralAndError(i,n_x+1, 1, j, prebk_err);
      }
      else if(fLessThan && !fPartnerLessThan && fAND){
        presig = SignalHistPartner->IntegralAndError(1,i, j, n_y+1, presig_err);
        prebk  = BackgroundHistPartner->IntegralAndError(1,i,j, n_y+1, prebk_err);
      }
      else if(!fLessThan && !fPartnerLessThan && fAND){
        presig = SignalHistPartner->IntegralAndError(i,n_x+1, j,n_y+1, presig_err);
        prebk  = BackgroundHistPartner->IntegralAndError(i,n_x+1,j,n_y+1, prebk_err);
      }
      else if(fLessThan && fPartnerLessThan && !fAND){
	double presig_err_two=0;
	presig = SignalHistPartner->IntegralAndError(1, i, j+1, n_y+1, presig_err) + SignalHistPartner->IntegralAndError(1, n_x+1, 1, j, presig_err_two);
	presig_err = TMath::Sqrt(presig_err*presig_err + presig_err_two*presig_err_two);
	double prebk_err_two=0;
	prebk  = BackgroundHistPartner->IntegralAndError(1, i, j+1, n_y+1, prebk_err) + BackgroundHistPartner->IntegralAndError(1, n_x+1, 1, j, prebk_err_two);
	prebk_err = TMath::Sqrt(prebk_err*prebk_err + prebk_err_two*prebk_err_two);
      }
      else if(!fLessThan && fPartnerLessThan && !fAND){
	double presig_err_two=0;
	presig = SignalHistPartner->IntegralAndError(i, n_x+1, j+1, n_y+1, presig_err) + SignalHistPartner->IntegralAndError(1, n_x+1, 1, j, presig_err_two);
	presig_err = TMath::Sqrt(presig_err*presig_err + presig_err_two*presig_err_two);
	double prebk_err_two=0;
	prebk  = BackgroundHistPartner->IntegralAndError(i, n_x+1, j+1, n_y+1, prebk_err) + BackgroundHistPartner->IntegralAndError(1, n_x+1, 1, j, prebk_err_two);
	prebk_err = TMath::Sqrt(prebk_err*prebk_err + prebk_err_two*prebk_err_two);
      }
      else if(fLessThan && !fPartnerLessThan && !fAND){
	double presig_err_two=0;
	presig = SignalHistPartner->IntegralAndError(1,i, 1, j-1, presig_err) + SignalHistPartner->IntegralAndError(1, n_x+1,j,n_y+1, presig_err_two);
	presig_err = TMath::Sqrt(presig_err*presig_err + presig_err_two*presig_err_two);
	double prebk_err_two=0;
	prebk  = BackgroundHistPartner->IntegralAndError(1,i,1, j-1, prebk_err) + BackgroundHistPartner->IntegralAndError(1, n_x+1,j,n_y+1, prebk_err_two);
	prebk_err = TMath::Sqrt(prebk_err*prebk_err + prebk_err_two*prebk_err_two);
      }
      else if(!fLessThan && !fPartnerLessThan && !fAND){
	double presig_err_two=0;
	presig = SignalHistPartner->IntegralAndError(i,n_x+1,1, j-1, presig_err) + SignalHistPartner->IntegralAndError(1, n_x+1,j,n_y+1, presig_err_two);
	presig_err = TMath::Sqrt(presig_err*presig_err + presig_err_two*presig_err_two);
	double prebk_err_two=0;
	prebk  = BackgroundHistPartner->IntegralAndError(i,n_x+1,1, j-1, prebk_err) + BackgroundHistPartner->IntegralAndError(1, n_x+1,j,n_y+1, prebk_err_two);
	prebk_err = TMath::Sqrt(prebk_err*prebk_err + prebk_err_two*prebk_err_two);
      }



      double efficencyErr = 0;
      double bkRejErr     = 0;
      double purityErr    = 0; 
      double effpurErr    = 0;
      double sigbkErr     = 0;

      efficiency[ipj] = presig / sigtot;
      efficencyErr = efficiency[ipj] * TMath::Sqrt(((presig_err / presig) * (presig_err / presig)) + (sigtot_err / sigtot) * (sigtot_err / sigtot));
      
      bkRej[ipj] = 1 - prebk / bktot;
      bkRejErr = bkRej[ipj] * TMath::Sqrt((prebk_err / prebk) * (prebk_err / prebk) + (bktot_err / bktot) * (bktot_err / bktot));
      
      if (presig + prebk != 0) {
	purity[ipj] = presig / (presig + prebk);
	double purityErr_dem = TMath::Sqrt((presig_err * presig_err) + (prebk_err * prebk_err));
	purityErr = purity[ipj] * TMath::Sqrt(((presig_err / presig) * (presig_err / presig)) + (purityErr_dem / (presig + prebk)) * (purityErr_dem / (presig + prebk)));
      } 
      else {
      purity[ipj] = 0;
      purityErr = 0;
      }
      
      effpur[ipj] = efficiency[ipj] * purity[ipj];
      if (efficiency[ipj] == 0 || purity[ipj] == 0) {
	effpurErr = 0;
      } 
      else {
	effpurErr = effpur[ipj] * TMath::Sqrt(((efficencyErr / efficiency[ipj]) * (efficencyErr / efficiency[ipj])) + (purityErr / purity[ipj]) * (purityErr / purity[ipj]));
      }
      
      sigbk[ipj] = efficiency[ipj] * bkRej[ipj];
      if (efficiency[ipj] == 0 || bkRej[ipj] == 0) {
	sigbkErr = 0;
      } 
      else {
	sigbkErr = sigbk[ipj] * TMath::Sqrt(((efficencyErr / efficiency[ipj]) * (efficencyErr / efficiency[ipj])) + (bkRejErr / bkRej[ipj]) * (bkRejErr / bkRej[ipj]));
      }

      xval[ipj] = ((TAxis*)SignalHistPartner->GetXaxis())->GetBinCenter(i);
      yval[ipj] = ((TAxis*)SignalHistPartner->GetYaxis())->GetBinCenter(j);

      //Keep the largest value.
      if (effpur[ipj] > maxeffpur) {
	max_ind = ipj;
	maxeffpur = effpur[ipj];
	maxeffpur_err = sigbkErr;
	maxefficiency_effpur_err = efficencyErr;
	maxpurity_err = purityErr;
	maxeffpur_i = i;
	maxeffpur_j = j;
      }
      if (sigbk[ipj] > maxsigbk) {
	maxsigbk = sigbk[ipj];
	maxsigbk_ind = ipj;
	maxsigbkErr = sigbkErr;
	maxefficiency_effbk_err = efficencyErr;
	maxbackgnd_err = bkRejErr;
	maxsigbk_i = i;
	maxsigbk_j = j;
      }
    }
  }

  //Draw the curves
  TGraph2D *efficency_graph  = new TGraph2D("eff2D","eff2D",n_x*n_y,xval,yval,efficiency);
  TGraph2D *background_graph = new TGraph2D("bkrej2D","bjreg2D",n_x*n_y,xval,yval,bkRej);
  TGraph2D *purity_graph     = new TGraph2D("purity2D","purity2D",n_x*n_y,xval,yval,purity);  
  TGraph2D *effbk_graph      = new TGraph2D("effbk2D","effbk2D",n_x*n_y,xval,yval,sigbk);
  TGraph2D *effpur_graph     = new TGraph2D("effpur2D","effpur2D",n_x*n_y,xval,yval,effpur);

  //Draw the histograms 
  SignalHistPartner->Write();
  BackgroundHistPartner->Write();

  auto signal_canvas = new TCanvas("signal2D_canvas", "signal2D_canvas", 600, 400);
  SignalHistPartner->Draw("colz");
  signal_canvas->Write();

  auto bkgrd_canvas = new TCanvas("background2D_canvas", "background2D_canvas", 600, 400);
  BackgroundHistPartner->Draw("colz");
  bkgrd_canvas->Write();

  //Draw the profiles
  TProfile* SignalProfileX = SignalHistPartner->ProfileX();
  TProfile* SignalProfileY = SignalHistPartner->ProfileY();
  TProfile* BackgroundProfileX = BackgroundHistPartner->ProfileX();
  TProfile* BackgroundProfileY = BackgroundHistPartner->ProfileY();

  SignalProfileX->Write();
  SignalProfileY->Write();
  BackgroundProfileX->Write();
  BackgroundProfileY->Write();

  efficency_graph->Write();
  background_graph->Write();
  purity_graph->Write();
  effbk_graph->Write();
  effpur_graph->Write();
  
  auto efficency_canvas = new TCanvas("eff2D_canvas", "eff2D_canvas", 600, 400);
  efficency_graph->Draw("colz");
  efficency_canvas->Write();

  auto background_canvas = new TCanvas("bk2D_canvas", "bk2D_canvas", 600, 400);
  background_graph->Draw("colz");
  background_canvas->Write();

  auto purity_canvas = new TCanvas("pur2D_canvas", "pur2D_canvas", 600, 400);
  purity_graph->Draw("colz");
  purity_canvas->Write();

  auto effbk_canvas = new TCanvas("effbk2D_canvas", "effbk2D_canvas", 600, 400);
  effbk_graph->Draw("colz");
  effbk_canvas->Write();

  auto effpur_canvas = new TCanvas("effpur2D_canvas", "effpur2D_canvas", 600, 400);
  effpur_graph->Draw("colz");
  effpur_canvas->Write();

  delete effpur_canvas;
  delete effbk_canvas;
  delete purity_canvas;
  delete background_canvas;
  delete efficency_canvas;
  delete efficency_graph;
  delete background_graph;
  delete purity_graph;
  delete effbk_graph;
  delete effpur_graph;

  //Make the 1D projections
  if(fPartnerLessThan){
    Perform1DCutFinder(SignalHistPartner->ProjectionX("_pxsig",1,maxsigbk_j,"e"), BackgroundHistPartner->ProjectionX("_pxbk",1,maxsigbk_j ,"e")); 
  }
  else{
    Perform1DCutFinder(SignalHistPartner->ProjectionX("_pxsig",maxsigbk_j,n_y+1,"e"), BackgroundHistPartner->ProjectionX("_pxbk",maxsigbk_j,n_y+1 ,"e"));
  }
  if(fLessThan){
    Perform1DCutFinder(SignalHistPartner->ProjectionY("_pysig",1,maxsigbk_i,"e"), BackgroundHistPartner->ProjectionY("_pjbk",1,maxsigbk_i,"e")); 
  }
  else{
    Perform1DCutFinder(SignalHistPartner->ProjectionY("_pysig",maxsigbk_i,n_x+1,"e"), BackgroundHistPartner->ProjectionY("_pjbk",maxsigbk_i,n_x+1,"e")); 
  }

  std::cout << "Best Cut is at x: " << xval[max_ind] << " y: " << yval[max_ind] << std::endl;
  std::cout << "Efficiency: " << efficiency[max_ind] << " +- " << maxefficiency_effpur_err << std::endl;
  std::cout << "Purity    : " << purity[max_ind] << " +- " << maxpurity_err << std::endl;
  std::cout << "Eff*Purity: " << effpur[max_ind] << " +- " << maxeffpur_err << std::endl;

  std::cout << "Best sig*bkrej Cut is at x: " << xval[maxsigbk_ind] << " y: " << yval[max_ind] << std::endl;
  std::cout << "Efficiency:     " << efficiency[maxsigbk_ind] << " +- " << maxefficiency_effbk_err << std::endl;
  std::cout << "Background Rej: " << bkRej[maxsigbk_ind] << " +- " << maxbackgnd_err << std::endl;
  std::cout << "Eff*Purity:     " << sigbk[maxsigbk_ind] << " +- " << maxsigbkErr << std::endl;


  return;
}

void optimiser::MetricHolder::Perform1DCutFinder(){
  Perform1DCutFinder(SignalHist,BackgroundHist);
}

void optimiser::MetricHolder::Perform1DCutFinder(TH1D* signalhist, TH1D* backgroundhist){

  // TH1::AddDirectory(kFALSE);
  //  signalhist->Scale(1/TotalSigPOT);
  //  backgroundhist->Scale(1/TotalBKPOT);
  
  //Find the min and max start positions
  float minsig = signalhist->GetBinCenter(0);
  float maxsig = signalhist->GetBinCenter(signalhist->GetNbinsX());
  float minbk = backgroundhist->GetBinCenter(0);
  float maxbk = backgroundhist->GetBinCenter(backgroundhist->GetNbinsX());

  if (minsig != minbk) {
    std::cout << "min point does not match" << std::endl;
    return;
  }
  if (maxbk != maxbk) {
    std::cout << "max point does not match" << std::endl;
    return;
  }
  if (signalhist->GetNbinsX() != backgroundhist->GetNbinsX()) {
    std::cout << "bin numbers sig " << signalhist->GetNbinsX() << " and background " << backgroundhist->GetNbinsX() << "do not match" << std::endl;
    return;
  }

  int max_ind = -999;
  int maxsigbk_ind = -999;
  int n = signalhist->GetNbinsX();
  double bktot_err = 0;
  double bktot = backgroundhist->IntegralAndError(1, n+1, bktot_err);
  double sigtot_err = 0;
  double sigtot = signalhist->IntegralAndError(1, n+1, sigtot_err);
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
    double presig = 0;
    double prebk  = 0;

    if(fLessThan){
      presig = signalhist->IntegralAndError(1, i, presig_err);
      prebk  = backgroundhist->IntegralAndError(1, i, prebk_err);
    }
    else{
      presig = signalhist->IntegralAndError(i, n+1, presig_err);
      prebk  = backgroundhist->IntegralAndError(i, n+1, prebk_err);
    }

    //    std::cout << "presig: " << presig << " prebk: " << prebk << std::endl;

    efficiency[i-1] = presig / sigtot;
    efficencyErr[i-1] = efficiency[i-1] * TMath::Sqrt(((presig_err / presig) * (presig_err / presig)) + (sigtot_err / sigtot) * (sigtot_err / sigtot));

    bkRej[i-1] = 1 - prebk / bktot;
    bkRejErr[i-1] = bkRej[i-1] * TMath::Sqrt((prebk_err / prebk) * (prebk_err / prebk) + (bktot_err / bktot) * (bktot_err / bktot));

    if (presig + prebk != 0) {
      purity[i-1] = presig / (presig + prebk);
      double purityErr_dem = TMath::Sqrt((presig_err * presig_err) + (prebk_err * prebk_err));
      purityErr[i-1] = purity[i-1] * TMath::Sqrt(((presig_err / presig) * (presig_err / presig)) + (purityErr_dem / (presig + prebk)) * (purityErr_dem / (presig + prebk)));
    } else {
      purity[i-1] = 0;
      purityErr[i-1] = 0;
    }

    effpur[i-1] = efficiency[i-1] * purity[i-1];
    if (efficiency[i-1] == 0 || purity[i-1] == 0) {
      effpurErr[i-1] = 0;
    } else {
      effpurErr[i-1] = effpur[i-1] * TMath::Sqrt(((efficencyErr[i-1] / efficiency[i-1]) * (efficencyErr[i-1] / efficiency[i-1])) + (purityErr[i-1] / purity[i-1]) * (purityErr[i-1] / purity[i-1]));
    }

    sigbk[i-1] = efficiency[i-1] * bkRej[i-1];
    if (efficiency[i-1] == 0 || bkRej[i-1] == 0) {
      sigbkErr[i-1] = 0;
    } else {
      sigbkErr[i-1] = sigbk[i-1] * TMath::Sqrt(((efficencyErr[i-1] / efficiency[i-1]) * (efficencyErr[i-1] / efficiency[i-1])) + (bkRejErr[i-1] / bkRej[i-1]) * (bkRejErr[i-1] / bkRej[i-1]));
    }

    xval[i-1] = signalhist->GetBinCenter(i);
    xvalErr[i-1] = signalhist->GetBinWidth(i)/2;

    //    std::cout << "xval " << xval[i-1] << " efficiency: " << efficiency[i-1] << " purity: " << purity[i-1] << " eff * pur "
    //          << effpur[i-1] << "background Rejection: " << bkRej[i-1] << " and sig*bkrej: " << sigbk[i-1] << std::endl;

    //Keep the largest value.
    if (effpur[i-1] > maxeffpur) {
      max_ind = i-1;
      maxeffpur = effpur[i-1];
    }
    if (sigbk[i-1] > maxsigbk) {
      maxsigbk = sigbk[i-1];
      maxsigbk_ind = i-1;
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
  effxpurPlot->SetLineColor(kBlack);
  effxpurPlot->Draw("a4");
  effxpurPlot->Write();

  TGraphErrors* sigbkPlot = new TGraphErrors(n, xval, sigbk, xvalErr, sigbkErr);
  sigbkPlot->GetYaxis()->SetTitle("Signal*Background Rejection");
  sigbkPlot->SetTitle("Signal*Background Rejection");
  sigbkPlot->SetLineColor(kBlack);
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

  double maxScale = TMath::Max(backgroundhist->GetBinContent(backgroundhist->GetMaximumBin()),
			       signalhist->GetBinContent(signalhist->GetMaximumBin()));

  signalhist->SetLineColor(kBlue);
  signalhist->SetFillStyle(3003);
  signalhist->SetFillColor(6);
  signalhist->Scale(1. / maxScale);
  signalhist->SetTitle("Signal");
  signalhist->Write();
  backgroundhist->SetLineColor(kRed);
  backgroundhist->SetFillStyle(3003);
  backgroundhist->SetFillColor(42);
  backgroundhist->Scale(1. / maxScale);
  backgroundhist->SetTitle("Background");
  backgroundhist->Write();
  
  auto sigbk_cut_canvas = new TCanvas("sigbk_cut_canvas", "sigbk_cut_canvas", 600, 400);
  sigbk_cut_canvas->cd();
  mg2->GetHistogram()->SetTitle(Form("Cut at %0.1f " + units + " with efficiency %0.2f and background rejection %0.2f", xval[maxsigbk_ind], efficiency[maxsigbk_ind], bkRej[maxsigbk_ind]));
  mg2->GetXaxis()->SetTitle(axisTitle);
  mg2->Draw("AP");
  signalhist->Draw("HIST SAME e");
  backgroundhist->Draw("HIST SAME e");
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
  signalhist->Draw("HIST SAME");
  backgroundhist->Draw("HIST SAME");
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


//This analysis the histograms are already integrated so no need for that. They also start off 2D.
void optimiser::MetricHolder::MakeStandardNumShowerAnalysisGraphs(float totalsig, float totalbk){

  //  SignalHistPartner->Scale(1/TotalSigPOT);
  //  BackgroundHistPartner->Scale(1/TotalBKPOT);

  TH2D* EfficiencyHist     = (TH2D*) SignalHistPartner->Clone();
  TH2D* PurityHist         = (TH2D*) SignalHistPartner->Clone(); 
  TH2D* BackgroundRejHist  = (TH2D*) BackgroundHistPartner->Clone();
  TH2D* EffBkRejHist       = (TH2D*) SignalHistPartner->Clone();
  TH2D* EffPurHist         = (TH2D*) SignalHistPartner->Clone();

  EfficiencyHist->SetName("Efficiency");
  PurityHist->SetName("Purity");
  BackgroundRejHist->SetName("Background Rejection");
  EffBkRejHist->SetName("Efficiency*Background Rejection");
  EffPurHist->SetName("Efficiency*Purity");
  EfficiencyHist->SetTitle("Efficiency");
  PurityHist->SetTitle("Purity");
  BackgroundRejHist->SetTitle("Background Rejection");
  EffBkRejHist->SetTitle("Efficiency*Background Rejection");
  EffPurHist->SetTitle("Efficiency*Purity");


  int n_x = SignalHistPartner->GetNbinsX();
  int n_y = SignalHistPartner->GetNbinsY();

  //  double totalsig = SignalHistPartner->GetEntries()/(SignalHistPartner->GetNbinsX()*SignalHistPartner->GetNbinsY());
  // double totalbk  = BackgroundHistPartner->GetEntries()/(SignalHistPartner->GetNbinsX()*SignalHistPartner->GetNbinsY());

  double max_effbk = -999;
  int max_effbk_i  = -999;
  int max_effbk_j  = -999;

  double max_effpur = -999;
  int max_effpur_i  = -999;
  int max_effpur_j  = -999;


  for(int i=1; i<n_x+1; ++i) {
    for(int j=1; j<n_y+1; ++j){ 

      int bin = EfficiencyHist->GetBin(i,j);

      EfficiencyHist->SetBinContent(i,j,SignalHistPartner->GetBinContent(i,j)/totalsig);
      EfficiencyHist->SetBinError(bin,EfficiencyHist->GetBinError(bin)/totalsig);

      if(SignalHistPartner->GetBinContent(i,j)+ BackgroundHistPartner->GetBinContent(i,j) > 0){
	PurityHist->SetBinContent(i,j,SignalHistPartner->GetBinContent(i,j)/(SignalHistPartner->GetBinContent(i,j)+ BackgroundHistPartner->GetBinContent(i,j)));
	double PurityDemErr = TMath::Sqrt(SignalHistPartner->GetBinError(bin)*SignalHistPartner->GetBinError(bin) + BackgroundHistPartner->GetBinError(bin)*BackgroundHistPartner->GetBinError(bin));
	double PurityError = TMath::Sqrt((SignalHistPartner->GetBinError(bin)/SignalHistPartner->GetBinContent(i,j))*(SignalHistPartner->GetBinError(bin)/SignalHistPartner->GetBinContent(i,j)) + (PurityDemErr/(SignalHistPartner->GetBinContent(i,j)+ BackgroundHistPartner->GetBinContent(i,j)))*(PurityDemErr/(SignalHistPartner->GetBinContent(i,j)+ BackgroundHistPartner->GetBinContent(i,j))));
	PurityHist->SetBinError(bin,PurityError); 
      }
      else{
	PurityHist->SetBinContent(i,j,0);
	PurityHist->SetBinError(bin,0);
      }

      BackgroundRejHist->SetBinContent(i,j,1-(BackgroundHistPartner->GetBinContent(i,j)/totalbk));
      BackgroundRejHist->SetBinError(bin,BackgroundRejHist->GetBinError(bin)/totalbk);

      EffBkRejHist->SetBinContent(i,j,EfficiencyHist->GetBinContent(i,j)*BackgroundRejHist->GetBinContent(i,j));
      double EffBkRejErr = TMath::Sqrt((EfficiencyHist->GetBinError(bin)/EfficiencyHist->GetBinContent(i,j))*(EfficiencyHist->GetBinError(bin)/EfficiencyHist->GetBinContent(i,j)) + (BackgroundRejHist->GetBinError(bin)/BackgroundRejHist->GetBinContent(i,j))*(BackgroundRejHist->GetBinError(bin)/BackgroundRejHist->GetBinContent(i,j)));
      EffBkRejHist->SetBinError(bin,EffBkRejErr);

      EffPurHist->SetBinContent(i,j,EfficiencyHist->GetBinContent(i,j)*PurityHist->GetBinContent(i,j));
      double EffPurErr = TMath::Sqrt((EfficiencyHist->GetBinError(bin)/EfficiencyHist->GetBinContent(i,j))*(EfficiencyHist->GetBinError(bin)/EfficiencyHist->GetBinContent(i,j)) + (PurityHist->GetBinError(bin)/PurityHist->GetBinContent(i,j))*(PurityHist->GetBinError(bin)/PurityHist->GetBinContent(i,j)));
      EffPurHist->SetBinError(bin,EffPurErr);

      //Store the max values 
      if(max_effbk < EffBkRejHist->GetBinContent(i,j)){
	max_effbk = EffBkRejHist->GetBinContent(i,j);
	max_effbk_i = i;
	max_effbk_j = j;
      } 
      if(max_effpur < EffPurHist->GetBinContent(i,j)){
	max_effpur = EffPurHist->GetBinContent(i,j);
	max_effpur_i = i;
	max_effpur_j = j;
      } 
    }
  }
  
  auto efficency_canvas = new TCanvas("eff2D_canvas", "eff2D_canvas", 600, 400);
  EfficiencyHist->Draw("colz");
  efficency_canvas->Write();

  auto background_canvas = new TCanvas("bk2D_canvas", "bk2D_canvas", 600, 400);
  BackgroundRejHist->Draw("colz");
  background_canvas->Write();

  auto purity_canvas = new TCanvas("pur2D_canvas", "pur2D_canvas", 600, 400);
  PurityHist->Draw("colz");
  purity_canvas->Write();

  auto effbk_canvas = new TCanvas("effbk2D_canvas", "effbk2D_canvas", 600, 400);
  EffBkRejHist->Draw("colz");
  effbk_canvas->Write();

  auto effpur_canvas = new TCanvas("effpur2D_canvas", "effpur2D_canvas", 600, 400);
  EffPurHist->Draw("colz");
  effpur_canvas->Write();

  //Make the 1D Graphs
  TH1D * Efficiency1D_sigbk    =  EfficiencyHist->ProjectionX("_pxsig",max_effbk_j,max_effbk_j,"e");
  Efficiency1D_sigbk->SetMarkerColor(kBlue);
  Efficiency1D_sigbk->SetLineColor(kBlue);
  Efficiency1D_sigbk->GetYaxis()->SetRangeUser(0, 1);
  
  TH1D * BackgroundRej1D_sigbk =  BackgroundRejHist->ProjectionX("_pxbk",max_effbk_j,max_effbk_j,"e");
  BackgroundRej1D_sigbk->SetMarkerColor(kRed);
  BackgroundRej1D_sigbk->SetLineColor(kRed);

  TH1D * EffBkRejHist1D_sigbk  =  EffBkRejHist->ProjectionX("_pxeffbk",max_effbk_j,max_effbk_j,"e");
  EffBkRejHist1D_sigbk->SetMarkerColor(kBlack);
  EffBkRejHist1D_sigbk->SetLineColor(kBlack);
  
  auto effbk1D_canvas = new TCanvas("effbk1D_canvas", "effbk1D_canvas", 600, 400);
  Efficiency1D_sigbk->Draw("ap");
  BackgroundRej1D_sigbk->Draw("hist same ap");
  EffBkRejHist1D_sigbk->Draw("hist same ap");
  effbk1D_canvas->BuildLegend();
  TLine* l1 = new TLine(((TAxis*)SignalHistPartner->GetYaxis())->GetBinCenter(max_effbk_j), 0, ((TAxis*)SignalHistPartner->GetYaxis())->GetBinCenter(max_effbk_j), 1);
  l1->SetLineColor(kBlack);
  l1->Draw();
  effbk1D_canvas->Write();
 

  std::cout << "Best Cut is at Energy: " << ((TAxis*)SignalHistPartner->GetXaxis())->GetBinCenter(max_effpur_i) << " Num Showers: " << ((TAxis*)SignalHistPartner->GetYaxis())->GetBinCenter(max_effpur_j) << std::endl;
  std::cout << "Efficiency: " << EfficiencyHist->GetBinContent(max_effpur_i,max_effpur_j) << " - " << EfficiencyHist->GetBinErrorLow(max_effpur_i,max_effpur_j) << " + " << EfficiencyHist->GetBinErrorUp(max_effpur_i,max_effpur_j)  << std::endl;
  std::cout << "Purity: " << PurityHist->GetBinContent(max_effpur_i,max_effpur_j) << " - " << PurityHist->GetBinErrorLow(max_effpur_i,max_effpur_j) << " + " << PurityHist->GetBinErrorUp(max_effpur_i,max_effpur_j)  << std::endl;
  std::cout << "Eff*Purity: " << EffPurHist->GetBinContent(max_effpur_i,max_effpur_j) << " - " << EffPurHist->GetBinErrorLow(max_effpur_i,max_effpur_j) << " + " << EffPurHist->GetBinErrorUp(max_effpur_i,max_effpur_j)  << std::endl;

  std::cout << "Best Cut Background Rej is at Energy: " << ((TAxis*)SignalHistPartner->GetXaxis())->GetBinCenter(max_effbk_i) << " Num Showers: " << ((TAxis*)SignalHistPartner->GetYaxis())->GetBinCenter(max_effbk_j) << std::endl;
  std::cout << "Efficiency: " << EfficiencyHist->GetBinContent(max_effbk_i,max_effbk_j) << " - " << EfficiencyHist->GetBinErrorLow(max_effbk_i,max_effbk_j) << " + " << EfficiencyHist->GetBinErrorUp(max_effbk_i,max_effbk_j)  << std::endl;
  std::cout << "Background Rejection: " << BackgroundRejHist->GetBinContent(max_effbk_i,max_effbk_j) << " - " << BackgroundRejHist->GetBinErrorLow(max_effbk_i,max_effbk_j) << " + " << BackgroundRejHist->GetBinErrorUp(max_effbk_i,max_effbk_j)  << std::endl;
  std::cout << "Eff*Purity: " << EffBkRejHist->GetBinContent(max_effbk_i,max_effbk_j) << " - " << EffBkRejHist->GetBinErrorLow(max_effbk_i,max_effbk_j) << " + " << EffBkRejHist->GetBinErrorUp(max_effbk_i,max_effbk_j)  << std::endl;

  delete efficency_canvas;
  delete background_canvas;
  delete purity_canvas;
  delete effbk_canvas;
  delete effpur_canvas;

  return;
}

void optimiser::MetricHolder::Plot1D(){
  Plot1D(SignalHist,BackgroundHist);
}


void optimiser::MetricHolder::Plot1D(TH1D* signalhist, TH1D* backgroundhist){

  double maxScale = TMath::Max(backgroundhist->GetBinContent(backgroundhist->GetMaximumBin()),
			       signalhist->GetBinContent(signalhist->GetMaximumBin()));

  signalhist->SetLineColor(kBlue);
  signalhist->SetFillStyle(3003);
  signalhist->SetFillColor(6);
  signalhist->Scale(1. / maxScale);
  signalhist->SetTitle("Signal");
  signalhist->Write();
  backgroundhist->SetLineColor(kRed);
  backgroundhist->SetFillStyle(3003);
  backgroundhist->SetFillColor(42);
  backgroundhist->Scale(1. / maxScale);
  backgroundhist->SetTitle("Background");
  backgroundhist->Write();
  
  auto single_canvas = new TCanvas("single_canvas", "single_canvas", 600, 400);
  single_canvas->cd();
  signalhist->Draw();
  backgroundhist->Draw("SAME");
  single_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
  single_canvas->Write();
}

void optimiser::MetricHolder::Plot2D(){
  Plot2D(SignalHistPartner,BackgroundHistPartner);
}

void optimiser::MetricHolder::Plot2D(TH2D* signalhist, TH2D* backgroundhist){

  auto signal_canvas = new TCanvas("signal_canvas", "signal_canvas", 600, 400);
  signalhist->Draw("colz");
  signal_canvas->Write();

  auto background_canvas = new TCanvas("background_canvas", "background_canvas", 600, 400);
  backgroundhist->Draw("colz");
  background_canvas->Write();
  

}

void optimiser::MetricHolder::MakeEfficiencyPlots(TH1D* postcut, TH1D* precut, TH1D* extra){

  TString eff_name = "eff_name";
  TString extraname = "extra_name";

  std::cout << "PostCut Size: " << postcut->GetEntries() << " PreCut Size; " << precut->GetEntries() << " Extra: " << extra->GetEntries() << std::endl;

  float maxSignal = precut->GetBinContent(precut->GetMaximumBin());
  
  if(TEfficiency::CheckConsistency(*postcut,*precut)){
    //Make the eff plot.
    TEfficiency* reco_eff = new TEfficiency(*postcut,*precut);
    //    reco_eff->SetTitle(title);
    //reco_eff->SetName(eff_name);
    reco_eff->SetTitle("title"); 
    reco_eff->SetName("eff_name"); 
    reco_eff->Write();

    postcut->Scale(1./maxSignal);
    precut->Scale(1./maxSignal);
      
    postcut->SetLineColor(kBlue);
    postcut->SetFillStyle(3003);
    postcut->SetFillColor(6);
    //    postcut->SetTitle(postcutname);
    postcut->SetTitle("postcutname"); 
    postcut->Write();
    precut->SetLineColor(kRed);
    precut->SetFillStyle(3003);
    precut->SetFillColor(42);
    precut->SetTitle("precutname");
    //    precut->SetTitle(precutname);
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

    std::cout << "Total Efficiency: " << (float) postcut->GetEntries() / (float) precut->GetEntries() << std::endl; 
    
  }

  return;

  if(TEfficiency::CheckConsistency(*extra,*precut)){
    //Make the eff plot.
    TEfficiency* reco_eff = new TEfficiency(*extra,*precut);
    //    reco_eff->SetTitle(title);
    //reco_eff->SetName(eff_name);
    reco_eff->SetTitle("reco eff");
    reco_eff->SetName("reco eff");
    
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
    //    precut->SetTitle(precutname);
    precut->SetTitle("precutname");
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


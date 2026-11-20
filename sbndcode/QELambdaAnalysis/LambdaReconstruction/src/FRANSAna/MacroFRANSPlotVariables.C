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


const int fSignalLS = kSolid;
const int fBackgroundLS = kDashed;
const int fSignalLC = kAzure+1;
const int fBackgroundLC = kOrange+1;

const std::vector<int> fColors = {kAzure-3, kOrange+7, kGreen+3, kBlack, kMagenta, kBlue, kRed, kOrange, kRed+2, kGreen+2, kMagenta+2, kBlack+2, kBlue+2, kOrange-3, kRed-3, kGreen-3, kMagenta-3, kBlack-3, kBlue-3, kOrange-6, kRed-6, kGreen-6, kMagenta-6, kBlack-6, kBlue-6, kOrange-9, kRed-9, kGreen-9, kMagenta-9, kBlack-9, kBlue-9};

const std::vector<int> fLines = {1, 2, 3, 4, 5, 6, 7, 8};

// function to draw ROC curve for one variable
TGraph GetROCCurve(TTree *tree, std::string varname, bool side, int nBins, double min=0, double max=1){
    TH1F* hSig = new TH1F(("hSig"+varname).c_str(), "Signal", nBins, min, max);
    TH1F* hBG = new TH1F(("hBG"+varname).c_str(), "Background", nBins, min, max);
    tree->Draw( (varname+">>hSig"+varname).c_str(), "IsSignal==1");
    tree->Draw( (varname+">>hBG"+varname).c_str(), "IsSignal!=1");
    
    // Normalize
    hSig->Scale(1./hSig->Integral());
    hBG->Scale(1./hBG->Integral());
    
    // Set colors
    hSig->SetLineColor(fSignalLC);
    hSig->SetLineStyle(fSignalLS);
    hBG->SetLineColor(fBackgroundLC);
    hBG->SetLineStyle(fBackgroundLS);

    double sigIntegral = hSig->Integral();
    double bgIntegral = hBG->Integral();

    // Calculate ROC curve
    TGraph ROC(nBins);
    // first point
    //ROC.SetPoint(0, 0, 1);
    for(int i=0; i<nBins; i++){
      double sig = hSig->Integral(0, i)/sigIntegral;
      double bg = hBG->Integral(0, i)/bgIntegral;
      if(!side) ROC.SetPoint(i, 1-sig, bg);
      else ROC.SetPoint(i, sig, 1-bg);
      
    }
    // add last point
    if(!side) ROC.SetPoint(nBins, 0, 1);
    else ROC.SetPoint(nBins, 1, 0);

    double integral = 0;
    double preX, preY;
    ROC.GetPoint(0, preX, preY);
    for(int i=1; i<ROC.GetN(); i++){
      double x, y;
      ROC.GetPoint(i, x, y);
      integral += (y+preY)/2 * (x-preX);
      preX = x;
      preY = y;
    }

    std::cout<<" Method "<<varname<<" Integral: "<<integral<<std::endl;

    return ROC;

}

int MacroFRANSPlotVariables(std::string fInputFileName="", int view = 2, std::string fTreeDirName = "FRANSCheatedVx/", std::string fTreeName = "FRANSTree")
{

  //--------- Configuration Parameters
  bool fUseMultiplicityVars = false; //fUseMultiplicityVars = true;
  std::string fBGLabel = "";  //fBGLabel = "QE"; //fBGLabel = "RES";
  std::string fVw = std::to_string(view);

  //--------- Input file
  TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
  fFile->ls();
  TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());

  SetFRANSStyle();
  gStyle->SetHistLineWidth(2);



  // Make ROC curve for each variable
  TCanvas* cROC = new TCanvas("cROC", "ROC", 800, 800);
  cROC->cd();
  // Left/bottom margins
  TPad *padROC = new TPad("padROC", "padROC", 0, 0, 1, 1);
  padROC->Draw();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  int fNBinsROC = 500;
  TH2F* hROCFrame = new TH2F("hROCFrame", ";Signal efficiency;Background rejection", 100, 0, 1, 100, 0, 1);
  hROCFrame->SetStats(0);
  TGraph ROCEta = GetROCCurve(fTree, "FRANSObj"+fVw+".fEta", 0, fNBinsROC, 0, 100);
  TGraph ROCDelta = GetROCCurve(fTree, "FRANSObj"+fVw+".fDelta", 0, fNBinsROC, 0, 1);
  TGraph ROCFitScore = GetROCCurve(fTree, "FRANSObj"+fVw+".fFitScore", 1, fNBinsROC, 0, 1);
  TGraph ROCIota = GetROCCurve(fTree, "FRANSObj"+fVw+".fIota", 1, fNBinsROC, 0, 10);

  // Colors and lines
  ROCEta.SetLineColor( fColors.at(1) );
  ROCEta.SetLineStyle( fLines.at(1) );
  ROCDelta.SetLineColor( fColors.at(2) );
  ROCDelta.SetLineStyle( fLines.at(2) );
  ROCFitScore.SetLineColor( fColors.at(3) );
  ROCFitScore.SetLineStyle( fLines.at(3) );
  ROCIota.SetLineColor( fColors.at(4) );
  ROCIota.SetLineStyle( fLines.at(5) );
  
  // Legend
  TLegend* ROClegend = new TLegend(0.2, 0.2, 0.5, 0.4);
  ROClegend->AddEntry(&ROCEta, "#eta", "l");
  ROClegend->AddEntry(&ROCDelta, "#Delta", "l");
  ROClegend->AddEntry(&ROCFitScore, "Fit score", "l");
  ROClegend->AddEntry(&ROCIota, "#iota", "l");

  // Draw
  hROCFrame->GetYaxis()->SetTitleOffset(1.3);
  hROCFrame->Draw();
  ROCEta.Draw("L same");
  ROCDelta.Draw("L same");
  ROCFitScore.Draw("L same");
  ROCIota.Draw("L same");
  
  // /Users/franciscojaviernicolas/Desktop/ThesisExamplesHypana/FRANSThesis/VarsROC/TMVAResultsBDT.root open this file and get ROC curve
  TFile* file1 = TFile::Open("/Users/franciscojaviernicolas/Desktop/ThesisExamplesHypana/FRANSThesis/VarsROC/TMVAResultsBDT.root");
  std::string fPlotName = "dataset/Method_BDT/BDT/MVA_BDT_rejBvsS";
  TH1F* roc = (TH1F*)file1->Get(fPlotName.c_str());
  roc->SetLineColor( fColors.at(0) );
  roc->SetLineStyle( fLines.at(0) );
  roc->Draw("L same");
  std::cout<<" Method BDT Integral: "<<roc->Integral()<<std::endl;
  ROClegend->AddEntry(roc, "BDT", "l");
  ROClegend->Draw("same");  

  cROC->cd();
  cROC->WaitPrimitive();
  cROC->Update();
  cROC->SaveAs("HistogramsFRANSPlotVariables/ROC.pdf");

  // Draw variables for signal and background
  TCanvas* c = new TCanvas("c", "Variables", 1000, 700);
  std::vector<TPad*> Tp = buildpadcanvas(3,2);
  Tp.at(1)->cd();

  TH1F* hAlphaFrame = new TH1F("hAlphaFrame", ";#alpha;AU", 50, 0, 10000);
  hAlphaFrame->SetStats(0);
  TH1F* hAlphaSig = new TH1F("hAlphaSig", "Signal;#alpha", hAlphaFrame->GetNbinsX(), hAlphaFrame->GetXaxis()->GetXmin(), hAlphaFrame->GetXaxis()->GetXmax());
  TH1F* hAlphaBG = new TH1F("hAlphaBG", "Background;#alpha", hAlphaFrame->GetNbinsX(), hAlphaFrame->GetXaxis()->GetXmin(), hAlphaFrame->GetXaxis()->GetXmax());
  fTree->Draw(("FRANSObj"+fVw+".fAlpha"+">>hAlphaSig").c_str(), "IsSignal==1");
  fTree->Draw(("FRANSObj"+fVw+".fAlpha"+">>hAlphaBG").c_str(), "IsSignal!=1");
  // Normalize
  hAlphaSig->Scale(1./hAlphaSig->Integral());
  hAlphaBG->Scale(1./hAlphaBG->Integral());
  // Set maximum
  hAlphaFrame->SetMaximum(std::max(hAlphaSig->GetMaximum(), hAlphaBG->GetMaximum())*1.2);
  // Set colors and line style
  hAlphaSig->SetLineColor(fSignalLC);
  hAlphaSig->SetLineStyle(fSignalLS);
  hAlphaBG->SetLineColor(fBackgroundLC);
  hAlphaBG->SetLineStyle(fBackgroundLS);

  hAlphaFrame->Draw();
  hAlphaSig->Draw("hist same");
  hAlphaBG->Draw("hist same");
  // Legend
  TLegend* legend = new TLegend(0.55, 0.7, 0.89, 0.89);
  legend->AddEntry(hAlphaSig, "Signal", "l");
  legend->AddEntry(hAlphaBG, "Background", "l");
  legend->SetBorderSize(0);
  legend->SetTextSize(0.05);
  legend->Draw("same");

  Tp.at(2)->cd();
  TH1F* hEtaFrame = new TH1F("hEtaFrame", ";#eta;AU", 50, 0, 20);
  hEtaFrame->SetStats(0);
  TH1F* hEtaSig = new TH1F("hEtaSig", "Signal;#eta", hEtaFrame->GetNbinsX(), hEtaFrame->GetXaxis()->GetXmin(), hEtaFrame->GetXaxis()->GetXmax());
  TH1F* hEtaBG = new TH1F("hEtaBG", "Background;#eta", hEtaFrame->GetNbinsX(), hEtaFrame->GetXaxis()->GetXmin(), hEtaFrame->GetXaxis()->GetXmax());
  fTree->Draw(("FRANSObj"+fVw+".fEta"+">>hEtaSig").c_str(), "IsSignal==1");
  fTree->Draw(("FRANSObj"+fVw+".fEta"+">>hEtaBG").c_str(), "IsSignal!=1");
  // Normalize
  hEtaSig->Scale(1./hEtaSig->Integral());
  hEtaBG->Scale(1./hEtaBG->Integral());
  // Set maximum
  hEtaFrame->SetMaximum(std::max(hEtaSig->GetMaximum(), hEtaBG->GetMaximum())*1.2);
  // Set colors
  hEtaSig->SetLineColor(fSignalLC);
  hEtaSig->SetLineStyle(fSignalLS);
  hEtaBG->SetLineColor(fBackgroundLC);
  hEtaBG->SetLineStyle(fBackgroundLS);

  hEtaFrame->Draw();
  hEtaSig->Draw("hist same");
  hEtaBG->Draw("hist same");
  legend->Draw("same");


  Tp.at(3)->cd();
  TH1F* hDeltaFrame = new TH1F("hDeltaFrame", ";#Delta;AU", 50, 0, 1);
  hDeltaFrame->SetStats(0);
  TH1F* hDeltaSig = new TH1F("hDeltaSig", "Signal;#Delta", hDeltaFrame->GetNbinsX(), hDeltaFrame->GetXaxis()->GetXmin(), hDeltaFrame->GetXaxis()->GetXmax());
  TH1F* hDeltaBG = new TH1F("hDeltaBG", "Background;#Delta", hDeltaFrame->GetNbinsX(), hDeltaFrame->GetXaxis()->GetXmin(), hDeltaFrame->GetXaxis()->GetXmax());
  fTree->Draw(("FRANSObj"+fVw+".fDelta"+">>hDeltaSig").c_str(), "IsSignal==1");
  fTree->Draw(("FRANSObj"+fVw+".fDelta"+">>hDeltaBG").c_str(), "IsSignal!=1");
  // Normalize
  hDeltaSig->Scale(1./hDeltaSig->Integral());
  hDeltaBG->Scale(1./hDeltaBG->Integral());
  // Set maximum
  hDeltaFrame->SetMaximum(std::max(hDeltaSig->GetMaximum(), hDeltaBG->GetMaximum())*1.2);
  // Set colors
  hDeltaSig->SetLineColor(fSignalLC);
  hDeltaSig->SetLineStyle(fSignalLS);
  hDeltaBG->SetLineColor(fBackgroundLC);
  hDeltaBG->SetLineStyle(fBackgroundLS);

  hDeltaFrame->Draw();
  hDeltaSig->Draw("hist same");
  hDeltaBG->Draw("hist same");
  legend->Draw("same");


  Tp.at(4)->cd();
  TH1F* hFitScoreFrame = new TH1F("hFitScoreFrame", ";Fit score;AU", 50, 0, 1.000001);
  hFitScoreFrame->SetStats(0);
  TH1F* hFitScoreSig = new TH1F("hFitScoreSig", "Signal;Fit score", hFitScoreFrame->GetNbinsX(), hFitScoreFrame->GetXaxis()->GetXmin(), hFitScoreFrame->GetXaxis()->GetXmax());
  TH1F* hFitScoreBG = new TH1F("hFitScoreBG", "Background;Fit score", hFitScoreFrame->GetNbinsX(), hFitScoreFrame->GetXaxis()->GetXmin(), hFitScoreFrame->GetXaxis()->GetXmax());
  fTree->Draw(("FRANSObj"+fVw+".fFitScore"+">>hFitScoreSig").c_str(), "IsSignal==1");
  fTree->Draw(("FRANSObj"+fVw+".fFitScore"+">>hFitScoreBG").c_str(), "IsSignal!=1");
  // Normalize
  hFitScoreSig->Scale(1./hFitScoreSig->Integral());
  hFitScoreBG->Scale(1./hFitScoreBG->Integral());
  // Set maximum
  hFitScoreFrame->SetMaximum(std::max(hFitScoreSig->GetMaximum(), hFitScoreBG->GetMaximum())*1.2);
  // Set colors
  hFitScoreSig->SetLineColor(fSignalLC);
  hFitScoreSig->SetLineStyle(fSignalLS);
  hFitScoreBG->SetLineColor(fBackgroundLC);
  hFitScoreBG->SetLineStyle(fBackgroundLS);

  hFitScoreFrame->Draw();
  hFitScoreSig->Draw("hist same");
  hFitScoreBG->Draw("hist same");
  legend->Draw("same");


  Tp.at(5)->cd();
  TH1F* hIotaFrame = new TH1F("hIotaFrame", ";#iota;AU", 20, 0, 10);
  hIotaFrame->SetStats(0);
  TH1F* hIotaSig = new TH1F("hIotaSig", "Signal;#iota", hIotaFrame->GetNbinsX(), hIotaFrame->GetXaxis()->GetXmin(), hIotaFrame->GetXaxis()->GetXmax());
  TH1F* hIotaBG = new TH1F("hIotaBG", "Background;#iota", hIotaFrame->GetNbinsX(), hIotaFrame->GetXaxis()->GetXmin(), hIotaFrame->GetXaxis()->GetXmax());
  fTree->Draw(("FRANSObj"+fVw+".fIota"+">>hIotaSig").c_str(), "IsSignal==1");
  fTree->Draw(("FRANSObj"+fVw+".fIota"+">>hIotaBG").c_str(), "IsSignal!=1");
  // Normalize
  hIotaSig->Scale(1./hIotaSig->Integral());
  hIotaBG->Scale(1./hIotaBG->Integral());
  // Set maximum
  hIotaFrame->SetMaximum(std::max(hIotaSig->GetMaximum(), hIotaBG->GetMaximum())*1.2);
  // Set colors
  hIotaSig->SetLineColor(fSignalLC);
  hIotaSig->SetLineStyle(fSignalLS);
  hIotaBG->SetLineColor(fBackgroundLC);
  hIotaBG->SetLineStyle(fBackgroundLS);

  hIotaFrame->Draw();
  hIotaSig->Draw("hist same");
  hIotaBG->Draw("hist same");
  legend->Draw("same");

  c->cd();
  c->Update();
  c->WaitPrimitive();
  // Remove and create output directory
  gSystem->Exec("rm -rf HistogramsFRANSPlotVariables");
  gSystem->Exec("mkdir HistogramsFRANSPlotVariables");
  // Save as pdf
  c->SaveAs(("HistogramsFRANSPlotVariables/FRANSPlotVariables_"+fVw+".pdf").c_str());
  // Save individual pads
  for(int i=1; i<Tp.size(); i++){
    Tp.at(i)->SaveAs(("HistogramsFRANSPlotVariables/FRANSPlotVariables_"+fVw+"_"+std::to_string(i)+".eps").c_str());
  }



  return 0;
}

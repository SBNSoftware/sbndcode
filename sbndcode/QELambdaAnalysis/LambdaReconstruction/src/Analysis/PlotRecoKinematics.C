#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>

#include "CutEfficienciesDefinitions.C"
#include "CutDefinitions.C"
#include "CutEfficienciesStyle.C"


//---------  Cuts
std::vector<PlotDef> fCutDefs = cutDefsPIDFull;

//--------- Signal and BG definitions
std::vector<SampleDef> sampleDefs = {
    {fTruthInAV+"IntOrigin==1 && IntDirt==0 && (IntNLambda>0 && IntMode==0 && abs(IntNuPDG)!=12)", "Signal", true, "Signal"}
    ,{fTruthInAV+"IntOrigin==1 &&  IntDirt==0 && !(IntNLambda>0 && IntMode==0 && abs(IntNuPDG)!=12)", "BG  #nu", false, "BG BNB"}
    ,{"IntOrigin==2 || IntDirt==1", "Dirt+Cosmic", false, "Dirt+Cosmic"}
};


//---------  Load function
void PlotRecoKinematics(){
    std::cout<<" Load LambdaLifetime.C"<<std::endl;
}


//--------- Function to plot the lifetime
void RunPlotLambdaLifetime(std::string fInputFileName="", std::string fTreeDirName = "originsAnaPost/", std::string fTreeName = "LambdaAnaTree")
{
    //---------  Settings
    double fMinTau = 0;
    double fMaxTau = 3;
    int fNBinsTau = 50;

    //---------  Remove all *.pdf with gSystem
    std::string fOutputDirName = "LambdaLifetimeOutputPlots";
    gSystem->Exec(("rm -rf "+fOutputDirName).c_str());
    gSystem->Exec(("mkdir "+fOutputDirName).c_str());

    //--------- Input TTree
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
    // Read TreeHeader 
    TTree *fTreeHeader = (TTree *)fFile->Get( (fTreeDirName+"TreeHeader").c_str() );
    
    //--------- Get the accumulated cut   
    TCut fAllCuts("");
    for (size_t i = 0; i < fCutDefs.size(); i++){
        if(fCutDefs[i].GetAccumulateCut()){
            fAllCuts+=TCut(fCutDefs[i].GetCut());
        }
    }
    std::cout<<"All cuts: "<<fAllCuts<<std::endl;
  
    //---------  Output canvas
    CutStyler *fStyler = new CutStyler(0);
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    TPad *pad1 = new TPad("pad1","pad1",0.01,0.01,0.49,0.49);
    TPad *pad2 = new TPad("pad2","pad2",0.51,0.01,0.99,0.49);
    TPad *pad3 = new TPad("pad3","pad3",0.01,0.51,0.49,0.99);
    TPad *pad4 = new TPad("pad4","pad4",0.51,0.51,0.99,0.99);
    pad1->Draw();
    pad2->Draw();
    pad3->Draw();
    pad4->Draw();
    // Set stats to 0
    gStyle->SetOptStat(0);

    //---------  Plot tau lab
    pad1->cd();
    TH1F *hTauLab = new TH1F("hTauLab","hTauLab",fNBinsTau,fMinTau,fMaxTau);
    fTree->Draw("10**9*RecoTauLab>>hTauLab",fAllCuts);
    hTauLab->SetLineColor(fStyler->GetColor(0));
    hTauLab->SetLineWidth(2);
    hTauLab->Draw();
    // Axis
    hTauLab->GetXaxis()->SetTitle("#tau^{#Lambda} [ns]");
    hTauLab->GetYaxis()->SetTitle("# events");


    //---------  Plot tau rest
    pad2->cd();
    TH1F *hTauRest = new TH1F("hTauRest","hTauRest",fNBinsTau,fMinTau,fMaxTau);
    fTree->Draw("10**9*RecoTau>>hTauRest",fAllCuts);
    hTauRest->SetLineColor(fStyler->GetColor(0));
    hTauRest->SetLineWidth(2);
    hTauRest->Draw();
    // Axis
    hTauRest->GetXaxis()->SetTitle("#tau_{0}^{#Lambda} [ns]");
    hTauRest->GetYaxis()->SetTitle("# events");
    hTauRest->SetLineColor(fStyler->GetColor(0));
    hTauRest->SetLineWidth(2);
    // Exponential fit
    // Define the function
    // Fit only from the bin with max number of entries
    double fMinTauFit = hTauRest->GetBinLowEdge(hTauRest->GetMaximumBin());
    TF1 *expoFitRest = new TF1("expoFitRest","[0]*exp(-x/[1])", fMinTauFit, fMaxTau);
    expoFitRest->SetParameter(0, 100);
    expoFitRest->SetParameter(1, 0.5);
    hTauRest->Fit("expoFitRest", "r");
    // Draw fit
    expoFitRest->SetLineColor(fStyler->GetColor(1));
    expoFitRest->SetLineWidth(2);
    expoFitRest->Draw("same");
    // Parameters
    double tauRest = expoFitRest->GetParameter(1);
    double tauRestErr = expoFitRest->GetParError(1);
    // Draw latex in the corner
    TLatex *latex = new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.08);
    latex->DrawLatex(0.4,0.5,Form("#tau_{0}^{#Lambda} = %.2f #pm %.2f ns",tauRest,tauRestErr));
    // Legend
    TLegend *leg1 = new TLegend(0.7,0.7,0.9,0.9);
    leg1->AddEntry(hTauRest,"Reco","l");
    leg1->AddEntry(expoFitRest,"Exponential fit","l");
    leg1->Draw();

    //---------  Plot reco gap
    pad3->cd();
    TH1F *hRecoGap = new TH1F("hRecoGap","hRecoGap",30,0,30);
    fTree->Draw("RecoGap>>hRecoGap",fAllCuts);
    TH1F *hTrueGap = new TH1F("hTrueGap","hTrueGap",30,0,30);
    fTree->Draw("Gap>>hTrueGap",fAllCuts);
    hRecoGap->SetLineColor(fStyler->GetColor(0));
    hRecoGap->SetLineWidth(2);
    hRecoGap->SetLineStyle(kDashed);
    hTrueGap->SetLineColor(fStyler->GetColor(1));
    hRecoGap->SetLineWidth(2);
    hRecoGap->Draw();
    hTrueGap->Draw("same");
    // Legend
    TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(hRecoGap,"Reco","l");
    leg->AddEntry(hTrueGap,"True","l");
    leg->Draw();
    // Axis
    hRecoGap->GetXaxis()->SetTitle("Gap [cm]");
    hRecoGap->GetYaxis()->SetTitle("# events");

    //---------  Plot reco gap resolution
    pad4->cd();
    TH1F *hRecoGapRes = new TH1F("hRecoGapRes","hRecoGapRes",30,-6,6);
    fTree->Draw("RecoGap-Gap>>hRecoGapRes",fAllCuts);
    hRecoGapRes->SetLineColor(fStyler->GetColor(0));
    hRecoGapRes->SetLineWidth(2);
    hRecoGapRes->Draw();
    // Axis
    hRecoGapRes->GetXaxis()->SetTitle("Reco Gap - True Gap [cm]");
    hRecoGapRes->GetYaxis()->SetTitle("# events");
    // Add stats latex
    latex->SetNDC();
    latex->DrawLatex(0.65,0.8,Form("Mean = %.2f",hRecoGapRes->GetMean()));
    latex->DrawLatex(0.65,0.7,Form("RMS = %.2f",hRecoGapRes->GetRMS()));

    
    

    c1->cd();
    c1->Update();
    c1->WaitPrimitive();
    c1->SaveAs((fOutputDirName+"/LambdaLifetime.pdf").c_str());
    c1->SaveAs((fOutputDirName+"/LambdaLifetime.png").c_str());
    c1->SaveAs((fOutputDirName+"/LambdaLifetime.root").c_str());



    return;
}


//--------- Function to plot the lifetime
void RunPlotArmenteros(std::string fInputFileName="", std::string fTreeDirName = "originsAnaPost/", std::string fTreeName = "LambdaAnaTree")
{
    //---------  Settings
    double fMinAlpha = -1;
    double fMaxAlpha = 1;
    int fNBinsAlpha = 100;
    double fMinPt = 0;
    double fMaxPt = 150;
    int fNBinsPt = 100;

    //---------  Remove all *.pdf with gSystem
    std::string fOutputDirName = "LambdaArmenterosOutputPlots";
    gSystem->Exec(("rm -rf "+fOutputDirName).c_str());
    gSystem->Exec(("mkdir "+fOutputDirName).c_str());

    //--------- Input TTree
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
    // Read TreeHeader 
    TTree *fTreeHeader = (TTree *)fFile->Get( (fTreeDirName+"TreeHeader").c_str() );
    
    //--------- Get the accumulated cut
    TCut fAllCuts("");
    for (size_t i = 0; i < fCutDefs.size(); i++){
        if(fCutDefs[i].GetAccumulateCut()){
            fAllCuts+=TCut(fCutDefs[i].GetCut());
        }
    }
    std::cout<<"All cuts: "<<fAllCuts<<std::endl;
  
    //---------  Output canvas
    CutStyler *fStyler = new CutStyler(0);
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    TPad *pad1 = new TPad("pad1","pad1",0.01,0.01,0.49,0.49);
    TPad *pad2 = new TPad("pad2","pad2",0.51,0.01,0.99,0.49);
    TPad *pad3 = new TPad("pad3","pad3",0.01,0.51,0.49,0.99);
    TPad *pad4 = new TPad("pad4","pad4",0.51,0.51,0.99,0.99);
    pad1->Draw();
    pad2->Draw();
    pad3->Draw();
    pad4->Draw();

    //---------  Plot the Armenteros-Podolanski plot (true level)
    pad1->cd();
    TH2F *hArmenterosTrue = new TH2F("hArmenterosTrue","hArmenterosTrue",fNBinsAlpha,fMinAlpha,fMaxAlpha,fNBinsPt,fMinPt,fMaxPt);
    hArmenterosTrue->SetStats(0);
    TH2F *hArmenterosTrueHighBeta = new TH2F("hArmenterosTrueHighBeta","hArmenterosTrueHighBeta",fNBinsAlpha,fMinAlpha,fMaxAlpha,fNBinsPt,fMinPt,fMaxPt);
    fTree->Draw("1000*LambdaDecayTransverseMomentum:LambdaDecayLongitudinalAsymmetry>>hArmenterosTrue", sampleDefs[0].GetVar());
    fTree->Draw("1000*LambdaDecayTransverseMomentum:LambdaDecayLongitudinalAsymmetry>>hArmenterosTrueHighBeta", sampleDefs[0].GetVar()+" && LambdaBeta>0.6");
    hArmenterosTrue->SetMarkerColor(fStyler->GetColor(0));
    hArmenterosTrue->SetMarkerStyle(20);
    hArmenterosTrue->SetMarkerSize(0.5);
    hArmenterosTrueHighBeta->SetMarkerColorAlpha(fStyler->GetColor(1), 0.5);
    hArmenterosTrueHighBeta->SetMarkerStyle(20);
    hArmenterosTrueHighBeta->SetMarkerSize(0.5);
    hArmenterosTrue->Draw();
    hArmenterosTrueHighBeta->Draw("same");
    // Axis
    hArmenterosTrue->GetXaxis()->SetTitle("#alpha");
    hArmenterosTrue->GetYaxis()->SetTitle("p_{T} [MeV/c]");
    // Legend
    TLegend *leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(hArmenterosTrue,"All","p");
    leg->AddEntry(hArmenterosTrueHighBeta,"#beta>0.6","p");
    leg->Draw();
    

    //---------  Plot the Armenteros-Podolanski plot (reco level)
    pad2->cd();
    TH2F *hArmenterosReco = new TH2F("hArmenterosReco","hArmenterosReco",fNBinsAlpha,fMinAlpha,fMaxAlpha,fNBinsPt,fMinPt,fMaxPt);
    hArmenterosReco->SetStats(0);
    TH2F *hArmenterosRecoHighBeta = new TH2F("hArmenterosRecoHighBeta","hArmenterosRecoHighBeta",fNBinsAlpha,fMinAlpha,fMaxAlpha,fNBinsPt,fMinPt,fMaxPt);
    fTree->Draw("RecoDecayTransverseMomentum:RecoDecayLongitudinalAsymmetry>>hArmenterosReco", fAllCuts);
    fTree->Draw("RecoDecayTransverseMomentum:RecoDecayLongitudinalAsymmetry>>hArmenterosRecoHighBeta",fAllCuts+TCut("RecoBeta>0.6") );
    hArmenterosReco->SetMarkerColor(fStyler->GetColor(0));
    hArmenterosReco->SetMarkerStyle(20);
    hArmenterosReco->SetMarkerSize(0.5);
    hArmenterosRecoHighBeta->SetMarkerColorAlpha(fStyler->GetColor(1), 0.5);
    hArmenterosRecoHighBeta->SetMarkerStyle(20);
    hArmenterosRecoHighBeta->SetMarkerSize(0.5);
    hArmenterosReco->Draw();
    hArmenterosRecoHighBeta->Draw("same");
    // Axis
    hArmenterosReco->GetXaxis()->SetTitle("#alpha");
    hArmenterosReco->GetYaxis()->SetTitle("p_{T} [MeV/c]");
    // Legend
    leg->Draw();


    //---------  Plot resolutions for alpha
    pad3->cd();
    TH1F *hAlphaRes = new TH1F("hAlphaRes","hAlphaRes",100,-1,1);
    fTree->Draw("RecoDecayLongitudinalAsymmetry-LambdaDecayLongitudinalAsymmetry>>hAlphaRes",fAllCuts);
    TH1F *hAlphaResFromGap = new TH1F("hAlphaResFromGap","hAlphaResFromGap",100,-1,1);
    fTree->Draw("RecoDecayLongitudinalAsymmetryFromGap-LambdaDecayLongitudinalAsymmetry>>hAlphaResFromGap",fAllCuts);
    hAlphaRes->SetLineColor(fStyler->GetColor(0));
    hAlphaRes->SetLineWidth(2);
    hAlphaRes->Draw();
    hAlphaResFromGap->SetLineColor(fStyler->GetColor(1));
    hAlphaResFromGap->SetLineWidth(2);
    hAlphaResFromGap->Draw("same");
    // Axis
    hAlphaRes->GetXaxis()->SetTitle("Reco #alpha - True #alpha");
    hAlphaRes->GetYaxis()->SetTitle("# events");
    
    //---------  Plot resolutions for pT
    pad4->cd();
    TH1F *hPtRes = new TH1F("hPtRes","hPtRes",100,-100,100);
    fTree->Draw("RecoDecayTransverseMomentum-1000*LambdaDecayTransverseMomentum>>hPtRes",fAllCuts);
    TH1F *hPtResFromGap = new TH1F("hPtResFromGap","hPtResFromGap",100,-100,100);
    fTree->Draw("RecoDecayTransverseMomentumFromGap-1000*LambdaDecayTransverseMomentum>>hPtResFromGap",fAllCuts);
    hPtRes->SetLineColor(fStyler->GetColor(0));
    hPtRes->SetLineWidth(2);
    hPtRes->Draw();
    hPtResFromGap->SetLineColor(fStyler->GetColor(1));
    hPtResFromGap->SetLineWidth(2);
    hPtResFromGap->Draw("same");
    // Axis
    hPtRes->GetXaxis()->SetTitle("Reco p_{T} - True p_{T} [MeV/c]");
    hPtRes->GetYaxis()->SetTitle("# events");


    

   
    c1->cd();
    c1->Update();
    c1->WaitPrimitive();
    c1->SaveAs((fOutputDirName+"/LambdaArmenteros.pdf").c_str());
    c1->SaveAs((fOutputDirName+"/LambdaArmenteros.png").c_str());
    c1->SaveAs((fOutputDirName+"/LambdaArmenteros.root").c_str());
    
    return;
}



//--------- Function to plot the lifetime
void RunPlotKinematicResolution(std::string fInputFileName="", std::string fTreeDirName = "originsAnaPost/", std::string fTreeName = "LambdaAnaTree")
{
    //---------  Settings
    double fMinKERes = -100;
    double fMaxKERes = 100;
    int fNBinsKERes = 100;

    //---------  Remove all *.pdf with gSystem
    std::string fOutputDirName = "LambdaKinematicResolutionOutputPlots";
    gSystem->Exec(("rm -rf "+fOutputDirName).c_str());
    gSystem->Exec(("mkdir "+fOutputDirName).c_str());

    //--------- Input TTree
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
    // Read TreeHeader 
    TTree *fTreeHeader = (TTree *)fFile->Get( (fTreeDirName+"TreeHeader").c_str() );
    
    //--------- Get the accumulated cut
    TCut fAllCuts("");
    for (size_t i = 0; i < fCutDefs.size(); i++){
        if(fCutDefs[i].GetAccumulateCut()){
            fAllCuts+=TCut(fCutDefs[i].GetCut());
        }
    }
    std::cout<<"All cuts: "<<fAllCuts<<std::endl;
  
    //---------  Output canvas
    CutStyler *fStyler = new CutStyler(0);
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    gStyle->SetOptStat(1100);
    // change position of stats box in relative coordinates
    gStyle->SetStatX(0.8);
    gStyle->SetStatY(0.8);
    gStyle->SetStatW(0.2);
    gStyle->SetStatH(0.2);
    gStyle->SetStatBorderSize(0);
    gStyle->SetStatFont(62);


    
    
    //---------  Plot the Proton KE resolution
    c1->cd();
    TH1F *hProtonKERes = new TH1F("hProtonKERes","Proton KE Resolution; KE_{reco} - KE_{MC} [MeV];# entries",fNBinsKERes,fMinKERes,fMaxKERes);
    fTree->Draw("KEHI-1000*ProtonKE>>hProtonKERes",fAllCuts);
    hProtonKERes->Draw();
    hProtonKERes->SetLineColor(fStyler->GetColor(0));
    hProtonKERes->SetLineWidth(2);
    c1->Update();
    c1->SaveAs((fOutputDirName+"/ProtonKERes.pdf").c_str());


    //---------  Plot the Pion KE resolution
    c1->cd();
    TH1F *hPionKERes = new TH1F("hPionKERes","Pion KE Resolution; KE_{reco} - KE_{MC} [MeV];# entries",fNBinsKERes,fMinKERes,fMaxKERes);
    fTree->Draw("KELI-1000*PionKE>>hPionKERes",fAllCuts);
    hPionKERes->Draw();
    hPionKERes->SetLineColor(fStyler->GetColor(1));
    hPionKERes->SetLineWidth(2);
    c1->Update();
    c1->SaveAs((fOutputDirName+"/PionKERes.pdf").c_str());


    //---------  Plot the Lambda KE resolution
    c1->cd();
    TH1F *hLambdaKERes = new TH1F("hLambdaKERes","Lambda KE Resolution; KE_{reco} - KE_{MC} [MeV];# entries",fNBinsKERes,fMinKERes,fMaxKERes);
    fTree->Draw("KEDecayedMother-1000*LambdaKE>>hLambdaKERes",fAllCuts);
    hLambdaKERes->Draw();
    hLambdaKERes->SetLineColor(fStyler->GetColor(2));
    hLambdaKERes->SetLineWidth(2);
    c1->Update();
    c1->SaveAs((fOutputDirName+"/LambdaKERes.pdf").c_str());


    //--------- Plot the Lambda Beta resolution
    c1->cd();
    TH1F *hLambdaBetaRes = new TH1F("hLambdaBetaRes","Lambda Beta Resolution; #beta_{reco} - #beta_{MC};# entries",100,-1,1);
    fTree->Draw("RecoBeta-LambdaBeta>>hLambdaBetaRes",fAllCuts);
    hLambdaBetaRes->Draw();
    hLambdaBetaRes->SetLineColor(fStyler->GetColor(0));
    hLambdaBetaRes->SetLineWidth(2);
    c1->Update();
    c1->SaveAs((fOutputDirName+"/LambdaBetaRes.pdf").c_str());


    //--------- Plot the Lambda gamma resolution
    c1->cd();
    TH1F *hLambdaGammaRes = new TH1F("hLambdaGammaRes","Lambda Gamma Resolution; #gamma_{reco} - #gamma_{MC};# entries",100,-1,1);
    fTree->Draw("RecoGamma-LambdaGamma>>hLambdaGammaRes",fAllCuts);
    hLambdaGammaRes->Draw();
    hLambdaGammaRes->SetLineColor(fStyler->GetColor(0));
    hLambdaGammaRes->SetLineWidth(2);
    c1->Update();
    c1->SaveAs((fOutputDirName+"/LambdaGammaRes.pdf").c_str());


    //--------- Plot the Lambda invariant mass
    c1->cd();
    TH1F *hLambdaInvMass = new TH1F("hLambdaInvMass","Lambda Invariant Mass; M_{#Lambda} [GeV];# entries",40,1.,1.5);
    fTree->Draw("0.001*InvariantMass>>hLambdaInvMass",fAllCuts);
    hLambdaInvMass->Draw();
    hLambdaInvMass->SetLineColor(fStyler->GetColor(0));
    hLambdaInvMass->SetLineWidth(2);
    // Draw marker for the true mass value
    TMarker *marker = new TMarker(1.115, 0, 33);
    marker->SetMarkerColor(fStyler->GetColor(1));
    marker->SetMarkerSize(3);
    marker->Draw();
    c1->Update();
    c1->SaveAs((fOutputDirName+"/LambdaInvMass.pdf").c_str());
    

    //--------- Plot the RecoDecayGapAngleDifference
    c1->cd();
    TH1F *hRecoDecayGapAngleDifference = new TH1F("hRecoDecayGapAngleDifference",";#phi;# entries",20,-20,100);
    fTree->Draw("180*RecoDecayGapAngleDifference/TMath::Pi()>>hRecoDecayGapAngleDifference",fAllCuts);
    hRecoDecayGapAngleDifference->Draw();
    c1->Update();
    c1->SaveAs((fOutputDirName+"/RecoDecayGapAngleDifference.pdf").c_str());



    return;
}
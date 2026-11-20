#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <string>
#include <vector>

#include "TMVA/Reader.h"

#include "../FRANS/FRANSObj.cpp"


void MacroCompareViews(std::string path="./",
                       const std::string& fInputFileName = "",
                       const std::string fTreeDirName = "FRANSCheatedVx/",
                       const std::string fTreeName = "FRANSTree")
{
    
    bool fTrainBDT = false; //fTrainBDT = true;
    bool fUseIota = false; fUseIota = true;
    bool fUseAlpha = false; //fUseAlpha = true;
    // map with label and file name
    std::map<std::string, std::string> fTMVAWeightsMap;
    // also for ROOT files
    std::map<std::string, std::string> fTMVARootFilesMap;

    // vector with view labels
    std::vector<std::string> viewLabels = {"0", "1", "2"};

    // initialize for the three views
    for (auto it = viewLabels.begin(); it != viewLabels.end(); ++it) {
        fTMVAWeightsMap[*it] = "";
        fTMVARootFilesMap[*it] = "";
    }

    // read the files in the directory
    TSystemDirectory dir(path.c_str(), path.c_str());
    TList* files = dir.GetListOfFiles();
    if (files) {
        TSystemFile* file;
        TString fname;
        TIter next(files);
        while ((file=(TSystemFile*)next())) {
            
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(".xml") && fname.BeginsWith("FRANSSelectionTMVA")) { 
               // check view label is in the name
                for (auto it = viewLabels.begin(); it != viewLabels.end(); ++it) {
                    std::string viewLabel = "_View"+*it;
                     if (fname.Contains(it->c_str())) {
                          fTMVAWeightsMap[*it] = path + fname.Data();
                     }
                }
            }

            if (!file->IsDirectory() && fname.EndsWith(".root") && fname.BeginsWith("TMVAResults")) { 
                // check view label is in the name
                for (auto it = viewLabels.begin(); it != viewLabels.end(); ++it) {
                    std::string viewLabel = "_View"+*it;
                    if (fname.Contains(it->c_str())) {
                        fTMVARootFilesMap[*it] = path + fname.Data();
                    }
                }
            }
        }
    }

    // cout files
    std::cout<<"Weight files to use:"<<std::endl;
    for (auto it = fTMVAWeightsMap.begin(); it != fTMVAWeightsMap.end(); ++it) {
        std::cout<<it->first<<" "<<it->second<<std::endl;
    }
    std::cout<<"ROOT files to use:"<<std::endl;
    for (auto it = fTMVARootFilesMap.begin(); it != fTMVARootFilesMap.end(); ++it) {
        std::cout<<it->first<<" "<<it->second<<std::endl;
    }
    
    // Load TMVA readers for each view
    std::map<std::string, TMVA::Reader*> fTMVAReaders;
    float alpha, eta, delta, fitScore, iota;
    float _gap, _protonKE, _pionKE;
    for (auto it = fTMVAWeightsMap.begin(); it != fTMVAWeightsMap.end(); ++it) {
        std::cout<<"Loading TMVA for "<<it->first<<" "<<it->second<<std::endl;
        TMVA::Reader* reader = new TMVA::Reader();

        if(fUseAlpha)
            reader->AddVariable("FRANSObj"+it->first+".fAlpha", &alpha);
        if(fUseIota)
            reader->AddVariable("FRANSObj"+it->first+".fIota", &iota);
        reader->AddVariable("FRANSObj"+it->first+".fEta", &eta);
        reader->AddVariable("FRANSObj"+it->first+".fDelta", &delta);
        reader->AddVariable("FRANSObj"+it->first+".fFitScore", &fitScore);

        reader->AddSpectator("Gap", &_gap);
        reader->AddSpectator("ProtonKE", &_protonKE);
        reader->AddSpectator("PionKE", &_pionKE);

        
        reader->BookMVA("FRANS BDT", it->second.c_str());
        fTMVAReaders[it->first] = reader;
        std::cout<<"  Loaded!"<<std::endl;
    }

    //--------- Input tree
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
    // Set Branch addresses for FRANSObj
    FRANSObj *fFRANSObj0 = new FRANSObj();
    FRANSObj *fFRANSObj1 = new FRANSObj();
    FRANSObj *fFRANSObj2 = new FRANSObj();
    fTree->SetBranchAddress("FRANSObj0", &fFRANSObj0);
    fTree->SetBranchAddress("FRANSObj1", &fFRANSObj1);
    fTree->SetBranchAddress("FRANSObj2", &fFRANSObj2);
    // Set Branch addresses for MC variables
    double gap, protonKE, pionKE;
    fTree->SetBranchAddress("Gap", &gap);
    fTree->SetBranchAddress("ProtonKE", &protonKE);
    fTree->SetBranchAddress("PionKE", &pionKE);
    // Set Branch addresses for true variables
    int isSignal;
    fTree->SetBranchAddress("IsSignal", &isSignal);

    //--------- Output file with combined TTree
    TFile* outputFile = TFile::Open( "CombinedViews.root", "RECREATE" );
    TTree* outputTree = new TTree("FRANSCombinedViews", "FRANS Combined Views Tree");
    // Set Branch addresses for scores
    double score0, score1, score2;
    outputTree->Branch("score0", &score0, "score0/D");
    outputTree->Branch("score1", &score1, "score1/D");
    outputTree->Branch("score2", &score2, "score2/D");
    // Set Branch addresses for true variables
    outputTree->Branch("isSignal", &isSignal, "isSignal/I");
    // Set Branch addresses for MC variables
    outputTree->Branch("gap", &gap, "gap/D");
    outputTree->Branch("protonKE", &protonKE, "protonKE/D");
    outputTree->Branch("pionKE", &pionKE, "pionKE/D");


    //--------- Loop over the tree
    for (Long64_t ievt = 0; ievt < fTree->GetEntries(); ievt++) {
        fTree->GetEntry(ievt);

        // Set variables for first view
        alpha = fFRANSObj0->GetAlpha();
        eta = fFRANSObj0->GetEta();
        delta = fFRANSObj0->GetDelta();
        fitScore = fFRANSObj0->GetFitScore();
        iota = fFRANSObj0->GetIota();
        // get the score
        score0 = fTMVAReaders["0"]->EvaluateMVA("FRANS BDT");

        // Set variables for second view
        alpha = fFRANSObj1->GetAlpha();
        eta = fFRANSObj1->GetEta();
        delta = fFRANSObj1->GetDelta();
        fitScore = fFRANSObj1->GetFitScore();
        iota = fFRANSObj1->GetIota();
        // get the score
        score1 = fTMVAReaders["1"]->EvaluateMVA("FRANS BDT");

        // Set variables for third view
        alpha = fFRANSObj2->GetAlpha();
        eta = fFRANSObj2->GetEta();
        delta = fFRANSObj2->GetDelta();
        fitScore = fFRANSObj2->GetFitScore();
        iota = fFRANSObj2->GetIota();
        // get the score
        score2 = fTMVAReaders["2"]->EvaluateMVA("FRANS BDT");

        // Print the scores
        //std::cout<<ievt<<" "<<score0<<" "<<score1<<" "<<score2<<std::endl;

        // Fill the output tree
        outputTree->Fill();
    }
    
    // Close
    outputFile->Write();
    outputFile->Close();

    // Read the tree back
    TFile* f = TFile::Open("CombinedViews.root");
    TTree* t = (TTree*)f->Get("FRANSCombinedViews");

    if(fTrainBDT){
        // Train a BDT with the three scores
        //--------- Output file
        gSystem->Exec( "rm -rf TMVAResults" );
        gSystem->Exec( "mkdir TMVAResults" );
        std::string fOutputTMVAROOtFileName = "TMVAResults/TMVACombinedViewsResults.root";
        TFile* outputFileTMVA = TFile::Open( fOutputTMVAROOtFileName.c_str(), "CREATE" );
        //---------  Load TMVA
        TMVA::Tools::Instance();
    
        TMVA::Factory *factory = new TMVA::Factory( "FRANSCombinedViews", outputFileTMVA,
                                                "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
        TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");

        // Add variables to the dataloader
        dataloader->AddVariable("score0", "Score 0", "", 'D');
        dataloader->AddVariable("score1", "Score 1", "", 'D');
        dataloader->AddVariable("score2", "Score 2", "", 'D');

        // Add signal and background trees
        Double_t signalWeight     = 1.0;
        Double_t backgroundWeight = 1.0;
        dataloader->AddSignalTree    ( t,     signalWeight );
        dataloader->AddBackgroundTree( t, backgroundWeight );

        //--------- Prepare training and test trees
        std::string tmva_options = ":SplitMode=Random:NormMode=NumEvents:!V";
        dataloader->PrepareTrainingAndTestTree( "isSignal==1", "isSignal==0", tmva_options );

        // Book the BDT
        factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );
        factory->BookMethod( dataloader, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );
        //factory->BookMethod( dataloader, TMVA::Types::kCuts, "Cuts","!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );
        //factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "Likelihood","H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );
        //factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
        //-------- Train MVAs using the set of training events
        factory->TrainAllMethods();

        //--------  Evaluate all MVAs using the set of test events
        factory->TestAllMethods();

        //--------  Evaluate and compare performance of all configured MVAs
        factory->EvaluateAllMethods();

        //--------  Save the output
        outputFileTMVA->Close();

        TMVA::TMVAGui( fOutputTMVAROOtFileName.c_str() );
    }

    TCut signalCut("isSignal==1");
    TCut bgCut("isSignal==0");



    // Output canvas
    gStyle->SetNdivisions(505, "X");
    gStyle->SetPadBottomMargin(0.15);
    TCanvas* c = new TCanvas("c", "Score comparison", 800, 800);
    std::vector<TPad*> Tp = buildpadcanvas(2,2);

    //--------- Score histograms for the different views
    double binLow = -0.5;
    double binHigh = 0.5;
    int nBins = 25;    
    TH1F hScore0Signal("hScore0Signal", "U plane;FRANS Score;AU", nBins, binLow, binHigh);
    TH1F hScore1Signal("hScore1Signal", "V plane;FRANS Score;AU", nBins, binLow, binHigh);
    TH1F hScore2Signal("hScore2Signal", "C pane;FRANS Score;AU", nBins, binLow, binHigh);
    TH1F hScore0Background("hScore0Background", "U plane;FRANS Score;AU", nBins, binLow, binHigh);
    TH1F hScore1Background("hScore1Background", "V plane;FRANS Score;AU", nBins, binLow, binHigh);
    TH1F hScore2Background("hScore2Background", "C plane;FRANS Score;AU", nBins, binLow, binHigh);
    // Set stats
    hScore0Signal.SetStats(0);
    hScore1Signal.SetStats(0);
    hScore2Signal.SetStats(0);
    hScore0Background.SetStats(0);
    hScore1Background.SetStats(0);
    hScore2Background.SetStats(0);
    // Line widths
    int fLineWidth = 3;

    int fSignalLS = kSolid;
    int fBackgroundLS = kDashed;
    int fSignalLC = kAzure+1;
    int fBackgroundLC = kOrange+1;

    hScore0Signal.SetLineWidth(fLineWidth);
    hScore1Signal.SetLineWidth(fLineWidth);
    hScore2Signal.SetLineWidth(fLineWidth);
    hScore0Background.SetLineWidth(fLineWidth);
    hScore1Background.SetLineWidth(fLineWidth);
    hScore2Background.SetLineWidth(fLineWidth);

    // Draw the histograms for signal and background (view 0)
    Tp.at(1)->cd();
    t->Draw("score0>>hScore0Signal", signalCut);
    t->Draw("score0>>hScore0Background", bgCut);
    // Area normalized
    hScore0Signal.Scale(1./hScore0Signal.Integral());
    hScore0Background.Scale(1./hScore0Background.Integral());
    hScore0Background.Draw("hist");
    hScore0Signal.Draw("hist same");
    //Colors
    hScore0Signal.SetLineColor(fSignalLC);
    hScore0Background.SetLineColor(fBackgroundLC);
    // Lines
    hScore0Signal.SetLineStyle(fSignalLS);
    hScore0Background.SetLineStyle(fBackgroundLS);
    // Offset
    hScore0Background.GetYaxis()->SetTitleOffset(1.3);
    // Legend
    TLegend* legend = new TLegend(0.55, 0.55, 0.8, 0.75);
    legend->AddEntry(&hScore0Signal, "Signal", "l");
    legend->AddEntry(&hScore0Background, "Background", "l");
    legend->SetBorderSize(0);
    legend->Draw("same");

    // Draw the histograms for signal and background (view 1)
    Tp.at(2)->cd();
    t->Draw("score1>>hScore1Background", bgCut);
    t->Draw("score1>>hScore1Signal", signalCut);
    // Area normalized
    hScore1Signal.Scale(1./hScore1Signal.Integral());
    hScore1Background.Scale(1./hScore1Background.Integral());
    hScore1Background.Draw("hist");
    hScore1Signal.Draw("hist same");
    // Colors
    hScore1Signal.SetLineColor(fSignalLC);
    hScore1Background.SetLineColor(fBackgroundLC);
    // Lines
    hScore1Signal.SetLineStyle(fSignalLS);
    hScore1Background.SetLineStyle(fBackgroundLS);
    // Offset
    hScore1Background.GetYaxis()->SetTitleOffset(1.3);
    // Legend
    legend->Draw("same");
    
    // Draw the histograms for signal and background (view 2)
    Tp.at(3)->cd();
    t->Draw("score2>>hScore2Signal", signalCut);
    t->Draw("score2>>hScore2Background", bgCut);
    // Area normalized
    hScore2Signal.Scale(1./hScore2Signal.Integral());
    hScore2Background.Scale(1./hScore2Background.Integral());
    hScore2Background.Draw("hist");
    hScore2Signal.Draw("hist same");
    //Colors
    hScore2Signal.SetLineColor(fSignalLC);
    hScore2Background.SetLineColor(fBackgroundLC);
    // Offset
    hScore2Background.GetYaxis()->SetTitleOffset(1.3);
    // Lines
    hScore2Signal.SetLineStyle(fSignalLS);
    hScore2Background.SetLineStyle(fBackgroundLS);
    // Legend
    legend->Draw("same");

    c->cd();
    c->Update();
    c->WaitPrimitive();
    c->Update();
    
    // Save
    c->SaveAs("ScoreComparison.pdf");
    // Save PADs
    Tp.at(1)->SaveAs("FRANSBDTScore0.eps");
    Tp.at(2)->SaveAs("FRANSBDTScore1.eps");
    Tp.at(3)->SaveAs("FRANSBDTScore2.eps");


    // Draw signal and background efficiency curves for each ROOT file
    TCanvas* c1 = new TCanvas("c1", "Efficiency comparison", 800, 600);
    c1->Divide(2, 2);
    // bottom margins
    
    // Store map ROC for each cut
    std::map<std::string, std::vector<double>> signalEfficiencies;
    std::map<std::string, std::vector<double>> backgroundRejections;
    std::map<std::string, std::vector<double>> significances;

    // Define the different cut combinations
    std::map<std::string, TString> cutsMap = {
        { "C", "score2>%.2f"}
        , { "V", "score1>%.2f"}
        , { "U", "score0>%.2f"}
        , { "C | (V & U)", "score2>%.2f || (score1>%.2f && score0>%.2f)"}
        , { "C & (V | U)", "score2>%.2f && (score1>%.2f || score0>%.2f || (score1>%.2f && score0>%.2f) )"}
        , { "C & V & U", "score2>%.2f && score1>%.2f && score0>%.2f"}
        , { "C | V | U", "score2>%.2f || score1>%.2f || score0>%.2f"}
    };

    std::vector<int> colors = {kRed, kBlue, kGreen+3, kViolet+2, kMagenta, kOrange+7, kBlack};

    // initialize vector maps
    for (std::pair<std::string, TString> cutIt : cutsMap) {
        std::string cutLabel = cutIt.first;
        signalEfficiencies[cutLabel] = {};
        backgroundRejections[cutLabel] = {};
    }

    // Score scan
    double scoreStep = 0.01;
    double scoreShift = 0.;
    double nSignal = t->GetEntries(signalCut);
    double nBackground = t->GetEntries(bgCut);
    // Vector to store the score
    std::vector<double> scores;
    
    // Loop over the scores
    for(double sc=-1; sc<1; sc+=scoreStep){
        scores.push_back(sc);

        std::cout<<"--- Score: "<<sc<<std::endl;
        // loop over the cuts map
        for (std::pair<std::string, TString> cutIt : cutsMap) {
            TString cut = cutIt.second;
            std::cout<<"  "<<cutIt.first<<"  Cut: "<<cut<<std::endl;
            cut = Form(cut.Data(), sc, sc+scoreShift, sc+scoreShift, sc+scoreShift, sc+scoreShift, sc+scoreShift);
            
            std::cout<<"Cut: "<<cut<<std::endl;
            // Get number od signal and background events
            double nSignalPass = t->GetEntries( signalCut && cut );
            double nBackgroundPass = t->GetEntries( bgCut && cut );

            // Calculate efficiencies
            double signalEfficiency = 100*nSignalPass/nSignal;
            double backgroundRejection = 100-100*nBackgroundPass/nBackground;

            //std::cout<<"  Signal efficiency: "<<signalEfficiency<<" "<<backgroundRejection<<std::endl;

            // Store
            signalEfficiencies[cutIt.first].push_back(signalEfficiency);
            backgroundRejections[cutIt.first].push_back(backgroundRejection);
            if(nBackgroundPass > 0)
                significances[cutIt.first].push_back(signalEfficiency/sqrt(nBackgroundPass));
            else
                significances[cutIt.first].push_back(0);
        }
    }

    // Draw the signal efficiency curves
    c1->cd(1);
    // aux histogram
    TH2F hSignalEffFrame("hSignalEffFrame", ";FRANS Score;Signal efficiency", 100, binLow, binHigh, 100, 0, 100);
    hSignalEffFrame.SetStats(0);
    hSignalEffFrame.Draw();
    // loop
    for(auto it = signalEfficiencies.begin(); it != signalEfficiencies.end(); ++it){
        TGraph* gSignalEff = new TGraph(scores.size(), &scores[0], &it->second[0]);
        gSignalEff->SetLineColor(colors.at(std::distance(signalEfficiencies.begin(), it)));
        gSignalEff->SetLineWidth(2);
        gSignalEff->Draw("L SAME");
    }
    // legend
    TLegend* legendROC = new TLegend(0.15, 0.15, 0.35, 0.45);
    for(auto it = signalEfficiencies.begin(); it != signalEfficiencies.end(); ++it){
        TGraph* gSignalEff = new TGraph(scores.size(), &scores[0], &it->second[0]);
        gSignalEff->SetLineColor(colors.at(std::distance(signalEfficiencies.begin(), it)));
        gSignalEff->SetLineWidth(2);
        legendROC->AddEntry(gSignalEff, (it->first).c_str(), "l");
    }
    legendROC->Draw("same");

    // Draw the background rejection curves
    c1->cd(2);
    // aux histogram
    TH2F hBackgroundRejFrame("hBackgroundRejFrame", ";FRANS Score;Background rejection", 100, -1, 1, 100, 0, 100);
    hBackgroundRejFrame.SetStats(0);
    hBackgroundRejFrame.Draw();
    // loop
    for(auto it = backgroundRejections.begin(); it != backgroundRejections.end(); ++it){
        TGraph* gBackgroundRej = new TGraph(scores.size(), &scores[0], &it->second[0]);
        gBackgroundRej->SetLineColor(colors.at(std::distance(backgroundRejections.begin(), it)));
        gBackgroundRej->SetLineWidth(2);
        gBackgroundRej->Draw("L SAME");
    }
    // legend
    legendROC->Draw("same");

    

    // Draw the ROC curve
    c1->cd(3);
    // aux histogram
    TH2F hROCFrame("hROCFrame", ";Signal efficiency;Background rejection", 500, 0, 100, 500, 0, 100);
    hROCFrame.SetStats(0);
    hROCFrame.Draw();
    // legend
    // Loop over the ROCs
    int i = 0;
    for (std::pair<std::string, TString> cutIt : cutsMap) {
        TGraph* gROC = new TGraph(signalEfficiencies[cutIt.first].size(), &signalEfficiencies[cutIt.first][0], &backgroundRejections[cutIt.first][0]);
        gROC->SetLineColor(colors.at(i));
        gROC->SetLineWidth(2);
        gROC->Draw("L SAME");
        i++;
        // cout ROC integral
        std::cout<<cutIt.first<<" "<<1e-4*gROC->Integral()<<std::endl;
    }
    legendROC->Draw("same");

    // Draw significance curve
    c1->cd(4);
    // aux histogram
    double maxSignificance = 0;
    for(auto it = significances.begin(); it != significances.end(); ++it){
        // get max
        double maxSignificanceCut = *std::max_element(it->second.begin(), it->second.end());
        if(maxSignificanceCut > maxSignificance) maxSignificance = maxSignificanceCut;
    }
    std::cout<<"Max significance: "<<maxSignificance<<std::endl;
    TH2F hSignificanceFrame("hSignificanceFrame", ";FRANS Score;S/#sqrt{BG}", 100, -1, 1, 100, 0, maxSignificance);
    hSignificanceFrame.SetStats(0);
    hSignificanceFrame.Draw();
    for(auto it = significances.begin(); it != significances.end(); ++it){
        TGraph* gSignificance = new TGraph(scores.size(), &scores[0], &it->second[0]);
        gSignificance->SetLineColor(colors.at(std::distance(significances.begin(), it)));
        gSignificance->SetLineWidth(2);
        gSignificance->Draw("L SAME");
    }
    // legend
    legendROC->Draw("same");
    
    
    


    c1->cd();
    c1->Update();


    // Make 2D scatter plots for different vies combinations, for signal and background
    TCanvas* c2 = new TCanvas("c2", "2D scatter plots", 800, 600);
    c2->Divide(2, 2);
    double fMarkerSize = 0.2;

    gStyle->SetOptStat(0);
    

    c2->cd(1);
    // Define the histograms
    TH2F hScore0Score1Signal("hScore0Score1Signal", "Signal;Score 0;Score 1", 100, -0.5, 0.5, 100, -0.5, 0.5);
    TH2F hScore0Score1Background("hScore0Score1Background", "Background;Score 0;Score 1", 100, -0.5, 0.5, 100, -0.5, 0.5);
    t->Draw("score0:score1>>hScore0Score1Signal", signalCut);
    t->Draw("score0:score1>>hScore0Score1Background", bgCut);
    hScore0Score1Background.Draw();
    hScore0Score1Signal.Draw("same");
    //Colors and markers
    hScore0Score1Signal.SetMarkerColorAlpha(kOrange-3, 0.5);
    hScore0Score1Background.SetMarkerColorAlpha(kAzure+2, 0.5);
    hScore0Score1Signal.SetMarkerStyle(20);
    hScore0Score1Background.SetMarkerStyle(29);
    hScore0Score1Signal.SetMarkerSize(fMarkerSize);
    hScore0Score1Background.SetMarkerSize(fMarkerSize);
    // Legend
    TLegend* legend2D = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend2D->AddEntry(&hScore0Score1Signal, "Signal", "p");
    legend2D->AddEntry(&hScore0Score1Background, "Background", "p");
    legend2D->Draw("same");

    c2->cd(3);
    // Define the histograms
    TH2F hScore2Score0Signal("hScore2Score0Signal", "Signal;Score 2;Score 0", 100, -0.5, 0.5, 100, -0.5, 0.5);
    TH2F hScore2Score0Background("hScore2Score0Background", "Background;Score 2;Score 0", 100, -0.5, 0.5, 100, -0.5, 0.5);
    t->Draw("score0:score2>>hScore2Score0Signal", signalCut);
    t->Draw("score0:score2>>hScore2Score0Background", bgCut);
    hScore2Score0Background.Draw();
    hScore2Score0Signal.Draw("same");
    //Colors and markers
    hScore2Score0Signal.SetMarkerColorAlpha(kOrange-3, 0.5);
    hScore2Score0Background.SetMarkerColorAlpha(kAzure+2, 0.5);
    hScore2Score0Signal.SetMarkerStyle(20);
    hScore2Score0Background.SetMarkerStyle(29);
    hScore2Score0Signal.SetMarkerSize(fMarkerSize);
    hScore2Score0Background.SetMarkerSize(fMarkerSize);
    legend2D->Draw("same");

    c2->cd(4);
    // Define the histograms
    TH2F hScore2Score1Signal("hScore2Score1Signal", "Signal;Score 2;Score 1", 100, -0.5, 0.5, 100, -0.5, 0.5);
    TH2F hScore2Score1Background("hScore2Score1Background", "Background;Score 2;Score 1", 100, -0.5, 0.5, 100, -0.5, 0.5);
    t->Draw("score1:score2>>hScore2Score1Signal", signalCut);
    t->Draw("score1:score2>>hScore2Score1Background", bgCut);
    hScore2Score1Background.Draw();
    hScore2Score1Signal.Draw("same");
    //Colors and markers
    hScore2Score1Signal.SetMarkerColorAlpha(kOrange-3, 0.5);
    hScore2Score1Background.SetMarkerColorAlpha(kAzure+2, 0.5);
    hScore2Score1Signal.SetMarkerStyle(20);
    hScore2Score1Background.SetMarkerStyle(29);
    hScore2Score1Signal.SetMarkerSize(fMarkerSize);
    hScore2Score1Background.SetMarkerSize(fMarkerSize);
    legend2D->Draw("same");

    c2->cd();
    c2->Update();


    c->WaitPrimitive();
    c->SaveAs("ScoreComparison.pdf");
    c1->WaitPrimitive();
    c1->SaveAs("EfficiencyComparison.pdf");
    c2->WaitPrimitive();
    c2->SaveAs("2DScatterPlots.pdf");





}
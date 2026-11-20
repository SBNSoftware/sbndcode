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

#include "CutEfficienciesDefinitions.C"

using namespace TMVA::Experimental;


// ---- Function to add variables to the dataloader and get the minimal cut
void GetMVAConfigurationVariables(TMVA::DataLoader *dataloader, TCut &minimalCut, std::string configFileMVA) {
    // Open config file
    std::ifstream configFile(configFileMVA);
    std::string line;

    // Minimal cut: first line of the file        
    std::getline(configFile, line);
    minimalCut = line.c_str();
    std::cout<<" Minimal cut: "<<minimalCut<<std::endl;

    // Loop over the rest of the lines and save to vector
    std::vector<std::string> variables;
    while (std::getline(configFile, line)) {
        variables.push_back(line);
        std::cout<<" Variable: "<<line<<std::endl;
    }

    // Function to check if a variable is in the vector
    auto is_in = [](const std::vector<std::string> &v, const std::string &value) {
        return std::find(v.begin(), v.end(), value) != v.end();
    };

    // Add variables to the dataloader
    // --- Angle variables
    if (is_in(variables, "NAngles")) dataloader->AddVariable( "NAngles", "N_{angles}", "", 'I' );
    if (is_in(variables, "AngleFRANSScore")) dataloader->AddVariable( "AngleFRANSScore", "FRANSScore Angle", "", 'D' );
    if (is_in(variables, "FRANSScorePANDORA")) dataloader->AddVariable( "FRANSScorePANDORA", "FRANSScore PANDORA", "", 'D' );
    if (is_in(variables, "AngleGap")) dataloader->AddVariable( "AngleGap", "Gap", "", 'D' );
    if (is_in(variables, "AngleNHitsMainTrack")) dataloader->AddVariable( "AngleNHitsMainTrack", "Main track # hits", "", 'I' );
    if (is_in(variables, "AngleNHitsTrack1")) dataloader->AddVariable( "AngleNHitsTrack1", "Track 1 # hits", "", 'I' );
    if (is_in(variables, "AngleNHitsTrack2")) dataloader->AddVariable( "AngleNHitsTrack2", "Track 2 # hits", "", 'I' );
    if (is_in(variables, "AngleMinNHits")) dataloader->AddVariable( "AngleMinNHits", "# hits min", "", 'I' );
    if (is_in(variables, "AngleLengthTrack1")) dataloader->AddVariable( "AngleLengthTrack1", "Track 1 Length [cm]", "", 'I' );
    if (is_in(variables, "AngleLengthTrack2")) dataloader->AddVariable( "AngleLengthTrack2", "Track 2 Length [cm]", "", 'I' );
    if (is_in(variables, "AngleLengthMainTrack")) dataloader->AddVariable( "AngleLengthMainTrack", "Main track Length [cm]", "", 'I' );
    if (is_in(variables, "AngleLongestIsMain")) dataloader->AddVariable( "AngleLongestIsMain", "LongestIsMain", "", 'I' );
    if (is_in(variables, "CRUMBSScore")) dataloader->AddVariable( "CRUMBSScore", "CRUMBS", "", 'D' );
    // --- Origin variables
    if (is_in(variables, "NUnOrigins")) dataloader->AddVariable( "NUnOrigins", "N", "", 'I' );
    if (is_in(variables, "NUnOriginsMultGT3")) dataloader->AddVariable( "NUnOriginsMultGT3", "N^{2}", "", 'I' );
    if (is_in(variables, "NOrigins")) dataloader->AddVariable( "NOrigins", "N_", "", 'I' );
    if (is_in(variables, "NOriginsMultGT3")) dataloader->AddVariable( "NOriginsMultGT3", "N^{2}", "", 'I' );
    // --- Cleaness
    if (is_in(variables, "AngleCoveredArea")) dataloader->AddVariable( "AngleCoveredArea", "Covered Area", "", 'D' );
    if (is_in(variables, "AngleDirtHits")) dataloader->AddVariable( "AngleDirtHits", "Dirt Hits", "", 'I' );
    if (is_in(variables, "NUnassociatedHits")) dataloader->AddVariable( "NUnassociatedHits", "# unassociated hits", "", 'I' );
    if (is_in(variables, "NMaxDirtUnassociatedHits")) dataloader->AddVariable( "NMaxDirtUnassociatedHits", "Max Dirt Unassociated Hits", "", 'I' );
    // --- Kinematics
    if (is_in(variables, "AngleDecayContainedDiff")) dataloader->AddVariable( "AngleDecayContainedDiff", "#alpha", "", 'D' );
    if (is_in(variables, "AngleOpeningAngle")) dataloader->AddVariable( "AngleOpeningAngle", "Opening Angle [º]", "", 'I' );
    if (is_in(variables, "AngleMainTrackOverlap")) dataloader->AddVariable( "AngleMainTrackOverlap", "Main Track Overlap", "", 'D' );
    if (is_in(variables, "AnglePzSign")) dataloader->AddVariable( "AnglePzSign", "AnglePzSign", "", 'D' );
    if (is_in(variables, "AngleGapOverlapWithAPAJuntion")) dataloader->AddVariable( "AngleGapOverlapWithAPAJuntion", "Gap Overlap", "", 'D' );
    // --- Calorimetry
    if (is_in(variables, "AngleTwoLinesChi2")) dataloader->AddVariable( "AngleTwoLinesChi2", "Two Lines Chi2", "", 'D' );
    if (is_in(variables, "AnglePassFit")) dataloader->AddVariable( "AnglePassFit", "PassFit", "", 'I' );
    if (is_in(variables, "AnglePassChargeFit")) dataloader->AddVariable( "AnglePassChargeFit", "Pass Charge Fit", "", 'I' );
    if (is_in(variables, "AngleBandOverlap")) dataloader->AddVariable( "AngleBandOverlap", "B and Overlap", "", 'D' );
    if (is_in(variables, "AngleBandCrossHits")) dataloader->AddVariable( "AngleBandCrossHits", "Band Overlap", "", 'D' );
    if (is_in(variables, "AngleChargeRatioFit")) dataloader->AddVariable( "AngleChargeRatioFit", "Charge Ratio Fit", "", 'D' );
    if (is_in(variables, "AngleChargeDifferenceFit")) dataloader->AddVariable( "AngleChargeDifferenceFit", "Charge Difference Fit", "", 'D' );
    if (is_in(variables, "AngleChargeRatioIntegral")) dataloader->AddVariable( "AngleChargeRatioIntegral", "Charge Ratio Integral", "", 'D' );
    if (is_in(variables, "AngleChargeDifferenceIntegral")) dataloader->AddVariable( "AngleChargeDifferenceIntegral", "Charge Difference", "", 'D' );
    if (is_in(variables, "AngleChargeRatioAverage")) dataloader->AddVariable( "AngleChargeRatioAverage", "Charge Ratio Average", "", 'D' );
    if (is_in(variables, "AngleVertexHitIntegralRatio")) dataloader->AddVariable( "AngleVertexHitIntegralRatio", "Vertex Hit Integral Ratio", "", 'D' );
    if (is_in(variables, "AngleTrackLengthRatio")) dataloader->AddVariable( "AngleTrackLengthRatio", "Track Length Ratio", "", 'D' );
    if (is_in(variables, "AngleResidualRange1RMS")) dataloader->AddVariable( "AngleResidualRange1RMS", "Residuals Range 1 RMS", "", 'D' );
    if (is_in(variables, "AngleResidualRange2RMS")) dataloader->AddVariable( "AngleResidualRange2RMS", "Residuals Range 2 RMS", "", 'D' );
    if (is_in(variables, "AngleResidualRangeMinRMS")) dataloader->AddVariable( "AngleResidualRangeMinRMS", "Residuals Range Min RMS", "", 'D' );
    if (is_in(variables, "AngleResidualRangeMaxAngleRMS")) dataloader->AddVariable( "AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS", "", 'D' );
    if (is_in(variables, "AngleNVertexHits")) dataloader->AddVariable( "AngleNVertexHits", "N Vertex Hits", "", 'I' );
    if (is_in(variables, "AngleNBulkHits")) dataloader->AddVariable( "AngleNBulkHits", "N Bulk Hits", "", 'I' );
    // --- Showers
    if (is_in(variables, "NShowers")) dataloader->AddVariable( "NShowers", "# showers", "", 'I' );
    if (is_in(variables, "NShowerHits")) dataloader->AddVariable( "NShowerHits", "# shower hits", "", 'I' );
    if (is_in(variables, "ShowerEnergy")) dataloader->AddVariable( "ShowerEnergy", "Shower Energy", "", 'D' );   
    // --- PID
    if (is_in(variables, "NTracksLI")) dataloader->AddVariable( "NTracksLI", "NTracksLI", "", 'I' );
    if (is_in(variables, "NTracksHI")) dataloader->AddVariable( "NTracksHI", "NTracksHI", "", 'I' );
}


//---------  Main function
void RunLambdaMVAAnalysis(std::string fInputFileName="", double nTrainFrac = -1, std::string configFile="configMVA.txt", std::string method="BDT", std::string fTreeDirName = "originsAna/", std::string fTreeName = "LambdaAnaTree")
{

    // ---- Output file
    gSystem->Exec( "rm -rf TMVAResultsTest" );
    gSystem->Exec( "mkdir TMVAResultsTest" );
    std::string fOutputTMVAROOtFileName = "TMVAResultsTest/TMVAResults.root";
    TFile* outputFile = TFile::Open( fOutputTMVAROOtFileName.c_str(), "CREATE" );

    // ---- Input file
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());

    // ---- Load the library
    TMVA::Tools::Instance();

    // ---- Define the MVA methods to be trained
    std::map<std::string,int> Use;
    // Rectangular cut optimisation
    if(method=="Cuts") Use["Cuts"] = 1;
    else Use["Cuts"] = 0;
    // 1-dimensional likelihood ("naive Bayes estimator")
    if(method=="Likelihood") Use["Likelihood"] = 1;
    else Use["Likelihood"] = 0;
    // BDT
    if(method=="BDT") Use["BDT"] = 1;
    else Use["BDT"] = 0;
    // Neural network (MLP)
    if(method=="MLP") Use["MLP"] = 1;
    else Use["MLP"] = 0;
    
    // ---- Define factory and data loader
    TMVA::Factory *factory = new TMVA::Factory( "FRANSSelectionTMVA", outputFile,
                                                "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");

    // ---- Add variables from configuration file
    TCut fMinimalCut;
    GetMVAConfigurationVariables(dataloader, fMinimalCut, configFile);
    
    Double_t signalWeight     = 1.0;
    Double_t backgroundWeight = 1.0;
    dataloader->AddSignalTree    ( fTree,     signalWeight );
    dataloader->AddBackgroundTree( fTree, backgroundWeight );
    dataloader->AddCut(fMinimalCut);

    // ----  Define signal and background cuts
    TCut fTruthInFV("TruthIsFiducial==1");
    TCut fTruthInAV("TruthIsAV==1");
    TCut signalCut = fTruthInFV + TCut("IntOrigin==1 && IntDirt==0 && (IntNLambda>0 && IntMode==0 && abs(IntNuPDG)!=12)");
    TCut bgCut = fTruthInFV && TCut("IntOrigin==1 && IntDirt==0 && !(IntNLambda>0 && IntMode==0 && abs(IntNuPDG)!=12)");
    int nMaxSignal = fTree->Draw(">>selectedEntries", signalCut + fMinimalCut, "entrylist")-1;
    int nMaxBg = fTree->Draw(">>selectedEntries", bgCut + fMinimalCut, "entrylist")-1;
    std::cout<<"nMaxSignal: "<<nMaxSignal<<std::endl;
    std::cout<<"nMaxBg: "<<nMaxBg<<std::endl;

    // ---- Prepare the training and test trees
    std::string tmva_options = "";
    if(nTrainFrac<0) tmva_options = "nTrain_Signal="+std::to_string(nMaxSignal)+":nTrain_Background="+to_string(nMaxBg);
    else tmva_options = "nTrain_Signal="+std::to_string((int)nTrainFrac*nMaxSignal)+":nTrain_Background="+to_string((int)nTrainFrac*nMaxBg);
    tmva_options+=":SplitMode=Random:NormMode=NumEvents:!V";
    dataloader->PrepareTrainingAndTestTree( signalCut, bgCut, tmva_options );

    // ---- Cut optimisation
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
    if (Use["RandomForest"])  // Adaptive Boost
        factory->BookMethod( dataloader, TMVA::Types::kBDT, "RandomForest",
                             "!H:!V:NTrees=2500:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

    // ---- Train MVAs using the set of training events
    factory->TrainAllMethods();

    // ---- Evaluate all MVAs using the set of test events
    factory->TestAllMethods();

    // ---- Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();


    // ---- Save the output
    outputFile->Close();

    delete factory;
    delete dataloader;

    // ---- Print cuts
    if(Use["Cuts"]) {
        std::cout<<"Cuts: "<<std::endl;
    }
    
    // ---- Launch the GUI for the root macros
    TMVA::TMVAGui( fOutputTMVAROOtFileName.c_str() );


}
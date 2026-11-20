#include <string>
#include <iostream>
#include <cmath>

#include "TString.h"
#include <TApplication.h>
#include "TImage.h"
#include "TF1.h"

#include "CommandLineParser.h"
#include "SEventHandle.h"
#include "LambdaAnaTree.h"
#include "SParserConfig.h"
#include "STPCAnalyzerTreeReader.h"
#include "TPCSimpleEvents.h"
#include "TPCLinesParameters.h"
#include "TPCLinesAlgo.h"

void RunAlgoTPCLines(const CommandLineParser& parser)
{

    // ------------------------------------------------------------------ 
    // Define the program control variables
    int Debug = parser.getDebug();
    std::string ConfPsetPath = parser.getPsetPath();
    if(ConfPsetPath==""){
        ConfPsetPath="config.txt";
    }
    int DebugMode = parser.getDebugMode();
    int n = parser.getN();
    int nskip = parser.getNskip();
    int event = parser.getEvent();
    int sr = parser.getSr();
    std::string file_name = parser.getFileName();
    std::string directory_path = parser.getDirectoryPath();
    std::string ext = parser.getExtension();

    // Get input files
    std::vector<TString> fFilePaths = GetInputFileList(file_name, ext, directory_path);

    // Set batch mode
    if(Debug==0) gROOT->SetBatch(true);

    int fNEv = n;
    int fEv = event;
    int fSubRun = sr;
    int fNEvSkip = nskip;
    int nEvents=0;

    double fFRANSScoreCut = 0.1;
    // ------------------------------------------------------------------ 


    // ---- TPCLines parameters ----------------------------------------
    CaloAlgorithmPsetType fPsetCalo = ReadCaloAlgorithmPset( FindFile("caloalg_config.fcl"), "CaloAlg:");
    fPsetCalo.Verbose = Debug;
    TrackFinderAlgorithmPsetType fPsetTrackFinder = ReadTrackFinderAlgorithmPset( FindFile("trackfinderalg_config.fcl"), "TrackFinderAlg:");
    fPsetTrackFinder.Verbose = Debug;
    HoughAlgorithmPsetType fPsetHough = ReadHoughAlgorithmPset( FindFile("houghalg_config.fcl"), "HoughAlg:");
    fPsetHough.Verbose = Debug;
    VertexFinderAlgorithmPsetType fPsetVertexFinder = ReadVertexFinderAlgorithmPset( FindFile("vertexfinderalg_config.fcl"), "VertexFinderAlg:");
    fPsetVertexFinder.Verbose = Debug;
    TPCLinesAlgoPsetType fPsetAnaView = ReadTPCLinesAlgoPset( FindFile("tpclinesalg_config.fcl"), "TPCLinesAlg:");
    fPsetAnaView.Verbose = Debug;
    fPsetAnaView.DebugMode = DebugMode;
    fPsetAnaView.HoughAlgorithmPset = fPsetHough;
    fPsetAnaView.TrackFinderAlgorithmPset = fPsetTrackFinder;
    fPsetAnaView.VertexFinderAlgorithmPset = fPsetVertexFinder;
    fPsetAnaView.CaloAlgorithmPset = fPsetCalo;
    //fPsetAnaView.View = parser.getView();
    fPsetAnaView.Print();
    // ---- FRANS parameters ----------------------------------------
    FRANSPsetType fPsetFRANS = ReadFRANSPset( FindFile("chargedensityalg_config.fcl"), "ChargeDensityAlg:");
    fPsetFRANS.Verbose = Debug;
    fPsetFRANS.TMVAFilename = FindFile(fPsetFRANS.TMVAFilename);
    if(parser.getView()==0)
        fPsetFRANS.TMVAFilename = FindFile("FRANSSelectionTMVA_BDT_ViewUWithWidth_Iota.weights.xml");
    else if(parser.getView()==1)
        fPsetFRANS.TMVAFilename = FindFile("FRANSSelectionTMVA_BDT_ViewVWithWidth_Iota.weights.xml");
    std::cout<<"FRANS TMVA file: "<<fPsetFRANS.TMVAFilename<<std::endl;
    fPsetFRANS.Print();

    // ---- Output files ----------------------------------------
    gSystem->Exec(("rm -rf "+fPsetAnaView.OutputPath).c_str());
    if (!gSystem->OpenDirectory(fPsetAnaView.OutputPath.c_str())) {    
        gSystem->Exec( ("mkdir "+fPsetAnaView.OutputPath).c_str());
        gSystem->Exec( ("mkdir "+fPsetAnaView.OutputPath+"/rootfiles").c_str());
    }
    TFile* anaOutputFile = new TFile("LambdaAnaOutput.root", "RECREATE");
    TDirectory *fAnaTreeHandleDirectory = anaOutputFile->mkdir("originsAna");
    fAnaTreeHandleDirectory->cd();
    // Output ana tree
    TTree *fTreeAna = new TTree("LambdaAnaTree", "LambdaAnaTree");
    LambdaAnaTree fAnaTreeHandle;
    fAnaTreeHandle.SetTree(fTreeAna);

    // Define TPC LINES ALGORITHM
    fPsetAnaView.View = parser.getView();
    TPCLinesAlgo _TPCLinesAlgo(fPsetAnaView, fPsetFRANS);
    fPsetAnaView.View = 0;
    fPsetFRANS.CalculateScore = false;
    TPCLinesAlgo _TPCLinesAlgoViewU(fPsetAnaView, fPsetFRANS);
    fPsetAnaView.View = 1;
    TPCLinesAlgo _TPCLinesAlgoViewV(fPsetAnaView, fPsetFRANS);
    
    // Effiency status
    EfficiencyCalculator _EfficiencyCalculator;

    // TTree loop
    std::vector<int> fNVertex;
    int nEntries = 0;
    int fileCounter = 0;
    for (const auto& filepath : fFilePaths) {
        fileCounter++;

        std::cout<<" ANALYZING THE FILE: "<<filepath<<std::endl;
        MyTPCTreeReader treeReader(filepath, (parser.getTreeName()+"/AnaTPCTree").c_str());
        
        for (int entry = 0; entry < treeReader.NEntries(); entry++) {
                
            // --- get the entry
            std::cout<<"Entry "<<entry<<std::endl;
            treeReader.GetEntry(entry);

            // --- Check program control variables
            SEventId ev(treeReader.runID, treeReader.subrunID, treeReader.eventID);
            if (fNEv > 0 && nEvents >= fNEv) continue;
            if (treeReader.eventID != fEv && fEv != -1) continue;
            if (treeReader.subrunID != fSubRun && fSubRun != -1) continue;
            nEntries++;
            if (nEntries <= fNEvSkip) continue; 
            nEvents++;
            std::cout << "\n\n ************** Analyzing: " << ev;
            std::cout << "   Interaction mode: "<<treeReader.intMode<<" NLambda: "<<treeReader.intNLambda<<" E = "<<treeReader.nuvE<<" GeV" << " T = " <<treeReader.nuvT<<" ns"<<std::endl;
            
            // --- True vertex
            std::vector<double> VertexXYZ = {treeReader.nuvX, treeReader.nuvY, treeReader.nuvZ};
            std::vector<int> VertexUVYT = {treeReader.nuvU, treeReader.nuvV, treeReader.nuvC, treeReader.nuvTimeTick};
            std::cout << "   - True vertex (X, Y, Z, T) " << VertexXYZ[0] << " " << VertexXYZ[1] << " " << VertexXYZ[2] << " " << treeReader.nuvT << " (U, V, C, TT): " << VertexUVYT[0] << " " << VertexUVYT[1] << " " << VertexUVYT[2] << " " << VertexUVYT[3] << std::endl;

            // --- Reco vertex
            std::vector<double> RecoVertexXYZ = {treeReader.recnuvX, treeReader.recnuvY, treeReader.recnuvZ};
            std::vector<int> RecoVertexUVYT = {treeReader.recnuvU, treeReader.recnuvV, treeReader.recnuvC, treeReader.recnuvTimeTick};
            std::cout << "  - Reco vertex (X, Y, Z) " << RecoVertexXYZ[0] << " " << RecoVertexXYZ[1] << " " << RecoVertexXYZ[2] << " (U, V, C, TT): " << RecoVertexUVYT[0] << " " << RecoVertexUVYT[1] << " " << RecoVertexUVYT[2] << " "<< RecoVertexUVYT[3]  << std::endl;

            // --- Set TPC using the drift coordinate of the reco vertex
            //TPC = (nuvX > 0) ? 1 : 0;
            //int TPC = (RecoVertexXYZ[0] > 0) ? 1 : 0;
            
            // --- We need a minimum number of hits to run the track finder
            size_t nhits = treeReader.hitsChannel->size();
            std::cout << "  - NHits: " << nhits << std::endl;
            if(nhits<=3 || treeReader.recnuvU==-1){
                std::cout<<"   SKIPPED NHits or RecoVertex \n";
                _EfficiencyCalculator.UpdateSkipped(ev);
                continue;
            }


            // --- Get and set the hits in the view
            int view = parser.getView();
            std::vector<SHit> hitList = treeReader.GetHitsInView(view);
            bool filled = _TPCLinesAlgo.SetHitList(view, RecoVertexUVYT, VertexUVYT, hitList);
            
            // --- Analyze
            SEvent recoEvent;
            if(filled){
                _TPCLinesAlgo.AnaView(ev.Label());
            }

            // --- Fill the TTree
            fAnaTreeHandle.ResetVars();
            // True vars
            fAnaTreeHandle.fEventID = treeReader.eventID;
            fAnaTreeHandle.fSubrunID = treeReader.subrunID;
            fAnaTreeHandle.fRunID = treeReader.runID;
            fAnaTreeHandle.fSliceID = 1;
            fAnaTreeHandle.fIntMode = treeReader.intMode;
            fAnaTreeHandle.fIntNLambda = treeReader.intNLambda;
            fAnaTreeHandle.fTruthIsFiducial = true;
            fAnaTreeHandle.fRecoIsFiducial = true;
            // Reco vars
            STriangle foundTriangle = _TPCLinesAlgo.FillLambdaAnaTree(fAnaTreeHandle);
            double fHyperonVertexX, fHyperonVertexY, fHyperonVertexZ;
            foundTriangle.GetVertexXYZ(fHyperonVertexX, fHyperonVertexY, fHyperonVertexZ);
            std::cout<<"  - Lambda reco vertex (X, Y, Z) " << fHyperonVertexX << " " << fHyperonVertexY << " " << fHyperonVertexZ << std::endl;
            // Fill the tree
            fAnaTreeHandleDirectory->cd();
            fAnaTreeHandle.FillTree();

            
            // --- Run for the other views
            if(parser.getThreeViews() && foundTriangle.GetMainVertex().X()!=-1 && foundTriangle.GetMainVertex().Y()!=-1){

                // check the V in the other views
                double triangleYlower = foundTriangle.GetMinY()+_TPCLinesAlgo.ShiftY();
                double triangleYupper = foundTriangle.GetMaxY()+_TPCLinesAlgo.ShiftY();
                // add porch
                double porch = 0.2*(triangleYupper-triangleYlower);
                triangleYlower -= porch;
                triangleYupper += porch;

                std::cout<<"  - Filling other views\n";
                // Get the hits in the view
                std::vector<SHit> hitListOther = treeReader.GetHitsInView(0);

                std::vector<SHit> hitListOtherFiltered;
                for(SHit &h:hitListOther){
                    if(h.Y()>triangleYlower && h.Y()<triangleYupper) hitListOtherFiltered.push_back(h);
                }
                std::cout<<" Filling second view 1\n";
                filled = _TPCLinesAlgoViewU.SetHitList(0, RecoVertexUVYT, VertexUVYT, hitListOtherFiltered);
                if(filled){
                    std::cout<<"Filled!\n";
                    _TPCLinesAlgoViewU.AnaView(ev.Label());
                    
                }

                hitListOther.clear();
                hitListOtherFiltered.clear();
                hitListOther = treeReader.GetHitsInView(1);

                
                for(SHit &h:hitListOther){
                    if(h.Y()>triangleYlower && h.Y()<triangleYupper) hitListOtherFiltered.push_back(h);
                }
                std::cout<<" Filling second view 2\n";
                filled = _TPCLinesAlgoViewV.SetHitList(1, RecoVertexUVYT, VertexUVYT, hitListOtherFiltered);
                if(filled){
                    std::cout<<"Filled!\n";
                    _TPCLinesAlgoViewV.AnaView(ev.Label());
                }
            }


            // --- Update the efficiency calculator
            bool accepted = fAnaTreeHandle.fNAngles>0 && fAnaTreeHandle.fAngleFRANSScore>fFRANSScoreCut;
            if(accepted)
                _EfficiencyCalculator.UpdateSelected(ev);
            else
                _EfficiencyCalculator.UpdateNotSelected(ev);
            std::cout<<_EfficiencyCalculator;

            
            // --- Displays
            std::string outNamePreffix = accepted? "Accepted":"Rejected";
            std::string outputLabel = "FinalReco" + outNamePreffix + ev.Label() + "_" + std::to_string(_EfficiencyCalculator.NEvents()); 
            TCanvas *cFRANS = new TCanvas(("canvasFRANS"+ev.Label()).c_str(),"CanvasFRANS", 1000, 0, 1000, 1000);
            TCanvas *cCalo = new TCanvas(("canvasCalo"+ev.Label()).c_str(),"CanvasCalorimetry", 1000, 0, 1000, 1000);
            TCanvas *cTPCDisplayViewU = new TCanvas( ("FinalRecoViewU"+ev.Label()).c_str(),  ("FinalRecoViewU"+ev.Label()).c_str(), 500, 500, 1000, 800);
            TCanvas *cTPCDisplayViewV = new TCanvas( ("FinalRecoViewV"+ev.Label()).c_str(),  ("FinalRecoViewV"+ev.Label()).c_str(), 500, 100, 1000, 800);
            TCanvas *cTPCDisplay = new TCanvas( ("FinalReco"+ev.Label()).c_str(),  ("FinalReco"+ev.Label()).c_str(), 0, 0, 1000, 800);
            _TPCLinesAlgo.Display("", cTPCDisplay, cCalo, cFRANS);
            cTPCDisplay->SaveAs( (fPsetAnaView.OutputPath+"/"+outputLabel+".pdf").c_str() );
            
            fAnaTreeHandle.PrintEventInfo();

            if(parser.getThreeViews()){
                _TPCLinesAlgoViewU.Display("", cTPCDisplayViewU);
                _TPCLinesAlgoViewV.Display("", cTPCDisplayViewV);
                cTPCDisplay->WaitPrimitive();
            }
            else{
                delete cTPCDisplayViewU;
                delete cTPCDisplayViewV;
            }

            cTPCDisplay->WaitPrimitive();
            delete cTPCDisplay;
            delete cCalo;
            delete cFRANS;
            if(parser.getThreeViews()){
                delete cTPCDisplayViewU;
                delete cTPCDisplayViewV;
            }

        }
        
    }

    // Origins Ana Results
    std::cout<<_EfficiencyCalculator;
    TDirectory *originsAnaDirectory = anaOutputFile->mkdir("originsAnaDirectory");
    TCanvas *cOriginsAna = new TCanvas("cOriginsAna", "cOriginsAna", 0, 0, 1400,900);
    originsAnaDirectory->cd();
    _EfficiencyCalculator.DrawHistograms(cOriginsAna);
    cOriginsAna->Write();
    anaOutputFile->Write();
    anaOutputFile->Close();

    std::cout<<" Processed files "<<fileCounter<<std::endl;
    return;

}

int main(int argc, char* argv[]){

    CommandLineParser parser(argc, argv);

    std::cout << "Debug: " << parser.getDebug() << std::endl;
    std::cout << "DebugMode: " << parser.getDebugMode() << std::endl;
    std::cout << "n: " << parser.getN() << std::endl;
    std::cout << "nskip: " << parser.getNskip() << std::endl;
    std::cout << "event: " << parser.getEvent() << std::endl;
    std::cout << "sr: " << parser.getSr() << std::endl;
    std::cout << "file_name: " << parser.getFileName() << std::endl;
    std::cout << "directory_path: " << parser.getDirectoryPath() << std::endl;
    std::cout << "ext: " << parser.getExtension() << std::endl;
    std::cout << "treeName: " << parser.getTreeName() << std::endl;
    std::cout << "plotFRANS: " << parser.getPlotFRANS() << std::endl;
    std::cout << "view: " << parser.getView() << std::endl;

    // Create a ROOT application object
    TApplication *myApp = new TApplication("myApp", &argc, argv);


    RunAlgoTPCLines(parser);


    // Run the ROOT event loop
    myApp->Run();    

    return 0;
}
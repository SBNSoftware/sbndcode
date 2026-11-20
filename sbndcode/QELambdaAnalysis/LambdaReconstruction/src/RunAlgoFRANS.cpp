#include <string>
#include <iostream>

#include "TString.h"
#include <TApplication.h>
#include "TSystem.h"
#include "TImage.h"

// Objects
#include "CommandLineParser.h"
#include "SEventHandle.h"
#include "SParserConfig.h"
#include "STPCAnalyzerTreeReader.h"
#include "TPCSimpleEvents.h"
#include "TPCSimpleHits.h"
#include "TPCLinesParameters.h"
#include "TPCLinesAlgo.h"
#include "ChargeDensity.h"
#include "ChargeDensityPset.h"
#include "FRANSTree.h"


void RunAlgoFRANS(const CommandLineParser& parser)
{

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
    int vertexOption = parser.getVertexOption();

    double fFRANSScoreCut = 0.15;

    // Input tree
    std::string fTreeName = "ana/AnaTPCTree"; 

    // Set batch mode
    if(Debug==0) gROOT->SetBatch(true);
    
    // Get input files
    std::vector<TString> fFilePaths = GetInputFileList(file_name, ext, directory_path);

    // ----------- ALGORITHM PARAMETERS ---------------------------------
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
    fPsetAnaView.View = parser.getView();
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

    
    // Define the program control variables
    int fNEv = n;
    int fEv = event;
    int fSubRun = sr;
    int fNEvSkip = nskip;
    int nEvents=0;


    // output ROOT files with analysis results
    TFile* anaOutputFile = new TFile("LambdaAnaOutput.root", "RECREATE");    
    TDirectory *fransPlotsDirectory = anaOutputFile->mkdir("FRANSPlots");

    // Output directory for analysis results
    std::string tree_dirname = "framsReco/";
    if(vertexOption==1) tree_dirname = "framsTrue/";
    if(vertexOption==2) tree_dirname = "framsMine/";
    std::string tree_name = "FRANSTree";
    TTree* fTree = new TTree((tree_name).c_str(), "FRANS Output Tree");
    //FRANSTTree myTree(fTree);
    FRANSTree myTree(fTree, false);

    // Define FRANS LINES ALGORITHM
    ChargeDensity _FRANSAlgo(fPsetFRANS, fPsetAnaView.View);
    ChargeDensity _FRANSAlgoPANDORA(fPsetFRANS, fPsetAnaView.View);    

    TPCLinesAlgo _TPCLinesAlgo(fPsetAnaView, fPsetFRANS);
    // Effiency status
    EfficiencyCalculator _EfficiencyCalculator;

    // TTree loop
    int nEntries = 0;
    for (const auto& filepath : fFilePaths) {

        std::cout<<" ANALYZING THE FILE: "<<filepath<<std::endl;
        MyTPCTreeReader treeReader(filepath, fTreeName);
        

        for (int entry = 0; entry < treeReader.NEntries(); entry++) {
            std::cout<<"Entry "<<entry<<std::endl;
            // get the entry
            treeReader.GetEntry(entry);

            // check program control variables
            SEventId ev(treeReader.runID, treeReader.subrunID, treeReader.eventID);
            if (fNEv > 0 && nEvents >= fNEv) continue;
            if (treeReader.eventID != fEv && fEv != -1) continue;
            if (treeReader.subrunID != fSubRun && fSubRun != -1) continue;
            nEntries++;
            if (nEntries <= fNEvSkip) continue;
        

            nEvents++;
            std::cout << "\n\n ************** Analyzing: " << ev;
            
            // True vertex
            std::vector<double> VertexXYZ = {treeReader.nuvX, treeReader.nuvY, treeReader.nuvZ};
            std::vector<int> VertexUVYT = {treeReader.nuvU, treeReader.nuvV, treeReader.nuvC, treeReader.nuvTimeTick};
            std::cout << "   - True vertex (X, Y, Z, T) " << VertexXYZ[0] << " " << VertexXYZ[1] << " " << VertexXYZ[2] << " " << treeReader.nuvT << " (U, V, C, TT): " << VertexUVYT[0] << " " << VertexUVYT[1] << " " << VertexUVYT[2] << " " << VertexUVYT[3] << std::endl;

            // Reco vertex
            std::vector<double> RecoVertexXYZ = {treeReader.recnuvX, treeReader.recnuvY, treeReader.recnuvZ};
            std::vector<int> RecoVertexUVYT = {treeReader.recnuvU, treeReader.recnuvV, treeReader.recnuvC, treeReader.recnuvTimeTick};
            std::cout << "  - Reco vertex (X, Y, Z) " << RecoVertexXYZ[0] << " " << RecoVertexXYZ[1] << " " << RecoVertexXYZ[2] << " (U, V, C, TT): " << RecoVertexUVYT[0] << " " << RecoVertexUVYT[1] << " " << RecoVertexUVYT[2] << " "<< RecoVertexUVYT[3]  << std::endl;

            // Set TPC using the drift coordinate of the reco vertex
            //TPC = (nuvX > 0) ? 1 : 0;
            int TPC = (RecoVertexXYZ[0] > 0) ? 1 : 0;
            
            // We need a minimum number of hits to run the track finder
            size_t nhits = treeReader.hitsChannel->size();
            std::cout << "  - NHits: " << nhits << std::endl;
            if(nhits<=3){
                std::cout<<"   SKIPPED NHits \n";
                _EfficiencyCalculator.UpdateSkipped(ev);
                continue;
            }
            if( (vertexOption==0 || vertexOption==2) and treeReader.recnuvC==-1){
                std::cout<<"   SKIPPED RecoVertex \n";
                _EfficiencyCalculator.UpdateSkipped(ev);
                continue;
            }
            if(vertexOption==1 and treeReader.nuvC==-1){
                std::cout<<"   SKIPPED TrueVertex \n";
                _EfficiencyCalculator.UpdateSkipped(ev);
                continue;
            }

            // Assing the view and TPC
            int view = fPsetAnaView.View;
            if(fPsetAnaView.View==-1){ // use best view
                int bestView=_TPCLinesAlgo.GetBestView(treeReader.hitsView, treeReader.hitsChi2);
                std::cout<<"  Using best view: "<<bestView<<" TPC="<<std::to_string(TPC)<<std::endl;
                view = bestView;
            }
            
            // Set the hits
            std::vector<SHit> hitList = treeReader.GetHitsInView(view);
            
            // Set the vertex
            // true
            double vertexXTrue = VertexUVYT[2];
            if ( view == 0 ) vertexXTrue = VertexUVYT[0];
            if ( view == 1 ) vertexXTrue = VertexUVYT[1];
            double vertexYTrue = VertexUVYT[3];
            SVertex fVertexTrue(SPoint(vertexXTrue, vertexYTrue), std::to_string(view));
            
            // reco
            double vertexXReco = RecoVertexUVYT[2];
            if ( view == 0 ) vertexXReco = RecoVertexUVYT[0];
            if ( view == 1 ) vertexXReco = RecoVertexUVYT[1];
            double vertexYReco = RecoVertexUVYT[3];
            SVertex fVertexReco(SPoint(vertexXReco, vertexYReco), std::to_string(view));

            SVertex fVertexMine;
            
            SEvent recoEvent;
            if(vertexOption==2){
                // Set the hits
                bool filled = _TPCLinesAlgo.SetHitList(view, RecoVertexUVYT, VertexUVYT, hitList);

                if(filled){
                    std::cout<<" Analyzing TPCLines\n";
                    // Analyze
                    _TPCLinesAlgo.AnaView(ev.Label());
                    recoEvent = _TPCLinesAlgo.GetRecoEvent();
                    fVertexMine = SVertex( SPoint((double)_TPCLinesAlgo.GetMainVertex().X()+_TPCLinesAlgo.ShiftX(),
                                                (double) _TPCLinesAlgo.GetMainVertex().Y()+_TPCLinesAlgo.ShiftY())
                                            , "");
                }
                else{
                    fVertexMine = fVertexReco;
                }
            }

            // set it
            SVertex fVertex = fVertexReco;
            if(vertexOption==1) fVertex = fVertexTrue;
            else if(vertexOption==2) fVertex = fVertexMine;

            
            _FRANSAlgoPANDORA.Fill(hitList, fVertexReco);
            _FRANSAlgo.Fill(hitList, fVertex);


            _FRANSAlgo.Score()>=fFRANSScoreCut ? _EfficiencyCalculator.UpdateSelected(ev):_EfficiencyCalculator.UpdateNotSelected(ev);
            std::string outputLabel = (_FRANSAlgo.Score()>=fFRANSScoreCut)? "plot_Accepted":"plot_Rejected";
            std::string outputLabel2 = (_FRANSAlgoPANDORA.Score()>=fFRANSScoreCut)? "Accepted":"Rejected";

                        
            if(Debug>=0){
                TCanvas *cDisplay = new TCanvas( (outputLabel+"_"+ev.Label()+"_vw"+std::to_string(view)).c_str(), outputLabel.c_str(), 600, 0, 800, 1200);
                _FRANSAlgo.Display(cDisplay);

                // Save TCanvas and pdf
                TImage *img = TImage::Create();
                img->FromPad(cDisplay);
                fransPlotsDirectory->cd();
                img->Write( ("image/"+outputLabel+"_"+ev.Label()+"_vw"+std::to_string(view)+".pdf").c_str() );
                cDisplay->Write();

                // Get TPads in the canvas
                TPad *pad1A = (TPad*)cDisplay->GetPrimitive("pad1A");
                pad1A->SaveAs("pad1A.pdf");
                TPad *pad1B = (TPad*)cDisplay->GetPrimitive("pad1B");
                pad1B->SaveAs("pad1B.pdf");
                TPad *pad2A = (TPad*)cDisplay->GetPrimitive("pad2A");
                pad2A->SaveAs("pad2A.pdf");
                TPad *pad2B = (TPad*)cDisplay->GetPrimitive("pad2B");
                pad2B->SaveAs("pad2B.pdf");
                TPad *pad3A = (TPad*)cDisplay->GetPrimitive("pad3A");
                pad3A->SaveAs("pad3A.pdf");
                TPad *pad3B = (TPad*)cDisplay->GetPrimitive("pad3B");
                pad3B->SaveAs("pad3B.pdf");

                
                cDisplay->WaitPrimitive();

                delete cDisplay;
                delete img;
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
        
    // Write the FRANS TTree
    anaOutputFile->cd();
    fTree->Write();
    delete fTree;

    anaOutputFile->Write();
    anaOutputFile->Close();

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

    RunAlgoFRANS(parser);


    // Run the ROOT event loop
    myApp->Run();    

    return 0;

}
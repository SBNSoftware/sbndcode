#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>

//---------  Load function
void PlotCalorimetryDisplay(){
    std::cout<<"LambdaAnalysis function loaded"<<std::endl;
    return;
}

//---------  Color list
std::vector<int> fColorList = { kAzure+7, kRed+2,  kAzure-5, kGreen+3, kOrange-3, kMagenta+1, kCyan- 3, kYellow+2, kViolet-1, kTeal-1};
int fColorLI = kAzure+7;
int fColorHI = kRed+2;

//-------- Draw all space points
bool fDrawAllSpacePoints = false;

//--------- Define a function to add points to the graph and create a TGraph2D
void addPointsAndCreateGraph(TGraph2D *graph, std::vector<double> x, std::vector<double> y, std::vector<double> z){

    for(size_t k = 0; k < x.size(); k++) {
        graph->SetPoint(k, x.at(k), y.at(k), z.at(k));
    }

    return;
}


// Define a function to draw the graph
void drawGraph(TGraph2D *graph, const char* title, Color_t color){
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(20);
    graph->SetTitle(title);
    graph->Draw("P same");
}


// Function to get minimum and maximum XYZ values of the vector of vectors
auto getMinMax = [](std::vector<double> v){
    double min = *std::min_element(v.begin(), v.end());
    double max = *std::max_element(v.begin(), v.end());
    return std::make_pair(min, max);
};

//---------  Main function
void RunCalorimetryDisplay(std::string fInputFileName="", int event=-1, std::string fTreeDirName = "lambdaPidPandoraAna/", std::string fTreeName = "CalorimetryTree")
{

    //--------- Read dEdx templates, path relative to LAMBDAANA_SRC environment variable
    std::string fTemplatesDir = std::getenv("LAMBDAANA_SRC");
    fTemplatesDir += "/Analysis/dEdxrestemplates.root";
    TFile *fFileTemplates = new TFile(fTemplatesDir.c_str(),"READ");
    TProfile *fProfiledEdxProton = (TProfile*)fFileTemplates->Get("dedx_range_pro");
    TProfile *fProfiledEMuon = (TProfile*)fFileTemplates->Get("dedx_range_mu");


    //--------- Input TTree
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());
    fTree->Print(); 

    //--------- Binning
    Int_t nBins = 200;
    Double_t xMin = 0;
    Double_t xMax = 20;
    Double_t yMin = 0;
    Double_t yMax = 40;
    Double_t yMaxCharge = 10000;;
    // Stats to 0
    gStyle->SetOptStat(0);
    
    //--------- Set Branches 
    Int_t RunID;
    Int_t SubRunID;
    Int_t EventID;
    fTree->SetBranchAddress("RunID", &RunID);
    fTree->SetBranchAddress("SubRunID", &SubRunID);
    fTree->SetBranchAddress("EventID", &EventID);
    Bool_t RecoStatus;
    fTree->SetBranchAddress("RecoStatus", &RecoStatus);
    Double_t InvariantMass, KELI, KEHI, OpeningAngle, LambdaEnergy, LambdaKE;
    fTree->SetBranchAddress("InvariantMass", &InvariantMass);
    fTree->SetBranchAddress("KELI", &KELI);
    fTree->SetBranchAddress("KEHI", &KEHI);
    fTree->SetBranchAddress("OpeningAngle", &OpeningAngle);
    fTree->SetBranchAddress("LambdaEnergy", &LambdaEnergy);
    fTree->SetBranchAddress("LambdaKE", &LambdaKE);

    // True variables
    Double_t TrueProtonKE, TruePionKE, TrueOpeningAngle, TrueLambdaE, TrueLambdaKE;
    fTree->SetBranchAddress("TrueProtonKE", &TrueProtonKE);
    fTree->SetBranchAddress("TruePionKE", &TruePionKE);
    fTree->SetBranchAddress("TrueOpeningAngle", &TrueOpeningAngle);
    fTree->SetBranchAddress("TrueLambdaE", &TrueLambdaE);
    fTree->SetBranchAddress("TrueLambdaKE", &TrueLambdaKE);

    // Space points
    std::vector<std::vector<double>> *SpacePointsX = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>> *SpacePointsY = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>> *SpacePointsZ = new std::vector<std::vector<double>>;
    fTree->SetBranchAddress("SpacePointsX", &SpacePointsX);
    fTree->SetBranchAddress("SpacePointsY", &SpacePointsY);
    fTree->SetBranchAddress("SpacePointsZ", &SpacePointsZ);

    // Space points (all)
    std::vector<double> *AllSpacePointsX = new std::vector<double>;
    std::vector<double> *AllSpacePointsY = new std::vector<double>;
    std::vector<double> *AllSpacePointsZ = new std::vector<double>;
    fTree->SetBranchAddress("AllSpacePointsX", &AllSpacePointsX);
    fTree->SetBranchAddress("AllSpacePointsY", &AllSpacePointsY);
    fTree->SetBranchAddress("AllSpacePointsZ", &AllSpacePointsZ);
    
    // Residual range LI
    std::vector<double> *ResidualRangeLI = new std::vector<double>;
    std::vector<double> *DepositedEnergyLI = new std::vector<double>;
    std::vector<double> *DepositedChargeLI = new std::vector<double>;
    fTree->SetBranchAddress("ResidualRangeLI", &ResidualRangeLI);
    fTree->SetBranchAddress("DepositedEnergyLI", &DepositedEnergyLI);
    fTree->SetBranchAddress("DepositedChargeLI", &DepositedChargeLI);
    // Residual range HI
    std::vector<double> *ResidualRangeHI = new std::vector<double>;
    std::vector<double> *DepositedEnergyHI = new std::vector<double>;
    std::vector<double> *DepositedChargeHI = new std::vector<double>;
    fTree->SetBranchAddress("ResidualRangeHI", &ResidualRangeHI);
    fTree->SetBranchAddress("DepositedEnergyHI", &DepositedEnergyHI);
    fTree->SetBranchAddress("DepositedChargeHI", &DepositedChargeHI);

    
    // deDx (all tracks)
    std::vector<std::vector<double>> *ResidualRange = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>> *DepositedEnergy = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>> *DepositedCharge = new std::vector<std::vector<double>>;
    std::vector<int> *BestPIDAll = new std::vector<int>;
    fTree->SetBranchAddress("ResidualRangeAll", &ResidualRange);
    fTree->SetBranchAddress("DepositedEnergyAll", &DepositedEnergy);
    fTree->SetBranchAddress("DepositedChargeAll", &DepositedCharge);
    fTree->SetBranchAddress("BestPIDAll", &BestPIDAll);


    //--------- Tree loop
    for(int i=0; i<fTree->GetEntries(); i++){
        
        fTree->GetEntry(i);
        std::cout<<i<<" RunID: "<<RunID<<", SubRunID: "<<SubRunID<<", EventID: "<<EventID<<std::endl;

        if( event>0 && EventID!=event) continue;
        if( event==-2 && RecoStatus!=1) continue;

        std::cout<<"InvariantMass: "<<InvariantMass<<std::endl;
        std::cout<<"KELI: "<<KELI<<" KEPion: "<<1000*TruePionKE<<std::endl;
        std::cout<<"KEHI: "<<KEHI<<" KEProton: "<<1000*TrueProtonKE<<std::endl;
        std::cout<<"OpeningAngle: "<<OpeningAngle<<" TrueOpeningAngle: "<<180*TrueOpeningAngle/TMath::Pi()<<std::endl;
        std::cout<<"TrueLambdaE: "<<1000*TrueLambdaE<<" TrueLambdaKE: "<<1000*TrueLambdaKE<<std::endl;
        std::cout<<"LambdaE: "<<LambdaEnergy<<" LambdaKE: "<<LambdaKE<<std::endl;



        //-------- Canvas
        std::string canvasTitle = "RunID: "+std::to_string(RunID)+", SubRunID: "+std::to_string(SubRunID)+", EventID: "+std::to_string(EventID);
        TCanvas *c1 = new TCanvas("c1",canvasTitle.c_str(),1200,650);
        gStyle->SetPadBottomMargin(0.2);
        gStyle->SetPadLeftMargin(0.15);
        gStyle->SetTitleYOffset(1.2);
        TPad *pad1 = new TPad("pad1","This is pad1",0.02,0.5,0.31,1);
        TPad *pad2 = new TPad("pad2","This is pad2",0.35,0.5,0.64,1);
        TPad *pad3 = new TPad("pad3","This is pad3",0.68,0.5,0.98,1);
        TPad *pad4 = new TPad("pad4","This is pad4",0.02,0.0,0.31,0.5);
        TPad *pad5 = new TPad("pad5","This is pad5",0.35,0.0,0.64,0.5);
        TPad *pad6 = new TPad("pad6","This is pad6",0.68,0.0,0.98,0.5);
        pad1->Draw();
        pad2->Draw();
        pad3->Draw();
        pad4->Draw();
        pad5->Draw();
        pad6->Draw();


        //-------- Plot dEdx from TTree
        std::cout<<" ResidualRangeLI Size: "<<ResidualRangeLI->size()<<std::endl;
        if(ResidualRangeLI->size()>0 && ResidualRangeHI->size()>0){
            pad1->cd();
            TH2F *h_dEdxLI = new TH2F("h_dEdxLI","h_dEdxLI",nBins,xMin,xMax, nBins,yMin,yMax);
            fTree->Draw("DepositedEnergyLI:ResidualRangeLI>>h_dEdxLI","","colz", 1, i);
            TH2F *h_dEdxHI = new TH2F("h_dEdxHI","h_dEdxHI",nBins,xMin,xMax, nBins,yMin,yMax);
            fTree->Draw("DepositedEnergyHI:ResidualRangeHI>>h_dEdxHI","","colz", 1, i);

            // Colors and markers
            h_dEdxLI->SetMarkerColor(fColorLI);
            h_dEdxHI->SetMarkerColor(fColorHI);
            h_dEdxLI->SetMarkerStyle(20);
            h_dEdxHI->SetMarkerStyle(20);
            
            //Draw histograms
            h_dEdxLI->GetYaxis()->SetTitle("dE/dx [MeV]");
            h_dEdxLI->GetXaxis()->SetTitle("Residual Range [cm]");
            h_dEdxLI->Draw("");
            h_dEdxHI->Draw("same");
            
            // Legend
            TLegend *leg = new TLegend(0.4,0.65,0.89,0.89);
            leg->AddEntry(h_dEdxLI,"Low Ionization","p");
            leg->AddEntry(h_dEdxHI,"High Ionization","p");
            leg->Draw("same");
            pad1->Update();

            //-------- Plot dQdx from TTree
            pad2->cd();
            TH2F *h_dQdxLI = new TH2F("h_dQdxLI","h_dQdxLI",nBins,xMin,xMax, nBins,yMin,yMaxCharge);
            fTree->Draw("DepositedChargeLI:ResidualRangeLI>>h_dQdxLI","","colz", 1, i);
            TH2F *h_dQdxHI = new TH2F("h_dQdxHI","h_dQdxHI",nBins,xMin,xMax, nBins,yMin,yMaxCharge);
            fTree->Draw("DepositedChargeHI:ResidualRangeHI>>h_dQdxHI","","colz", 1, i);

            // Colors and markers
            h_dQdxLI->SetMarkerColor(fColorLI);
            h_dQdxHI->SetMarkerColor(fColorHI);
            h_dQdxLI->SetMarkerStyle(20);
            h_dQdxHI->SetMarkerStyle(20);

            //Draw histograms
            h_dQdxLI->GetYaxis()->SetTitle("dQ/dx [ADC/cm]");
            h_dQdxLI->GetXaxis()->SetTitle("Residual range [cm]");
            h_dQdxLI->Draw("");
            h_dQdxHI->Draw("same");

            // Legend
            leg->Draw("same");

        }

        //-------- Plot dEdx (all tracks)
        pad3->cd();
        std::vector<TH1F*> h_dEdxAll;
        THStack *hs_dEdx = new THStack("hs_dEdx",";Residual range [cm];dE/dx [MeV]");
        for(int j=0; j<ResidualRange->size(); j++){
            std::string histName = "h_dEdx_"+std::to_string(j);
            TH1F *h_dEdx = new TH1F(histName.c_str(),";Residual range [cm];dE/dx [MeV];", nBins,xMin,xMax);
            for(int k=0; k<ResidualRange->at(j).size(); k++){
                h_dEdx->Fill(ResidualRange->at(j).at(k), DepositedEnergy->at(j).at(k));
            }
            h_dEdx->SetMarkerColor(fColorList.at(j));
            h_dEdx->SetMarkerStyle(20);
            h_dEdx->SetLineColor(fColorList.at(j));
            hs_dEdx->Add(h_dEdx);
            h_dEdxAll.push_back(h_dEdx);
        }

        // TH2 Frame with minimum and maximum values
        TH2F *hFrame2 = new TH2F("hFrame",";Residual range [cm]; dE/dx [MeV]", nBins,xMin,xMax, nBins,yMin,yMax);
        hFrame2->Draw();
        hs_dEdx->Draw("nostack same  PLC HIST P");
        
    
        // Draw templates for proton and muon
        fProfiledEdxProton->SetLineColor(kOrange);
        fProfiledEdxProton->SetFillColorAlpha(kOrange, 0.5);
        fProfiledEMuon->SetLineColor(kAzure);
        fProfiledEMuon->SetFillColorAlpha(kAzure, 0.5);
        fProfiledEdxProton->Draw("E4 same");
        fProfiledEMuon->Draw("E4 same");

        TLegend *leg = new TLegend(0.4,0.65,0.89,0.89);
        for(int j=0; j<h_dEdxAll.size(); j++){
            std::string pidLabel = "Track "+std::to_string(j)+" (best PDG="+std::to_string(BestPIDAll->at(j)) + ")";
            leg->AddEntry(h_dEdxAll.at(j), pidLabel.c_str(),"p");
        }
        leg->AddEntry(fProfiledEdxProton, "Proton template","l");
        leg->AddEntry(fProfiledEMuon, "Pion template","l");
        leg->Draw("same");
        leg->SetTextSize(0.04);
        



        //-------- Draw SpacePoints
        pad4->cd();

        std::vector<TGraph2D*> graphList;
        std::vector<int> graphColorList;

        // Global min and max values
        std::pair<double, double> xMinMaxGlobal = {200, -200};
        std::pair<double, double> yMinMaxGlobal = {200, -200};
        std::pair<double, double> zMinMaxGlobal = {500, 0};
        
        // Draw AllSpacePoints
        if(AllSpacePointsX->size()>0 && fDrawAllSpacePoints){
            TGraph2D *graphAll = new TGraph2D();
            addPointsAndCreateGraph(graphAll, *AllSpacePointsZ, *AllSpacePointsX, *AllSpacePointsY);
            graphAll->SetMarkerColorAlpha(kBlack, 0.5);
            graphAll->SetMarkerStyle(20);
            graphAll->Draw("P");
            graphList.push_back(graphAll);
            graphColorList.push_back(kBlack);

            // Get min and max values of X, Y and Z
            std::pair<double, double> xMinMax = getMinMax(*AllSpacePointsX);
            std::pair<double, double> yMinMax = getMinMax(*AllSpacePointsY);
            std::pair<double, double> zMinMax = getMinMax(*AllSpacePointsZ);
            // Update global min and max values
            xMinMaxGlobal = std::make_pair(std::min(xMinMax.first, xMinMaxGlobal.first), std::max(xMinMax.second, xMinMaxGlobal.second));
            yMinMaxGlobal = std::make_pair(std::min(yMinMax.first, yMinMaxGlobal.first), std::max(yMinMax.second, yMinMaxGlobal.second));
            zMinMaxGlobal = std::make_pair(std::min(zMinMax.first, zMinMaxGlobal.first), std::max(zMinMax.second, zMinMaxGlobal.second));
        }

        // Draw Calo SpacePoints
        std::cout<<" SpacePoints Size: "<<SpacePointsX->size()<<std::endl;
        if(SpacePointsX->size()>0){
            // Graphs
            int graphCounter = 0;
            for(int j=0; j<SpacePointsX->size(); j++){
                std::cout<<"SpacePointsX size: "<<SpacePointsX->at(j).size()<<std::endl;
                if(SpacePointsX->at(j).size()==0) continue;
                
                // Get min and max values of X, Y and Z
                std::pair<double, double> xMinMax = getMinMax(SpacePointsX->at(j));
                std::pair<double, double> yMinMax = getMinMax(SpacePointsY->at(j));
                std::pair<double, double> zMinMax = getMinMax(SpacePointsZ->at(j));
                // Update global min and max values
                xMinMaxGlobal = std::make_pair(std::min(xMinMax.first, xMinMaxGlobal.first), std::max(xMinMax.second, xMinMaxGlobal.second));
                yMinMaxGlobal = std::make_pair(std::min(yMinMax.first, yMinMaxGlobal.first), std::max(yMinMax.second, yMinMaxGlobal.second));
                zMinMaxGlobal = std::make_pair(std::min(zMinMax.first, zMinMaxGlobal.first), std::max(zMinMax.second, zMinMaxGlobal.second));
            
                TGraph2D *graph = new TGraph2D();
                addPointsAndCreateGraph(graph, SpacePointsZ->at(j), SpacePointsX->at(j), SpacePointsY->at(j));    
                graphList.push_back(graph);
                graphColorList.push_back(fColorList.at(j));            
            }
        }


        TLegend *legSP = new TLegend(0.75,0.8,0.95,0.95);
        TH3F *hFrame = new TH3F("frame",";Z [cm]; X [cm]; Y [cm]", nBins,zMinMaxGlobal.first,zMinMaxGlobal.second, nBins, xMinMaxGlobal.first,xMinMaxGlobal.second, nBins,yMinMaxGlobal.first,yMinMaxGlobal.second);
        hFrame->Draw();
        for(int j=0; j<graphList.size(); j++){
            drawGraph(graphList.at(j), "SpacePoints", graphColorList.at(j));
            legSP->AddEntry(graphList[j], ("Track "+std::to_string(j)).c_str(), "p");
        }
        // Title offsets
        hFrame->GetXaxis()->SetTitleOffset(1.3);
        hFrame->GetYaxis()->SetTitleOffset(1.3);
        hFrame->GetZaxis()->SetTitleOffset(1.3);
        // Left margin
        gPad->SetLeftMargin(0.15);
        // NDivisions
        hFrame->GetXaxis()->SetNdivisions(505);
        hFrame->GetYaxis()->SetNdivisions(505);
        hFrame->GetZaxis()->SetNdivisions(505);
        // Draw legend
        legSP->Draw("same");

        pad4->Update();

        
        //-------- Draw ZX projections
        pad5->cd();
        TH2F *hZX = new TH2F("hZX",";Z [cm]; X [cm]", nBins,zMinMaxGlobal.first,zMinMaxGlobal.second, nBins,xMinMaxGlobal.first,xMinMaxGlobal.second);
        hZX->Draw();
        if(AllSpacePointsX->size()>0){
            TGraph *graphZX = new TGraph(AllSpacePointsZ->size(), &AllSpacePointsZ->at(0), &AllSpacePointsX->at(0));
            graphZX->SetMarkerColor(kBlack);
            graphZX->SetMarkerStyle(20);
            graphZX->Draw("P same");
        }

        //-------- Draw ZY projections
        pad6->cd();
        TH2F *hZY = new TH2F("hZY",";Z [cm]; Y [cm]", nBins,zMinMaxGlobal.first,zMinMaxGlobal.second, nBins,yMinMaxGlobal.first,yMinMaxGlobal.second);
        hZY->Draw();
        if(AllSpacePointsZ->size()>0){
            TGraph *graphZY = new TGraph(AllSpacePointsZ->size(), &AllSpacePointsZ->at(0), &AllSpacePointsY->at(0));
            graphZY->SetMarkerColor(kBlack);
            graphZY->SetMarkerStyle(20);
            graphZY->Draw("P same");
        }




        c1->cd();
        c1->Update();
        c1->WaitPrimitive();
        
    }
    

}



//---------  Main function
void RunPlotAveragedEdx(std::string fInputFileName="", std::string fTreeDirName = "lambdaPidPandoraAna/", std::string fTreeName = "CalorimetryTree")
{

    //--------- Input TTree
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());

    //--------- Binning
    Int_t nBinsX = 100;
    Double_t xMin = 0;
    Double_t xMax = 20;
    Int_t nBinsY = 50;
    Double_t yMin = 0;
    Double_t yMax = 40;

    // TH2F
    TCanvas *c1 = new TCanvas("c1","c1",1200,650);
    TH2F *h_dEdxLI = new TH2F("h_dEdxLI",";Residual range [cm]; dE/dx [MeV]",nBinsX,xMin,xMax, nBinsY,yMin,yMax);
    TH2F *h_dEdxHI = new TH2F("h_dEdxHI",";Residual range [cm]; dE/dx [MeV]",nBinsX,xMin,xMax, nBinsY,yMin,yMax);

    // Draw
    fTree->Draw("DepositedEnergyLI:ResidualRangeLI>>h_dEdxLI","","colz");
    fTree->Draw("DepositedEnergyHI:ResidualRangeHI>>h_dEdxHI","","colz");

    c1->cd();
    h_dEdxLI->Draw("colz");
    h_dEdxHI->Draw("colz same");
    c1->Update();
    c1->WaitPrimitive();

}

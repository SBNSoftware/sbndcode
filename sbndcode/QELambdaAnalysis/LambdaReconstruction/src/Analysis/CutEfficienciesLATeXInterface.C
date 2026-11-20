#ifndef CUTEFFICIENCIES_LATEXINTERFACE_H
#define CUTEFFICIENCIES_LATEXINTERFACE_H

#define WriteRelativeEfficiencies 0


//--------- Read the POT normalization
void ReadPOT(TFile *fFile, double totalPOTNorm, double& potScaling, double& potScalingSignal, std::string fTreeName = "originsAna/pottree"){
    
    // Read TreePOT
    std::cout<<"Reading POT tree "<<fTreeName<<std::endl;
    TTree *fTreePOT = (TTree *)fFile->Get( fTreeName.c_str() );
    if (!fTreePOT) {
        std::cerr << "Failed to read the POT tree." << std::endl;
        return;
    }
    
    // -------- Get the accumulated POT
    double pot = 0;
    double averageintmode;
    bool inclusive;
    int nnuinteractions;
    fTreePOT->SetBranchAddress("pot",&pot);
    fTreePOT->SetBranchAddress("averageintmode",&averageintmode);
    fTreePOT->SetBranchAddress("inclusive",&inclusive);
    fTreePOT->SetBranchAddress("nnuinteractions",&nnuinteractions);
    double totalPOT = 0;
    double totalPOTSignal = 0;
    for (int i = 0; i < fTreePOT->GetEntries(); ++i) {
        fTreePOT->GetEntry(i);
        if(averageintmode==0 && inclusive==0) totalPOTSignal+=pot;
        else if(pot<1e17) totalPOT+=pot;
    }
    
    potScaling = totalPOTNorm/totalPOT;
    potScalingSignal = totalPOTNorm/totalPOTSignal;
    std::cout<<" Accumulated POT: "<<totalPOT<<std::endl;
    std::cout<<" Accumulated POT Signal: "<<totalPOTSignal<<std::endl;
    std::cout<<" POT scaling factors: "<<potScaling<<" "<<potScalingSignal<<std::endl;

    return;   
}


// --------- Create output list with events that pass the cuts
void CreateHandScanList(TTree *fTree, TTree *fTreeHeader, TCut cut, std::vector<SampleDef> sampleDefs, bool bgOnly = false){

    //--------- Loop over the TreeHeader
    unsigned int fSubRunId;
    fTreeHeader->SetBranchAddress("SubRunID", &fSubRunId);
    unsigned int fRunId;
    fTreeHeader->SetBranchAddress("RunID", &fRunId);
    std::string *fLArInputFileName = new std::string;
    fTreeHeader->SetBranchAddress("InputFileName", &fLArInputFileName);
    // Create map of run-subrun filename
    std::map<std::string, std::vector<std::string>> subrunFilenameMap;

    std::map<std::string, int> fileCounterMap;

    for(size_t i=0; i<fTreeHeader->GetEntries(); ++i){
        fTreeHeader->GetEntry(i);
        
        // check the string includes the substring "V0Lambda"
        //if(fLArInputFileName->find("V0Lambda") != std::string::npos || fLArInputFileName->find("V0Overlay") != std::string::npos) continue;
        

        if(bgOnly && (fLArInputFileName->find("Inclusive") == std::string::npos
                        && fLArInputFileName->find("cosmic_rockbox_sbnd") == std::string::npos)
                    ) continue;
        

        
        std::string runSubrunLabel = std::to_string(fRunId) + ":" + std::to_string(fSubRunId);

        // if the subrun is not in the map, add it
        if(subrunFilenameMap.find(runSubrunLabel) == subrunFilenameMap.end()){
            subrunFilenameMap[runSubrunLabel] = std::vector<std::string>();
        }

        // add the filename to the vector
        subrunFilenameMap[runSubrunLabel].push_back(*fLArInputFileName);


        // if the file is not in the map, add it
        if(fileCounterMap.find(*fLArInputFileName) == fileCounterMap.end()){
            fileCounterMap[*fLArInputFileName] = 1;
        }
        else{
            fileCounterMap[*fLArInputFileName]++;
        }
          
    }

    // print counter map
    for(auto& fileCounter : fileCounterMap){
        std::cout << "File: " << fileCounter.first << " Count: " << fileCounter.second << std::endl;
    }
    // cout map size
    std::cout << "Map size: " << fileCounterMap.size() << std::endl;

    // Output 
    std::ofstream handScanEvents;
    handScanEvents.open("OutputPlots/HandScanEvents.txt");

    // loop over signal types
    for (size_t k = 0; k < sampleDefs.size(); k++) {
        // Skip if it is the signal
        if (sampleDefs[k].IsSignal()) continue;

        std::cout << " Sample: " << sampleDefs[k].GetLabel() << std::endl;

        // Create the histogram
        TH3D* h3_all = new TH3D("h3_all", "3D Histogram", 2, 1, 3, 10000, 1, 10001, 2100, 1, 2101);

        // Create the cut
        TCut debugCut = cut && (TCut(sampleDefs[k].GetVar()));

        // Draw the histogram
        // Z axis is the run number, Y axis is event ID, and X axis is subrun ID
        fTree->Draw("EventID:SubrunID:RunID>>h3_all", debugCut);

        // Get the number of bins in the histogram
        int numBinsX = h3_all->GetNbinsX();
        int numBinsY = h3_all->GetNbinsY();
        int numBinsZ = h3_all->GetNbinsZ();

        for (int run = 1; run <= numBinsX; run++) {
            for (int sr = 1; sr <= numBinsY; sr++) {
                for (int e = 1; e <= numBinsZ; e++) {
                    double binValue = h3_all->GetBinContent(run, sr, e);
                    if (binValue > 0) {
                        std::cout << "RunNumber: " << run << " SubRunID: " << sr << " EventID: " << e << " BinValue: " << binValue << std::endl;
                        std::string runSubrunLabel = std::to_string(run) + ":" + std::to_string(sr);
                        for (size_t i = 0; i < subrunFilenameMap[runSubrunLabel].size(); i++){
                            std::cout << "   Filename: " << subrunFilenameMap[runSubrunLabel][i] << std::endl;
                            handScanEvents << run << ":" << sr << ":" << e << " " << subrunFilenameMap[runSubrunLabel][i] << std::endl;
                        }
                    }
                }
            }
        }
    }


    handScanEvents.close();

    return;

}

//--------- Generate and compile the LaTeX table
void GenerateAndCompileTeXTable(
    const std::vector<PlotDef> plotDefs,
    const std::vector<SampleDef> sampleDefs,
    const std::vector<std::vector<int>> matrixCounter,
    const std::string fileName,
    const std::string tableCaption,
    const std::string dirLabel,
    double potNorm=1,
    double potNormSignal=1,
    double totalPOTNorm=1,
    int denominatorIndex=0
) {
    //--- Open a file stream for writing the .tex file
    std::ofstream texFile(fileName+".tex");

    //--- Check if the file is open
    if (!texFile.is_open()) {
        std::cerr << "Failed to open the output file." << std::endl;
        return;
    }

    //--- Denominator weights
    std::vector<int> counts0 = matrixCounter[denominatorIndex];

    //--- POT scaling map weights
    std::vector<double> potScalingMap;
    for(auto& sampleName : sampleDefs){
        if(sampleName.GetLabelS()=="Signal")
            potScalingMap.push_back( potNormSignal );
        else
            potScalingMap.push_back( potNorm );
    }

    //--- Write the LaTeX document preamble
    texFile << "\\documentclass{article}" << std::endl;
    texFile << "\\usepackage{graphicx}" << std::endl;
    texFile << "\\usepackage{pdflscape}"<<std::endl;
    texFile << "\\begin{document}" << std::endl;
    texFile << "\\begin{landscape}" << std::endl;


    //--- Write the LaTeX code for the table
    texFile << "\\begin{table}[h]" << std::endl;
    texFile << "\\centering" << std::endl;
    texFile << "\\begin{tabular}{|c|";
    for(size_t k = 0; k<sampleDefs.size(); k++)
        texFile<<"c|";
    texFile << "}" << std::endl;
    texFile << "\\hline" << std::endl;

    

    //--- Write the header row
    texFile << " Cut & ";
    int cont=0;  
    for(auto& sampleName : sampleDefs){
        if( cont == sampleDefs.size()-1 )
            texFile << sampleName.GetLatexLabel();
        else
            texFile << sampleName.GetLatexLabel() << " & ";
        
        cont++;
    }
    texFile << "\\\\ \\hline" << std::endl;


    // Write the data rows
    // last stored index
    size_t lastStoredIndex = 0;
    for (size_t i = 0; i < matrixCounter.size(); ++i) {

        texFile << "$ {\\rm " << plotDefs[i].GetCutLabelS() << "}$" << " & ";

        std::vector<int> counts = matrixCounter[i];
        std::vector<int> countsPre = matrixCounter[lastStoredIndex];

        cont=0;
        for (size_t j = 0; j < matrixCounter[i].size(); ++j) {

            double potScaling = potScalingMap[j];
            std::string sampleName = sampleDefs[j].GetLabelS();
            int counts = matrixCounter[i][j];
            
            std::ostringstream streamObjEff;
            streamObjEff << " (" << std::fixed << std::setprecision(3) << 100.*counts/counts0[j];
            if(WriteRelativeEfficiencies==1){
                streamObjEff << " \\%) \\ $\\epsilon_{r}$= " << std::setprecision(0) << 100.*counts/countsPre[j] << "\\%";
            }
            else{
                streamObjEff << " \\%)";
            }

            // Get string from out
            if(  cont == sampleDefs.size()-1 )
                texFile << (int)(potScaling*counts) << streamObjEff.str();
            else
                texFile << (int)(potScaling*counts) << streamObjEff.str() << " & ";
            cont++;
        }

        texFile << "\\\\ \\hline" << std::endl;

        lastStoredIndex = i;
    }
    


    /*
    texFile << "$ {\\rm " << "Final eff" << "}$" << " & ";
    
    cont=0;
    std::map<std::string, int> histogramCounts = anaPlots[ lastStoredIndex].GetCountsV();
    std::map<std::string, int> initialCounts;
    for(auto& sample:sampleDefs){
        initialCounts[sample.GetLabelS()] = sample.GetNEvents();
    }

    for(auto& sampleName : histogramCounts){

        double potScaling = potScalingMap[sampleName.first];

        std::ostringstream streamObjEff;
        streamObjEff <<  std::fixed << std::setprecision(2) << 100.*sampleName.second/initialCounts[sampleName.first];
        // Get string from out
        if(  cont == sampleDefs.size()-1 )
            texFile << potScaling*sampleName.second <<" ("<< streamObjEff.str() <<" \\%)";
        else
            texFile << potScaling*sampleName.second <<" ("<< streamObjEff.str() <<" \\%)"  << " & ";
        cont++;
    }


    texFile << "\\\\ \\hline" << std::endl;
    */

    // --- Write the caption and end the table    
    std::string newTableCaption = tableCaption;
    if(totalPOTNorm!=1){
        // total pot in scientific notation
        std::ostringstream streamObj;
        streamObj <<  std::scientific << std::setprecision(2) << totalPOTNorm;
        // Get string from out
        std::string totalPOTNormStr = streamObj.str();
        newTableCaption += " (Normalized to "+totalPOTNormStr+" POT)";
    }
    texFile << "\\end{tabular}" << std::endl;
    texFile << "\\caption{" << newTableCaption << "}" << std::endl;
    texFile << "\\end{table}" << std::endl;
    texFile << "\\end{landscape}" << std::endl;


    //--- Close the LaTeX document
    texFile << "\\end{document}" << std::endl;

    //--- Close the file
    texFile.close();

    //--- Compile the .tex file
    std::string compileCommand = "pdflatex " + fileName;
    int compileResult = std::system(compileCommand.c_str());
    if (compileResult != 0) {
        std::cerr << "Failed to compile the .tex file." << std::endl;
    }
    gSystem->Exec( ("mv "+fileName+".pdf "+dirLabel).c_str() );
    gSystem->Exec( ("rm "+fileName+".{aux,log,tex,pdf}").c_str() );

    return;

}

void SetCustomYAxisLabelsLog(TCanvas *c, double minBinContentDecade, double yLabelXPos=0.5){
    c->cd();
    int nDecades = static_cast<int>(floor(log10(100 / minBinContentDecade)));

    
    for (int j = 0; j <= nDecades; ++j) {
        double decade = minBinContentDecade * std::pow(10, j);
        int exponent = static_cast<int>(log10(decade));
        std::string label = "10^{" + std::to_string(exponent) + "}";
        TLatex *latex = new TLatex(yLabelXPos+0.25, decade, label.c_str());
        latex->SetTextSize(0.04);
        latex->SetTextAlign(12);
        latex->SetTextAngle(90);
        latex->Draw();
    }

}


void SetCustomYAxisLabels(TCanvas *c, int nLabels, double yLabelXPos=0.5){
    c->cd();
   
    double labelStep = 100./(nLabels-1);
    
    for (int j = 0; j <= nLabels; ++j) {
        
        // label to integer
        int labelInt = (int)j*labelStep;
        std::string label = std::to_string(labelInt);
        TLatex *latex = new TLatex(yLabelXPos+0.25, j*labelStep, label.c_str());
        latex->SetTextSize(0.04);
        latex->SetTextAlign(12);
        latex->SetTextAngle(90);
        latex->Draw();
    }

}

//--------- Make the cut flow plot
void MakeCutFlowPlot(const std::vector<PlotDef> plotDefs,
                    const std::vector<SampleDef> sampleDefs,
                    const std::vector<std::vector<int>> matrixCounter,
                    const std::string dirLabel,
                    double potNorm=1,
                    double potNormSignal=1,
                    double totalPOTNorm=1
)
{

    //--- POT scaling map weights
    std::vector<double> potScalingMap;
    for(auto& sampleName : sampleDefs){
        if(sampleName.GetLabelS()=="Signal")
            potScalingMap.push_back( potNormSignal );
        else
            potScalingMap.push_back( potNorm );
    }

    //--- Styler
    CutStyler *fStyler(new CutStyler(0));
    double yLabelXPos = plotDefs.size();
    double fLabelTextSize=0.04;
    double fEffLabelSize=0.04;
    
    // Vector of TH1F, one per sample, for efficiency
    std::vector<TH1F*> hEfficiencyV;
    // Vector of TH1F, one per sample, for fraction of events
    std::vector<TH1F*> hFractionV;
    int nCutBins = plotDefs.size();
    for (size_t i = 0; i < sampleDefs.size(); ++i) {
        hEfficiencyV.push_back(new TH1F(("hEfficiency"+sampleDefs[i].GetLabelS()).c_str(), ";;Efficiency [%]", plotDefs.size(), 0, plotDefs.size()));
        hFractionV.push_back(new TH1F(("hFraction"+sampleDefs[i].GetLabelS()).c_str(), ("Fraction "+sampleDefs[i].GetLabelS()+";;Fraction [%]").c_str(), plotDefs.size(), 0, plotDefs.size()));
        // Set the labels
        for (size_t j = 0; j <nCutBins; ++j) {
            hEfficiencyV[i]->GetXaxis()->SetBinLabel(j+1, plotDefs[j].GetVarLabel());
            hFractionV[i]->GetXaxis()->SetBinLabel(j+1, plotDefs[j].GetVarLabel());
        }

        // center bins
        hEfficiencyV[i]->GetXaxis()->SetLabelSize(0.05);
        hFractionV[i]->GetXaxis()->SetLabelSize(0.05);

        // Hide Y axis ticks
        hEfficiencyV[i]->GetXaxis()->SetTickLength(0);
        hFractionV[i]->GetXaxis()->SetTickLength(0);
    }
    // For purity
    TH1F* hPurity = new TH1F("hPurity", "Purity;Cut;Purity [%]", plotDefs.size(), 0, plotDefs.size());
    TH1F* hEfficiencySignal = new TH1F("hEfficiencySignal", "Signal efficiency;Cut;Efficiency [%]", plotDefs.size(), 0, plotDefs.size());
    for (size_t j = 0; j <nCutBins; ++j) {
        hPurity->GetXaxis()->SetBinLabel(j+1, plotDefs[j].GetVarLabel());
        hEfficiencySignal->GetXaxis()->SetBinLabel(j+1, plotDefs[j].GetVarLabel());
    }
    // center bins
    hPurity->GetXaxis()->SetLabelSize(0.05);
    hEfficiencySignal->GetXaxis()->SetLabelSize(0.05);
    // Hide Y axis ticks
    hPurity->GetXaxis()->SetTickLength(0);
    hEfficiencySignal->GetXaxis()->SetTickLength(0);
    
    size_t lastStoredIndex = 0;
    size_t cont = 0;
    size_t signalIndex = 0;
    for (size_t i = 0; i < matrixCounter.size(); ++i) {

        std::vector<int> counts = matrixCounter[i];
        cont=0;
        std::vector<double> eventsFraction;
        
        for (size_t j = 0; j < counts.size(); ++j) {
            cont++;
            double eff = 100. * counts[j]/(double)matrixCounter[0][j];
            hEfficiencyV[j]->SetBinContent(i+1, eff);

            if(sampleDefs[j].IsSignal()){
                signalIndex = j;
                hEfficiencySignal->SetBinContent(i+1, eff);
            }

            eventsFraction.push_back( potScalingMap[j]*counts[j] );
            
        }

        // Fraction of events
        double totalEvents = std::accumulate(eventsFraction.begin(), eventsFraction.end(), 0.);
        for (size_t j = 0; j < eventsFraction.size(); ++j) {
            hFractionV[j]->SetBinContent(i+1, 100. * eventsFraction[j]/totalEvents);
            std::cout << "Sample: " << sampleDefs[j].GetLabelS() << " Events: " << eventsFraction[j] << " Total: " << totalEvents << " Fraction: " << 100. * eventsFraction[j]/totalEvents << std::endl;
        }
        double purity = 100. * eventsFraction[signalIndex]/totalEvents;
        hPurity->SetBinContent(i+1, purity);

    }

    // Canvas per sample
    double minDecade = 1;
    for(int i=0; i<sampleDefs.size(); i++){
        TCanvas *c = new TCanvas(("cEfficiency"+sampleDefs[i].GetLabelS()).c_str(), ("Efficiency "+sampleDefs[i].GetLabelS()).c_str(), 600, 800);
        c->SetFrameLineColor(0);    
        gPad->SetLeftMargin(0.1);
        gPad->SetRightMargin(0.15);
        gPad->SetBottomMargin(0.25);
        gPad->SetTopMargin(0.06);  
        gPad->SetLogy(1); 
        gPad->Update();
    
        hEfficiencyV[i]->SetStats(0);
        hEfficiencyV[i]->SetFillColor(fStyler->GetColor( sampleDefs[i].Style() ));
        hEfficiencyV[i]->SetBarWidth(0.45);
        hEfficiencyV[i]->SetBarOffset(0.1);
       
        hEfficiencyV[i]->GetXaxis()->LabelsOption("v"); 
        hEfficiencyV[i]->Draw("bar Y+");

        // Custom Y axis labels
        hEfficiencyV[i]->GetYaxis()->SetLabelSize(0);
        double minBinContent = hEfficiencyV[i]->GetMinimum(1e-5);
        double minBinContentDecade = std::pow(10, floor(log10(minBinContent)));
        hEfficiencyV[i]->GetYaxis()->SetRangeUser(minBinContentDecade, 100);
        SetCustomYAxisLabelsLog(c, minBinContentDecade, yLabelXPos);
        if(minBinContentDecade < minDecade) minDecade = minBinContentDecade;

        // Add label with the efficiency
        for (size_t j = 1; j <=nCutBins; ++j) {
            double eff = hEfficiencyV[i]->GetBinContent(j);
            double xpos = j-1+0.25;
            const char *formatEff = (eff < 1) ? "%.1g" : "%.1f";
            // Add the label
            TLatex *latex = new TLatex(xpos, eff, Form(formatEff, eff));
            latex->SetTextSize(fLabelTextSize);
            latex->SetTextAlign(12);
            latex->SetTextAngle(90);
            latex->Draw();
        }

        // Save
        c->SaveAs((dirLabel+"/Efficiency"+sampleDefs[i].GetLabelS()+".pdf").c_str());
    }


    // All in the same canvas
    double boxesMinY = 0.2;
    double boxesMaxY = 0.9;
    double sampleLegSpace = (boxesMaxY-boxesMinY)/sampleDefs.size();
    double boxRelSize = 0.2 * sampleLegSpace;


    TCanvas *cAll = new TCanvas("cEfficiencyAll", "Efficiency All",  600, 800);
    cAll->SetFrameLineColor(0);    
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.25);
    gPad->SetTopMargin(0.05);
    gPad->SetLogy(1);
    gPad->Update();
    THStack *hsAll = new THStack("hsAll", "");
    for(int i=0; i<sampleDefs.size(); i++){
        hEfficiencyV[i]->SetStats(0);
        hEfficiencyV[i]->SetFillColor(fStyler->GetColor( sampleDefs[i].Style() ));
        hEfficiencyV[i]->SetBarWidth(0.45);
        hEfficiencyV[i]->SetBarOffset(0.1);
        hsAll->Add(hEfficiencyV[i]);        
    }
    hsAll->Draw("nostackb Y+");
    // Title and labels
    hsAll->GetYaxis()->SetTitle("Efficiency [%]");
    hsAll->GetXaxis()->LabelsOption("v"); 
    

    // Custom Y axis labels
    hsAll->GetYaxis()->SetLabelSize(0);
    double minBinContentDecade = std::pow(10, floor(log10(minDecade)));
    hsAll->GetYaxis()->SetRangeUser(minBinContentDecade, 100);
    SetCustomYAxisLabelsLog(cAll, minBinContentDecade, yLabelXPos);

    // Draw Text equispaced for each sample
    for(int i=0; i<sampleDefs.size(); i++){
        TLatex *latex = new TLatex(0.05, boxesMinY+sampleLegSpace*i+boxRelSize, sampleDefs[i].GetLabelS().c_str());
        latex->SetTextSize(fLabelTextSize);
        latex->SetTextAngle(90);
        latex->SetNDC(kTRUE);
        latex->Draw();
    }

    // Draw TLatex with efficiency numbers for each cut and sample
    double sampleStep = 1./sampleDefs.size();
    for (size_t j = 0; j <nCutBins; ++j) {
        for(int i=0; i<sampleDefs.size(); i++){
            double eff = hEfficiencyV[i]->GetBinContent(j+1);
            const char *formatEff = (eff < 1) ? "%.1g" : "%.1f";
            TLatex *latexEff = new TLatex(j+sampleStep*(i+0.5), eff, Form(formatEff, eff));
            latexEff->SetTextSize(fEffLabelSize);
            latexEff->SetTextAlign(12);
            latexEff->SetTextAngle(90);
            latexEff->SetTextColor(fStyler->GetColor( sampleDefs[i].Style() ));
            latexEff->Draw();
        }
    }

    TPad *pLeg = new TPad("p","p",0., 0., 0.1, 1.);
    pLeg->SetFillStyle(0);
    pLeg->Draw();
    pLeg->cd();
    // Draw boxes equispaced for each sample
    for(int i=0; i<sampleDefs.size(); i++){
        TBox *box = new TBox(0.25, boxesMinY+sampleLegSpace*i, 0.75, boxesMinY+sampleLegSpace*i+boxRelSize);
        box->SetFillColor(fStyler->GetColor( sampleDefs[i].Style() ));
        box->Draw();
    }


    cAll->SaveAs((dirLabel+"/EfficiencyAll.pdf").c_str());

    
    // Draw signal efficiency and purity together
    int fEffColor = kAzure+7;
    int fPurColor = kOrange+7;
    int fPurHatched = 3144;
    int fLineWidth = 10;
    TCanvas *cSignal = new TCanvas("cSignal", "Signal efficiency and purity",  600, 800);
    cSignal->SetFrameLineColor(0);    
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.25);
    gPad->SetTopMargin(0.05);
    gPad->Update();

    THStack *hsSignal = new THStack("hsSignal", "");
    hEfficiencySignal->SetFillColor(fEffColor);
    hPurity->SetFillColor(fPurColor);
    hPurity->SetFillStyle(fPurHatched);

    hsSignal->Add(hEfficiencySignal);
    hsSignal->Add(hPurity);

    hsSignal->Draw("nostackb Y+");
    hsSignal->GetXaxis()->LabelsOption("v");
    hsSignal->GetYaxis()->SetRangeUser(0, 100);
    hsSignal->GetYaxis()->SetLabelSize(0);
    SetCustomYAxisLabels(cSignal, 6, yLabelXPos);



    // Add latex with the purity and efficiency
    for (size_t j = 0; j <nCutBins; ++j) {
        double eff = hEfficiencySignal->GetBinContent(j+1);
        double pur = hPurity->GetBinContent(j+1);
        TLatex *latexEff = new TLatex(j+0.25,hEfficiencySignal->GetBinContent(j+1), Form("%.1f", eff));
        latexEff->SetTextSize(fEffLabelSize);
        latexEff->SetTextAlign(12);
        latexEff->SetTextAngle(90);
        latexEff->SetTextColor(fEffColor);
        latexEff->Draw();
        TLatex *latexPur = new TLatex(j+0.75,hPurity->GetBinContent(j+1), Form("%.1f", pur));
        latexPur->SetTextSize(fEffLabelSize);
        latexPur->SetTextAlign(12);
        latexPur->SetTextAngle(90);
        latexPur->SetTextColor(fPurColor);
        latexPur->Draw();
    }

    TLatex tEff(.075,.32,"Efficiency");  
    tEff.SetTextAngle(90);
    tEff.SetNDC(kTRUE);
    tEff.Draw();
    TLatex tPur(.075,.67,"Purity");
    tPur.SetTextAngle(90);
    tPur.SetNDC(kTRUE);
    tPur.Draw();

    // Draw TBox
    TPad *p = new TPad("p","p",0., 0., 0.1, 1.);
    p->SetFillStyle(0);
    p->Draw();
    p->cd();
    TBox *boxEff = new TBox(0.25, 0.25, 0.75, 0.3);
    boxEff->SetFillColor(fEffColor);
    boxEff ->Draw();
    TBox *boxPur = new TBox(0.25, 0.6, 0.75, 0.65);
    boxPur->SetFillColor(fPurColor);
    boxPur->SetFillStyle(fPurHatched);
    boxPur ->Draw();

    cSignal->SaveAs((dirLabel+"/EfficiencyAndPuritySignal.pdf").c_str());


    // Fraction of events, in THStack
    TCanvas *cFraction = new TCanvas("cFraction", "Fraction of events", 600, 800);
    cFraction->SetFrameLineColor(0);
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.25);
    gPad->SetTopMargin(0.05);
    gPad->Update();
    THStack *hs = new THStack("hs", "");
    for(int i=0; i<sampleDefs.size(); i++){
        hFractionV[i]->SetStats(0);
        hFractionV[i]->SetFillColor(fStyler->GetColor( sampleDefs[i].Style() ));
        hFractionV[i]->SetBarWidth(0.45);
        hFractionV[i]->SetBarOffset(0.1);
        hs->Add(hFractionV[i]);
    }
    hs->Draw("");
    hs->GetXaxis()->LabelsOption("v"); 
    hs->GetYaxis()->SetRangeUser(0, 100);
    hs->GetYaxis()->SetLabelSize(0);
    hs->GetYaxis()->SetTickLength(0);

    TGaxis *axis = new TGaxis(yLabelXPos, 0, yLabelXPos, 100, 0, 100, 510, "+L");
    axis->SetLabelSize(0.);
    //Colot Y axis
    axis->SetLineColor(kBlack);
    axis->Draw();
    axis->SetTitle("Fraction of events [%]");
    axis->SetTitleSize(0.05);
    axis->SetTitleOffset(0.85);
    SetCustomYAxisLabels(cFraction, 6, yLabelXPos);

    // Draw Text equispaced for each sample
    for(int i=0; i<sampleDefs.size(); i++){
        TLatex *latex = new TLatex(0.075, boxesMinY+sampleLegSpace*i+boxRelSize, sampleDefs[i].GetLabelS().c_str());
        latex->SetTextSize(fLabelTextSize);
        latex->SetTextAngle(90);
        latex->SetNDC(kTRUE);
        latex->Draw();
    }

    pLeg->Draw();
    pLeg->cd();
    // Draw boxes equispaced for each sample
    for(int i=0; i<sampleDefs.size(); i++){
        TBox *box = new TBox(0.25, boxesMinY+sampleLegSpace*i, 0.75, boxesMinY+sampleLegSpace*i+boxRelSize);
        box->SetFillColor(fStyler->GetColor( sampleDefs[i].Style() ));
        box->Draw();
    }

    cFraction->SaveAs((dirLabel+"/FractionOfEvents.pdf").c_str());


    




}

#endif
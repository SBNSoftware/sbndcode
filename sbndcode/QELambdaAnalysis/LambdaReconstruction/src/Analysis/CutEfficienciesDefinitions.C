#ifndef CUTEFFICIENCIES_DEFINITIONS_H
#define CUTEFFICIENCIES_DEFINITIONS_H

#include "CutEfficienciesStyle.C"


// ------- Class for the definition of the samples -------
class SampleDef {
public:
    SampleDef(const std::string& var = "", const std::string& label = "", bool isSignal = false, const std::string& latexLabel = "", int styleIndex=-1,  const std::string& weight = "1")
        : fVar(var),
        fLabel(label),
        fIsSignal(isSignal),
        fLatexLabel(latexLabel),
        fStyleIndex(styleIndex),
        fWeight(weight),
        fNEvents(0) {
    }

    // Getter methods
    TString GetVar() const {
        return fVar;
    }

    TString GetLabel() const {
        return fLabel;
    }

    std::string GetLatexLabel() const {
        return fLatexLabel;
    }

    std::string GetVarS() const {
        return fVar.Data();
    }

    std::string GetLabelS() const {
        return fLabel.Data();
    }

    TString GetWeight() const {
        return fWeight;
    }

    std::string GetWeightS() const {
        return fWeight.Data();
    }

    bool IsSignal() const {
        return fIsSignal;
    }

    // Setter methods
    void SetVar(const TString& var) {
        fVar = var;
    }

    void SetLabel(const TString& label) {
        fLabel = label;
    }

    void SetIsSignal(bool isSignal) {
        fIsSignal = isSignal;
    }

    void SetWeight(const TString& weight) {
        fWeight = weight;
    }

    void SetNEvents(int nEvents) {
        fNEvents = nEvents;
    }

    int GetNEvents() const {
        return fNEvents;
    }

    int Style() const {
        return fStyleIndex;
    }
private:
    TString fVar;
    TString fLabel;
    std::string fLatexLabel;
    bool fIsSignal;
    int fStyleIndex;
    TString fWeight;
    int fNEvents;
};


// ------- Class for the definition of the binning -------
class Binning {
public:
    Binning(double x1 = 0, double x2 = 1, int nBins = 2)
        : fX1(x1), fX2(x2), fNBins(nBins) {
    }

    // Getter methods
    double GetX1() const {
        return fX1;
    }

    double GetX2() const {
        return fX2;
    }

    int GetNBins() const {
        return fNBins;
    }

    // Setter methods
    void SetX1(double x1) {
        fX1 = x1;
    }

    void SetX2(double x2) {
        fX2 = x2;
    }

    void SetNBins(int nBins) {
        fNBins = nBins;
    }

private:
    double fX1;
    double fX2;
    int fNBins;
};


// ------- Class for the definition of the cuts -------
enum class CutType {
    kNone=-1,
    kRight=0,
    kLeft,
    kRightInt,
    kLeftInt,
    kCenter,
    k2D
};


// ------- Class for the definition of the plots -------
class PlotDef {
public:

    // --- Constructor
    PlotDef(const TString& var = "",
            const TString& cut = "",
            const CutType cutType = CutType::kRight,
            double cutValue = 0,
            const Binning& bins = Binning(0, 1, 2),
            bool accumulateCut = true,
            const TString& varLabel = "",
            const TString& cutLabel = "",
            bool norm = false,
            bool log = false,
            const TString& suffix = "")
        : fVar(var),
          fCut(cut),
          fCutType(cutType),
          fCutValue(cutValue),
          fBins(bins),
          fAccumulateCut(accumulateCut),
          fVarLabel(varLabel),
          fCutLabel(cutLabel),
          fNormalize(norm),
          fLog(log),
          fSuffix(suffix)
    {
        if(cutType == CutType::kNone){
            fCut = "0==0";
        }
        else if(cutType == CutType::kCenter){
            fCut =  fCut + ("==" + std::to_string(cutValue)).c_str();
        }
        else if(cutType == CutType::kRight){
            fCut =  fCut + (">=" + std::to_string(cutValue)).c_str();
            std::ostringstream streamObj;
            streamObj <<  std::fixed << std::setprecision(2) << cutValue;
            fCutLabel = fCutLabel + ("\\ >= " + streamObj.str()).c_str();
        }
        else if(cutType == CutType::kLeft || cutType == CutType::kLeftInt){
            fCut =  fCut + ("<=" + std::to_string(cutValue)).c_str();
            std::ostringstream streamObj;
            streamObj <<  std::fixed << std::setprecision(2) << cutValue;
            fCutLabel = fCutLabel + ("\\ <=" + streamObj.str()).c_str();
        }
    }

    // --- Getter methods
    TString GetSuffix() const {
        return fSuffix;
    }

    TString GetVar() const {
        return fVar;
    }

    std::string GetVarS() const {
        return fVar.Data();
    }
    
    TString GetCut() const {
        return fCut;
    }

    double GetCutValue() const {
        return fCutValue;
    }

    TString GetVarLabel(bool useStdString=0) const {
        return fVarLabel;
    }

    std::string GetVarLabelS(bool useStdString=0) const {
        return fVarLabel.Data();
    }

    TString GetCutLabel() const {
        return fCutLabel;
    }

    TString GetCutLabelS() const {
        return fCutLabel.Data();
    }

    CutType GetCutType() const {
        return fCutType;
    }

    Binning GetBins() const {
        return fBins;
    }

    bool GetAccumulateCut() const {
        return fAccumulateCut;
    }

    bool GetLog() const {
        return fLog;
    }

    bool GetNormalize() const {
        return fNormalize;
    }

    // --- Setter methods
    void SetSuffix(const TString& suffix) {
        fSuffix = suffix;
    }

    void SetVar(const TString& var) {
        fVar = var;
    }

    void SetCut(const TString& cut) {
        fCut = cut;
    }

    void SetVarLabel(const TString& varLabel) {
        fVarLabel = varLabel;
    }

    void SetCutLabel(const TString& cutLabel) {
        fCutLabel = cutLabel;
    }

    void SetBins(const Binning& bins) {
        fBins = bins;
    }

    void SetAccumulateCut(bool accumulateCut) {
        fAccumulateCut = accumulateCut;
    }

    void SetLog(bool log) {
        fLog = log;
    }

    bool operator<(const PlotDef& other) const {
        return this->GetVarS() < other.GetVarS();
    }

private:
    TString fVar;
    TString fCut;
    CutType fCutType;
    double fCutValue;
    Binning fBins;
    bool fAccumulateCut;
    TString fVarLabel;
    TString fCutLabel;
    bool fNormalize;
    bool fLog;
    TString fSuffix;
};


// ------- Class to make the actual plotting -------
class AnaPlot{
    public:
        AnaPlot(int plotIndex, PlotDef plotDef, std::vector<SampleDef> intTypes, std::vector<PlotDef> phaseSpaceVaribles={}, std::vector<PlotDef> otherCuts={});

        void DrawHistograms(TTree* fTree, TCut currentCut, bool afterCut=0, double signalScale=1, double backgroundScale=1);

        PlotDef GetPlotDef() const {
            return fPlotDef;
        };

        std::map<std::string, int> GetCountsV() const {
            return fCountsV;
        };



    private:
        
        int fPlotIndex;
        PlotDef fPlotDef;
        
        int fNTypes;
        std::vector<SampleDef> fIntTypes;

        // Main histograms
        std::map<std::string, TH1F*> fHistV;
        std::map<std::string, TH1F*> fHistCumulativeV;
        std::map<std::string, int> fCountsV;

        // 2D histograms
        std::map<std::string, TH2F*> fHist2DV;
        std::map<std::string, TH2F*> fHist2DCumulativeV;
        
        // Phase space histograms
        std::vector<PlotDef> fPhaseSpaceVars;
        std::vector<TH1F*> fHistVPhaseSpace0;
        std::vector<TH1F*> fHistVPhaseSpace1;

        // Other distributions
        std::vector<PlotDef> fOtherCuts;
        std::vector< std::map<std::string, TH1F*> > fHistVOther;

        // Canvases
        TCanvas *fCanvas;
        TCanvas *fCanvasPhaseSpace;
        TCanvas *fCanvasOther;
        CutStyler *fStyler;

        bool fDrawPhaseSpace;
        bool fDrawOtherDistributions;

        void DrawVerticalLineWithArrow(double x, double x1, double x2, double y1, double y2, CutType cutType);
};

// --- Constructor
AnaPlot::AnaPlot(int plotIndex, PlotDef plotDef, std::vector<SampleDef> intTypes, std::vector<PlotDef> phaseSpaceVaribles, std::vector<PlotDef> otherCuts)
    : fPlotIndex(plotIndex),
    fPlotDef(plotDef),
    fIntTypes(intTypes),
    fNTypes(intTypes.size()),
    fPhaseSpaceVars(phaseSpaceVaribles),
    fOtherCuts(otherCuts),
    fCanvas(new TCanvas( ("c_"+std::to_string(fPlotIndex)).c_str(), "Selection", 0, 0, 1500, 1000)),
    fCanvasPhaseSpace( new TCanvas( ("cPS_"+std::to_string(fPlotIndex)).c_str(),"Phase Space",0, 0, 1200, 850)),
    fCanvasOther( new TCanvas( ("cOther_"+std::to_string(fPlotIndex)).c_str(),"Other", 0, 0, 1500, 1000)),
    fDrawPhaseSpace(phaseSpaceVaribles.size()>0),
    fDrawOtherDistributions(plotDef.GetAccumulateCut() && otherCuts.size()>0),
    fStyler(new CutStyler(0))
{

    //--------- Hisogram initializations
    for (size_t j = 0; j < fIntTypes.size(); ++j){
    
        // Sample type labels
        std::string intTypeLabel = fIntTypes[j].GetLabelS();
        std::string plotLabel = std::to_string(fPlotIndex)+"_"+std::to_string(j);

        // Create histograms
        TH1F *hAux = new TH1F(plotLabel.c_str(), plotLabel.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());
        TH1F *hAuxCumulative = new TH1F( (plotLabel+"C").c_str(), plotLabel.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());

        // Initialize histograms
        fCountsV[intTypeLabel] = 0;
        fHistV[intTypeLabel] = hAux;
        fHistCumulativeV[intTypeLabel] = hAuxCumulative;

        // 2D cut option 
        if(fPlotDef.GetCutType()==CutType::k2D){
            TH2F *hAux2D = new TH2F(plotLabel.c_str(), plotLabel.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());
            TH2F *hAux2DCumulative = new TH2F((plotLabel+"C").c_str(), plotLabel.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());
            
            fHist2DV[intTypeLabel] = hAux2D;
            fHist2DCumulativeV[intTypeLabel] = hAux2DCumulative;
        }
        
    }

    //--------- Phase space initializations
    for (size_t j = 0; j < fPhaseSpaceVars.size(); ++j){
        std::string plotLabel0 = std::to_string(fPlotIndex)+"_"+std::to_string(j)+"_PS0";
        TH1F *hAux0 = new TH1F(plotLabel0.c_str(), (";"+fPhaseSpaceVars[j].GetVarLabelS()+";# events").c_str(), fPhaseSpaceVars[j].GetBins().GetNBins(), fPhaseSpaceVars[j].GetBins().GetX1(), fPhaseSpaceVars[j].GetBins().GetX2());
        
        std::string plotLabel1 = std::to_string(fPlotIndex)+"_"+std::to_string(j)+"_PS1";
        TH1F *hAux1 = new TH1F(plotLabel1.c_str(), (";"+fPhaseSpaceVars[j].GetVarLabelS()+";# events").c_str(), fPhaseSpaceVars[j].GetBins().GetNBins(), fPhaseSpaceVars[j].GetBins().GetX1(), fPhaseSpaceVars[j].GetBins().GetX2());
        
        hAux0->SetStats(0);
        hAux1->SetStats(0);
        hAux0->SetLineColor(kOrange-3);
        hAux1->SetLineColor(kAzure+2);
        hAux0->SetLineStyle(kSolid);
        hAux1->SetLineStyle(kDashed);

        fHistVPhaseSpace0.push_back(hAux0);
        fHistVPhaseSpace1.push_back(hAux1);    
    }

    //--------- Other distributions initializations
    for (size_t j = 0; j < fOtherCuts.size(); ++j){
        std::map<std::string, TH1F*> histMap;
        for (size_t k = 0; k < fIntTypes.size(); ++k){
            std::string plotLabel = std::to_string(fPlotIndex)+"_"+std::to_string(j)+"_"+std::to_string(k);
            TH1F *hAux = new TH1F(plotLabel.c_str(), (";"+fOtherCuts[j].GetVarLabelS()+";# events").c_str(), fOtherCuts[j].GetBins().GetNBins(), fOtherCuts[j].GetBins().GetX1(), fOtherCuts[j].GetBins().GetX2());
            hAux->SetStats(0);
            histMap[fIntTypes[k].GetLabelS()] = hAux;
        }
        fHistVOther.push_back(histMap);
    }


}

// --- Util to draw vertical line
void AnaPlot::DrawVerticalLineWithArrow(double x, double x1, double x2, double y1, double y2, CutType cutType) {

    if(cutType == CutType::kNone) return;
    // Get Y range value
    double minY = std::min(y1, y2);
    double maxY = std::max(y1, y2);

    // Get X range value
    double x_min = std::min(x1, x2);
    double x_max = std::max(x1, x2);    

    
    double minX = x;
    double maxX = x;
    if(cutType == CutType::kRight ){
        maxX = x_max;
    }
    else if(cutType == CutType::kLeft || cutType == CutType::kLeftInt ){
        minX = x_min;
    }
    
    if(cutType == CutType::kLeftInt || cutType == CutType::kCenter ){
        maxX = maxX+1;
    }


    if(cutType == CutType::kCenter ){
        // Create a TLine object for the vertical line
        TLine *line1 = new TLine(minX, minY, minX, maxY);
        line1->SetLineColor(40);
        line1->SetLineWidth(2);
        line1->Draw("same");
        TLine *line2 = new TLine(maxX, minY, maxX, maxY);
        line2->SetLineColor(40);
        line2->SetLineWidth(2);
        line2->Draw("same");

        // Calculate the position for the arrowhead (in the middle of the line)
        double yMid = (line1->GetY1() + line1->GetY2()) / 2;
        std::string marker = "|>";
        TArrow *arrow1 = new TArrow(minX - 0.02, yMid, minX + 0.02, yMid, 0.02, marker.c_str());
        arrow1->SetLineColor(40);
        arrow1->SetFillColor(40);
        arrow1->Draw();

        // Calculate the position for the arrowhead (in the middle of the line)
        yMid = (line2->GetY1() + line2->GetY2()) / 2;
        marker = "<|";
        TArrow *arrow2 = new TArrow(maxX - 0.02, yMid, maxX + 0.02, yMid, 0.02, marker.c_str());
        arrow2->SetLineColor(40);
        arrow2->SetFillColor(40);
        arrow2->Draw();
    }
    else{
        // Create a TLine object for the vertical line
        if(cutType==CutType::kLeftInt) x+=1;
        TLine *line = new TLine(x, minY, x, maxY);
        line->SetLineColor(40);
        line->SetLineWidth(2);
        line->Draw("same");

        // Calculate the position for the arrowhead (in the middle of the line)
        double yMid = (line->GetY1() + line->GetY2()) / 2;
        std::string marker = (cutType == CutType::kRight)? "|>":"<|";
        TArrow *arrow = new TArrow(x - 0.02, yMid, x + 0.02, yMid, 0.02, marker.c_str());
        arrow->SetLineColor(40);
        arrow->SetFillColor(40);
        arrow->Draw();
    }
    
    // Create a TBox object for the highlight box
    std::cout<<"minX "<<minX<<" maxX "<<maxX<<" minY "<<minY<<" maxY "<<maxY<<std::endl;
    TBox *highlightBox = new TBox(minX, minY, maxX, maxY);
    highlightBox->SetFillColorAlpha(40, 0.2); // Set transparent green color
    highlightBox->SetFillStyle(1001); // Set the fill style
    highlightBox->Draw("same");

    return;
}

// --- Drawing function
void AnaPlot::DrawHistograms(TTree* fTree, TCut currentCut, bool afterCut, double signalScale, double backgroundScale){

    // --- Frame histograms and legend
    std::string plotIndex = std::to_string(fPlotIndex);
    std::cout<<"Making plot "<<fPlotIndex<<std::endl;
    

    // --- 2D cut option (not well impplemented yet)
    if(fPlotDef.GetCutType()==CutType::k2D){
        std::string plotLabel = "hAux2D"+plotIndex;
        std::cout<<"2D plots not implemented yet"<<std::endl;
        TCanvas *c3 = new TCanvas(("c3"+plotIndex).c_str(),"2DSelection",0, 0, 1000, 750);
        c3->cd();

        TPad *pad1 = new TPad("pad1","pad1",0,0.5,0.5,1);
        pad1->Draw();
        TPad *pad2 = new TPad("pad2","pad2",0.5,0.5,1,1);
        pad2->Draw();
        TPad *pad3 = new TPad("pad3","pad3",0,0,0.5,0.5);
        pad3->Draw();
        TPad *pad4 = new TPad("pad4","pad4",0.5,0,1,0.5);
        pad4->Draw();
        
         
        std::string plotAxisLabels="1;# hits track 1;# hits track 2";

        plotLabel = "hAux2D";
        TH2F *hAux2D = new TH2F(plotLabel.c_str(), plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());

        std::cout<<" Cut name " << fPlotDef.GetVarS() << std::endl; 
        for (size_t j = 0; j < fIntTypes.size(); ++j) {
            
            TCut sampelCut(fIntTypes[j].GetVar());
            std::string plotLabel = "hAux2D";

            // --- Draw distribution
            pad1->cd();
            fTree->Draw( (fPlotDef.GetVarS()+">>"+std::to_string(fPlotIndex)+"_"+std::to_string(j)).c_str(), sampelCut && currentCut );
        }

        hAux2D->Draw();
        
        for (size_t j = 0; j < fIntTypes.size(); ++j) {
            std::string intTypeLabel = fIntTypes[j].GetLabelS();
            fHist2DV[intTypeLabel]->Draw("same");
            // set marker styles and colors
            fHist2DV[intTypeLabel]->SetLineColor(fStyler->GetColor(j));
            fHist2DV[intTypeLabel]->SetLineWidth(2);
            fHist2DV[intTypeLabel]->SetMarkerStyle(fStyler->GetMarkerStyle(j));
            fHist2DV[intTypeLabel]->SetMarkerSize(1);
            fHist2DV[intTypeLabel]->SetMarkerColorAlpha(fStyler->GetColor(j), 0.6);
        }

        c3->Update();
        c3->SaveAs( ("OutputPlots/test2D_"+std::to_string(fPlotIndex)+".pdf").c_str());
        return;
    }

   
    fCanvas->cd();
    std::vector<TPad*> padV = buildpadcanvas(3, 2);

    // Top margins
    for (size_t i = 0; i < padV.size(); ++i) {
        padV[i]->SetTopMargin(0.15);
    }
    std::string plotAxisLabels=";"+fPlotDef.GetVarLabelS()+"; # events";

    std::string plotLabel = "hAux"+std::to_string(fPlotIndex);
    TH1F *hAux = new TH1F( plotLabel.c_str(), plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());

    plotLabel = "hAuxEff"+std::to_string(fPlotIndex);
    TH1F *hAuxEff = new TH1F(plotLabel.c_str(), plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());

    plotLabel = "hAuxSignificance"+std::to_string(fPlotIndex);
    TH1F *hAuxSignificance = new TH1F(plotLabel.c_str(), plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());

    plotLabel = "hCounter"+std::to_string(fPlotIndex);
    TH1F *hCounter = new TH1F( plotLabel.c_str(), plotAxisLabels.c_str(), 1, 0, 1);

    double legX1 = 0.6;
    double legX2 = 0.95;
    double legY1 = 0.6;
    double legY2 = 0.95;
    legX1 = 0.15;
    legX2 = .95;
    legY1 = 0.85;
    legY2 = 1.;
    TLegend *legend = new TLegend(legX1, legY1, legX2, legY2);
    legend->SetNColumns(fIntTypes.size());
    legend->SetLineColorAlpha(0, 0);
    legend->SetFillColorAlpha(0, 0);
       
    padV[1]->cd();
    hAux->Draw();
    double maxVal=0;
    double maxValInt=0;
    double minVal=1000000;
    double minValInt=0;

    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        
        // --- Interaction type
        fHistV[fIntTypes[j].GetLabelS()]->Draw();
        TCut sampelCut(fIntTypes[j].GetVar());
       
        std::string intTypeLabel = fIntTypes[j].GetLabelS();
        std::string plotLabel = std::to_string(fPlotIndex)+"_"+std::to_string(j);

        // --- Draw distribution
        padV[1]->cd();
        fTree->Draw( (fPlotDef.GetVarS()+">>"+plotLabel).c_str(), sampelCut && currentCut );
        
        // Scale distribution
        if(fIntTypes[j].IsSignal())
            fHistV[intTypeLabel]->Scale(signalScale);
        else
            fHistV[intTypeLabel]->Scale(backgroundScale);
 
        // --- Draw cumulative distribution
        padV[3]->cd();

        if(fPlotDef.GetCutType()==CutType::kCenter){
            fHistCumulativeV[intTypeLabel] = (TH1F*)fHistV[intTypeLabel]->Clone();
        }
        else if(fPlotDef.GetCutType()==CutType::kRight){
            if(fIntTypes[j].IsSignal()){
                fHistCumulativeV[intTypeLabel] = (TH1F*)fHistV[intTypeLabel]->GetCumulative(false);
            }
            else{
                fHistCumulativeV[intTypeLabel] = (TH1F*)fHistV[intTypeLabel]->GetCumulative(false);
            }
        }
        else{
            if(fIntTypes[j].IsSignal()){
                fHistCumulativeV[intTypeLabel] = (TH1F*)fHistV[intTypeLabel]->GetCumulative();
            }
            else{
                fHistCumulativeV[intTypeLabel] = (TH1F*)fHistV[intTypeLabel]->GetCumulative();
            }
        }


        // --- Get counts
        padV[1]->cd();
        plotLabel = "hCounter"+std::to_string(fPlotIndex);
        fTree->Draw( ("0>>"+plotLabel).c_str(), fIntTypes[j].GetWeight()*TCut(sampelCut && currentCut && fPlotDef.GetCut()), "hist");
        fCountsV[intTypeLabel] = hCounter->Integral();


        // --- Set styles
        fHistV[intTypeLabel]->SetLineColor(fStyler->GetColor( fIntTypes[j].Style() )); 
        
        fHistV[intTypeLabel]->SetLineStyle(fStyler->GetLineStyle( fIntTypes[j].Style() ));
        fHistV[intTypeLabel]->SetMarkerStyle(fStyler->GetMarkerStyle( fIntTypes[j].Style() ));
        fHistV[intTypeLabel]->SetMarkerColor(fStyler->GetColor( fIntTypes[j].Style() ));
        fHistV[intTypeLabel]->SetStats(0);
        legend->AddEntry(fHistV[intTypeLabel], fIntTypes[j].GetLabelS().c_str(), "lp");
        

        fHistCumulativeV[intTypeLabel]->SetLineColor(fStyler->GetColor( fIntTypes[j].Style() ));
        fHistCumulativeV[intTypeLabel]->SetMarkerStyle(fStyler->GetMarkerStyle( fIntTypes[j].Style() ));
        fHistCumulativeV[intTypeLabel]->SetMarkerColor(fStyler->GetColor( fIntTypes[j].Style() ));
        fHistCumulativeV[intTypeLabel]->Scale(100./fHistV[intTypeLabel]->Integral());
        fHistCumulativeV[intTypeLabel]->SetLineStyle(fStyler->GetLineStyle( fIntTypes[j].Style() ));
        fHistCumulativeV[intTypeLabel]->SetStats(0);
        

        // --- Set maximum
        if(fHistV[intTypeLabel]->GetMaximum()>maxVal){
            maxVal = fHistV[intTypeLabel]->GetMaximum();
            maxValInt = fHistV[intTypeLabel]->Integral();
        }

        if(fHistV[intTypeLabel]->GetMinimum(0)<minVal){
            minVal = fHistV[intTypeLabel]->GetMinimum(0);
            minValInt = fHistV[intTypeLabel]->Integral();
        }

        // --- Draw phase space
        if(fDrawPhaseSpace){
            // --- Draw phase space
            if(fIntTypes[j].IsSignal()){
                for(int i=0; i<fPhaseSpaceVars.size(); i++){
                    std::string plotLabel0 = std::to_string(fPlotIndex)+"_"+std::to_string(i)+"_PS0";
                    fTree->Draw( (fPhaseSpaceVars[i].GetVarS()+">>"+plotLabel0).c_str(), sampelCut, "hist");
                    std::string plotLabel1 = std::to_string(fPlotIndex)+"_"+std::to_string(i)+"_PS1";
                    fTree->Draw( (fPhaseSpaceVars[i].GetVarS()+">>"+plotLabel1).c_str(), sampelCut && currentCut && fPlotDef.GetCut(), "hist");
                }
            }
        }

        // --- Draw other distributions
        if(fDrawOtherDistributions && afterCut){
            for(int i=0; i<fOtherCuts.size(); i++){
                std::string plotLabel = std::to_string(fPlotIndex)+"_"+std::to_string(i)+"_"+std::to_string(j);
                for(int k=0; k<fIntTypes.size(); k++){
                    fTree->Draw( (fOtherCuts[i].GetVarS()+">>"+plotLabel).c_str(), fIntTypes[k].GetWeight()*TCut(sampelCut && currentCut && fPlotDef.GetCut()), "hist");
                }
            }
        }


    }

    // ------ Draw histograms ------ 
    // Frame histogram
    bool fNormalize = fPlotDef.GetNormalize();
    plotLabel = "hAuxFrame"+std::to_string(fPlotIndex);
    double frameYMin = fNormalize? 0.8*(1.*minVal/minValInt):0.9*minVal;
    double frameYMax = fNormalize? 1.2*(1.*maxVal/maxValInt):maxVal*1.1;

    padV[1]->cd();
    TH2F *hAuxFrame = new TH2F( plotLabel.c_str(), plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2(), 100, frameYMin, frameYMax);
    hAuxFrame->SetStats(0);
    if(fNormalize) hAuxFrame->GetYaxis()->SetTitle("AU");
    hAuxFrame->Draw();
    
    // Draw histograms
    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        if(fNormalize)
            fHistV[fIntTypes[j].GetLabelS()]->Scale(1./fHistV[fIntTypes[j].GetLabelS()]->Integral());
        fHistV[fIntTypes[j].GetLabelS()]->Draw("same hist");
    }
    
    // Legend and vertical lines
    DrawVerticalLineWithArrow( fPlotDef.GetCutValue(), hAux->GetXaxis()->GetXmin(),  hAux->GetXaxis()->GetXmax(), 0, maxVal*1.1, fPlotDef.GetCutType() );
    legend->SetFillColorAlpha(0, 0);
    legend->Draw("same");

    // ------ Draw histograms (log) ------ 
    padV[2]->cd();
    padV[2]->SetLogy();
    TH2F *hAuxFrame2 = new TH2F( plotLabel.c_str(), plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2(), 100, frameYMin, frameYMax);
    hAuxFrame2->SetStats(0);
    if(fNormalize) hAuxFrame2->GetYaxis()->SetTitle("AU");
    hAuxFrame2->Draw();

    // Draw histograms
    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        fHistV[fIntTypes[j].GetLabelS()]->Draw("same hist");
    }

    
    
    // Legend and vertical lines
    DrawVerticalLineWithArrow( fPlotDef.GetCutValue(), hAux->GetXaxis()->GetXmin(),  hAux->GetXaxis()->GetXmax(), 0, maxVal*1.1, fPlotDef.GetCutType() );
    legend->SetFillColorAlpha(0, 0);
    legend->Draw("same"); 
    

    // ------ Draw efficiency curves ------ 
    padV[3]->cd();
    hAuxEff->GetYaxis()->SetRangeUser(-0.1, 101);
    hAuxEff->GetYaxis()->SetTitle("#varepsilon [%]");
    hAuxEff->SetStats(0);
    hAuxEff->Draw();


    gStyle->SetOptLogy(0);
    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        fHistCumulativeV[fIntTypes[j].GetLabelS()]->Draw("same p lhist");
    }

    legend->Draw("same");
    padV[3]->Update();


    // ------ Draw Significance ------ 
    padV[4]->cd();
    hAuxSignificance->SetStats(0);
    

    hAuxSignificance->SetStats(0);
    hAuxSignificance->GetYaxis()->SetTitle("S/#sqrt{BG}");
    hAuxSignificance->Draw();
    // new Th1 to sotre the signal
    TH1F *hAuxSigSignal = new TH1F("hAuxSigSignal", plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());
    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        if(fIntTypes[j].IsSignal()){
            hAuxSigSignal->Add(fHistCumulativeV[fIntTypes[j].GetLabelS()]);
        }
    }

    // Legend for significance
    TLegend *legendSig = new TLegend(legX1, legY1, legX2, legY2);
    legendSig->SetNColumns(fIntTypes.size()-1);
    legendSig->SetLineColorAlpha(0, 0);
    legendSig->SetFillColorAlpha(0, 0);

    // now draw all the other backgrounds
    // vector to store TGrapsh
    std::vector<TGraph*> grV;
    // minmium and maximum Y values
    double minY = 1000000;
    double maxY = 0;
    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        if(!fIntTypes[j].IsSignal()){
           
            //create x and y vectors and loop over bins, add signal/sqrt(bg)
            std::vector<double> x;
            std::vector<double> y;
            for (int i = 1; i <= fPlotDef.GetBins().GetNBins(); i++) {
                double binValueBG = (fHistCumulativeV[fIntTypes[j].GetLabelS()]->GetBinContent(i));
                double binValueSignal = hAuxSigSignal->GetBinContent(i);
                double binCenter = fHistV[fIntTypes[j].GetLabelS()]->GetBinCenter(i);
                if(binValueBG>0){
                    x.push_back(binCenter);
                    y.push_back(binValueSignal/std::sqrt(binValueBG));
                }
            }

            // create TGraph
            TGraph *gr = new TGraph(x.size(), &x[0], &y[0]);
            gr->SetLineColor(fStyler->GetColor( fIntTypes[j].Style() ));
            gr->SetMarkerStyle(fStyler->GetMarkerStyle( fIntTypes[j].Style() ));
            gr->SetMarkerColor(fStyler->GetColor( fIntTypes[j].Style() ));
            gr->SetLineStyle(fStyler->GetLineStyle( fIntTypes[j].Style() ));

            grV.push_back(gr);

            // add to legend
            legendSig->AddEntry(gr, fIntTypes[j].GetLabelS().c_str(), "lp");

            double min = gr->GetYaxis()->GetXmin();
            double max = gr->GetYaxis()->GetXmax();
            if(min<minY) minY = min;
            if(max>maxY) maxY = max;
        }
    }

    hAuxSignificance->GetYaxis()->SetRangeUser(minY*0.9, maxY*1.1);
    for (size_t j = 0; j < grV.size(); ++j) {
        grV[j]->Draw("same lp");
    }
    
    if(grV.size()>1){
        legendSig->Draw("same");
    }
    padV[4]->Update();

    // ------ Draw ROC curve ------ 
    bool fDrawROC = false;
    if(fDrawROC){
        padV[5]->cd();

        // TH1 to sotre the signal
        TH1F *hAuxEffSignal = new TH1F("hAuxEffSignal", plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());
        for (size_t j = 0; j < fIntTypes.size(); ++j) {
            if(fIntTypes[j].IsSignal()){
                hAuxEffSignal->Add(fHistCumulativeV[fIntTypes[j].GetLabelS()]);
            }
        }

        TH1F *hAuxROC = new TH1F("hAuxROC", "ROC;;", 200, 0, 100);
        hAuxROC->SetStats(0);
        hAuxROC->GetXaxis()->SetTitle("Signal efficiency [%]");
        hAuxROC->GetYaxis()->SetTitle("BG rejection [%]");
        hAuxROC->GetYaxis()->SetRangeUser(-0.1, 101);
        hAuxROC->Draw();

        std::vector<TGraph*> grVROC;
        for (size_t j = 0; j < fIntTypes.size(); ++j) {

            if(!fIntTypes[j].IsSignal()){
                //create x and y vectors and loop over bins, add signal/sqrt(bg)
                std::vector<double> x;
                std::vector<double> y;
                for (int i = 1; i <= fPlotDef.GetBins().GetNBins(); i++) {
                    double binValueBG = (fHistCumulativeV[fIntTypes[j].GetLabelS()]->GetBinContent(i));
                    double binValueSignal = hAuxEffSignal->GetBinContent(i);
                    
                    if(binValueBG>0){
                        x.push_back(binValueSignal);
                        y.push_back(binValueBG);
                    }
                }

                // create TGraph
                TGraph *gr = new TGraph(x.size(), &x[0], &y[0]);
                gr->SetLineColor(fStyler->GetColor( fIntTypes[j].Style() ));
                gr->SetMarkerStyle(fStyler->GetMarkerStyle( fIntTypes[j].Style() ));
                gr->SetMarkerColor(fStyler->GetColor( fIntTypes[j].Style() ));
                gr->SetLineStyle(fStyler->GetLineStyle( fIntTypes[j].Style() ));

                grVROC.push_back(gr);   
            }
        }
        hAuxROC->Draw();
        legend->Draw("same");

        for (size_t j = 0; j < grVROC.size(); ++j) {
            grVROC[j]->Draw("same l");
        }
        padV[5]->Update();
    }

    // ------ Draw purity curve ------
    padV[6]->cd();
    TH1F *hAuxPurity = new TH1F("hAuxPurity", plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());
    hAuxPurity->SetStats(0);
    hAuxPurity->GetYaxis()->SetTitle("#rho [%]");
    hAuxPurity->Draw();

    // Map to store the events from hits integral
    std::map<std::string, double> events;
    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        events[fIntTypes[j].GetLabelS()] = fHistV[fIntTypes[j].GetLabelS()]->Integral();
    }

    std::vector<TGraph*> grVPurity;
    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        if(!fIntTypes[j].IsSignal()){
            //create x and y vectors and loop over bins, add signal/sqrt(bg)
            std::vector<double> x;
            std::vector<double> y;
            double bgIntegral = fHistV[fIntTypes[j].GetLabelS()]->Integral();
            double signalIntegral = fHistV["Signal"]->Integral();
            for (int i = 0; i <= fPlotDef.GetBins().GetNBins(); i++) {
                
                double binValueBG = 0.01*fHistCumulativeV[fIntTypes[j].GetLabelS()]->GetBinContent(i) * bgIntegral;
                double binValueSignal = 0.01*fHistCumulativeV["Signal"]->GetBinContent(i) * signalIntegral;
                double binCenter = fHistV[fIntTypes[j].GetLabelS()]->GetBinCenter(i);
                std::cout<<i<<" "<<binCenter<<" binValueBG "<<binValueBG<<" binValueSignal "<<binValueSignal<<std::endl;
                if(binValueBG+binValueSignal>0){
                    double purity = 100*binValueSignal/(binValueSignal+binValueBG);
                    x.push_back( binCenter );
                    y.push_back( purity );
                    std::cout<<"  purity="<<purity<<std::endl;
                }
            }

            // create TGraph
            TGraph *gr = new TGraph(x.size(), &x[0], &y[0]);
            gr->SetLineColor(fStyler->GetColor( fIntTypes[j].Style() ));
            gr->SetMarkerStyle(fStyler->GetMarkerStyle( fIntTypes[j].Style() ));
            gr->SetMarkerColor(fStyler->GetColor( fIntTypes[j].Style() ));
            gr->SetLineStyle(fStyler->GetLineStyle( fIntTypes[j].Style() ));

            grVPurity.push_back(gr);
        }
    }

    hAuxPurity->Draw();
    legendSig->Draw("same");
    
    double minYPurity = 1000000;
    double maxYPurity = 0;
    for (size_t j = 0; j < grVPurity.size(); ++j) {
        grVPurity[j]->Draw("same l");
        double min = grVPurity[j]->GetYaxis()->GetXmin();
        double max = grVPurity[j]->GetYaxis()->GetXmax();
        if(min<minYPurity) minYPurity = min;
        if(max>maxYPurity) maxYPurity = max;
    }
    hAuxPurity->GetYaxis()->SetRangeUser(minYPurity*0.9, maxYPurity*1.1);

    padV[6]->Update();


    // ------ Draw efficiency times purity curve ------
    padV[5]->cd();
    TH1F *hAuxEffPurity = new TH1F("hAuxEffPurity", plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());
    hAuxEffPurity->SetStats(0);
    hAuxEffPurity->GetYaxis()->SetTitle("#varepsilon#rho");
    hAuxEffPurity->Draw();

    std::vector<TGraph*> grVEffPurity;
    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        if(!fIntTypes[j].IsSignal()){
            //create x and y vectors and loop over bins, add signal/sqrt(bg)
            std::vector<double> x;
            std::vector<double> y;
            double bgIntegral = fHistV[fIntTypes[j].GetLabelS()]->Integral();
            double signalIntegral = fHistV["Signal"]->Integral();

            for (int i = 1; i <= fPlotDef.GetBins().GetNBins(); i++) {
                
                double binValueBG = 0.01*fHistCumulativeV[fIntTypes[j].GetLabelS()]->GetBinContent(i) * bgIntegral;
                double binValueSignal = 0.01*fHistCumulativeV["Signal"]->GetBinContent(i) * signalIntegral;
                double binCenter = fHistV[fIntTypes[j].GetLabelS()]->GetBinCenter(i);

                if(binValueBG+binValueSignal>0){
                    double purity = binValueSignal/(binValueSignal+binValueBG);
                    double eff = binValueSignal/signalIntegral;
                    x.push_back( binCenter );
                    y.push_back( eff*purity );
                }
            }

            // create TGraph
            TGraph *gr = new TGraph(x.size(), &x[0], &y[0]);
            gr->SetLineColor(fStyler->GetColor( fIntTypes[j].Style() ));
            gr->SetMarkerStyle(fStyler->GetMarkerStyle( fIntTypes[j].Style() ));
            gr->SetMarkerColor(fStyler->GetColor( fIntTypes[j].Style() ));
            gr->SetLineStyle(fStyler->GetLineStyle( fIntTypes[j].Style() ));

            grVEffPurity.push_back(gr);
        }
    }

    
   

    hAuxEffPurity->Draw();
    legendSig->Draw("same");
    // Range to maximum of gr
    double minYEffPurity = 1000000;
    double maxYEffPurity = 0;
    for (size_t j = 0; j < grVEffPurity.size(); ++j) {
        grVEffPurity[j]->Draw("same l");
        double min = grVEffPurity[j]->GetYaxis()->GetXmin();
        double max = grVEffPurity[j]->GetYaxis()->GetXmax();
        if(min<minYEffPurity) minYEffPurity = min;
        if(max>maxYEffPurity) maxYEffPurity = max;
    }

    hAuxEffPurity->GetYaxis()->SetRangeUser(minYEffPurity*0.9, maxYEffPurity*1.1);

    padV[5]->Update();


    // ------ Draw phase space ------ 
    if(fDrawPhaseSpace){
        fCanvasPhaseSpace->cd();
        fCanvasPhaseSpace->SetLogy();
        // Create the pads
        int nRows = std::ceil(std::sqrt( fPhaseSpaceVars.size() ));
        std::vector<TPad*> padVPhaseSpace = buildpadcanvas( nRows, nRows);

        TLegend *legendRP = new TLegend(0.7, 0.7, 0.9, 0.9);
        legendRP->SetLineColorAlpha(0, 0);
        legendRP->SetFillColorAlpha(0, 0);


        // Draw each plot
        for(int i=0; i<fPhaseSpaceVars.size(); i++){
            padVPhaseSpace[i+1]->cd();
            //fHistVPhaseSpace1[i]->GetXaxis()->SetTitle("x");
            fHistVPhaseSpace1[i]->GetYaxis()->SetTitle("y");
            fHistVPhaseSpace0[i]->GetYaxis()->SetTitle("y");
            // POT scaling
            fHistVPhaseSpace1[i]->Scale(signalScale);
            fHistVPhaseSpace0[i]->Scale(signalScale);
            TRatioPlot * rp = new TRatioPlot(fHistVPhaseSpace1[i], fHistVPhaseSpace0[i]);
            padVPhaseSpace[i+1]->SetLogy();
            rp->SetH1DrawOpt("HIST");
            rp->SetH2DrawOpt("HIST SAME");
            rp->SetGraphDrawOpt("lp");
            
            rp->Draw();
            rp->SetLowBottomMargin(0.3);
            rp->SetLeftMargin(0.15);
            rp->SetSplitFraction(0.35);
            rp->SetSeparationMargin(0.);
            rp->GetLowerRefYaxis()->SetTitle("#varepsilon");
            rp->GetUpperRefYaxis()->SetTitle("# entries");
            rp->GetLowerRefYaxis()->SetTitleOffset(1.);
            rp->GetUpperRefYaxis()->SetTitleOffset(1.);
            rp->GetLowYaxis()->SetNdivisions(1005);
            rp->GetUpperRefYaxis()->SetRangeUser(0.1, fHistVPhaseSpace0[i]->GetMaximum()*1.1);
            rp->GetLowerRefGraph()->SetMarkerStyle(20);
            rp->GetLowerRefGraph()->SetMarkerSize(1);
            rp->GetLowerRefGraph()->SetLineColor(kBlack);
            if(i==0){
                legendRP->AddEntry(fHistVPhaseSpace0[i],"All","l");
                legendRP->AddEntry(fHistVPhaseSpace1[i],"After cut","l");
            }
        
            legendRP->Draw("same");
        }
        
        fCanvasPhaseSpace->WaitPrimitive();
    }


    // ------ Draw other distributions ------
    std::cout<<"fDrawOtherDistributions "<<fDrawOtherDistributions<<std::endl;
    if(fDrawOtherDistributions && afterCut){
        fCanvasOther->cd();
        int nRows = std::ceil(std::sqrt( fOtherCuts.size() ));
        std::vector<TPad*> padVOther = buildpadcanvas( nRows, nRows);
        for(int i=0; i<fOtherCuts.size(); i++){
            padVOther[i+1]->cd();
            
            // log scale
            padVOther[i+1]->SetLogy();
            fHistVOther[i][fIntTypes[0].GetLabelS()]->Draw();
            for (size_t j = 0; j < fIntTypes.size(); ++j) {
                fHistVOther[i][fIntTypes[j].GetLabelS()]->Draw("same hist");
                // --- Set styles
                fHistVOther[i][fIntTypes[j].GetLabelS()]->SetLineColor(fStyler->GetColor( fIntTypes[j].Style() ));
                fHistVOther[i][fIntTypes[j].GetLabelS()]->SetMarkerStyle(fStyler->GetMarkerStyle( fIntTypes[j].Style() ));
                fHistVOther[i][fIntTypes[j].GetLabelS()]->SetMarkerColorAlpha(fStyler->GetColor( fIntTypes[j].Style() ), 0.6);
            }

            // Set limits with 0 and maxVal
            fHistVOther[i][fIntTypes[0].GetLabelS()]->GetYaxis()->SetRangeUser(0.5, maxVal*1.1);

            // Legend
            legend->Draw("same");

        }
    }
        

    // --- Save canvas
    fCanvas->WaitPrimitive();
    fCanvas->Update();
    fCanvasPhaseSpace->Update();
    fCanvasOther->Update();
    std::string stageLabel = afterCut ? "1after" : "0before";
    fCanvas->SaveAs(("OutputPlots/plot"+std::to_string(fPlotIndex)+stageLabel+fPlotDef.GetVarLabelS()+".pdf").c_str());
    if(fDrawPhaseSpace) fCanvasPhaseSpace->SaveAs(("OutputPlots/PhaseSpace/zphaseSpacePlot"+std::to_string(fPlotIndex)+stageLabel+fPlotDef.GetVarLabelS()+"PhaseSpace.pdf").c_str());
    if(fDrawOtherDistributions && afterCut) fCanvasOther->SaveAs(("OutputPlots/OtherDistributions/otherDistributionsPlot"+std::to_string(fPlotIndex)+stageLabel+fPlotDef.GetVarS()+"OtherDistributions.pdf").c_str());
    TFile *fFile = new TFile("OutputPlots/OutputPlots.root","UPDATE");
    fFile->cd();
    fCanvas->Write();
    fCanvasPhaseSpace->Write();
    if(fDrawOtherDistributions && afterCut) fCanvasOther->Write();
    fFile->Close();
    
    return;
    
}


#endif
////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleObjects.h
//
// \brief Definition of SimpleTPCObjects
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCSimpleCalo.h"

#include "TCanvas.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMath.h"
#include "TArc.h"
#include "Fit/Fitter.h"
#include <Math/Functor.h>
#include "TH3F.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TPaveText.h"
#include "TStyle.h"
#include <Math/WrappedMultiTF1.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>

// Calorimetry conversion factor (ADC to charge)
// (1 / 0.0201293) e-/ADC*time_ticks x  23.6e-6 MeV/e-  x)



// Constructor to initialize the collection
SCalo::SCalo(const std::vector<SHit>& hits, double trackAngle, double yAngle)
    : fHitIntegralToEnergy( (1 / 0.0201293) *  23.6e-6 ),
    fStepXToLength(0.3),
    fStepYToLength(0.075),
    fTrackLength(0),
    fCosTrackAngle( std::abs(std::cos(trackAngle)) ),
    fCosY( std::abs(std::cos(yAngle)) ),
    fTrackAngle(trackAngle),
    fHitList(hits)
{
    CalculateResidualRange();
}


// Function to calculate the distance between two hit points
double SCalo::CalculateDistance(const SHit& p1, const SHit& p2) {
    double dx = fStepXToLength * ( p2.X() - p1.X() );
    double dy = fStepYToLength * ( p2.Y() - p1.Y() );
    return std::sqrt(dx * dx + dy * dy);
}

// Method to calculate and store the path lengths
void SCalo::CalculatePathLengths() {
    fPathLengths.clear();

    if(fHitList.size()>1){
        // first hit
        fPathLengths.push_back( CalculateDistance(fHitList[0], fHitList[1]) );

        for (size_t i = 1; i < fHitList.size()-1; ++i) {
            double distance1 = CalculateDistance(fHitList[i - 1], fHitList[i]);
            double distance2 = CalculateDistance(fHitList[i], fHitList[i + 1]);
            double halfDistance = 0.5 * (distance1 + distance2);
            
            fPathLengths.push_back(halfDistance);
        }

        // last hit
        fPathLengths.push_back( CalculateDistance(fHitList[fHitList.size()-1], fHitList[fHitList.size()-2]) );
    }
    else{
        fPathLengths.push_back(0);
    }

}

void SCalo::CalculateResidualRange() {
    
    CalculatePathLengths();
    
    fDepositedEnergy.clear();
    fResidualRange.clear();

    fTrackLength = 0;
    if(fPathLengths.size()>1){
        fTrackLength = std::accumulate(fPathLengths.begin(), fPathLengths.end(), fTrackLength);

        double pathLength = 0;
        for (size_t i = 0; i < fHitList.size(); ++i) {
            fResidualRange.push_back( fTrackLength - pathLength);
            fDepositedEnergy.push_back( (fHitList[i].Integral() * fHitIntegralToEnergy) / fPathLengths[i]);
            //fDepositedEnergy.push_back( (fHitList[i].Integral() * fHitIntegralToEnergy) * fCosTrackAngle );
            pathLength+=fPathLengths[i];
        }

        // Check last points
        double minE = *std::min_element(fDepositedEnergy.begin(), fDepositedEnergy.end());

        if(fResidualRange.size()>3){
            while (!fDepositedEnergy.empty() && fDepositedEnergy.back() <= minE) {
                fDepositedEnergy.pop_back();
                fResidualRange.pop_back();
                fPathLengths.pop_back();
                fHitList.pop_back();
                minE = *std::min_element(fDepositedEnergy.begin(), fDepositedEnergy.end());
            }
        }

        double lastResidualRange = fResidualRange.back();
        // reshift 
        for(size_t i=0; i<fResidualRange.size(); i++){
            fResidualRange[i] = fResidualRange[i] - lastResidualRange;
        }
        fTrackLength = fTrackLength - lastResidualRange;
    }

    if(fResidualRange.size()>=2){
        
        std::vector<double> rangeDiff;
        for(size_t i=0; i<fResidualRange.size()-1; i++){
            rangeDiff.push_back( fResidualRange[i+1] - fResidualRange[i] );
        }
        double sq_sum = std::inner_product(rangeDiff.begin(), rangeDiff.end(), rangeDiff.begin(), 0.0);
        double mean = std::accumulate(rangeDiff.begin(), rangeDiff.end(), 0.0) / rangeDiff.size();
        fResidualRangeRMS = std::sqrt(sq_sum / rangeDiff.size() - mean * mean);

    }
    else{
        fResidualRangeRMS = -1;
    }
    
}


void SCalo::Display(){
    std::cout<<"Displaying calorimetry...track length: "<<fTrackLength<<"\n";

    for (size_t i = 0; i < fPathLengths.size(); ++i) {
        std::cout << " Residual range / dEdx: " << fResidualRange[i] << " / " << fDepositedEnergy[i] << std::endl;
    }

}


// ------ STriangleCalo class ------

double CalculateMean( std::vector<double>& values) {

    double sum = 0;
    for ( double& value : values) {
        sum += value;
    }
    return sum / values.size();
}

double CalculateStdDev( std::vector<double>& values) {

    if(values.size()==0) return -1;
    double mean = CalculateMean(values);
    double sumSqDiff = 0;
    for ( double& value : values) {
        double diff = value - mean;
        sumSqDiff += diff * diff;
    }
    return std::sqrt(sumSqDiff / values.size());
}

// --- Broken line fit ---
Double_t BrokenLine(Double_t *x,Double_t *par) {
  
  Double_t arg = 0.;
  if(x[0]<par[0])
    arg = par[1]+par[2]*(x[0]-par[0]);
  else
    arg = par[1]+par[3]*(x[0]-par[0]);
  
  return arg; 
}


// --- Fisher discriminant 2D ---
TVectorD FisherDiscriminant2D(const std::vector<SHit>& sample1, const std::vector<SHit>& sample2) {
    // Calculate means for each sample
    double mean1X = 0.0, mean1Y = 0.0;
    double mean2X = 0.0, mean2Y = 0.0;

    for (const auto& point : sample1) {
        mean1X += point.X();
        mean1Y += point.Y();
    }

    for (const auto& point : sample2) {
        mean2X += point.X();
        mean2Y += point.Y();
    }

    mean1X /= sample1.size();
    mean1Y /= sample1.size();

    mean2X /= sample2.size();
    mean2Y /= sample2.size();

    std::cout<<mean1X<<" "<<mean1Y<<" "<<mean2X<<" "<<mean2Y<<std::endl;

    // Calculate covariance matrices for each sample
    TMatrixD cov1(2, 2), cov2(2, 2);

    for (const auto& point : sample1) {
        cov1(0, 0) += (point.X() - mean1X) * (point.X() - mean1X);
        cov1(0, 1) += (point.X() - mean1X) * (point.Y() - mean1Y);
        cov1(1, 0) += (point.Y() - mean1Y) * (point.X() - mean1X);
        cov1(1, 1) += (point.Y() - mean1Y) * (point.Y() - mean1Y);
    }

    for (const auto& point : sample2) {
        cov2(0, 0) += (point.X() - mean2X) * (point.X() - mean2X);
        cov2(0, 1) += (point.X() - mean2X) * (point.Y() - mean2Y);
        cov2(1, 0) += (point.Y() - mean2Y) * (point.X() - mean2X);
        cov2(1, 1) += (point.Y() - mean2Y) * (point.Y() - mean2Y);
    }

    cov1 *= 1.0 / (sample1.size() - 1);
    cov2 *= 1.0 / (sample2.size() - 1);

    std::cout<<" Cov "<<cov1(0,0)<<" "<<cov1(0,1)<<" "<<cov1(1,0)<<" "<<cov1(1,1)<<std::endl;

    TMatrixD cov1Inv (cov1.Invert());

    std::cout<<" CovInv "<<cov1Inv(0,0)<<" "<<cov1Inv(0,1)<<" "<<cov1Inv(1,0)<<" "<<cov1Inv(1,1)<<std::endl;
    TVectorD menDiffVector(0, 1, mean1X - mean2X, mean1Y - mean2Y, "END");
    std::cout<<" menDiffVector "<<menDiffVector(0)<<" "<<menDiffVector(1)<<std::endl;
    
    // Calculate the Fisher's discriminant direction
    TVectorD w = cov1Inv * menDiffVector;

    // Print the discriminant direction
    std::cout << "Fisher's Discriminant Direction: (" << w[0] << ", " << w[1] << ")" << std::endl;

    return w;
}


// --- Constructor ---
STriangleCalo::STriangleCalo(CaloAlgorithmPsetType caloPset) : fCaloPset(caloPset)
{
    fRatioUpperLimit = 15.;
}

void STriangleCalo::SetTriangle(STriangle triangle){

    fAllHits.clear();
    fHitsTrack1.clear();
    fHitsTrack2.clear();
    fVertexHits.clear();
    fResidualHits.clear();
    fFitHits.clear();
    
    fTriangle = triangle;
   
   /*float yCorrFactor = 0.075/0.3;

    SPoint mainVertex(triangle.GetMainVertex().X(), triangle.GetMainVertex().Y()*yCorrFactor);
    SPoint vertexB(triangle.GetVertexB().X(), triangle.GetVertexB().Y()*yCorrFactor);
    SPoint vertexC(triangle.GetVertexC().X(), triangle.GetVertexC().Y()*yCorrFactor);
    SHit mainHit(mainVertex.X(), mainVertex.Y());
    std::vector<SHit> track1Hits = triangle.GetTrack1().GetHits();
    for(SHit &h:track1Hits){
        h.SetY( h.Y()*yCorrFactor );
        h.SetWidth( h.Width()*yCorrFactor );
    }
    SLinearCluster newTrack1(track1Hits);
    std::vector<SHit> track2Hits = triangle.GetTrack2().GetHits();
    for(SHit &h:track2Hits){
        h.SetY( h.Y()*yCorrFactor );
        h.SetWidth( h.Width()*yCorrFactor );
    }
    
   
    SLinearCluster newTrack2(track2Hits);

    STriangle newTriangle(mainVertex, vertexB, vertexC, mainHit, newTrack1, newTrack2, triangle.GetMainTrack(), 1, 1);
    
    std::cout<<triangle.GetTrack1().GetTrackEquation().Slope()<<" "<<triangle.GetTrack2().GetTrackEquation().Slope()<<std::endl;
    fTriangle = newTriangle;
    std::cout<<fTriangle.GetTrack1().GetTrackEquation().Slope()<<" "<<fTriangle.GetTrack2().GetTrackEquation().Slope()<<std::endl;*/
    
}


// --- JointFit function ---
void STriangleCalo::JointFitAnalysisFisher(){

    TCanvas *c1 = new TCanvas("c1","A Simple Graph Example", 800, 0, 800,1200);
    c1->SetGrid();
    // Create TPads
    c1->cd();
    TPad *pad1 = new TPad("pad1","This is pad1",0.02,0.52,0.98,0.98);
    TPad *pad2 = new TPad("pad2","This is pad2",0.02,0.02,0.98,0.48);
    pad1->Draw();
    pad2->Draw();

    size_t nhits1 = fTriangle.GetNHitsTrack1();
    size_t nhits2 = fTriangle.GetNHitsTrack2();
    size_t maxNHits = std::min(nhits1, nhits2);
    maxNHits = std::min(nhits1, nhits2);
   

    // x-y values from hits
    size_t hitCont = 0;
    std::vector<double> xV1, yV1, xV2, yV2;
    for(SHit &h:fTriangle.GetTrack1().GetHits()){
        if(hitCont>maxNHits) continue;
        xV1.push_back(h.X());
        yV1.push_back(h.Y());
        hitCont++;
        
    }
    hitCont=0;
    for(SHit &h:fTriangle.GetTrack2().GetHits()){
        if(hitCont>maxNHits) continue;
        xV2.push_back(h.X());
        yV2.push_back(h.Y());
        hitCont++;
    }

    double xVertex = fTriangle.GetMainVertex().X();
    double yVertex = fTriangle.GetMainVertex().Y();

    // Fisher direcrtion
    TVectorD fWVector = FisherDiscriminant2D(fTriangle.GetTrack1().GetHits(), fTriangle.GetTrack2().GetHits());

    // opposite side
    double x1 = fTriangle.GetVertexB().X();
    double y1 = fTriangle.GetVertexB().Y();
    double x2 = fTriangle.GetVertexC().X();
    double y2 = fTriangle.GetVertexC().Y();

    double maxD = 0;
    for(size_t k=0; k<xV1.size(); k++){
        double d = std::hypot( xV1[k]-xVertex, yV1[k]-yVertex);
        if(d>maxD){
            x1 = xV1[k];
            y1 = yV1[k];
            maxD = d;
        }
    }

    maxD = 0;
    for(size_t k=0; k<xV2.size(); k++){
        double d = std::hypot( xV2[k]-xVertex, yV2[k]-yVertex);
        if(d>maxD){
            x2 = xV2[k];
            y2 = yV2[k];
            maxD = d;
        }
    }

    double m = (y2-y1)/(x2-x1);

    double rawSlope1 = fTriangle.GetTrack1().GetTrackEquation().Slope();
    double rawSlope2 = fTriangle.GetTrack2().GetTrackEquation().Slope();

    std::vector<double> x, y;
    x.insert(x.end(), xV1.begin(), xV1.end());
    x.insert(x.end(), xV2.begin(), xV2.end());
    y.insert(y.end(), yV1.begin(), yV1.end());
    y.insert(y.end(), yV2.begin(), yV2.end());

    // project the points
    for(size_t i=0; i<x.size(); i++){
        x[i] = x[i]-xVertex;
        y[i] = y[i]-yVertex;
    }

    size_t n = x.size();


    // project into the line with points (17.3, 12.4) and (13.7, 6.65)
    // get the line equation
    // y = m*x + q
    pad1->cd();
    // create a graph
    TGraph *grRaw = new TGraph(n, &x[0], &y[0]);
    grRaw->SetMarkerColor(4);
    grRaw->SetMarkerStyle(21);
    grRaw->Draw("AP");
    // axis labels
    grRaw->GetXaxis()->SetTitle("Wire ID");
    grRaw->GetYaxis()->SetTitle("Time Tick [0.5#mu s]");
    // axis label offset
    grRaw->GetYaxis()->SetTitleOffset(.7);
    grRaw->GetXaxis()->SetTitleOffset(.83);
    grRaw->SetTitle("View (raw)");
    
    // angle
    double angle = atan(m);
    std::cout << "m: " << m << std::endl;
    std::cout << "angle: " << 180 * angle / TMath::Pi() << std::endl;
    std::cout << "opening angle: " << fTriangle.GetOpeningAngle() << std::endl;
    double isocAngle = ( 180 - fTriangle.GetOpeningAngle() ) / 2;
    std::cout << "isocAngle: " << isocAngle << std::endl;

    // convert the ray slopes
    rawSlope1 = sin(-angle)+rawSlope1*cos(-angle) / (cos(-angle)-rawSlope1*sin(-angle));
    rawSlope2 = sin(-angle)+rawSlope2*cos(-angle) / (cos(-angle)-rawSlope2*sin(-angle));
   
    // distances to the main vertex
    double d1 = std::hypot( x1, y1);
    double d2 = std::hypot( x2, y2);

    // min track angle
    double minTrackAngle1 = std::acos( x1/std::hypot(x1, y1) );
    double minTrackAngle2 = std::acos( x2/std::hypot(x2, y2) );
    std::cout<<" Side track angle: "<<180*minTrackAngle1/TMath::Pi()<<" "<<180*minTrackAngle2/TMath::Pi()<<std::endl;

    // Set angle to the isoceles triangle
    double minTrackAngle = (d1>d2)? minTrackAngle1:minTrackAngle2;
    double newSlope =  180 - 180*minTrackAngle/TMath::Pi() - isocAngle;
    newSlope = TMath::Pi() * newSlope / 180;
    angle = newSlope;
    std::cout<<" New slope: "<< newSlope << "angle "<<angle<<std::endl;

    // Set angle to the Fisher Direction
    angle  = fWVector(1)/fWVector(0);

    // project the points
    for(size_t i=0; i<x.size(); i++){
        double x0 = x[i];
        double y0 = y[i];
        double xproj = x0*cos(angle) + y0*sin(angle);
        double yproj = -x0*sin(angle) + y0*cos(angle);
        x[i] = xproj;
        y[i] = yproj;
    }
    // get min and max x values
    double xmin = *std::min_element(x.begin(), x.end());
    double xmax = *std::max_element(x.begin(), x.end());
    xmin -= 0.1 * (xmax-xmin);
    xmax += 0.1 * (xmax-xmin);
    std::cout << "xmin: " << xmin << std::endl;
    std::cout << "xmax: " << xmax << std::endl;


    // create a graph
    pad2->cd();
    
    TGraph *gr = new TGraph(n, &x[0], &y[0]);
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->GetXaxis()->SetTitle("Wire ID");
    gr->GetYaxis()->SetTitle("Time Tick [0.5#mu s]");
    // axis label offset
    gr->GetYaxis()->SetTitleOffset(.7);
    gr->GetXaxis()->SetTitleOffset(.83);

    gr->SetTitle("Rotated view");
    gr->Draw("AP");
    c1->Update();


    // Fit to the brokenLine function
    TF1 *brokenLine = new TF1("brokenLine",BrokenLine, xmin, xmax, (int)(xmax-xmin));
    brokenLine->SetParameters(0, 0, rawSlope1, rawSlope2);
    brokenLine->SetParNames("breakPoint","y0","slope1","slope2");
    brokenLine->SetLineColor(kRed);
    brokenLine->SetLineWidth(2);
    gr->Fit(brokenLine,"RQ");

    // Draw the fitted function
    brokenLine->Draw("same");

    pad1->cd();
    // Undo the transformation
    double breakPoint = brokenLine->GetParameter(0);
    double y0 = brokenLine->GetParameter(1);
    double slope1 = brokenLine->GetParameter(2);
    double slope2 = brokenLine->GetParameter(3);
    
    // Convert the parameters
    breakPoint = breakPoint*cos(angle) - y0*sin(angle);
    y0 = breakPoint*sin(angle) + y0*cos(angle);
    
    slope1 = ( sin(angle)+slope1*cos(angle) ) / (cos(angle)-slope1*sin(angle));

    slope2 = ( sin(angle)+slope2*cos(angle) )/ (cos(angle)-slope2*sin(angle));

    std::cout<<" Break point: "<<breakPoint<<" y0: "<<y0<<" slope1: "<<slope1<<" slope2: "<<slope2<<std::endl;


    // Draw the fitted function
    xmin = *std::min_element(x.begin(), x.end());
    xmax = *std::max_element(x.begin(), x.end());
    TF1 *line1 = new TF1("line1","[0]+[1]*x", xmin, xmax);
    double q1 = y0 - slope1*breakPoint;
    line1->SetParameters(q1, slope1);
    line1->SetLineColor(kOrange+7);
    line1->SetLineWidth(2);
    line1->Draw("same");
    TF1 *line2 = new TF1("line2","[0]+[1]*x", xmin, xmax);
    double q2 = y0 - slope2*breakPoint;
    line2->SetParameters(q2, slope2);
    line2->SetLineColor(kGreen+2);
    line2->SetLineWidth(2);
    line2->Draw("same");

    // Draw m/q line
    TF1 *line3 = new TF1("line3","[0]+[1]*x", xmin, xmax);
    double q = y1-yVertex - m*(x1-xVertex);
    line3->SetParameters(q, m);
    line3->SetLineColor(kPink+2);
    line3->SetLineWidth(2);
    line3->Draw("same");

    
    c1->Update();

    // Print the results
    std::cout << "Fit results:" << std::endl;
    std::cout << "breakPoint: " << brokenLine->GetParameter(0) << std::endl;

}

double CheckCrossHits(TH1F *h1, TH1F *h2, TGraph *gr1, TGraph *gr2){

    // Total number of points
    int nPoints = gr1->GetN();
    nPoints+=gr2->GetN();

    int nPointsCross = 0;
    // Check points in graph 1 compatible with band in h1
    for(int i=0; i<gr1->GetN(); i++){
        double x = gr1->GetX()[i];
        double y = gr1->GetY()[i];
        
        double yHist = h2->GetBinContent(h2->FindBin(x));
        double yHistError = h2->GetBinError(h2->FindBin(x));
        double yBandLow = yHist-yHistError;
        double yBandHigh = yHist+yHistError;
        
        if(y>yBandLow && y<yBandHigh){
            nPointsCross++;
        }
    }

    // Check points in graph 2 compatible with band in h2
    for(int i=0; i<gr2->GetN(); i++){
        double x = gr2->GetX()[i];
        double y = gr2->GetY()[i];
        
        double yHist = h1->GetBinContent(h1->FindBin(x));
        double yHistError = h1->GetBinError(h1->FindBin(x));
        double yBandLow = yHist-yHistError;
        double yBandHigh = yHist+yHistError;
        
        if(y>yBandLow && y<yBandHigh){
            nPointsCross++;
        }
    }

    return (double)nPointsCross/(double)nPoints;
}

// --- Function to check overlap between two contours plots ---
double CheckOverlapRegion(TH1F *h1, TH1F *h2, double x1Max, double x2Max){
    double overlapArea = 0;

    double x1Min = h1->GetXaxis()->GetXmin();
    double x2Min = h2->GetXaxis()->GetXmin();
    
    double x1BinWidth = h1->GetXaxis()->GetBinWidth(1);
    double x2BinWidth = h2->GetXaxis()->GetBinWidth(1);

    double xMax = std::min(x1Max, x2Max);
    // First get total area of the two histograms in the range min-max
    // Now get the overlap
    double x1 = x1Min;
    double x2 = x2Min;    

    double area1 = 0;
    double area2 = 0;
    while(x1<xMax){
        area1+=2*h1->GetBinError(h1->FindBin(x1))*x1BinWidth;
        x1+=x1BinWidth;
    }
    while(x2<xMax){
        area2+=2*h2->GetBinError(h2->FindBin(x2))*x2BinWidth;
        x2+=x2BinWidth;
    }

    // Now get the overlap
    x1 = x1Min;
    x2 = x2Min;    
    while(x1<xMax && x2<xMax){
        // Get y values with errors
        double y1Error = h1->GetBinError(h1->FindBin(x1));
        double y1 = h1->GetBinContent(h1->FindBin(x1));
        double y2Error = h2->GetBinError(h2->FindBin(x2));
        double y2 = h2->GetBinContent(h2->FindBin(x2));

        double y1Low = y1-y1Error;
        double y1High = y1+y1Error;
        double y2Low = y2-y2Error;
        double y2High = y2+y2Error;

        // Calculate the overlap integral
        double overlap = ( std::min(y1High, y2High) - std::max(y1Low, y2Low) ) * x1BinWidth;
        if (overlap < 0.0) overlap = 0.0;

        x1+=x1BinWidth;
        x2+=x2BinWidth;
        overlapArea+=2*overlap; // overlap counts to the two ranges
    }
    
    return overlapArea/(area1+area2);
}


// --- Function to get confidence interval ---
bool GetFitConfidenceInterval(TH1F *h, TF1 *f, TFitResultPtr fitResult, double confidenceLevel, double maxX){

    if(fitResult->Status()!=0) return false;
    else{
        // Get range from histogram
        int nPoints = h->GetNbinsX();
        
        std::vector<double> xPoints;
        for(int i=1; i<=nPoints; i++){
            double x = h->GetBinCenter(i);
            xPoints.push_back(x);
        }
        std::vector<double> confInterval(nPoints);
        
        fitResult->GetConfidenceIntervals(nPoints, 1, 1, &xPoints[0], &confInterval[0], confidenceLevel);

        // Fill the histogram
        for(int i=1; i<=nPoints; i++){
            if(xPoints[i]<maxX){
                h->SetBinContent(i, f->Eval(xPoints[i-1]));
                h->SetBinError(i, confInterval[i]);
            }
            else{
                h->SetBinContent(i, 0);
                h->SetBinError(i, 0);
            }
        }

        return true;
    }
   
}

// --- Function to display dEdx ---
void STriangleCalo::MakeEnergyLossVsResidualRangePlot(SCalo calo1, SCalo calo2, TPad *pad) {
    
    if(pad!=nullptr){
        pad->cd();
        gStyle->SetOptStat(0);
    }
    double fCondifenceLevel = 0.95;
    
    // Get vectors
    const std::vector<double> residualRange1 = calo1.GetResidualRange();
    const std::vector<double> depositedEnergy1 = calo1.GetDepositedEnergy();
    const std::vector<double> residualRange2 = calo2.GetResidualRange();
    const std::vector<double> depositedEnergy2 = calo2.GetDepositedEnergy();

    // Define graph 1
    TGraph* graph1 = new TGraph(residualRange1.size(), &residualRange1[0], &depositedEnergy1[0]);
    graph1->SetTitle("Track 1");
    graph1->SetLineColor(fColor1);
    graph1->SetLineStyle(kDashed);
    graph1->SetMarkerColor(fColor1);
    graph1->SetMarkerStyle(20);

    // Define graph 2
    TGraph* graph2 = new TGraph(residualRange2.size(), &residualRange2[0], &depositedEnergy2[0]);
    graph2->SetTitle("Track 2");
    graph2->SetLineColor(fColor2);
    graph2->SetMarkerColor(fColor2);
    graph2->SetLineStyle(kDashed);
    graph2->SetMarkerStyle(20);

    // Frame
    double maxLength = std::max(calo1.GetTrackLength(), calo2.GetTrackLength());   
    double maxY1 = 0, maxY2 = 0;
    if(depositedEnergy1.size()>1)
        maxY1 = *std::max_element(depositedEnergy1.begin(), depositedEnergy1.end());
    if(depositedEnergy2.size()>1)
        maxY2 = *std::max_element(depositedEnergy2.begin(), depositedEnergy2.end()); 
    double maxY = std::max(maxY1, maxY2);
    TH2F *hFrame = new TH2F("hFrame", ";Residual range [cm];dQ/dx [AU]", 200, 0, 1.1*maxLength, 100, 0, 1.1*maxY);
    hFrame->GetYaxis()->SetTitleOffset(.8);
    hFrame->GetXaxis()->SetTitleOffset(.8);
    hFrame->SetStats(0);
    hFrame->Draw();
    graph1->Draw("lp same");
    graph2->Draw("lp same");

    // TLegend
    TLegend *leg = new TLegend(0.7, 0.5, 0.85, 0.65);
    leg->AddEntry(graph1, "Track1", "lp");
    leg->AddEntry(graph2, "Track2", "lp");
    leg->Draw("same");
    
    // Fit results
    bool passFit1 = false;
    bool passFit1Exp = false;
    bool passFit2 = false;
    bool passFit2Exp = false;
    // Confidence interval histograms
    TH1F *hint1 = new TH1F("hint1", "Fitted Gaussian with .95 conf.band", 100, 0, maxLength);
    TH1F *hint1Exp = new TH1F("hint1Exp", "Fitted Gaussian with .95 conf.band", 100, 0, maxLength);
    TH1F *hint2 = new TH1F("hint2", "Fitted Gaussian with .95 conf.band", 100, 0, maxLength);
    TH1F *hint2Exp = new TH1F("hint2Exp", "Fitted Gaussian with .95 conf.band", 100, 0, maxLength);
    hint1->SetFillColorAlpha(fColor1,0.5);
    hint1Exp->SetFillColorAlpha(fColor1,0.5);
    hint2->SetFillColorAlpha(fColor2,0.5);
    hint2Exp->SetFillColorAlpha(fColor2,0.5);
    TH1F *h1 = new TH1F("h2", "h2", 100, 0, maxLength);
    TH1F *h2 = new TH1F("h1", "h1", 100, 0, maxLength);
    // Fit functions
    TF1 *fit1 = new TF1("fit1", "[0]+[1]*x", 0, calo1.GetTrackLength());
    TF1 *fit1Exp = new TF1("fit1Exp", "[0]*exp(x*[1])+[2]", 0, calo1.GetTrackLength());
    TF1 *fit2 = new TF1("fit2", "[0]+[1]*x", 0, calo2.GetTrackLength());
    TF1 *fit2Exp = new TF1("fit2Exp", "[0]*exp(x*[1])+[2]", 0, calo2.GetTrackLength());

    // --- First track fitting
    if(graph1->GetN()>0){
        // Constant
        fit1->SetLineColor(fColor1);
        fit1->SetLineWidth(2);
        TFitResultPtr fitResult1 = graph1->Fit(fit1, "RSQ");
        passFit1 = GetFitConfidenceInterval(hint1, fit1, fitResult1, fCondifenceLevel, calo1.GetTrackLength());
        
        // Exponential
        fit1Exp->SetLineColor(fColor1);
        fit1Exp->SetLineWidth(2);
        fit1Exp->SetParameters(maxY, -1, fit1->GetParameter(0));
        fit1Exp->SetParLimits(1, -1e6, 0);
        fit1Exp->SetParLimits(2, 0, 1e6);
        TFitResultPtr fitResult1Exp = graph1->Fit(fit1Exp, "BRSQ");    
        passFit1Exp = GetFitConfidenceInterval(hint1Exp, fit1Exp, fitResult1Exp, fCondifenceLevel, calo1.GetTrackLength());
    }    

    // --- Second track fitting
    if(graph2->GetN()>0){ 
        // Constant
        fit2->SetLineColor(fColor2);
        fit2->SetLineWidth(2);
        TFitResultPtr fitResult2 = graph2->Fit(fit2, "RSQ");
        passFit2 = GetFitConfidenceInterval(hint2, fit2, fitResult2, fCondifenceLevel, calo2.GetTrackLength());

        // Exponential
        fit2Exp->SetLineColor(fColor2);
        fit2Exp->SetLineWidth(2);
        fit2Exp->SetParameters(maxY, -1, fit2->GetParameter(0));
        fit2Exp->SetParLimits(1, -1e6, 0);
        fit2Exp->SetParLimits(2, 0, 1e6);
        TFitResultPtr fitResult2Exp = graph2->Fit(fit2Exp, "BRSQ");
        passFit2Exp = GetFitConfidenceInterval(hint2Exp, fit2Exp, fitResult2Exp, fCondifenceLevel, calo2.GetTrackLength());
    
    }

    // --- Draw the fitted functions
    // std::cout<<" *** Fits status: "<<passFit1<<" "<<passFit1Exp<<" "<<passFit2<<" "<<passFit2Exp<<std::endl;
    bool pass1 = passFit1 || passFit1Exp;
    bool pass2 = passFit2 || passFit2Exp;

    double p0Fit_1 = -1;
    double p0Fit_2 = -1;

    if(pass1){
        fit1->Draw("same");
        fit1Exp->Draw("same");

        //std::cout<<" *** "<<fit1Exp->GetParameter(0)<<" "<<fit1Exp->GetParameter(1)<<" "<<hint1Exp->Integral()<<std::endl;
        if(fit1Exp->GetParameter(0)>0 && fit1Exp->GetParameter(1)<0 && hint1Exp->Integral()>0){
            hint1Exp->Draw("e3 same");
            p0Fit_1 = fit1Exp->GetParameter(2);
            // copy to h1
            for(int i=1; i<=hint1Exp->GetNbinsX(); i++){
                h1->SetBinContent(i, hint1Exp->GetBinContent(i));
                h1->SetBinError(i, hint1Exp->GetBinError(i));
            }
        }
        else{
            hint1->Draw("e3 same");
            p0Fit_1 = fit1->GetParameter(0);
            // copy to h1
            for(int i=1; i<=hint1->GetNbinsX(); i++){
                h1->SetBinContent(i, hint1->GetBinContent(i));
                h1->SetBinError(i, hint1->GetBinError(i));
            }
            if( fit1->GetParameter(1)-fit1->GetParError(1) > 0 ) pass1 = false;
        }

    }
    if(pass2){
        fit2->Draw("same");
        fit2Exp->Draw("same");

        if(fit2Exp->GetParameter(0)>0 && fit2Exp->GetParameter(1)<0 && hint2Exp->Integral()>0){
            hint2Exp->Draw("e3 same");
            p0Fit_2 = fit2Exp->GetParameter(2);
            // copy to h2
            for(int i=1; i<=hint2Exp->GetNbinsX(); i++){
                h2->SetBinContent(i, hint2Exp->GetBinContent(i));
                h2->SetBinError(i, hint2Exp->GetBinError(i));
            }
        }
        else{
            hint2->Draw("e3 same");
            p0Fit_2 = fit2->GetParameter(0);
            // copy to h2
            for(int i=1; i<=hint2->GetNbinsX(); i++){
                h2->SetBinContent(i, hint2->GetBinContent(i));
                h2->SetBinError(i, hint2->GetBinError(i));
            }
            if( fit2->GetParameter(1)-fit2->GetParError(1) > 0 ) pass2 = false;
        }
    }

    double charge1=0, charge2=0;
    if(depositedEnergy1.size()>1)
        charge1 = std::accumulate(depositedEnergy1.begin(), depositedEnergy1.end(), 0.0) / depositedEnergy1.size();
    if(depositedEnergy2.size()>1)
        charge2 = std::accumulate(depositedEnergy2.begin(), depositedEnergy2.end(), 0.0) / depositedEnergy2.size();
    double charge1Average = std::max(charge1, charge2);
    double charge2Average = std::min(charge1, charge2);
    if(charge2Average>0 && charge1Average>0)
        fChargeRatioAverage = charge1Average / charge2Average;
    else
        fChargeRatioAverage = 0;
    if(fChargeRatioAverage>fRatioUpperLimit) fChargeRatioAverage = fRatioUpperLimit;
    fChargeDifferenceAverage = charge2Average - charge1Average;
    fChargeRelativeDifferenceAverage = fChargeDifferenceAverage / charge1Average;


    // Integral method near the end track
    double fEndTrackLength = 4.;
    int bin1_End = h1->FindBin(3.);
    int bin2_End = h2->FindBin(3.);
    double length1 = fEndTrackLength;
    double length2 = fEndTrackLength;
    if(calo1.GetTrackLength()<fEndTrackLength) length1 = calo1.GetTrackLength();
    if(calo2.GetTrackLength()<fEndTrackLength) length2 = calo2.GetTrackLength();
    
    double charge1Integral = 0;
    double charge2Integral = 0;
    if(pass1) charge1Integral = h1->Integral(1, bin1_End, "width");
    if(pass2) charge2Integral = h2->Integral(1, bin2_End, "width");
    
    if(length1>0)
        charge1Integral/=length1;
    else
        charge1Integral = 0;
    if(length2>0)
        charge2Integral/=length2;
    else
        charge2Integral = 0;
    double charge1IntegralMax = std::max(charge1Integral, charge2Integral);
    double charge2IntegralMin = std::min(charge1Integral, charge2Integral);
    if(charge2IntegralMin>0 && charge1IntegralMax>0)
        fChargeRatioIntegral = charge1IntegralMax / charge2IntegralMin;
    else
        fChargeRatioIntegral = 0;
    if(fChargeRatioIntegral>fRatioUpperLimit) fChargeRatioIntegral = fRatioUpperLimit;
    fChargeDifferenceIntegral = charge2IntegralMin - charge1IntegralMax;
    fChargeRelativeDifferenceIntegral = fChargeDifferenceIntegral / charge1IntegralMax;


    fPassChargeFit = pass1 && pass2;
    if(fPassChargeFit){
        double charge1Fit = std::max(p0Fit_1, p0Fit_2);
        double charge2Fit = std::min(p0Fit_1, p0Fit_2);

        if(charge2Fit!=0)
            fChargeRatioFit = charge1Fit / charge2Fit;
        else
            fChargeRatioFit = 0;
        if(fChargeRatioFit>fRatioUpperLimit) fChargeRatioFit = fRatioUpperLimit;
        fChargeDifferenceFit = charge2Fit - charge1Fit;
        fChargeRelativeDifferenceFit = fChargeDifferenceFit / charge1Fit;

        // Check the overlap
        double overlap = CheckOverlapRegion(h1, h2, calo1.GetTrackLength(), calo2.GetTrackLength());
        // check overlap is not a nan
        if(isnan(std::abs(overlap))) overlap = 1;

        fBandOverlap = overlap;

        double crossHitsFraction = CheckCrossHits(h1, h2, graph1, graph2);
        // check overlap is not a nan
        if(isnan(std::abs(crossHitsFraction))) crossHitsFraction = 1;

        fBandCrossHits = crossHitsFraction;

    }
    else{
        fChargeRatioFit = 1;
        fChargeDifferenceFit = 0;
        fChargeRelativeDifferenceFit = 0;
        fBandOverlap = 1;
        fBandCrossHits = 1;
    }

}


void STriangleCalo::CreateEnergyLossVsResidualRangePlot() {
    TCanvas *c1 = new TCanvas("canvasCalo","Calorimetry", 600,1200);
    c1->cd();
    TPad *pad1 = new TPad("pad1","This is pad1",0.02,0.02,0.98,0.98);
    MakeEnergyLossVsResidualRangePlot(fCalo1, fCalo2, pad1);
    delete c1;
    delete pad1;
}


// --- Sort hits by distance ---
void SortSHitVectorByDistance(std::vector<SHit>& hits, double externalX, double externalY) {
    // Lambda expression for comparing SHit objects based on distance
    auto distanceComparator = [externalX, externalY](const SHit& a, const SHit& b) {
        // Function to calculate the distance between two points
        auto calculateDistance = [](const SHit& point, double externalX, double externalY) {
            double deltaX = point.X() - externalX;
            double deltaY = 0.25*(point.Y() - externalY);
            return std::sqrt(deltaX * deltaX + deltaY * deltaY);
        };

        double distanceA = calculateDistance(a, externalX, externalY);
        double distanceB = calculateDistance(b, externalX, externalY);
        return distanceA < distanceB;
    };

    // Sort the vector using the distanceComparator lambda
    std::sort(hits.begin(), hits.end(), distanceComparator);
}

// --- Get mean and std deviation of vector of hits ---
void GetStatisticsHitVector(std::vector<SHit> hits, double & mean, double & stddev){
    double sum = 0;
    double sum2 = 0;
    for(SHit & h:hits){
        sum+=h.Integral();
        sum2+=h.Integral() * h.Integral();
    }
    sum/=hits.size();
    sum2/=hits.size();
    mean = sum;
    stddev = std::sqrt(sum2-sum*sum);

    return;

}

// --- Right/left hits with respect to vertex ---
std::vector<SHit> GetSideHits(std::vector<SHit> hitV, double xVertex, std::vector<SHit> & residualHits){

    std::vector<SHit> vLeft, vRight;
    for(SHit & h:hitV){
        if(h.X()>=xVertex) vRight.push_back(h);
        else  vLeft.push_back(h);
    }

    if(vLeft.size()>vRight.size()) {
        residualHits.clear();
        residualHits = vRight;
        return vLeft;
    }
    else{
        residualHits.clear();
        residualHits = vLeft;
        return vRight;
    }
}

// --- Print results function ---
void STriangleCalo::PrintResults(){
    std::cout << "PassFit: " << fPassFit <<std::endl;
    std::cout << "NVertexHits: " << fNVertexHits <<std::endl;
    std::cout << "NBulkHits: " << fNBulkHits <<std::endl;
    std::cout << "VertexHitIntegralRatio: "<<fVertexHitIntegralRatio<<std::endl;
    std::cout << "VertexHitIntegralDifference: "<<fVertexHitIntegralDifference<<std::endl;
    std::cout << "VertexHitIntegralRelDifference: "<<fVertexHitIntegralRelativeDifference<<std::endl;
    std::cout << "--- TrackLengthVariables: "<<std::endl;
    std::cout << "TrackLengthRatio: "<<fTrackLengthRatio<<std::endl;
    std::cout << "Track 1 length: "<<fTrackLength1<<std::endl;
    std::cout << "Track 2 length: "<<fTrackLength2<<std::endl;
    std::cout << "Track 1 residual range RMS: "<<fResidualRange1RMS<<std::endl;
    std::cout << "Track 2 residual range RMS: "<<fResidualRange2RMS<<std::endl;
    std::cout << "Min track residual range RMS: "<<fResidualRangeMinRMS<<std::endl;
    std::cout << "--- Charge variables: "<<std::endl;
    std::cout << "Charge ratio (average): " << fChargeRatioAverage << std::endl;
    std::cout << "Charge difference (average): " << fChargeDifferenceAverage << std::endl;
    std::cout << "Charge relative difference (average): " << fChargeRelativeDifferenceAverage << std::endl;
    std::cout << "Charge ratio (fit) variables: "<<std::endl;
    std::cout << "Pass charge fit: " << fPassChargeFit << std::endl;
    std::cout << "Band overlap: " << fBandOverlap << std::endl;
    std::cout << "Band cross hits: " << fBandCrossHits << std::endl;
    std::cout << "Charge ratio (fit): " << fChargeRatioFit << std::endl;
    std::cout << "Charge difference (fit): " << fChargeDifferenceFit << std::endl;
    std::cout << "Charge relative difference (fit): " << fChargeRelativeDifferenceFit << std::endl;
    std::cout << "Charge ratio (integral): " << fChargeRatioIntegral << std::endl;
    std::cout << "Charge difference (integral): " << fChargeDifferenceIntegral << std::endl;
    std::cout << "Charge relative difference (integral): " << fChargeRelativeDifferenceIntegral << std::endl;

    return;
}

double lineDistance(double m, double n, double x, double y){
    return std::abs( -m*x+y-n )/std::sqrt(m*m+1);        
}


// --- JointFit analysis---
void STriangleCalo::JointFitAnalysis(std::vector<double> *pProton, std::vector<double> *pPion){

    unsigned int maxHits = fCaloPset.MaxTrackHits;
    // --- Get the hits and the vertex ---x
    // Vertex
    fXVertex = fTriangle.GetMainVertex().X();
    fYVertex = fTriangle.GetMainVertex().Y();

    // Hit coordinates for the first maxHits hits
    std::vector<double> xV, yV, widthV;
    size_t n = 0;

    // Get tracks start/end slope based on distance to main vertex
    SLinearCluster track1 = fTriangle.GetTrack1();
    SLinearCluster track2 = fTriangle.GetTrack2();
    SPoint mainVertex = fTriangle.GetMainVertex();
    double m1, n1, m2, n2;

    
    if (std::abs(mainVertex.X() - track1.GetMinX()) < std::abs(mainVertex.X()  - track1.GetMaxX())) {
        m1 = track1.GetTrackEquationStart().Slope();
        n1 = track1.GetTrackEquationStart().Intercept();
    }
    else{
        m1 = track1.GetTrackEquationEnd().Slope();
        n1 = track1.GetTrackEquationEnd().Intercept();
    }
    if (std::abs(mainVertex.X() - track2.GetMinX()) < std::abs(mainVertex.X()  - track2.GetMaxX())) {
        m2 = track2.GetTrackEquationStart().Slope();
        n2 = track2.GetTrackEquationStart().Intercept();
    }
    else{
        m2 = track2.GetTrackEquationEnd().Slope();
        n2 = track2.GetTrackEquationEnd().Intercept();
    }


    for(SHit &h:fTriangle.GetTrack1().GetHits()){
        if(n>=maxHits) break;
        xV.push_back(h.X()-fXVertex);
        yV.push_back(h.Y()-fYVertex);
        widthV.push_back(h.Width());
        fFitHits.push_back(h);
        // if compatible with the two tracks, discard
        if( std::abs( h.Y() - m2*h.X() - n2 ) > h.Width() )
            n++;
    }
    n=0;
    for(SHit &h:fTriangle.GetTrack2().GetHits()){
        if(n>=maxHits) break;
        xV.push_back(h.X()-fXVertex);
        yV.push_back(h.Y()-fYVertex);
        widthV.push_back(h.Width());
        fFitHits.push_back(h);
        if( std::abs( h.Y() - m1*h.X() - n1 ) > h.Width() )
            n++;
    }

    // --- Get frame values ---
    fMinX = *std::min_element(xV.begin(), xV.end());
    fMaxX = *std::max_element(xV.begin(), xV.end());
    fMinY = *std::min_element(yV.begin(), yV.end());
    fMaxY = *std::max_element(yV.begin(), yV.end());

    // --- Create the TGraph ---
    n = xV.size();
    TGraphErrors* gr = new TGraphErrors(n, &xV[0], &yV[0], 0, &widthV[0]);
    gr->SetMarkerStyle(3);

    if(fCaloPset.MakeFit){
        // --- Chi2 minimizer function for two lines with common vertex ---
        auto Chi2TwoLines = [&](const double *par) {
            // minimisation function computing the sum of squares of residuals
            int np = gr->GetN();
            double f = 0;
            double *x = gr->GetX();
            double *y = gr->GetY();
            double *error = gr->GetEY();

            // loop over the data points
            for (int i=0;i<np;i++) {
                // get values
                double xH = x[i];
                double yH = y[i];
                double errorH =  fCaloPset.UseHitRMSInFit? error[i]:1;

                // get residuals to the two lines
                //double dr1 = ( yH - ( par[0]*xH + par[1] ) ) / errorH ;
                //double dr2 = ( yH - ( par[2]*xH + (par[0]-par[2])*par[3] + par[1] ) ) / errorH;
                double dr1 = lineDistance(par[0], par[1], xH, yH) / errorH;
                double dr2 = lineDistance(par[2], (par[0]-par[2])*par[3], xH, yH) / errorH;

                //dr1 = dr1*dr1;
                //dr2 = dr2*dr2;
                // add contribution to the chi2 (minimum residual to the two lines)
                double dr = std::min(dr1,dr2);
                f += dr;
            }

            return f;
        };

        // --- Fit the graph to the function ---
        ROOT::Math::Functor fcn(Chi2TwoLines, 4);
        ROOT::Fit::Fitter  fitter;
        // Initializations
        double pStart[4] = {fTriangle.GetTrack1().GetTrackEquation().Slope(),
                            0,
                            fTriangle.GetTrack2().GetTrackEquation().Slope(),
                            0};
        fitter.SetFCN(fcn, pStart);
        fitter.Config().ParSettings(0).SetName("m1");
        fitter.Config().ParSettings(1).SetName("n1");
        fitter.Config().ParSettings(2).SetName("m2");
        fitter.Config().ParSettings(3).SetName("xB");

        // Do the fit
        bool ok = fitter.FitFCN();
        if (!ok) {
            fPassFit = false;
        }
        else{
            fPassFit = true;
        }

        // --- Fit results ---
        const ROOT::Fit::FitResult & result = fitter.Result();
        fM1 = result.Parameter(0);
        fN1 = result.Parameter(1);
        fM2 = result.Parameter(2);
        fXB = result.Parameter(3);
        fN2 = (fM1-fM2)*fXB+fN1;
        fXFitVtx = fXB;
        fYFitVtx = fM1*fXFitVtx+fN1;
        fTwoLinesChi2 = result.MinFcnValue();
        result.Print(std::cout);
    }
    else{
        fPassFit = true;
        fM1 = fTriangle.GetTrack1().GetTrackEquation().Slope();
        fN1 = 0;
        fM2 = fTriangle.GetTrack2().GetTrackEquation().Slope();
        fN2 = 0;
        fXB = (fN2-fN1)/(fM1-fM2);
        fXFitVtx = fXB;
        fYFitVtx = fM1*fXFitVtx+fN1;
        fTwoLinesChi2 = 0;
    }


    // --- Assing each hit to a track ---
    fMaxIntegral = 0;
    fHitsTrack1.clear();
    fHitsTrack2.clear();
    fVertexHits.clear();
    std::vector<SHit> hitList = fTriangle.GetAllHits();
    fAllHits = hitList;
    double fYCorrFactor = 0.075/0.3;
    for(SHit &h:hitList){

        // shifted values 
        double x = h.X()-fXVertex;
        double y = h.Y()-fYVertex;

        if(x<fMinX) fMinX = x;
        if(y<fMinY) fMinY = y;
        if(x>fMaxX) fMaxX = x;
        if(y>fMaxY) fMaxY = y;

        // values from fitted slopes
        double yH1 = fM1*x+fN1;
        double yH2 = fM2*x+fN2;

        // closest tracks
        bool isTrack1 = std::abs( y-yH1 )<fCaloPset.HitWidthTolInFit*h.Width();
        bool isTrack2 = std::abs( y-yH2 )<fCaloPset.HitWidthTolInFit*h.Width();

        // if shared hit, split the integral in two hits
        if(isTrack1 && isTrack2){
            //SHit sharedHit(-1, h.X(), h.Y(), h.Width(), 0.5*h.Integral());
            //fHitsTrack1.push_back(sharedHit);
            //fHitsTrack2.push_back(sharedHit);
            fVertexHits.push_back(h);
        }
        else if(isTrack1){
            fHitsTrack1.push_back(h);
        }
        else if(isTrack2){
            fHitsTrack2.push_back(h);
        }
        else{
            // if not width compatible, add to the closest
            double d1 = lineDistance(fYCorrFactor*fM1, fYCorrFactor*fN1, x, fYCorrFactor*y);
            double d2 = lineDistance(fYCorrFactor*fM2, fYCorrFactor*fN2, x, fYCorrFactor*y);
            d1<d2? fHitsTrack1.push_back(h):fHitsTrack2.push_back(h);
        }

        // update maximum integral
        if(h.Integral()>fMaxIntegral) fMaxIntegral = h.Integral();
    }

    // --- Remove outliers hits ---
    std::vector<SHit> resHits;
    fHitsTrack1 = GetSideHits(fHitsTrack1, fXVertex+fXFitVtx, resHits);
    if(resHits.size()>0){
        fResidualHits.insert(fResidualHits.end(), resHits.begin(), resHits.end());
        resHits.clear();
    }
    fHitsTrack2 = GetSideHits(fHitsTrack2, fXVertex+fXFitVtx, resHits);
    if(resHits.size()>0){
        fResidualHits.insert(fResidualHits.end(), resHits.begin(), resHits.end());
        resHits.clear();
    }
    // Average track 1 hit intetgral
    double meanIntegral1 = 0, stdDev1 = 0;
    double meanIntegral2 = 0, stdDev2 = 0;
    GetStatisticsHitVector(fHitsTrack1, meanIntegral1, stdDev1);
    GetStatisticsHitVector(fHitsTrack2, meanIntegral2, stdDev2);

    double meanIntegralShared = 0, stdDevShared = 0;
    GetStatisticsHitVector(fVertexHits, meanIntegralShared, stdDevShared);
    std::vector<SHit> auxVertexHits = fVertexHits;
    fVertexHits.clear();
    for(SHit & h:auxVertexHits){
        fVertexHits.push_back(h);
    }

    // --- Calorimetry objects ---
    // Sort by distance to fitted vertex ---
    SortSHitVectorByDistance(fHitsTrack1, fXVertex+fXFitVtx, fYVertex+fYFitVtx);
    SortSHitVectorByDistance(fHitsTrack2, fXVertex+fXFitVtx, fYVertex+fYFitVtx);
    SortSHitVectorByDistance(fVertexHits, fXVertex+fXFitVtx, fYVertex+fYFitVtx);

    //fHitsTrack1.pop_back(); fHitsTrack2.pop_back();
    // Create the objects
    double yAngle1 = 0;
    double yAngle2 = 0;
    if(pProton!=nullptr){
        if(pProton->size()>0)
            yAngle1 = std::acos( pProton->at(1) );
    }
    if(pPion!=nullptr){
        if(pPion->size()>0)
            yAngle2 = std::acos( pPion->at(1) );
    }
    SCalo calo1(fHitsTrack1, std::atan(fM1), yAngle1);
    SCalo calo2(fHitsTrack2, std::atan(fM2), yAngle2);
    fCalo1 = calo1;
    fCalo2 = calo2;

    // Angle deviation
    std::vector<double> anglesTrack1, anglesTrack2;
    if(fHitsTrack1.size()>=2){
        for(size_t k=0; k<fHitsTrack1.size()-1; k++){
            SHit h1 = fHitsTrack1[k];
            SHit h2 = fHitsTrack1[k+1];
            double angle = 180 * std::atan2((h2.Y()-h1.Y()), (h2.X()-h1.X())) / TMath::Pi();
            anglesTrack1.push_back(angle);
        }
    }
    
    if(fHitsTrack2.size()>=2){
        for(size_t k=0; k<fHitsTrack2.size()-1; k++){
            SHit h1 = fHitsTrack2[k];
            SHit h2 = fHitsTrack2[k+1];
            double angle = 180 * std::atan2((h2.Y()-h1.Y()), (h2.X()-h1.X())) / TMath::Pi();
            anglesTrack2.push_back(angle);
        }
    }

    std::vector<double> anglesDiffTrack1;
    std::vector<double> anglesDiffTrack2;
    if(anglesTrack1.size()>=2){
        for(size_t k=0; k<anglesTrack1.size()-1; k++){
            anglesDiffTrack1.push_back( anglesTrack1[k+1]-anglesTrack1[k] );
        }
    }
    if(anglesTrack2.size()>=2){
        for(size_t k=0; k<anglesTrack2.size()-1; k++){
            anglesDiffTrack2.push_back( anglesTrack2[k+1]-anglesTrack2[k] );
        }
    }

    fResidualRange1AngleRMS = CalculateStdDev(anglesDiffTrack1);
    fResidualRange2AngleRMS = CalculateStdDev(anglesDiffTrack2);
    fResidualRangeMinAngleRMS = std::min(fResidualRange1AngleRMS, fResidualRange2AngleRMS);
    fResidualRangeMaxAngleRMS = std::max(fResidualRange1AngleRMS, fResidualRange2AngleRMS);

    
    
    MakeEnergyLossVsResidualRangePlot(fCalo1, fCalo2, nullptr);
   
    // --- Return values ---
    fMinX-=0.1*(fMaxX-fMinX);
    fMaxX+=0.1*(fMaxX-fMinX);
    fMinY-=0.1*(fMaxY-fMinY);
    fMaxY+=0.1*(fMaxY-fMinY);

    // vertex charge vs track charge
    SortSHitVectorByDistance(hitList, fXVertex+fXFitVtx, fYVertex+fYFitVtx);
    double vertexHitsIntegral = 0;
    double bulkHitsIntegral = 0;
    double maxVertexCharge = 0;
    int nVertexHits = 0;
    int nBulkHits = 0;
    // loop over the hits
    for(SHit &h:fVertexHits){
        vertexHitsIntegral+=h.Integral();
        if(h.Integral()>maxVertexCharge) maxVertexCharge = h.Integral();
    }
    nVertexHits = fVertexHits.size();
    int nCounter=0;
    for(SHit &h:fHitsTrack1){
        if(nCounter>(int)fHitsTrack1.size()/2) continue;
        bulkHitsIntegral+=h.Integral();
        nBulkHits++;
        nCounter++;
    }
    nCounter=0;
    for(SHit &h:fHitsTrack2){
        if(nCounter>(int)fHitsTrack2.size()/2) continue;
        bulkHitsIntegral+=h.Integral();
        nBulkHits++;
        nCounter++;
    }

    //double vertexCharge = vertexHitsIntegral/nVertexHits;
    double vertexCharge = maxVertexCharge;
    double bulkCharge = bulkHitsIntegral/nBulkHits;
    
    fNVertexHits = nVertexHits;
    fNBulkHits = nBulkHits;
    if(nVertexHits>0 and nBulkHits>0){
        fVertexHitIntegralRatio = vertexCharge/bulkCharge;
        fVertexHitIntegralDifference = vertexCharge-bulkCharge;
        fVertexHitIntegralRelativeDifference = fVertexHitIntegralDifference/bulkCharge;
        if(fVertexHitIntegralRatio>fRatioUpperLimit) fVertexHitIntegralRatio = fRatioUpperLimit;
    }
    else{
        fVertexHitIntegralRatio = 0;
        fVertexHitIntegralDifference = 0;
        fVertexHitIntegralRelativeDifference = 0;
    }

    double trackLength1 = calo1.GetTrackLength();
    double trackLength2 = calo2.GetTrackLength();
    fTrackLength1 = std::max(trackLength1, trackLength2);
    fTrackLength2 = std::min(trackLength1, trackLength2);
    if(fTrackLength2!=0)
        fTrackLengthRatio = fTrackLength1/fTrackLength2;
    else
        fTrackLengthRatio = 0;
    fResidualRange1RMS = calo1.GetResidualRangeRMS();
    fResidualRange2RMS = calo2.GetResidualRangeRMS();
    fResidualRangeMinRMS = std::min(fResidualRange1RMS, fResidualRange2RMS);
    fResidualRangeMaxRMS = std::max(fResidualRange1RMS, fResidualRange2RMS);
    
    return;

}

// Define a function to add points to the graph and create a TGraph2D
void addPointsAndCreateGraph(TGraph2D *graph, const std::vector<SHit>& hits) {

    int n=0;
    for(size_t k = 0; k < hits.size(); k++) {
        if(hits[k].SPZ() < -100 || hits[k].SPY() == 0) continue;
        graph->SetPoint(n, hits[k].X(), hits[k].Y(), hits[k].SPY());
        n++;
    }

    //if(n==0) graph->SetPoint(0, 0, 0, 0);

    return;
}


// Define a function to draw the graph
void drawGraph(TGraph2D *graph, const char* title, Color_t color){
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(20);
    graph->SetTitle(title);
    graph->Draw("P same");
}


// Define a function to create a TH3 histogram with limits set by the limits of three graphs
TH3* createTH3WithGraphLimits(const TGraph2D* graph1, const TGraph2D* graph2, const TGraph2D* graph3) {
    // Get the minimum and maximum values for X, Y, and Z from the three graphs

    // get minimum for each graph, if points are available
    double min1X = (graph1->GetN()>=1)? graph1->GetXmin():1e3;    
    double min2X = (graph2->GetN()>=1)? graph2->GetXmin():1e3;   
    double min3X = (graph3->GetN()>=1)? graph3->GetXmin():1e3;
    double min1Y = (graph1->GetN()>=1)? graph1->GetYmin():1e3;
    double min2Y = (graph2->GetN()>=1)? graph2->GetYmin():1e3;
    double min3Y = (graph3->GetN()>=1)? graph3->GetYmin():1e3;
    double min1Z = (graph1->GetN()>=1)? graph1->GetZmin():1e3;
    double min2Z = (graph2->GetN()>=1)? graph2->GetZmin():1e3;
    double min3Z = (graph3->GetN()>=1)? graph3->GetZmin():1e3;

    double xmin = std::min({min1X, min2X, min3X});
    double ymin = std::min({min1Y, min2Y, min3Y});
    double zmin = std::min({min1Z, min2Z, min3Z});

    // get maximum for each graph, if points are available
    double max1X = (graph1->GetN()>=1)? graph1->GetXmax():-1e3;
    double max2X = (graph2->GetN()>=1)? graph2->GetXmax():-1e3;
    double max3X = (graph3->GetN()>=1)? graph3->GetXmax():-1e3;
    double max1Y = (graph1->GetN()>=1)? graph1->GetYmax():-1e3;
    double max2Y = (graph2->GetN()>=1)? graph2->GetYmax():-1e3;
    double max3Y = (graph3->GetN()>=1)? graph3->GetYmax():-1e3;
    double max1Z = (graph1->GetN()>=1)? graph1->GetZmax():-1e3;
    double max2Z = (graph2->GetN()>=1)? graph2->GetZmax():-1e3;
    double max3Z = (graph3->GetN()>=1)? graph3->GetZmax():-1e3;

    double xmax = std::max({max1X, max2X, max3X});
    double ymax = std::max({max1Y, max2Y, max3Y});
    double zmax = std::max({max1Z, max2Z, max3Z});

    // Create the TH3 histogram with limits set by the minimum and maximum values
    TH3* hist = new TH3D("hist", ";Z [cm]; X [cm]; Y [cm]", 100, xmin, xmax, 100, ymin, ymax, 100, zmin, zmax);

    return hist;
}

// --- Display function ---
void STriangleCalo::Display(TCanvas *c1){
    // --- Create the canvas ---
    if(c1==nullptr) {
        c1 = new TCanvas("canvasCalo","Calorimetry", 600,1200);
    }
    c1->cd();
    // 4 pads
    TPad *pad1 = new TPad("pad1","This is pad1",0.02,0.52,0.48,0.98);
    TPad *pad2 = new TPad("pad2","This is pad2",0.02,0.02,0.48,0.48);
    TPad *pad3 = new TPad("pad3","This is pad3",0.52,0.52,0.98,0.98);
    TPad *pad4 = new TPad("pad4","This is pad4",0.52,0.02,0.98,0.48);
    pad1->Draw();
    pad2->Draw();
    pad3->Draw();
    pad4->Draw();

    // --- Draw the graph and the fit ---
    pad1->cd();
    TH2F * hFrame = new TH2F("hFrame", ";Wire ID;Time Tick [0.5 #mu s]", 200, fMinX, fMaxX, 200, fMinY, fMaxY);
    hFrame->GetYaxis()->SetTitleOffset(.8);
    hFrame->GetXaxis()->SetTitleOffset(.8);
    hFrame->SetStats(0);

    // Line1
    double y1 = fM1*fMinX+fN1;
    double y2 = fM1*fMaxX+fN1;
    TLine *line1 = new TLine(fMinX,y1,fMaxX,y2);
    line1->SetLineColor(fColor1);
    line1->SetLineWidth(4); 

    // Line2
    y1 = fM2*fMinX+fN2;
    y2 = fM2*fMaxX+fN2;
    TLine *line2 = new TLine(fMinX,y1,fMaxX,y2);
    line2->SetLineColor(fColor2);
    line2->SetLineWidth(4);

    // Fitted vertex
    TMarker *marker = new TMarker(fXB, fM1*fXB+fN1, 29);
    marker->SetMarkerColor(kRed);
    marker->SetMarkerSize(2);

    // TLegend
    TLegend *leg = new TLegend(0.7, 0.4, 0.85, 0.65);
    leg->AddEntry(line1, "Fitted line 1", "l");
    leg->AddEntry(line2, "Fitted line 2", "l");
    leg->AddEntry(marker, "Fitted vertex", "p");

    // Create text with the fit results
    TPaveText *pt = new TPaveText(0.7, 0.7, 0.85, 0.85, "NDC");
    pt->AddText(Form("Fit result #chi^{2} = %.2f\n",fTwoLinesChi2));
    pt->AddText(Form("m1 = %.2f #pm %.2f",fM1, 0.));
    pt->AddText(Form("n1 = %.2f #pm %.2f",fN1, 0.));
    pt->AddText(Form("m2 = %.2f #pm %.2f",fM2, 0.));
    pt->AddText(Form("xB = %.2f #pm %.2f",fXB, 0.));

    // Draw
    hFrame->Draw();
    line1->Draw("same");
    line2->Draw("same");
    marker->Draw("same");
    pt->Draw("same");
    
    // --- Draw xV and yV markers after assignment ---

    // TH1F histogram around the vertex
    TH1F *h1 = new TH1F("h1", ";Residual range [cm];Energy loss [MeV/cm]", (int)(fMaxX-fMinX), fMinX, fMaxX);

    for(SHit &h:fAllHits){
        double size = (h.Integral()/fMaxIntegral) * fMaxMarkerSize;
        TMarker *marker = new TMarker(h.X()-fXVertex, h.Y()-fYVertex, 20);
        marker->SetMarkerColorAlpha(kGray, fAlpha);
        marker->SetMarkerSize(size);
        marker->Draw("same");
        h1->Fill(h.X()-fXVertex, h.Integral());
    }

    for(SHit &h:fResidualHits){
        double size = (h.Integral()/fMaxIntegral) * fMaxMarkerSize;
        TMarker *marker = new TMarker(h.X()-fXVertex, h.Y()-fYVertex, 20);
        marker->SetMarkerColorAlpha(kGreen, fAlpha);
        marker->SetMarkerSize(size);
        marker->Draw("same");
    }

    for(SHit &h:fHitsTrack1){
        double size = (h.Integral()/fMaxIntegral) * fMaxMarkerSize;
        TMarker *marker = new TMarker(h.X()-fXVertex, h.Y()-fYVertex, 22);
        marker->SetMarkerColorAlpha(fColor1, fAlpha);
        marker->SetMarkerSize(size);
        marker->Draw("same");
        if (&h == &fHitsTrack1.front()) leg->AddEntry(marker, "Track 1 hits", "p");
    }

    for(SHit &h:fHitsTrack2){
        double size = (h.Integral()/fMaxIntegral) * fMaxMarkerSize;
        TMarker *marker = new TMarker(h.X()-fXVertex, h.Y()-fYVertex, 23);
        marker->SetMarkerColorAlpha(fColor2, 1);
        marker->SetMarkerSize(size);
        marker->Draw("same");
        if (&h == &fHitsTrack2.front()) leg->AddEntry(marker, "Track 2 hits",  "p");
    }

    for(SHit &h:fVertexHits){
        double size = (h.Integral()/fMaxIntegral) * fMaxMarkerSize;
        TMarker *marker = new TMarker(h.X()-fXVertex, h.Y()-fYVertex, 20);
        marker->SetMarkerColorAlpha(kViolet-6, 1);
        marker->SetMarkerSize(size);
        marker->Draw("same");
        if (&h == &fVertexHits.front()) leg->AddEntry(marker, "Vertex hits",  "p");
    }

    for(SHit &h:fFitHits){
        double size = (h.Integral()/fMaxIntegral) * fMaxMarkerSize;
        TMarker *marker = new TMarker(h.X()-fXVertex, h.Y()-fYVertex, 20);
        marker->SetMarkerColorAlpha(kBlack, fAlpha);
        marker->SetMarkerSize(size);
        marker->SetMarkerStyle(24);
        marker->Draw("same");
    }
    
    // Draw arrows joining hits for each track
    if(fHitsTrack1.size()>=2){
        for(size_t k=0; k<fHitsTrack1.size()-1; k++){
            SHit h1 = fHitsTrack1[k];
            SHit h2 = fHitsTrack1[k+1];
            TArrow *arrow = new TArrow(h1.X()-fXVertex, h1.Y()-fYVertex, h2.X()-fXVertex, h2.Y()-fYVertex, 0.01, "|>");
            arrow->Draw();
            
        }
    }
    if(fHitsTrack2.size()>=2){
        for(size_t k=0; k<fHitsTrack2.size()-1; k++){
            SHit h1 = fHitsTrack2[k];
            SHit h2 = fHitsTrack2[k+1];
            TArrow *arrow = new TArrow(h1.X()-fXVertex, h1.Y()-fYVertex, h2.X()-fXVertex, h2.Y()-fYVertex, 0.01, "|>");
            arrow->Draw();
        }
    }

    leg->Draw("same");


    // ---- Plot 3D points
    bool make3DPlot = true;
    if(make3DPlot){
        pad3->cd();

        // Frame TH3

        TGraph2D *graphVertex = new TGraph2D();
        TGraph2D *graph1 = new TGraph2D();
        TGraph2D *graph2 = new TGraph2D();
        addPointsAndCreateGraph(graphVertex, fVertexHits);
        addPointsAndCreateGraph(graph1, fHitsTrack1);
        addPointsAndCreateGraph(graph2, fHitsTrack2);

        // Draw the frame
        TH3* h3Frame = createTH3WithGraphLimits(graphVertex, graph1, graph2);
        h3Frame->Draw();
    
        // Set marker style and color for graphs
        drawGraph(graphVertex, "Graph title", kViolet-6);
        drawGraph(graph1, "Track1", fColor1);
        drawGraph(graph2, "Track2", fColor2);
        
    }
    


    pad4->cd();
    h1->Draw();    

    MakeEnergyLossVsResidualRangePlot(fCalo1, fCalo2, pad2); 

    PrintResults();

    c1->Update();
    c1->cd();

    return;
}
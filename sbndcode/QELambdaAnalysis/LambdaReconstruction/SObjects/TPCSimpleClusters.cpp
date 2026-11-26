////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleClusters.cpp
//
// \brief Definition of SimpleTPCClusters
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCSimpleClusters.h"

SCluster::SCluster(std::vector<SHit> hitList){

    fNHits = hitList.size();

    // sort hits by X coordinate
    std::sort(hitList.begin(), hitList.end(), [](SHit& a, SHit& b) {
        return a.X() < b.X();
    });
    if (fNHits > 0) {

        for(auto &hit: hitList){
            double dminComp = DefaultMaxDistance;
            double dminConn = DefaultMaxDistance;
            double dminX = DefaultMaxDistance;
            
            for(auto &hit2: hitList){
                if (&hit == &hit2) continue;

                // distance between centers               
                double d = TPCLinesDistanceUtils::GetHitDistance(hit, hit2);
                if(d<dminComp) dminComp=d;

                // calculate connectedeness
                d = TPCLinesDistanceUtils::GetHitDistanceW(hit, hit2);
                //d = GetHitDistanceOverlap(hit, hit2)
                if(d<dminConn) dminConn=d;

                // calculate connectedness 1D
                d = TPCLinesDistanceUtils::GetHitDistance1D(hit, hit2);
                if(d<dminX) dminX=d;
            }

            SHit newHit = hit;
            newHit.SetHitConnectivity(dminComp, dminConn, dminX);

            fCompactnessV.push_back(dminComp);
            fConnectednessV.push_back(dminConn);
            fConnectedness1DV.push_back(dminX);
            fWidths.push_back(newHit.Width());
            fHitList.push_back(newHit);
        }
            
        fCompactness = TPCLinesDistanceUtils::CalculateMean(fCompactnessV);
        fCompactnessRMS = TPCLinesDistanceUtils::CalculateStdDev(fCompactnessV);
        fConnectedness = TPCLinesDistanceUtils::CalculateMean(fConnectednessV);
        fConnectednessRMS = TPCLinesDistanceUtils::CalculateStdDev(fConnectednessV);
        fConnectedness1D = TPCLinesDistanceUtils::CalculateMean(fConnectedness1DV);
        fConnectedness1DRMS = TPCLinesDistanceUtils::CalculateStdDev(fConnectedness1DV);
        fAverageWidth = TPCLinesDistanceUtils::CalculateMean(fWidths);

        CalculateMainClusterID();
    }
    else {
        fCompactness = 0;
        fCompactnessRMS = 0;
        fConnectedness = 0;
        fConnectednessRMS = 0;
        fConnectedness1D = 0;
        fConnectedness1DRMS = 0;
        fAverageWidth = 0;
        fMainClusterId = -1;
        fClusterIDs.clear();
    }
}

void SCluster::AddHit(SHit& hit) {
    fHitList.push_back(hit);
    fNHits++;
}

void SCluster::CalculateMainClusterID(){

    // map to count entries per cluster id
    std::map<int, int> clusterIdCounter;
    for (auto& hit : fHitList) {
        clusterIdCounter[hit.ClusterId()] = 0;
    }

    // count entries per cluster id
    for (auto& hit : fHitList) {
        clusterIdCounter[hit.ClusterId()]++;
    }

    // get cluster id with max entries
    int mainClusterID = -1;
    int maxEntries = 0;
    fClusterIDs.clear();
    for (auto& entry : clusterIdCounter) {
        fClusterIDs.push_back(entry.first);
        if (entry.second > maxEntries) {
            maxEntries = entry.second;
            mainClusterID = entry.first;
        }
    }

    fMainClusterId = mainClusterID;

    return;
}

std::ostream& operator<<(std::ostream& out, SCluster & c)
{
    out <<  " --- Created Simple Cluster:\n"
    << "  Comp=" << c.fCompactness << " pm " << c.fCompactnessRMS
    << "  Conn=" << c.fConnectedness << " pm " << c.fConnectednessRMS
    << "  Conn1D=" <<c.fConnectedness1D << " pm " << c.fConnectedness1DRMS;

    return out;
}

template <typename T>
double SCluster::GetMinDistanceToCluster(const T& h) {
    double minD = 1e6;
    for (auto& hit : fHitList) {
        double d =TPCLinesDistanceUtils::GetHitDistance(hit, h);
        if (d < minD) minD = d;
    }
    return minD;
}

template <typename T>
double SCluster::GetMinDistanceToClusterW(const T& h) {
    double minD = 1e4;
    for (auto& hit : fHitList) {
        double d = TPCLinesDistanceUtils::GetHitDistanceW(hit, h);
        if (d < minD) minD = d;
    }
    return minD;
}

template <typename T>
double SCluster::GetMinDistanceToClusterOverlap(const T& h) {
    double minD = 1e4;
    for (auto& hit : fHitList) {
        double d = TPCLinesDistanceUtils::GetHitDistanceOverlap(hit, h);
        if (d < minD) minD = d;
    }
    return minD;
}

template <typename T>
double SCluster::GetMinDistanceToCluster1D(const T& h) {
    double minD = 1e4;
    for (auto& hit : fHitList) {
        double d = TPCLinesDistanceUtils::GetHitDistance1D(hit, h);
        if (d < minD) minD = d;
    }
    return minD;
}

template <typename T>
std::pair<SHit, double> SCluster::GetClosestHitToPoint(const T& h) {
    double minD = 1e4;
    size_t closestHitIx = 0;

    for (size_t ix =0; ix<fHitList.size(); ix++){
        SHit hit = fHitList[ix];
        double d = TPCLinesDistanceUtils::GetHitDistance(hit, h);
        if (d < minD){
            minD = d;
            closestHitIx = ix;
        }
    }

    std::pair<SHit, double> result(fHitList[closestHitIx], minD);
    return result;
}


// Function to calculate the slope between two data points
double SLinearCluster::CalculateSlope(const SHit& p1, const SHit& p2) {
    return (p2.Y() - p1.Y()) / (p2.X() - p1.X());
}

    
// Function to find the data point with the maximum slope variation using N neighbors
SHit SLinearCluster::FindMaxVariationPointSlidingWindow(std::vector<SHit> data, size_t N) {
    if (data.size() < N + 2 + 1) {
        // Not enough data points for N neighbors; return an invalid point.
        return SHit(-1, -1, -1);
    }

    std::sort(data.begin(), data.end(), [&](SHit& a, SHit& b) {return a.X() < b.X();});
    TPCLinesPCA pcaAlgo;
    double maxVariation = 0.0;
    SHit maxVariationPoint = SHit(-1, -1, -1);

    std::vector<double> slopesVector;
    for (size_t i = N; i < data.size() - N; ++i) {

        std::vector<SHit> subset(data.begin() + (i - N), data.begin() + (i + N + 1));
        double subSlope = pcaAlgo.PerformPCA2D(subset).Slope();

        slopesVector.push_back(subSlope);
    }

    for (size_t j = 1; j < slopesVector.size(); j++) {

        double slopeVariation = std::abs(slopesVector[j]-slopesVector[j-1]);
        
        if (slopeVariation > maxVariation) {
            maxVariation = slopeVariation;
            maxVariationPoint = data[j+N];
        }
    }

    return maxVariationPoint;
}


//------------------------------  SLinearCluster
SLinearCluster::SLinearCluster(std::vector<SHit> hitList){

    fId = -1;
    fHitCluster = SCluster(hitList);

    fMinX = DefaultMax;
    fMinY = DefaultMax;
    fMaxX = 0;
    fMaxY = 0;
    fYAtMaxX = 0;
    fYAtMinX = 0;
    fMeanX = 0;
    fMeanY = 0;

    std::vector<double> yAtMinX;
    std::vector<double> yAtMaxX;

    for (auto &hit : hitList) {
        if (hit.X() < fMinX) {
            fMinX = hit.X();
        }
        if (hit.X() > fMaxX) {
            fMaxX = hit.X();
        }
        if (hit.Y() < fMinY) fMinY = hit.Y();
        if (hit.Y() > fMaxY) fMaxY = hit.Y();
        fMeanX += hit.X();
        fMeanY += hit.Y();
    }

    if (!hitList.empty()) {
        fMeanX /= hitList.size();
        fMeanY /= hitList.size();
        fCoMPoint = SPoint(fMeanX, fMeanY);
    }


    for (auto &hit : hitList) {
        if (hit.X() == fMinX) {
            yAtMinX.push_back(hit.Y());
        }
        if (hit.X() == fMaxX) {
            yAtMaxX.push_back(hit.Y());
        }
    }

    if(yAtMinX.size()>0 && yAtMaxX.size()>0){

        double meanYAtMinX = 0;
        double meanYAtMaxX = 0;
        meanYAtMinX = std::accumulate(yAtMinX.begin(), yAtMinX.end(), 0) / yAtMinX.size();
        meanYAtMaxX = std::accumulate(yAtMaxX.begin(), yAtMaxX.end(), 0) / yAtMaxX.size();

        if(meanYAtMaxX>meanYAtMinX){
            fYAtMaxX = *std::max_element(yAtMaxX.begin(), yAtMaxX.end());
            fYAtMinX = *std::min_element(yAtMinX.begin(), yAtMinX.end());
        }
        else{
            fYAtMaxX = *std::min_element(yAtMaxX.begin(), yAtMaxX.end());
            fYAtMinX = *std::max_element(yAtMinX.begin(), yAtMinX.end());
        }
    }
    

    // fill track equations
    if(hitList.size()>2){
        TPCLinesPCA pcaAlgo;
        fTrackEquation = pcaAlgo.PerformPCA2D(hitList);
        if(hitList.size()>6){
            double trackPercentageCut = 0.5;

            fTrackEquationStart = pcaAlgo.PerformPCA2DThreshold(hitList, trackPercentageCut );
            fTrackEquationEnd = pcaAlgo.PerformPCA2DThreshold(hitList, 1-trackPercentageCut, true);
        }
        else{
            fTrackEquationStart = fTrackEquation;
            fTrackEquationEnd = fTrackEquation;
        }
    }

    double yHypoAtMinX = fTrackEquationStart.EvaluateX(fMinX);
    double yHypoAtMaxX = fTrackEquationEnd.EvaluateX(fMaxX);

    float closestYAtMinX = fYAtMinX;
    float closestYAtMaxX = fYAtMaxX;
    float closestXAtMinX = fMinX;
    float closestXAtMaxX = fMaxX;

    // get closest hit to each edge
    double dMinAtMinX = DefaultMax;
    double dMinAtMaxX = DefaultMax;
    for(SHit &h:GetHits()){
        double d = std::hypot(fMinX-h.X(), 0.25*(yHypoAtMinX-h.Y()) );
        if(d<dMinAtMinX){
            dMinAtMinX = d;
            closestYAtMinX = h.Y();
            closestXAtMinX = h.X();
        }
        d = std::hypot(fMaxX-h.X(), 0.25*(yHypoAtMaxX-h.Y()) );
        if(d<dMinAtMaxX){
            dMinAtMaxX = d;
            closestYAtMaxX = h.Y();
            closestXAtMaxX = h.X();
        }
    }
    
    fStartPoint = SPoint(closestXAtMinX, closestYAtMinX);
    fEndPoint = SPoint(closestXAtMaxX, closestYAtMaxX);

    // members ro fill later
    fHasResidualHits = false;
    fHasStartEndPoints = false;
    
}


float SLinearCluster::GetHitDensity(){
    return NHits()/(fMaxX-fMinX);
}

double SLinearCluster::GetIntegral(){
    double w = 0;
    for (SHit hit : GetHits()) {
        w += hit.Integral();
    }
    return w;
}


std::vector<int> SLinearCluster::detect_outliers_iqr(std::vector<float> data, float threshold) {
    std::vector<int> outlier_indices;

    // Calculate the first and third quartiles (25th and 75th percentiles)
    size_t n = data.size();
    std::vector<float> sorted_data = data;
    std::sort(sorted_data.begin(), sorted_data.end());

    size_t q1_index = static_cast<size_t>(n * 0.25);
    size_t q3_index = static_cast<size_t>(n * 0.75);
    float q1 = sorted_data[q1_index];
    float q3 = sorted_data[q3_index];

    // Calculate the IQR
    float iqr = q3 - q1;

    // Calculate the lower and upper bounds for outliers
    float lower_bound = q1 - threshold * iqr;
    float upper_bound = q3 + threshold * iqr;

    // Find the outliers in the data
    for (size_t i = 0; i < n; ++i) {
        if (data[i] < lower_bound || data[i] > upper_bound) {
            outlier_indices.push_back(i);
        }
    }

    return outlier_indices;
}


// Function to detect outliers using IQR method
std::vector<int> SLinearCluster::detect_outliers_iqr2(std::vector<double> data, double threshold) {
    // Calculate the first and third quartiles (25th and 75th percentiles)
    std::sort(data.begin(), data.end());
    double q1 = data[data.size() / 4];
    double q3 = data[data.size() * 3 / 4];

    // Calculate the IQR
    double iqr = q3 - q1;

    // Calculate the lower and upper bounds for outliers
    double lower_bound = q1 - threshold * iqr;
    double upper_bound = q3 + threshold * iqr;

    // Find the outliers in the data
    std::vector<int> outlier_indices;
    for (size_t i = 0; i < data.size(); i++) {
        if (data[i] < lower_bound || data[i] > upper_bound) {
            outlier_indices.push_back(i);
        }
    }

    return outlier_indices;
}


void SLinearCluster::FillResidualHits(bool customKinkPoint) {
    fHasResidualHits = true;
    //std::cout << "\n\n +-+-+-+-+-+-+- Filling residual hits +-+-+-+-+-+-+-" << std::endl;

    if (NHits() > 0) {
        // Get average residual
        std::vector<float> dV;
        for (SHit hit : GetHits()) {
            float d = fTrackEquation.GetDistance(SPoint(hit.X(), hit.Y()));
            dV.push_back(d);
        }

        std::vector<int> outliersIx = detect_outliers_iqr(dV, 1.5);
        // std::vector<int> outliersIx = detect_outliers_zscore(dV, 3);

        std::vector<SHit> resHitList;
        std::vector<SHit> mainHitList;
        for (int ix = 0; ix < NHits(); ++ix) {
            if (std::find(outliersIx.begin(), outliersIx.end(), ix) != outliersIx.end()) {
                resHitList.push_back(fHitCluster.GetHits()[ix]);
            } else {
                mainHitList.push_back(fHitCluster.GetHits()[ix]);
            }
        }

        fResidualHitCluster = SCluster(resHitList);
        fMainHitCluster = SCluster(mainHitList);

        TPCLinesPCA pcaAlgo;
        double trackPercentageCut = 0.5;

        fTrackEquation= pcaAlgo.PerformPCA2D(mainHitList);
        if (mainHitList.size() > 6) {
            if(customKinkPoint){
                SHit kinkHit = FindMaxVariationPointSlidingWindow(mainHitList, 3);
                if(kinkHit.X()!=-1){
                    trackPercentageCut = (kinkHit.X()-fMinX) / (fMaxX-fMinX);
                }
            }
            
            fTrackEquationStart = pcaAlgo.PerformPCA2DThreshold(mainHitList, trackPercentageCut );
            fTrackEquationEnd = pcaAlgo.PerformPCA2DThreshold(mainHitList, 1-trackPercentageCut, true);
        } else {
            fTrackEquationStart = fTrackEquation;
            fTrackEquationEnd = fTrackEquation;
        }
    }
}


double SLinearCluster::GetOccupancy(){
    int numBins = NHits();
    std::vector<int> bins(numBins, 0);

    double trackLenght = std::hypot( fMaxX-fMinX, fYAtMaxX-fYAtMinX );
    double step = trackLenght/numBins;

    for (const SHit& h : GetHits()){
        SPoint cloPoint = fTrackEquation.GetLineClosestPoint(SPoint(h.X(), h.Y()));

        double d = std::hypot( cloPoint.X()-fMinX, cloPoint.Y()-fYAtMinX);

        unsigned int binIndex = static_cast<unsigned int> (d / step );
        
        if(binIndex<bins.size()) bins[binIndex]++;
    }

    int slotsFilled = 0;
    for (int count : bins) {
        if (count > 0) {
            slotsFilled++;
        }
    }

    return (1.*slotsFilled)/numBins;
}

double SLinearCluster::GetOccupancy1D(){
    int nX = fMaxX-fMinX;    

    if(nX<1) return 0;

    std::vector<int> bins(nX, 0);

    for (const SHit& h : GetHits()){
        unsigned int binIndex = static_cast<unsigned int> (h.X()-fMinX);
        if(binIndex<bins.size()) bins[binIndex]++;
    }

    int slotsFilled = 0;
    for (int count : bins) {
        slotsFilled+=count;
    }

    return (1.*slotsFilled)/(nX+1);
}

SPoint SLinearCluster::GetEdgeHit(SPoint p){
    
    std::vector<SPoint> projectedHits;
    double minX=1e4;
    double maxX = -1e4;
    SHit hitAtMinX;
    SHit hitAtMaxX;
   
    for(SHit &h:GetHits()){
        SPoint projP = fTrackEquation.GetLineClosestPoint(h);
        if(projP.X()>maxX) {
            maxX = projP.X();
            hitAtMaxX = h;
        }
        if(projP.X()<minX){
            minX = projP.X();
            hitAtMinX = h;
        }

        projectedHits.push_back(projP);
    }

    SHit edgeHit = ( std::abs(hitAtMinX.X()-p.X()) < std::abs(hitAtMaxX.X()-p.X()) )? hitAtMinX:hitAtMaxX;

    return SPoint(edgeHit.X(), edgeHit.Y());
}

SHit SLinearCluster::GetOppositeEdgeHit(SPoint p){
    
    double minX=1e4;
    double maxX = -1e4;
    SHit hitAtMinX;
    SHit hitAtMaxX;
   
    for(SHit &h:GetHits()){
        SPoint projP = fTrackEquation.GetLineClosestPoint(h);
        if(projP.X()>maxX) {
            maxX = projP.X();
            hitAtMaxX = h;
        }
        if(projP.X()<minX){
            minX = projP.X();
            hitAtMinX = h;
        }

    }

    if( std::abs(p.X()-minX) < std::abs(p.X()-maxX) )
        return hitAtMaxX;
    else
        return hitAtMinX;

}

double SLinearCluster::GetLength(){
    return std::hypot( 0.3*(fMaxX-fMinX), 0.075*(fYAtMaxX-fYAtMinX) );
}

// Function to get slope i¡of the track using Sliding window algorithm
void SLinearCluster::FillSlidingWindowLineEquations(size_t N){
    
    // reset
    fTrackEquationSlidingWindowV.clear();

    // round N to be next even number
    if(N%2!=0) N++;
    size_t halfStep = N/2;

    // sort hits in a map
    std::vector<SHit> hits = GetHits();
    std::map<int, std::vector<SHit>> hitsMap;

    for(SHit & h:hits){
        hitsMap[h.X()]={};
    }
    for(SHit & h:hits){
        hitsMap[h.X()].push_back(h);
    }

    std::vector<int> xValues;
    for(auto & pos:hitsMap){
        xValues.push_back(pos.first);
    }


    if (xValues.size() >= N + 2 + 1) {

       
        TPCLinesPCA pcaAlgo;

        std::vector<double> slopesVector;
        for (size_t i = halfStep; i < xValues.size() - halfStep; ++i) {
            int xCenter = xValues[i];
            std::vector<int> xSubset(xValues.begin() + (i - halfStep), xValues.begin() + (i + halfStep));
            std::vector<SHit> hitsSubset;
            for(int &x:xSubset){
                for(SHit &h:hitsMap[x]){
                    hitsSubset.push_back(h);
                }
            }

            LineEquation pcaDir = pcaAlgo.PerformPCA2D(hitsSubset);

            fTrackEquationSlidingWindowV[xCenter] = pcaDir;
        }
    
    }
    else{
        fTrackEquationSlidingWindowV[fMinX] = fTrackEquationStart;
        fTrackEquationSlidingWindowV[fMaxX] = fTrackEquationEnd;
    }

    
    return;
}

LineEquation SLinearCluster::GetLineEquationAtX(int x){
    // return the sliding window line equation at X or closest to x
    if(fTrackEquationSlidingWindowV.find(x)!=fTrackEquationSlidingWindowV.end()){
        return fTrackEquationSlidingWindowV[x];
    }
    else{//get closest X
        int minD = 1e4;
        int xIt = -1;
        for(auto &lineEqs:fTrackEquationSlidingWindowV){
            if( std::abs(lineEqs.first-x)<minD ){
                minD = std::abs(lineEqs.first-x);
                xIt = lineEqs.first;
            }
        }
        if(xIt!=-1){
            return fTrackEquationSlidingWindowV[xIt];
        }
        else{
            return fTrackEquation;
        }
    }
}

// instances
template double SCluster::GetMinDistanceToCluster(const SHit& h);
template double SCluster::GetMinDistanceToCluster(const SPoint& h);
template double SCluster::GetMinDistanceToClusterW(const SHit& h);
template double SCluster::GetMinDistanceToClusterOverlap(const SHit& h);
template double SCluster::GetMinDistanceToCluster1D(const SHit& h);
template double SCluster::GetMinDistanceToCluster1D(const SPoint& h);
template std::pair<SHit, double>  SCluster::GetClosestHitToPoint(const SHit& h);
template std::pair<SHit, double>  SCluster::GetClosestHitToPoint(const SPoint& h);
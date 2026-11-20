////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesAlgo.h
//
// \brief Definition of TPCLinesAlgo
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCLinesAlgo.h"

//----------------------------------------------------------------------
// Constructor
TPCLinesAlgo::TPCLinesAlgo(TPCLinesAlgoPsetType tpcLinesAlgoPset, FRANSPsetType fransPset):
    fTPCLinesPset(tpcLinesAlgoPset),
    fNTotalHits(0),
    fHitList({}),
    fHitListOutROI({}),
    fUnmatchedHits({}),
    fVertex(SPoint(-1, -1), ""),
    fHoughAlgo(tpcLinesAlgoPset.HoughAlgorithmPset),
    fTrackFinder(tpcLinesAlgoPset.TrackFinderAlgorithmPset),
    fVertexFinder(tpcLinesAlgoPset.VertexFinderAlgorithmPset),
    fFRANSAlgo(fransPset, tpcLinesAlgoPset.View),
    fTriangleCalo(tpcLinesAlgoPset.CaloAlgorithmPset),
    fHasTriangle(false),
    fDisplay(TPCLinesDisplay(tpcLinesAlgoPset.Verbose>0, tpcLinesAlgoPset.OutputPath)),
    fBadChannelsTPC0({4802, 4803, 4804, 4805, 4806, 4807}),
    fBadChannelsTPC1({10440, 10441, 10442, 10443, 10444, 10445})
{
}


// ---------------------------------------------------------------------
// Check overlap with APA
int TPCLinesAlgo::GetNOverlappedChannelsWithAPABadChannels(int ch1, int ch2){
    int nOverlappedChannels = 0;
    int minCh = std::min(ch1, ch2);
    int maxCh = std::max(ch1, ch2);
    for(int ch=minCh; ch<=maxCh; ch++){
        if(ch>=fBadChannelsTPC0.front() && ch<=fBadChannelsTPC0.back()) nOverlappedChannels++;
        if(ch>=fBadChannelsTPC1.front() && ch<=fBadChannelsTPC1.back()) nOverlappedChannels++;
    }

    return nOverlappedChannels;
}


//----------------------------------------------------------------------
// Display function
void TPCLinesAlgo::Display(std::string name, TCanvas *canvas, TCanvas *canvasCalo, TCanvas *canvasFRANS){

    std::vector<SHit> hitList = fHitList;
    hitList.insert(hitList.end(), fHitListOutROI.begin(), fHitListOutROI.end());

    fDisplay.Show(false, name, hitList, LineEquation(0, 0), fUnmatchedHits, fFinalTrackList, fMainDirection, fAngleList, fVertex, fVertexTrue, fOrigins, canvas);

    if(fHasTriangle){
        //fTriangleCalo.Display(canvasCalo);
        fFRANSAlgo.Display(canvasFRANS);
    }

    return;
}


//----------------------------------------------------------------------
// Get the view with the less Chi2
int TPCLinesAlgo::GetBestView(std::vector<int> *_View, std::vector<double> *_Chi2){

        std::map<int, double> meanChi2Map;

        std::map<int, double> sumChi2Map;
        std::map<int, int> countMap;

        for (size_t i = 0; i < _View->size(); ++i) {
            int viewValue = (*_View)[i];
            double chi2Value = (*_Chi2)[i];

            sumChi2Map[viewValue] += chi2Value;
            countMap[viewValue]++;
        }

        int minChi2View = -1;
        double minChi2=1e4;
        for (const auto& entry : sumChi2Map) {
            int viewValue = entry.first;
            double sumChi2 = entry.second;
            int count = countMap[viewValue];

            
            if (count != 0) {
                double meanChi2 = sumChi2/count; 
                if(meanChi2<minChi2){
                    minChi2=meanChi2;
                    minChi2View=viewValue;
                }
            }

        }

        return minChi2View;
}


//----------------------------------------------------------------------
// Function to get the average hit density
double TPCLinesAlgo::GetAverageHitDensity(){
    // initialize hit map
    std::map<int, int> fHitMap;
    for(int i=0; i<=fMaxX; i++){
        fHitMap[i]=0;
    }

    for(SHit & h:fHitList){
        fHitMap[h.X()]+=1;
    }

    double sum=0;
    int nwires=0;
    for(int i=0; i<=fMaxX; i++){
        if(fHitMap[i]>0 and 0.3*std::abs(i-fVertex.X())<40 ){
            sum+=fHitMap[i];
            nwires++;
        }
        
    }

    return sum/nwires;
}


//----------------------------------------------------------------------
// Remove isolated hits
std::vector<SHit> TPCLinesAlgo::RemoveIsolatedHits(std::vector<SHit> hitListForHough, std::vector<SHit>& discardedHits, double maxD, int minNeighbours) {

    if(fTPCLinesPset.Verbose>=2) std::cout << "\n\n +-+-+-+-+-+-+- Removing isolated hits +-+-+-+-+-+-+-" << std::endl;

    std::vector<SHit> hitListForHoughWithoutIsolatedHits;

    for (SHit& hit : hitListForHough) {
        int NNeighboursHits = 0;
        for (SHit& hit2 : hitListForHough) {
            if (hit.Id() == hit2.Id()) {
                continue;
            }
            //double d = GetHitDistance1D(hit, hit2);
            double d = std::abs(hit.X()-hit2.X());
            if (d <= maxD) {
                NNeighboursHits++;
            }

        }
        
        if (NNeighboursHits >= minNeighbours) {
            hitListForHoughWithoutIsolatedHits.push_back(hit);
        }
        else{
            discardedHits.push_back(hit);
        }
    }

    return hitListForHoughWithoutIsolatedHits;
}


//----------------------------------------------------------------------
// Merge isolated hits
std::vector<SLinearCluster> TPCLinesAlgo::MergeIsolatedHits(std::vector<SLinearCluster> recoTrackList, std::vector<SHit> hitList, double dCleaning1D, std::vector<SHit> & discardedHits, double dTh)
{
    if(fTPCLinesPset.Verbose>=2) std::cout << "\n\n +-+-+-+-+-+-+- Isolated hit merger +-+-+-+-+-+-+-" << std::endl;
    std::vector<SLinearCluster> newRecoTrackList;

    // initializations
    std::vector<std::vector<int>> hitToTrackDict(hitList.size());
    std::vector<bool> mergedHitsCounter(hitList.size(), false);
    for (size_t hitix = 0; hitix < hitList.size(); ++hitix) {
        hitToTrackDict[hitix] = {};
    }

    // first collect potential hits to merge
    for (size_t hitix = 0; hitix < hitList.size(); ++hitix) {
        SHit hit = hitList[hitix];
        double minD = 1e3;
        for (size_t trkix = 0; trkix < recoTrackList.size(); ++trkix) {
            SLinearCluster track = recoTrackList[trkix];
            double trackComp = track.GetCompactness();

            double hitTrackDist = track.GetHitCluster().GetMinDistanceToCluster1D(hit);
            if (hitTrackDist < dCleaning1D) {
                double d = track.GetTrackEquation().GetDistance(SPoint(hit.X(), hit.Y()));
                //double hypoY = track.GetTrackEquation().Slope() * hit.X() + track.GetTrackEquation().Intercept();
                if (d < dTh * trackComp && d < minD) {
                    minD = d;
                    hitToTrackDict[hitix].push_back(trkix);
                }
            }
        }
    }

    for (size_t trkix = 0; trkix < recoTrackList.size(); ++trkix) {
        SLinearCluster track = recoTrackList[trkix];
        double trackComp = track.GetCompactness();
        double trackConn = track.GetConnectedness();

        // Get average residual
        double trackAvResidual=0.;
        LineEquation trk_eq = track.GetTrackEquation();
        for (SHit hit : track.GetHits()) {
            trackAvResidual = trackAvResidual +  trk_eq.GetDistance(SPoint(hit.X(), hit.Y()));
        }
        trackAvResidual/=track.NHits();
        if(fTPCLinesPset.Verbose>=2) std::cout << "\n Merging track " << trkix << " " <<track.GetMinX()<<" "<<track.GetMaxX()<< " COMP " << trackComp << " CONN " << trackConn << " AvResidual: "<<trackAvResidual<<std::endl;

        std::vector<SHit> candidateHits;
        std::vector<int> candidateHitsIx;
        for (size_t hitix = 0; hitix < hitList.size(); ++hitix) {
            if (!mergedHitsCounter[hitix]) {
                if (std::find(hitToTrackDict[hitix].begin(), hitToTrackDict[hitix].end(), trkix) != hitToTrackDict[hitix].end()) {
                    candidateHits.push_back(hitList[hitix]);
                    candidateHitsIx.push_back(hitix);
                }
            }
        }

        std::vector<double> candidateHitsDist;
        for (SHit hit : candidateHits) {
            double distToCluster = track.GetHitCluster().GetMinDistanceToCluster(hit);
            candidateHitsDist.push_back(distToCluster);
        }

        // sort hits by distance to the cluster
        std::vector<int> sortedCandidateHits(candidateHits.size());
        std::iota(sortedCandidateHits.begin(), sortedCandidateHits.end(), 0);
        std::sort(sortedCandidateHits.begin(), sortedCandidateHits.end(), [&](int i, int j)->bool {
            return candidateHitsDist[i] < candidateHitsDist[j];
        });

        std::vector<int> sortedCandidateHitsIx(candidateHitsIx.size());
        std::iota(sortedCandidateHitsIx.begin(), sortedCandidateHitsIx.end(), 0);
        std::sort(sortedCandidateHitsIx.begin(), sortedCandidateHitsIx.end(), [&
        ](int i, int j) {
            return candidateHitsDist[i] < candidateHitsDist[j];
        });

        std::vector<SHit> newHitList = track.GetHits();
        SCluster currentCluster = track.GetHitCluster();

        if(fTPCLinesPset.Verbose>=2) std::cout << "  candidates hits: " << std::endl;
        for (size_t ix = 0; ix < candidateHits.size(); ++ix) {
            SHit hit = candidateHits[sortedCandidateHits[ix]];

            LineEquation trk_eqClosest = (std::abs(hit.X()-track.GetMinX())<std::abs(hit.X()-track.GetMaxX()))? track.GetTrackEquationStart():track.GetTrackEquationEnd();

            double minD = currentCluster.GetMinDistanceToCluster(hit);
            double minDConn = currentCluster.GetMinDistanceToClusterOverlap(hit);

            double hypoY = track.GetTrackEquation().Slope() * hit.X() + track.GetTrackEquation().Intercept();
            double minDHypo = std::abs(hypoY-hit.Y());
            
            double residual =  trk_eqClosest.GetDistance(SPoint(hit.X(), hit.Y()));

            if(fTPCLinesPset.Verbose>=2) std::cout << "       " << hit.Id() << "  X=" << hit.X() << " Y=" << hit.Y() << " " << " d " << minD << " dConn " << minDConn << " dCompHypo" << minDHypo << " residual: "<<residual << std::endl;
            if (minDConn < dTh * trackConn && minDHypo < dTh * trackComp && minD < dTh * trackComp){
                newHitList.push_back(hit);
                mergedHitsCounter[candidateHitsIx[sortedCandidateHits[ix]]] = true;
                currentCluster = SCluster(newHitList);
            }
        }

        newRecoTrackList.push_back(SLinearCluster(newHitList));
    
    }

    // add non merged hits to the discarded vector
    for (size_t hitix = 0; hitix < hitList.size(); ++hitix) {
        if (!mergedHitsCounter[hitix]) {
            discardedHits.push_back(hitList[hitix]);
        }
    }

    return newRecoTrackList;
}


//----------------------------------------------------------------------
// Wire distance subtracting the APA gap
int TPCLinesAlgo::GetWireDistance(int ch1, int ch2){
    ch1+=fMinX;
    ch2+=fMinX;
    
    // overlapped channels
    int nOverlappedChannels = GetNOverlappedChannelsWithAPABadChannels(ch1, ch2);
    
    return std::abs(ch2-ch1)-nOverlappedChannels;
}

//----------------------------------------------------------------------
// Merge out of ROI hits
std::vector<SLinearCluster> TPCLinesAlgo::MergeOutOfROIHits(std::vector<SLinearCluster> recoTrackList, std::vector<SHit> hitList, std::vector<SHit> & remainingHitList, int xStepTol)
{
    if(fTPCLinesPset.Verbose>=1) std::cout<<" --- Ina MergeOutOfROIHits function ---\n";
    remainingHitList.clear();

    // map of vector of hits
    std::map<int, std::vector<SHit>> newHitsMap;
    std::map<int, bool> usedHitsMap;
    for(size_t k=0; k<recoTrackList.size(); k++){
        newHitsMap[k] = {};
    }

    // loop over the free hits
    for(SHit & h:hitList){
        
        usedHitsMap[h.Id()]=false;

        // get the potential tracks to be added
        double minD = 1e4;
        int matchedTrackId = -1;
        for(size_t k=0; k<recoTrackList.size(); k++){

            if( h.ClusterId()==recoTrackList[k].GetHitCluster().GetMainClusterID() ){
                
                // calculate the min distance to the track
                double d = recoTrackList[k].GetHitCluster().GetMinDistanceToCluster(h);
                if(d<minD){
                    minD = d;
                    matchedTrackId = k;
                }
            }
        }
        if(matchedTrackId!=-1)
            newHitsMap[matchedTrackId].push_back(h);
    }


    // Create the new clusters with the new added hits
    std::vector<SLinearCluster> newRecoTrackList;
    if(newHitsMap.size()==0) return newRecoTrackList;

    for(size_t k=0; k<recoTrackList.size(); k++){
        std::vector<SHit> newHitList = recoTrackList[k].GetHits();    
        SLinearCluster updatedTrack = recoTrackList[k];

        std::map<int, std::vector<SHit>> hitsToAddMap;
        for(SHit & h:newHitsMap[k]){
            if(hitsToAddMap.find(h.X())==hitsToAddMap.end()){
                hitsToAddMap[h.X()]={h};
            }
            else{
                hitsToAddMap[h.X()].push_back(h);
            }
        }
        
        // loop over the hits
        int lastMinX = updatedTrack.GetMinX();
        int lastMaxX = updatedTrack.GetMaxX();
        for(auto & hitMap:hitsToAddMap){
            int x = hitMap.first;

            if(x>updatedTrack.GetMinX() && x<updatedTrack.GetMaxX()) continue;
            
            //LineEquation trackSlope = updatedTrack.GetTrackEquation();
            LineEquation trackSlope = updatedTrack.GetLineEquationAtX( x );

            bool tryToAdd = false;
            if(x<updatedTrack.GetMinX() && GetWireDistance(x, lastMinX)<=xStepTol && x!=lastMinX){ // 6 is maximum APA junction
                lastMinX = x;
                tryToAdd = true;        
            }
            else if(x>updatedTrack.GetMaxX() && GetWireDistance(x, lastMaxX)<=xStepTol && x!=lastMaxX){
                lastMaxX = x;
                tryToAdd = true;
            }

            if(tryToAdd){
                double minD = 1e4;
                std::vector<SHit> hits = hitMap.second;
                SHit selHit = hits[0];
                for(SHit &h:hits){
                    double d = trackSlope.GetDistance(h.GetPoint());
                    if(d<minD){
                        minD = d;
                        selHit = h;
                    }
                }
                
                //std::cout<<x<<" "<<minD<<" COMP "<<updatedTrack.GetCompactness()<<std::endl;
                if(minD<5*updatedTrack.GetCompactness()){
                    newHitList.push_back(selHit);
                    updatedTrack.AddHit(selHit);
                    updatedTrack.FillSlidingWindowLineEquations(10);
                    usedHitsMap[selHit.Id()]=true;
                }
            }
            
        }

        // redefine the cluster
        SLinearCluster newCluster(newHitList);
        newRecoTrackList.push_back(newCluster);
    }

    // fill unused hits
    for(SHit &h:hitList){
        if(usedHitsMap[h.Id()]==false){
            remainingHitList.push_back(h);
        }
    }
        
    return newRecoTrackList;
}


//----------------------------------------------------------------------
// Create out of ROI clusters
std::vector<SLinearCluster> TPCLinesAlgo::CreateOutOfROIClusters(std::vector<SLinearCluster> recoTrackList)
{

    if(fTPCLinesPset.Verbose>=1)  std::cout<<"--- In CreateOutOfROIClusters function ---\n";

    // get the used clusters IDs
    std::unordered_set<int> usedClusterID;
    for(SLinearCluster &cluster:recoTrackList){
        std::vector<int> clusterIDs = cluster.GetHitCluster().GetClusterIDs();
        for(int &id:clusterIDs){
            usedClusterID.insert(id);
        }
    }

    // loop over the out of ROI hits
    std::vector<SHit> hitsToAdd;
    std::unordered_set<int> freeClusterIDs;
    for(SHit &h:fHitListOutROI){
        // check the hit cluster ID is not in the unordered set
        if( usedClusterID.find(h.ClusterId()) == usedClusterID.end() ){
            hitsToAdd.push_back(h);
            freeClusterIDs.insert(h.ClusterId());
        }
    }


    std::map<int, std::vector<SHit>> newHitsMap;
    for (const auto& id:freeClusterIDs){
        newHitsMap[id] = {};
    }

    for(SHit &h:hitsToAdd){
        newHitsMap[h.ClusterId()].push_back(h);
    }

    std::vector<SLinearCluster> newRecoTrackList = recoTrackList;
    for(auto& clusterPair:newHitsMap){
        SLinearCluster newCluster(clusterPair.second);
        newRecoTrackList.push_back(newCluster);
    }


    return newRecoTrackList;
}


//----------------------------------------------------------------------
// Set the input hits
bool TPCLinesAlgo::SetHitList(int view,
                            std::vector<int>& vertex,
                            std::vector<int>& vertexTrue,
                            std::vector<SHit> hits)
{   
    // reset the class variables
    fHitList.clear();
    fHitListOutROI.clear();
    fNTotalHits=0;
    fClusterIdCounter.clear();

    // Define the reco vertex
    double vertexX = vertex[2];
    if (view == 0) vertexX = vertex[0];
    if (view == 1) vertexX = vertex[1];
    double vertexY = vertex[3] / fTPCLinesPset.DriftConversion;
    
    // Define the true vertex
    double vertexXTrue = vertexTrue[2];
    if (view == 0) vertexXTrue = vertexTrue[0];
    if (view == 1) vertexXTrue = vertexTrue[1];
    double vertexYTrue = vertexTrue[3] / fTPCLinesPset.DriftConversion;
    if(fTPCLinesPset.Verbose>=2) std::cout << "  ** Vertex XY: " << vertexX << " " << vertexY << std::endl;
    
    if (vertexX != -1){

        // Get hits in the ROI circle around the vertex
        // Hits are already filtered by plane type
        for (size_t i = 0; i < hits.size(); i++) { 
            
            double x = hits[i].X();
            double y = hits[i].Y();
            
            double d = std::sqrt(std::pow(x - vertexX, 2) + std::pow(y - vertexY, 2));
            
            if (d < fTPCLinesPset.MaxRadius) {
                fHitList.push_back(hits[i]);
            }
            else if(std::abs(x-vertexX)<1500){ // make sure is the same TPC
                fHitListOutROI.push_back(hits[i]);
            }
        }

        if ( (int)fHitList.size() >= fTPCLinesPset.MinTrackHits) {
            
            // Shift hits to have origin in (0,0), leave 3 ticks under/overflow
            double minX = std::min_element(fHitList.begin(), fHitList.end(), [](const SHit& h1, const SHit& h2) {return h1.X()<h2.X();})->X() - 3;
            double minY = std::min_element(fHitList.begin(), fHitList.end(), [](const SHit& h1, const SHit& h2) {return h1.Y()<h2.Y();})->Y() - 3;
            double maxX = std::max_element(fHitList.begin(), fHitList.end(), [](const SHit& h1, const SHit& h2) {return h1.X()>h2.X();})->X() - 3;        
            // Fill cluster counter map
            int maxClusterId = std::max_element(fHitList.begin(), fHitList.end(), [](const SHit& h1, const SHit& h2) {return h1.ClusterId()<h2.ClusterId();} )->ClusterId();        

            // Set cluster id counter
            for (int cID = 0; cID <= maxClusterId; cID++) {
                fClusterIdCounter[cID] = 0;
            }

            for (size_t i = 0; i < fHitList.size(); i++) {
                SHit h = fHitList[i];
                fHitList[i] = SHit(i, h.X() - minX, h.Y() - minY,
                                    h.Width(), h.Integral(),
                                    h.StartT() - minY, h.EndT() - minY,
                                    h.Chi2(), h.ClusterId(),
                                    h.SPX(), h.SPY(), h.SPZ()
                                    );
                fClusterIdCounter[h.ClusterId()]++;
            }

            for (size_t i = 0; i < fHitListOutROI.size(); i++) {
                SHit h = fHitListOutROI[i];
                fHitListOutROI[i] = SHit(fHitList.size()+i, h.X() - minX, h.Y() - minY,
                                        h.Width(), h.Integral(),
                                        h.StartT() - minY, h.EndT() - minY,
                                        h.Chi2(), h.ClusterId()
                                    );
            }

            fNTotalHits = fHitList.size();
    
            // create vertex with common origin
            vertexX = vertexX - minX;
            vertexY = vertexY - minY;
            vertexXTrue = vertexXTrue - minX;
            vertexYTrue = vertexYTrue - minY;
            fVertex = SVertex(SPoint(vertexX, vertexY), std::to_string(view));
            fVertexTrue = SVertex(SPoint(vertexXTrue, vertexYTrue), std::to_string(view));

            // set edge variables
            fMaxX = maxX - minX;
            fMinX = minX;
            fMinY = minY;
            if(fTPCLinesPset.Verbose>=1) std::cout << "  ** Origin vertex: " << fVertex << "  NInputHits:"<<fNTotalHits<<std::endl;

            return true;
        }
        else {
            if(fTPCLinesPset.Verbose>=1) std::cout << "   SKIPPED NHits selected near the vertex " << fHitList.size() << std::endl;
            return false;
        }

    }

    return false;
}


//----------------------------------------------------------------------
// Get the hits in a given cluster
std::vector<SHit> TPCLinesAlgo::GetHitsInCluster(int clusterId){

    std::vector<SHit> hits;
    for(SHit &h:fHitList){
        if(h.ClusterId()==clusterId) hits.push_back(h);
    }

    return hits;
}

double GetClusterHitDensity(std::vector<SHit> hits){

    std::map<int, int> occuppiedIDs;
    for(SHit &h:hits){
        if(occuppiedIDs.find(h.X())==occuppiedIDs.end()){
            occuppiedIDs[h.X()]=1;
        }
        else{
            occuppiedIDs[h.X()]++;
        }
    }

    double meanOccupation = 0;
    for(auto &occ:occuppiedIDs){
        meanOccupation+=occ.second;
    }
    meanOccupation/=occuppiedIDs.size();

    return meanOccupation;
}


//----------------------------------------------------------------------
// Main function
void TPCLinesAlgo::AnaView(std::string eventLabel)
{
    // --- Reconstructed objects
    std::vector<SHit> discardedHits;
    std::vector<SLinearCluster> finalLinearClusterV;


    if(fTPCLinesPset.Verbose>=2){
        fDisplay.Show(true, "Input hits", fHitList, LineEquation(-1, -1));
    }

    for(std::pair<int, int> clusterPair: fClusterIdCounter){
            
        //std::vector<SHit> hitListForHough = fHitList;
        std::vector<SHit> hitListForHough = GetHitsInCluster(clusterPair.first);
        double clusterMeanOccupation = GetClusterHitDensity(hitListForHough);
        if(fTPCLinesPset.Verbose>=1){
            std::cout<<"   Analyzing cluster "<<clusterPair.first<<" with "<<clusterPair.second<<" hits \n";
            std::cout<<"    MeanOccupation: "<<clusterMeanOccupation<<std::endl;
        }

        if(clusterMeanOccupation>3){
            SLinearCluster track(hitListForHough);
            finalLinearClusterV.push_back(track);
            continue;
        }
    
        int trkIterCluster = (int)hitListForHough.size()>=fTPCLinesPset.MinTrackHits ? 0 : fTPCLinesPset.MaxHoughTracks;
 
        // --- Loop over the hough tracks
        while(trkIterCluster<fTPCLinesPset.MaxHoughTracks){

            if(fTPCLinesPset.Verbose>=2) std::cout<<" **** Track finder iteration: "<<trkIterCluster<< " # hough hits:"<<hitListForHough.size()<<std::endl;

            // -- Get the best Hough line       
            HoughLine houghLine = fHoughAlgo.GetBestHoughLine(hitListForHough, fVertex);
            if(fTPCLinesPset.Verbose>=3)
                if(fTPCLinesPset.Verbose>=2) std::cout<<"    Hough line results: "<<houghLine.NHits()<<" Score: "<<houghLine.Score()<<std::endl;

            // -- Skip if not enough hits
            if(houghLine.NHits() < fTPCLinesPset.MinTrackHits){
                if(fTPCLinesPset.Verbose>=3)
                    if(fTPCLinesPset.Verbose>=2) std::cout<<"   Hough lines has <"<<fTPCLinesPset.MinTrackHits<<", ending the search\n";
                trkIterCluster = fTPCLinesPset.MaxHoughTracks;
            }

            // -- Call the track finder otherwise
            else{

                // -- Reconstruct the tracks from the Hough hits
                std::vector<SLinearCluster> linearClusterV = fTrackFinder.ReconstructTracksFromHoughDirection(hitListForHough, houghLine.GetLineEquation(), trkIterCluster);

                // -- check the found tracks has enough hits
                if(linearClusterV.size()!=0){
                    
                    finalLinearClusterV.insert(finalLinearClusterV.begin(), linearClusterV.begin(), linearClusterV.end());
                    
                    if(fTPCLinesPset.Verbose>=2){
                        std::cout<<"   Found "<<linearClusterV.size()<<" new tracks! Current track list: \n";
                        for(SLinearCluster &trk:finalLinearClusterV){
                            std::cout<<"   "<<trk.GetId()<<" "<<trk.GetMinX()<<" "<<trk.GetMaxX()<<std::endl;
                        }
                    }

                    // -- check there's enough hits for the next Hough iteration
                    if((int)hitListForHough.size() < fTPCLinesPset.MinTrackHits){
                        if(fTPCLinesPset.Verbose>=3)
                            if(fTPCLinesPset.Verbose>=2) std::cout<<"   Remaining hits are "<<hitListForHough.size()<<", ending the search\n";
                        trkIterCluster = fTPCLinesPset.MaxHoughTracks;
                    }
                }
                // end Hough loop otherwise
                else{
                    if(fTPCLinesPset.Verbose>=3)
                        if(fTPCLinesPset.Verbose>=2) std::cout<<"   Not found new tracks\n";
                    trkIterCluster = fTPCLinesPset.MaxHoughTracks;
                }
                
                // -- Remove isolated hits
                if (fTPCLinesPset.RemoveIsolatedHits) {
                    hitListForHough = RemoveIsolatedHits(hitListForHough, discardedHits, fTPCLinesPset.MaxNeighbourDistance, fTPCLinesPset.MinNeighboursHits);
                }

            }

            trkIterCluster++;
        } // -- end Hough loop

        discardedHits.insert(discardedHits.end(), hitListForHough.begin(), hitListForHough.end());

    } // -- end cluster loop

    
    
    // -- Slope track merger
    // sort by minX
    std::sort(finalLinearClusterV.begin(), finalLinearClusterV.end(), [&](SLinearCluster& l1, SLinearCluster& l2) {return l1.GetMinX() < l2.GetMinX();} );    
    finalLinearClusterV = TPCLinesDirectionUtils::SlopeTrackMerger(finalLinearClusterV, 2, 15, fTPCLinesPset.Verbose); 

    // -- Characterize the tracks
    for(size_t ix = 0; ix<finalLinearClusterV.size(); ix++){
        finalLinearClusterV[ix].FillResidualHits(fTPCLinesPset.CustomKinkPoint);
        finalLinearClusterV[ix].AssignId(ix);
    }
        

    // -- Get the short/vertex tracks and remove them
    std::map<int, std::vector<int>> shortToLongTrackDict;
    std::map<int, int> shortConnectionTrackDict;
    std::vector<SLinearCluster> shortTrackList = TPCLinesDirectionUtils::GetVertexTracks(finalLinearClusterV, shortToLongTrackDict, shortConnectionTrackDict, 6, 3, 5, fTPCLinesPset.Verbose); 
    // remove the short tracks
    std::vector<SLinearCluster> newTrackList;
    for (SLinearCluster& track : finalLinearClusterV) {
    
        if (shortToLongTrackDict.find(track.GetId()) == shortToLongTrackDict.end()) {

            SLinearCluster newTrack = track;
            
            for(auto & sTrackId:shortToLongTrackDict){
                if( sTrackId.second[0]==track.GetId() || sTrackId.second[1]==track.GetId() ) {
                    std::vector<SHit> auxHits = track.GetHits();
                    std::vector<SHit> auxHitsShort;
                    for(SLinearCluster sTrack:shortTrackList){
                        
                        if(sTrack.GetId()==sTrackId.first){
                            auxHitsShort = sTrack.GetHits(); 
                        }
                    }
                    
                    auxHits.insert(auxHits.begin(), auxHitsShort.begin(), auxHitsShort.end());
                    newTrack = SLinearCluster(auxHits);
                }
                
            }
            
            newTrack.FillResidualHits(fTPCLinesPset.CustomKinkPoint);
            newTrack.AssignId(track.GetId());
            newTrackList.push_back(newTrack);
        }
    }

    finalLinearClusterV.clear();
    for(SLinearCluster & trk: newTrackList){
        finalLinearClusterV.push_back(trk);
    }

    std::vector<SLinearCluster> NewTrackList;
    NewTrackList.clear();
    SLinearCluster mainDirection;
    std::vector<SOrigin> intersectionsInBall;
    intersectionsInBall.clear();
    std::vector<STriangle> vertexList;
    vertexList.clear();
    std::vector<SOrigin> associatedOrigins;
    associatedOrigins.clear();
    
    if(finalLinearClusterV.size()>0){
        
        // -- Get the parallel tracks
        std::vector<std::vector<SLinearCluster>> parallelTracks = TPCLinesDirectionUtils::GetParallelTracks(finalLinearClusterV, -2, 15, fTPCLinesPset.ParallelTracksMaxD, fTPCLinesPset.Verbose);
        NewTrackList.clear();

        for(size_t ix = 0; ix<parallelTracks.size(); ix++){ 
            std::vector<SHit> hitList = parallelTracks[ix][0].GetHits();
            for(size_t jx = 1; jx<parallelTracks[ix].size(); jx++){
                std::vector<SHit> moreHits = parallelTracks[ix][jx].GetHits();
                hitList.insert(hitList.end(), moreHits.begin(), moreHits.end());
            }
            SLinearCluster newTrack(hitList);
            newTrack.FillResidualHits(fTPCLinesPset.CustomKinkPoint);
            newTrack.AssignId(ix);
            NewTrackList.push_back( newTrack );    
        }

        // -- Isolated hit merger
        std::vector<SHit> remainingHits;// = hitListForHough;
        remainingHits.insert(remainingHits.end(), discardedHits.begin(), discardedHits.end());
        discardedHits.clear();
        NewTrackList = MergeIsolatedHits(NewTrackList, remainingHits, 10, discardedHits);
        // Characterize the tracks
        for(size_t ix = 0; ix<NewTrackList.size(); ix++){
            NewTrackList[ix].FillResidualHits(fTPCLinesPset.CustomKinkPoint);
            NewTrackList[ix].AssignId(ix);
        }

        // -- Out of ROI hits merger
        // Fill sliding window intersections
        for(SLinearCluster &trk:NewTrackList) trk.FillSlidingWindowLineEquations(10);
        
        std::vector<SHit> missingHits = discardedHits;
        missingHits.insert(missingHits.end(), fHitListOutROI.begin(), fHitListOutROI.end());
        std::vector<SHit> remainingHitsAfterROIMerge;
        NewTrackList = MergeOutOfROIHits(NewTrackList, missingHits, remainingHitsAfterROIMerge, 3);    
        // Characterize the tracks
        for(size_t ix = 0; ix<NewTrackList.size(); ix++){
            NewTrackList[ix].FillResidualHits(fTPCLinesPset.CustomKinkPoint);
            NewTrackList[ix].AssignId(ix);
        }
        discardedHits.clear();
        discardedHits = remainingHitsAfterROIMerge;

        // --- Create Out of ROI clusters ---
        NewTrackList = CreateOutOfROIClusters(NewTrackList);
        // Characterize the tracks
        for(size_t ix = 0; ix<NewTrackList.size(); ix++){
            NewTrackList[ix].FillResidualHits(fTPCLinesPset.CustomKinkPoint);
            NewTrackList[ix].AssignId(ix);
        }

        // Display
        if(fTPCLinesPset.Verbose>=2)
            fDisplay.Show(true, "Reconstructed line clusters", fHitList, LineEquation(-1, -1), {}, NewTrackList);
        

        // -- Find secondary vertices
        if(NewTrackList.size()>0){
            intersectionsInBall = fVertexFinder.GetAngleVertices(NewTrackList, fVertex.Point(), vertexList, associatedOrigins, mainDirection);
        }

        // Set the main vertex
        fMainVertex = mainDirection.GetStartPoint();
        float d1 = std::hypot( 0.3*(mainDirection.GetStartPoint().X() - fVertex.Point().X()), 0.075*(mainDirection.GetStartPoint().Y() - fVertex.Point().Y()) );
        float d2 = std::hypot( 0.3*(mainDirection.GetEndPoint().X() - fVertex.Point().X()), 0.075*(mainDirection.GetEndPoint().Y() - fVertex.Point().Y()) );
        fMainVertex = (d1<d2)? mainDirection.GetStartPoint() : mainDirection.GetEndPoint();

        // Save final results
        fFinalTrackList = NewTrackList;
        fMainDirection = mainDirection;
        fAngleList = vertexList;
        fOrigins = intersectionsInBall;

        for(SOrigin & ori: intersectionsInBall){
            std::cout<<ori;
        }

    }

    double hitDensity = GetAverageHitDensity();

    // Save the reco event
    SEvent recoEvent(NewTrackList, intersectionsInBall, vertexList, associatedOrigins, hitDensity, discardedHits);
    fRecoEvent = recoEvent;


    return;
}



STriangle TPCLinesAlgo::FillLambdaAnaTree(LambdaAnaTree &lambdaAnaTree){

    // --- Hits for FRANS ---
    std::vector<SHit> allHits = fHitList;
    if(fHitListOutROI.size()>0)
        allHits.insert(allHits.end(), fHitListOutROI.begin(), fHitListOutROI.end());
    
    // --- FRANS with PANDORA vertex ---
    std::cout << "** Computing FRANS with Pandora vertex" << std::endl;
    fFRANSAlgo.Fill(allHits, fVertex);
    double FRANSScorePANDORA = fFRANSAlgo.Score();
    lambdaAnaTree.fFRANSScorePANDORA = FRANSScorePANDORA;
    
    // --- FRANS with custom origin: get the best triangle ---
    std::cout << "** Computing FRANS with reconstructed origins" << std::endl;
    std::vector<STriangle> angleList = fRecoEvent.GetAngleList();
    std::vector<SOrigin> associatedOrigins = fRecoEvent.GetAssociatedOrigins();
    
    int bestTriangleIx = -1;
    double bestFRANSScore = -1000;
    double bestFRANSEta = -1000;
    double bestFRANSIota = -1000;
    double bestFRANSDelta = -1000;
    double bestFRANSFitScore = -1000;

    fHasTriangle = false;
    for(size_t orix=0; orix<angleList.size(); orix++){
        SVertex fVertexMine(associatedOrigins[orix].GetPoint(), "");

        fFRANSAlgo.Fill(allHits, fVertexMine);
        double score = fFRANSAlgo.Score();

        if(score>bestFRANSScore){
            bestFRANSScore = score;
            bestTriangleIx = orix;
            bestFRANSEta = fFRANSAlgo.Eta();
            bestFRANSIota = fFRANSAlgo.Iota();
            bestFRANSDelta = fFRANSAlgo.Delta();
            bestFRANSFitScore = fFRANSAlgo.FitScore();
        }
    }

    std::cout << "Best FRANS score: " << bestFRANSScore << " with bestTriangleIx " << bestTriangleIx << std::endl;

    // --- Number of origins variables ---
    lambdaAnaTree.fNOrigins = fRecoEvent.GetNOrigins();
    lambdaAnaTree.fNOriginsMult1 = fRecoEvent.GetNOriginsMult(1);
    lambdaAnaTree.fNOriginsMult2 = fRecoEvent.GetNOriginsMult(2);
    lambdaAnaTree.fNOriginsMultGT3 = fRecoEvent.GetNOriginsMultGt(3);
    lambdaAnaTree.fNOriginsPairOneTwo = fRecoEvent.GetNOriginsMult(2) * fRecoEvent.GetNOriginsMult(1);
    lambdaAnaTree.fNAngles = fRecoEvent.GetNAngles();

    // --- Best triangle variables ---
    if(bestTriangleIx!=-1){
        fHasTriangle = true;
        STriangle bestTriangle = angleList[bestTriangleIx];
        SOrigin bestOrigin = associatedOrigins[bestTriangleIx];
        std::cout<<" --- Best triangle with tracks "<<bestTriangle.GetTrack1().GetId()<<" and "<<bestTriangle.GetTrack2().GetId()<<" with associated origin track "<<bestOrigin.GetTrack(0).GetId()<<std::endl;


        // --- Triangle variables ---
        lambdaAnaTree.fAngleFRANSScore = bestFRANSScore;
        lambdaAnaTree.fAngleFRANSEta = bestFRANSEta;
        lambdaAnaTree.fAngleFRANSIota = bestFRANSIota;
        lambdaAnaTree.fAngleFRANSDelta = bestFRANSDelta;
        lambdaAnaTree.fAngleFRANSFitScore = bestFRANSFitScore;

        lambdaAnaTree.fAngleGap = bestTriangle.GetGap();
        lambdaAnaTree.fAngleDecayContainedDiff = bestTriangle.GetDecayAngleDifference();
        lambdaAnaTree.fAngleNHits = bestTriangle.GetNHitsTriangle();
        lambdaAnaTree.fAngleNHitsTrack1 = bestTriangle.GetNHitsTrack1();
        lambdaAnaTree.fAngleNHitsTrack2 = bestTriangle.GetNHitsTrack2();
        lambdaAnaTree.fAngleNHitsMainTrack = bestTriangle.GetNHitsMainTrack();
        lambdaAnaTree.fAngleLengthTrack1 = bestTriangle.GetLengthTrack1();
        lambdaAnaTree.fAngleLengthTrack2 = bestTriangle.GetLengthTrack2();
        lambdaAnaTree.fAngleLengthMainTrack = bestTriangle.GetLengthMainTrack();
        lambdaAnaTree.fAngleMinNHits = std::min(bestTriangle.GetNHitsTrack1(), bestTriangle.GetNHitsTrack2());
        lambdaAnaTree.fAngleOpeningAngle = bestTriangle.GetOpeningAngle();
        
        // --- Unassociated origins ---
        SEvent notAssociatedRecoEvent = fRecoEvent.UnassociatedOrigins(bestTriangle);
        lambdaAnaTree.fNUnOrigins = notAssociatedRecoEvent.GetNOrigins();
        lambdaAnaTree.fNUnOriginsMult1 = notAssociatedRecoEvent.GetNOriginsMult(1);
        lambdaAnaTree.fNUnOriginsMult2 = notAssociatedRecoEvent.GetNOriginsMult(2);
        lambdaAnaTree.fNUnOriginsMultGT3 = notAssociatedRecoEvent.GetNOriginsMultGt(3);
        if( bestTriangle.GetLengthMainTrack()> bestTriangle.GetLengthTrack1() && bestTriangle.GetLengthMainTrack()> bestTriangle.GetLengthTrack2() ){
            lambdaAnaTree.fAngleLongestIsMain = true;
        }
        else{
            lambdaAnaTree.fAngleLongestIsMain = false;
        }
        lambdaAnaTree.fAngleCoveredArea = bestTriangle.ComputeCoveredArea();
        // Print track connections
        fRecoEvent.PrintTrackConnections();

        // Angle-main track overlap
        lambdaAnaTree.fAngleMainTrackOverlap = bestTriangle.GetOverlapWithMainTrack();
        // Kinematics
        int PzLambda = bestTriangle.GetMainVertex().X()-bestOrigin.GetPoint().X() > 0 ? 1 : -1;
        SLinearCluster muonTrack = bestOrigin.GetTrack(0);
        int PzMuon = muonTrack.GetMeanX() - bestOrigin.GetPoint().X() > 0 ? 1 : -1;
        lambdaAnaTree.fAnglePzSignLambda = PzLambda;
        lambdaAnaTree.fAnglePzSignMuon = PzMuon;
        lambdaAnaTree.fAnglePzSign = PzLambda + PzMuon;
        // APA gap
        int lambdaVertexCh = bestTriangle.GetMainVertex().X() + ShiftX();
        int muonVertexCh = bestOrigin.GetPoint().X() + ShiftX();
        int gapMinCh = std::min(lambdaVertexCh, muonVertexCh);
        int gapMaxCh = std::max(lambdaVertexCh, muonVertexCh);
        int nOverlappedChannelsWithAPA = GetNOverlappedChannelsWithAPABadChannels(gapMinCh, gapMaxCh);
        double gapSize = gapMaxCh - gapMinCh;
        if(gapSize>0){
            lambdaAnaTree.fAngleGapOverlapWithAPAJuntion = 1.*nOverlappedChannelsWithAPA/gapSize;
        }
        else{
            lambdaAnaTree.fAngleGapOverlapWithAPAJuntion = 0;
        }
        
        // --- Triangle cleaness ---
        int nDirtHitsInTriangle = 0;
        double nFractionDirtHitsInTriangle = 0;
        int nDirtHitsInTriangleWires = 0;
        double nFractionDirtHitsInTriangleWires = 0;
        fRecoEvent.FreeHitsAroundTriangle(bestTriangle,
                                        nDirtHitsInTriangle,
                                        nFractionDirtHitsInTriangle,
                                        nDirtHitsInTriangleWires,
                                        nFractionDirtHitsInTriangleWires);
        
        lambdaAnaTree.fAngleDirtHits = nDirtHitsInTriangle;
        lambdaAnaTree.fAngleDirtHitsWires = nDirtHitsInTriangleWires;
        lambdaAnaTree.fAngleDirtHitsRatio = nFractionDirtHitsInTriangle;
        lambdaAnaTree.fAngleDirtHitsWiresRatio = nFractionDirtHitsInTriangleWires;

        // --- Unassociated hits ---
        int nFreeHits = 0;
        int nUnassociatedHits = 0;
        fRecoEvent.GetUnassociatedHits(bestTriangle, nFreeHits, nUnassociatedHits);               
        lambdaAnaTree.fNUnassociatedHits = nUnassociatedHits;
        lambdaAnaTree.fNFreeHits = nFreeHits;

        lambdaAnaTree.fNMaxDirtUnassociatedHits = std::max(nUnassociatedHits, nDirtHitsInTriangle);


        // // --- Triangle calorimetry ---
        // fTriangleCalo.SetTriangle(bestTriangle);
        // fTriangleCalo.JointFitAnalysis();

        // lambdaAnaTree.fAnglePassFit = fTriangleCalo.PassFit();
        // lambdaAnaTree.fAngleTwoLinesChi2 = fTriangleCalo.TwoLinesChi2();
        // lambdaAnaTree.fAngleNVertexHits = fTriangleCalo.NVertexHits();
        // lambdaAnaTree.fAngleNBulkHits = fTriangleCalo.NBulkHits();
        // lambdaAnaTree.fAngleVertexHitIntegralRatio = fTriangleCalo.VertexHitIntegralRatio();
        // lambdaAnaTree.fAngleVertexHitIntegralDifference = fTriangleCalo.VertexHitIntegralDifference();
        // lambdaAnaTree.fAngleVertexHitIntegralRelativeDifference = fTriangleCalo.VertexHitIntegralRelativeDifference();
        // lambdaAnaTree.fAngleTrackLength1 = fTriangleCalo.TrackLength1();
        // lambdaAnaTree.fAngleTrackLength2 = fTriangleCalo.TrackLength2();
        // lambdaAnaTree.fAngleTrackLengthRatio = fTriangleCalo.TrackLengthRatio();
        // lambdaAnaTree.fAngleResidualRange1RMS = fTriangleCalo.ResidualRange1RMS();
        // lambdaAnaTree.fAngleResidualRange2RMS = fTriangleCalo.ResidualRange2RMS();
        // lambdaAnaTree.fAngleResidualRangeMinRMS = fTriangleCalo.ResidualRangeMinRMS();
        // lambdaAnaTree.fAngleResidualRangeMaxRMS = fTriangleCalo.ResidualRangeMaxRMS();
        // lambdaAnaTree.fAngleResidualRange1AngleRMS = fTriangleCalo.ResidualRange1AngleRMS();
        // lambdaAnaTree.fAngleResidualRange2AngleRMS = fTriangleCalo.ResidualRange2AngleRMS();
        // lambdaAnaTree.fAngleResidualRangeMinAngleRMS = fTriangleCalo.ResidualRangeMinAngleRMS();
        // lambdaAnaTree.fAngleResidualRangeMaxAngleRMS = fTriangleCalo.ResidualRangeMaxAngleRMS();
        // lambdaAnaTree.fAngleChargeRatioAverage = fTriangleCalo.ChargeRatioAverage();
        // lambdaAnaTree.fAngleChargeDifferenceAverage = fTriangleCalo.ChargeDifferenceAverage();
        // lambdaAnaTree.fAngleChargeRelativeDifferenceAverage = fTriangleCalo.ChargeRelativeDifferenceAverage();
        // lambdaAnaTree.fAnglePassChargeFit = fTriangleCalo.PassChargeFit();
        // lambdaAnaTree.fAngleBandOverlap = fTriangleCalo.BandOverlap();
        // lambdaAnaTree.fAngleBandCrossHits = fTriangleCalo.BandCrossHits();
        // lambdaAnaTree.fAngleChargeRatioFit = fTriangleCalo.ChargeRatioFit();
        // lambdaAnaTree.fAngleChargeDifferenceFit = fTriangleCalo.ChargeDifferenceFit();
        // lambdaAnaTree.fAngleChargeRelativeDifferenceFit = fTriangleCalo.ChargeRelativeDifferenceFit();
        // lambdaAnaTree.fAngleChargeRatioIntegral = fTriangleCalo.ChargeRatioIntegral();
        // lambdaAnaTree.fAngleChargeDifferenceIntegral = fTriangleCalo.ChargeDifferenceIntegral();
        // lambdaAnaTree.fAngleChargeRelativeDifferenceIntegral = fTriangleCalo.ChargeRelativeDifferenceIntegral();

        return bestTriangle;
    }

    return STriangle();
}
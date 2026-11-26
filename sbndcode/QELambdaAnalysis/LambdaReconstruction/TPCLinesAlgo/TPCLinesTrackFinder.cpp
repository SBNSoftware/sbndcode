////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesTrackFinder.h
//
// \brief Definition of TPCLinesTrackFinder
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCLinesTrackFinder.h"

//----------------------------------------------------------------------
// Constructor
TPCLinesTrackFinder::TPCLinesTrackFinder(TrackFinderAlgorithmPsetType tpcLinesTrackFinderPset):
    fTPCLinesTrackFinderPset(tpcLinesTrackFinderPset),
    fDisplay(TPCLinesDisplay(tpcLinesTrackFinderPset.Verbose>0))
{}


//----------------------------------------------------------------------
// Get hii density
double TPCLinesTrackFinder::GetHitLinearDensity(std::vector<SHit> data){
    int minX = 1e6;
    int maxX = -1e6;
    for (SHit& h : data) {
        if(h.X()<minX) minX = h.X();
        if(h.X()>maxX) maxX = h.X();
    }
    
    return 1.*data.size()/(maxX+1-minX);
}
    

// Get Hits in tube
void TPCLinesTrackFinder::GetHitsInTube(std::vector<SHit> & hitsIn, std::vector<SHit> & hitsOut, LineEquation& line, std::vector<SHit>& hitList, double maxD) {
    hitsIn.clear();
    hitsOut.clear();

    for (auto& hit : hitList) {
        SPoint pProj = line.GetLineClosestPoint( SPoint(hit.X(), hit.Y()) );
        double d = std::sqrt(std::pow(hit.X() - pProj.X(), 2) + std::pow(hit.Y() - pProj.Y(), 2));
        if (d < maxD) {
            hitsIn.push_back(hit);
        } else {
            hitsOut.push_back(hit);
        }
    }

    return;
}


//----------------------------------------------------------------------
// Get Hits in tube, single wire mode
void TPCLinesTrackFinder::GetHitsInTubeSingleWire(std::vector<SHit> & hitsIn, std::vector<SHit> & hitsOut, LineEquation& line, std::vector<SHit>& hitList, double maxD){
    hitsIn.clear();
    hitsOut.clear();

    std::unordered_map<int, int> hitIndexDict;
    std::unordered_map<int, double> hitDistDict;
    std::unordered_map<int, std::vector<SHit>> discardedhitsDict;

    for (size_t hitix = 0; hitix < hitList.size(); hitix++) {
        SHit hit = hitList[hitix];
        SPoint pProj = line.GetLineClosestPoint( SPoint(hit.X(), hit.Y()) );
        double d = std::sqrt(std::pow(hit.X() - pProj.X(), 2) + std::pow(hit.Y() - pProj.Y(), 2));
        if (d < maxD) {
            int xIndex = static_cast<int>(hit.X());
            //if(fTPCLinesTrackFinderPset.Verbose>=2)  std::cout<<hit<<std::endl;
            if (hitIndexDict.count(xIndex) > 0) {
                if (d < hitDistDict[xIndex]) {
                    int refHitIndex = hitIndexDict[static_cast<int>(hit.X())];
                    SHit oldHit = hitsIn[refHitIndex];
                    hitsIn[refHitIndex] = hit;
                    hitDistDict[xIndex] = d;
                    discardedhitsDict[xIndex].push_back(oldHit);
                    //if(fTPCLinesTrackFinderPset.Verbose>=2)  std::cout<<" Found index, Replacing\n";
                } else {
                    hitsOut.push_back(hit);
                    //if(fTPCLinesTrackFinderPset.Verbose>=2)  std::cout<<"Found index, Not replacing\n";
                }
            } else {
                //if(fTPCLinesTrackFinderPset.Verbose>=2)  std::cout<<" Not found index, adding\n";
                hitsIn.push_back(hit);
                hitIndexDict[xIndex] = hitsIn.size() - 1;
                hitDistDict[xIndex] = d;
            }
        }
        else {
            hitsOut.push_back(hit);
        }
    }

    for (auto& entry : discardedhitsDict) {
        int wireix = entry.first;
        auto& discardedHitList = entry.second;
        int refHitIndex = hitIndexDict[wireix];
        std::vector<SHit> refHitV = { hitsIn[refHitIndex] };
        std::sort(discardedHitList.begin(), discardedHitList.end(), [&](SHit& a, SHit& b) {
            return std::abs(a.Y() - refHitV[0].Y()) < std::abs(b.Y() - refHitV[0].Y());
        });

        for (auto& hit : discardedHitList) {
            bool doOverlap = false;
            for (auto& refHit : refHitV) {
                if (TPCLinesDistanceUtils::HitWidthOverlap(hit, refHit)) {
                    doOverlap = true;
                }
            }
            if (doOverlap) {
                refHitV.push_back(hit);
                hitsIn.push_back(hit);
            } else {
                hitsOut.push_back(hit);
            }
        }
    }

    return;
}


//----------------------------------------------------------------------
// Vector to map of frequencies
std::vector<std::pair<int, int>> TPCLinesTrackFinder::GetListFrequency(std::vector<int> vList) {

    // Count the frequency of each element in the list
    std::unordered_map<int, int> frequency_dict;
    for (int element : vList) {
        // skip -1 (no associated cluster)
        if (element == -1) {
            continue;
        }
        if (frequency_dict.find(element) != frequency_dict.end()) {
            frequency_dict[element]++;
        } else {
            frequency_dict[element] = 1;
        }
    }

    // Sort the dictionary by frequency in descending order
    std::vector<std::pair<int, int>> ordered_values(frequency_dict.begin(), frequency_dict.end());
    std::sort(ordered_values.begin(), ordered_values.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
        return a.second > b.second;});

    return ordered_values;

}


//----------------------------------------------------------------------
// DBSCAN 1D clusters
std::pair<double, double> TPCLinesTrackFinder::ComputeConnectivityMode(std::vector<double> v, int minClusterHits, double weightWidth) {
    // Sort connV
    std::vector<double> connV(v);
    std::sort(connV.begin(), connV.end());


    // Get data ready for DBSCAN
    //NOTE: Using SHit format as it's the one I have implemented, should improve this
    std::vector<SHit> data;
    for (const auto& x : connV) {
        data.push_back(SHit(0, x, 0, 0, 0, 0, 0));
    }

    // Initialize DBSCAN
    DBSCAN dbscan(weightWidth, minClusterHits);
    dbscan.setDistanceFunction(DBSCANHitEuclidianDistance);
    
    // Fit the points
    dbscan.fit(data);
    std::vector<int>& clusterAssignment = dbscan.getClusterAssignment();

    // Get the cluster frequencies maps
    std::vector<std::pair<int, int>> clusterFrequencies = GetListFrequency(clusterAssignment);

    // Prints
    if(fTPCLinesTrackFinderPset.Verbose>=2){

        // Hit assigments
        for (size_t i = 0; i < data.size(); ++i) {
            std::cout << data[i].X() << " " << data[i].X() << " Cluster=" << clusterAssignment[i] << std::endl;
        }

        // Cluster frequencies
        for(auto & x:clusterFrequencies){
            if(fTPCLinesTrackFinderPset.Verbose>=2)  std::cout<<" Cluster "<<x.first<<" with frequency "<<x.second<<std::endl;
        }
    }


    // Get the connectivity for the cluster with the largest frequency
    if(clusterFrequencies.size()>0){
        std::vector<double> mainClusterConnV;
        int mainClusterIx = clusterFrequencies[0].first;
        for(size_t ix=0; ix<connV.size(); ix++){
            if(clusterAssignment[ix]==mainClusterIx){
                mainClusterConnV.push_back(connV[ix]);
            }
        }

        double epsilon = std::accumulate(mainClusterConnV.begin(), mainClusterConnV.end(), 0.0) / mainClusterConnV.size();
        double epsilonRMS = 0.0;
        for (size_t i = 0; i < mainClusterConnV.size(); ++i) {
            epsilonRMS += (mainClusterConnV[i] - epsilon)*(mainClusterConnV[i] - epsilon);
        }
        epsilonRMS = std::sqrt(epsilonRMS / mainClusterConnV.size());

        return std::make_pair(epsilon, epsilonRMS);
    }
    else {
        return std::make_pair(-1.0, -1.0);
    }
    
}


//----------------------------------------------------------------------
// DBSCAN 2D clusters
std::vector<SLinearCluster> TPCLinesTrackFinder::Get2DClusters(std::vector<SHit> & hits, double epsilon, int minPts, std::string option){

    // Initialize DBSCAN
    DBSCAN dbscan(epsilon, minPts);
    //dbscan.setDistanceFunction(manhattanDistance);  // You can change this to euclideanDistance or any other custom function
    
    std::vector<SHit> & data = hits;
    // sort hits by distance to the center
    double meanX = 0;
    double meanY = .0;
    for(SHit & h:data){
        meanX+=h.X();
        meanY+=h.Y();

    }
    meanX/=data.size();
    meanY/=data.size();
    SPoint CoM(meanX, meanY);
    std::sort( data.begin(), data.end(), [&](SHit& h1, SHit& h2) {return h1.X() < h2.X();} );
    //std::sort( data.begin(), data.end(), [&](SHit& h1, SHit& h2) {return GetPointDistance(CoM, SPoint(h1.X(), h1.Y()))<GetPointDistance(CoM, SPoint(h2.X(), h2.Y()));} );
    //std::sort( data.begin(), data.end(), [&](SHit& h1, SHit& h2) {return std::abs(h1.X()-CoM.X()) < std::abs(h2.X()-CoM.X());} );

    // Set the distance function
    if(option=="DistanceWidth")
        dbscan.setDistanceFunction(DBSCANHitWidthDistance);
    else if(option=="DistanceWidth2")
        dbscan.setDistanceFunction(DBSCANHitWidthDistance2);
    else if(option=="DistanceOverlap")
        dbscan.setDistanceFunction(DBSCANHitOverlapDistance);
    else if(option=="DistanceEuclidianDriftConversion")
        dbscan.setDistanceFunction(DBSCANHitEuclidianDistanceDriftConversion);
    else if(option=="DistanceEuclidian")
        dbscan.setDistanceFunction(DBSCANHitEuclidianDistance);
    
    // Fit the points
    dbscan.fit(data);
    
    // Print results
    std::vector<int>& clusterAssignment = dbscan.getClusterAssignment();
    for (size_t i = 0; i < data.size(); ++i) {
        if (clusterAssignment[i] == -1) {
            if(fTPCLinesTrackFinderPset.Verbose>=2)  std::cout << "SPoint (" << data[i].X() << ", " << data[i].Y() << ") is noise" << std::endl;
        } else {
            if(fTPCLinesTrackFinderPset.Verbose>=2)  std::cout << "SPoint (" << data[i].X() << ", " << data[i].Y() << ") belongs to cluster " << clusterAssignment[i] << std::endl;
        }
    }

    std::vector<std::pair<int, int>> clusterFrequencies = GetListFrequency(clusterAssignment);
    for(auto & x:clusterFrequencies){
        if(fTPCLinesTrackFinderPset.Verbose>=2)  std::cout<<" Cluster "<<x.first<<" with frequency "<<x.second<<std::endl;
    }
    

    // create the line cluster
    std::vector<SLinearCluster> linearClusterV;
    for(size_t cIx=0; cIx<clusterFrequencies.size(); cIx++){
        int cID = (int)clusterFrequencies[cIx].first;
        int nhitscluster = clusterFrequencies[cIx].second;
        if(nhitscluster>=minPts){

            std::vector<SHit> clusterHits;
            for(size_t k=0; k<data.size(); k++){
                if(clusterAssignment[k]==cID)
                    clusterHits.push_back(data[k]);
            }

            linearClusterV.push_back(SLinearCluster(clusterHits));

        }
    }
        
    return linearClusterV;
}


//----------------------------------------------------------------------
// Capture missing hits
// Define the C++ equivalent function
std::vector<SLinearCluster> TPCLinesTrackFinder::CaptureMissingHits(std::vector<SLinearCluster> trackV, std::vector<SHit> hitList) {
    
    // vector to store the udpated tracks
    std::vector<SLinearCluster> newRecoTrackV;

    // hit list to be added
    std::vector<SHit> freeHitsList = hitList;

    for(size_t trkix=0; trkix<trackV.size(); trkix++){

        SLinearCluster newRecoTrack = trackV[trkix];
        std::vector<SHit> candidateHitList;
        std::vector<SHit> unmatchedHitList;

        LineEquation trackLine = newRecoTrack.GetTrackEquation();

        for (SHit& hit : freeHitsList) {
            // Get track line hypothesis for the hit
            double pEvalY = trackLine.Evaluate(hit.X());

            // If hypothesis is within the hit RMS, keep
            if (std::abs(hit.Y() - pEvalY)< hit.Width() ){
                candidateHitList.push_back(hit);
                if(fTPCLinesTrackFinderPset.Verbose>=2)  std::cout<<"MERGING "<<hit.X()<<std::endl;
            }
            else {
                unmatchedHitList.push_back(hit);
            }
        }

        if (!candidateHitList.empty()) {
            
            // Sort by distance to track
            SCluster hitCluster = newRecoTrack.GetHitCluster();
            
            std::sort(candidateHitList.begin(), candidateHitList.end(), [&](SHit& h1, SHit& h2) {
                return hitCluster.GetMinDistanceToCluster(h1) < hitCluster.GetMinDistanceToCluster(h2);
            });

            for (SHit& hit : candidateHitList) {
                if(fTPCLinesTrackFinderPset.Verbose>=2)  std::cout<<"   Check conn "<<hit.X()<<std::endl;
                double epsilon = hitCluster.GetConnectedness() + 3 * hitCluster.GetConnectednessRMS();
                double conn = hitCluster.GetMinDistanceToClusterW(hit);

                // if connected, update the track
                if (conn < epsilon) {
                    std::vector<SHit> newHitList = newRecoTrack.GetHits();
                    newHitList.push_back(hit);
                    newRecoTrack = SLinearCluster(newHitList);
                    hitCluster = newRecoTrack.GetHitCluster();
                }
                // else, add to the unmatched hits
                else {
                    unmatchedHitList.push_back(hit);
                }
            }

        }

        // add the updated track
        newRecoTrackV.push_back(newRecoTrack);

        // Return the updated newRecoTrack and unmatchedHitList and update the free hit list
        freeHitsList.clear();
        for (SHit& hit : unmatchedHitList){
            freeHitsList.push_back(hit);
        }
    }

    return newRecoTrackV;
}

//----------------------------------------------------------------------
// Make tracks function
std::vector<SLinearCluster> TPCLinesTrackFinder::MakeTrack(std::vector<SLinearCluster> linearClusterList) {
    
    // Order cluster list by the number of hits
    std::sort(linearClusterList.begin(), linearClusterList.end(), [](SLinearCluster& lCluster1, SLinearCluster& lCluster2) {
        return lCluster1.NHits() > lCluster2.NHits();
    });

    std::vector<SLinearCluster> finalTrackClusterList;

    for (SLinearCluster & lCluster : linearClusterList) {
        LineEquation trackEquation = lCluster.GetTrackEquation();

        std::vector<SHit> trackHitList;

        for (SHit& hit : lCluster.GetHits()) {
            double d = trackEquation.GetDistance(SPoint(hit.X(), hit.Y()));

            if (d < 2 * hit.Width()) {
                trackHitList.push_back(hit);
            }
        }

        SLinearCluster cluster(trackHitList);
        finalTrackClusterList.emplace_back(cluster);
    }

    finalTrackClusterList = TPCLinesDirectionUtils::SlopeTrackMerger(finalTrackClusterList, 2, 5, fTPCLinesTrackFinderPset.Verbose);


    std::sort(finalTrackClusterList.begin(), finalTrackClusterList.end(), [](SLinearCluster& lCluster1, SLinearCluster& lCluster2) {
        return lCluster1.GetHits().size() > lCluster2.GetHits().size();
    });

    return finalTrackClusterList;


}



//----------------------------------------------------------------------
// Main function
// Input: vector of SHits
// Returns vector of SLinearClusters
// Algorithm: Find PCA line, get hits in tube around PCA line, make connectedness DBSCAN clusters,
// make compactness DBSCAN clusters, and creat SClusters
std::vector<SLinearCluster> TPCLinesTrackFinder::ReconstructTracksFromHoughDirection(std::vector<SHit> & hitList, LineEquation houghLine, int trkIter) {


    //---------- Hough tube block
    // Get hits in the Hough tube
    std::vector<SHit> hitHoughTubeList, hitOutHoughTubeList;
    GetHitsInTube(hitHoughTubeList, hitOutHoughTubeList, houghLine, hitList,  fTPCLinesTrackFinderPset.MaxDTube);
    
    // Compute hit density fo the Hough hits
    double hitDensity = GetHitLinearDensity(hitHoughTubeList);
    if(fTPCLinesTrackFinderPset.Verbose>=2)  std::cout<<"    Hit density:"<<hitDensity<<std::endl;
    
    

    //::DISPLAY
    if(fTPCLinesTrackFinderPset.Verbose>=2){
        fDisplay.Show(true, "Pre-hough direction", hitList, LineEquation(0, -1));
        fDisplay.Show(true, "Hough direction", hitList, houghLine, hitHoughTubeList);
    }
    
    //---------- PCA tube block
    // Get hits in the refined PCA direction
    TPCLinesPCA pcaAlgo;
    LineEquation pcaLine = pcaAlgo.PerformPCA2D(hitHoughTubeList);
    std::vector<SHit> hitPCATubeList, hitPCAOutTubeList;
    if (!fTPCLinesTrackFinderPset.SingleWireMode && hitDensity >= fTPCLinesTrackFinderPset.HitDensityThreshold) {
        GetHitsInTube(hitPCATubeList, hitPCAOutTubeList, pcaLine, hitList,  fTPCLinesTrackFinderPset.MaxDCluster);
    }
    else {
        if(fTPCLinesTrackFinderPset.Verbose>=2) std::cout<<"Single wire mode\n";
        GetHitsInTubeSingleWire(hitPCATubeList, hitPCAOutTubeList, pcaLine, hitList,  fTPCLinesTrackFinderPset.MaxDCluster);
    }

    // Create the PCA cluster, return if not enough hits
    SCluster pcaCluster(hitPCATubeList);
    if(pcaCluster.NHits()< fTPCLinesTrackFinderPset.MinTrackHits){
        return {};
    }

    //::DISPLAY
    if(fTPCLinesTrackFinderPset.Verbose>=2){
        fDisplay.Show(true, "PCA direction", hitList, pcaLine, hitPCATubeList);
    }


    //---------- Connectedness clusters block
    if(fTPCLinesTrackFinderPset.Verbose>=2) {
        std::cout<<"\n******** Making connectedness clusters\n";
        std::cout<<pcaCluster<<std::endl;
    }

    // Get the epsilon parameter for the 2D clustering
    std::vector<double> clusterConnV = pcaCluster.GetConnectednessV();
    
    double weightWidth = fTPCLinesTrackFinderPset.ConnectednessWidthTol * pcaCluster.GetAverageWidth();
    std::pair<double, double> connectivityResult = ComputeConnectivityMode(clusterConnV, fTPCLinesTrackFinderPset.MinClusterHits, weightWidth);
    double epsilon=connectivityResult.first + fTPCLinesTrackFinderPset.ConnectednessTol * connectivityResult.second;
    if(fTPCLinesTrackFinderPset.Verbose>=2) std::cout << "EpsilonMean: " << connectivityResult.first << ", EpsilonRMS: " << connectivityResult.second << " Epsilon: "<< epsilon << std::endl;
    
    // Perform PCA 2D clustering
    std::vector<SLinearCluster> connectedLinearClustersV;
    if(epsilon>0){
        std::vector<SHit> pcaClusterHits = pcaCluster.GetHits();
        connectedLinearClustersV = Get2DClusters(pcaClusterHits, epsilon, fTPCLinesTrackFinderPset.MinClusterHits, "DistanceWidth");
    }

    //::DISPLAY
    if(fTPCLinesTrackFinderPset.Verbose>=2){
        fDisplay.Show(true, "DBSCAN clusters: connectedness", hitList, pcaLine, hitPCATubeList, connectedLinearClustersV);
    }

    // Capture missing hits by the original PCA direction
    if(fTPCLinesTrackFinderPset.CaptureMissingHits==true and hitPCAOutTubeList.size()>1){
        if(fTPCLinesTrackFinderPset.Verbose>=2)  std::cout<<"Trying to merge additional hits out of the main tube\n";
        std::vector<SHit> unmatchedHits = hitPCAOutTubeList;
        connectedLinearClustersV = CaptureMissingHits(connectedLinearClustersV, unmatchedHits);
    }
    
    //::DISPLAY
    if(fTPCLinesTrackFinderPset.Verbose>=2){
        fDisplay.Show(true, "DBSCAN clusters: connectedness", hitList, houghLine, hitPCATubeList, connectedLinearClustersV);
    }
    
    //---------- Compactness clusters block
    // Vector to store the compact clusters
    std::vector<SLinearCluster> finalClustersV = connectedLinearClustersV;
    
    if(fTPCLinesTrackFinderPset.UseCompactness == true){
        if(fTPCLinesTrackFinderPset.Verbose>=2)  std::cout<<"\n******** Making compactness clusters\n";

        // Vector to store the compact clusters
        std::vector<SLinearCluster> compactLinearClustersV;

        for(SLinearCluster &lCluster:connectedLinearClustersV){

            double clusterCompactness = lCluster.GetCompactness();
            double clusterConnectedness = lCluster.GetConnectedness();
            if(fTPCLinesTrackFinderPset.Verbose>=2)  std::cout<<"\n  Cluster compactness "<<clusterCompactness<<"  Connectedness: "<<clusterConnectedness<<std::endl;

            std::vector<double> compactnessV = lCluster.GetCompactnessV();
            std::pair<double, double> compactnessResult = ComputeConnectivityMode(compactnessV, fTPCLinesTrackFinderPset.MinClusterHits, clusterConnectedness);

            double epsilon;
            //epsilon=compactnessResult.first+5*compactnessResult.second;
            epsilon=fTPCLinesTrackFinderPset.CompactnessTol*compactnessResult.first;

            if(fTPCLinesTrackFinderPset.Verbose>=2)  std::cout<<"   --- Compacness Eps="<<compactnessResult.first<< " pm "<<compactnessResult.second<<" Epsilon="<<epsilon<<std::endl;
            if(epsilon>0){
                std::vector<SHit> hits = lCluster.GetHits();
                std::vector<SLinearCluster> compactClusterList = Get2DClusters(hits, epsilon, fTPCLinesTrackFinderPset.MinClusterHits, "DistanceEuclidianDriftConversion");
                //std::vector<SLinearCluster> compactClusterList = Get2DClusters(hits, epsilon, fTPCLinesTrackFinderPset.MinClusterHits, "DistanceEuclidian");
                if(compactClusterList.size()>0)
                    compactLinearClustersV.insert(compactLinearClustersV.end(),compactClusterList.begin(),compactClusterList.end());
                else
                    compactLinearClustersV.push_back(lCluster);
            }
            else
                compactLinearClustersV.push_back(lCluster);
        }
        
        // Check kinks in the tracks
        if(compactLinearClustersV.size()>0){

            std::vector<SLinearCluster> finalCompactLinearClustersV;
            finalCompactLinearClustersV.clear();

            for(SLinearCluster &lCluster:compactLinearClustersV){

                double angle1 = std::atan(lCluster.GetTrackEquationStart().Slope()) * 180.0 / M_PI;
                double angle2 = std::atan(lCluster.GetTrackEquationEnd().Slope()) * 180.0 / M_PI;
                double slope_angle_difference = std::abs(angle1 - angle2);
                double hitDensity = lCluster.GetHitDensity();
                if(slope_angle_difference>fTPCLinesTrackFinderPset.ClusterAngleCut && hitDensity<fTPCLinesTrackFinderPset.HitDensityThreshold){
                    if(fTPCLinesTrackFinderPset.Verbose>=2)  std::cout<<"Slope difference: "<<slope_angle_difference<<std::endl;
                    SHit kinkHit = lCluster.FindMaxVariationPointSlidingWindow(lCluster.GetHits(), 3);
                    if(kinkHit.X()!=-1){
                        // split the cluster
                        std::vector<SHit> startHits, endHits;
                        for(SHit &h:lCluster.GetHits()){
                            if(h.X()<kinkHit.X()) startHits.push_back(h);
                            else if(h.X()>kinkHit.X()) endHits.push_back(h);
                        }

                        SLinearCluster startCluster(startHits);
                        SLinearCluster endCluster(endHits);

                        // assign the kink hit to the closest line
                        double dStart = startCluster.GetTrackEquation().GetDistance(SPoint(kinkHit.X(), kinkHit.Y()));
                        double dEnd = endCluster.GetTrackEquation().GetDistance(SPoint(kinkHit.X(), kinkHit.Y()));

                        if(dStart<dEnd){
                            startHits.push_back(kinkHit);
                            startCluster = SLinearCluster(startHits);
                        }
                        else{
                            endHits.push_back(kinkHit);
                            endCluster = SLinearCluster(endHits);
                        }

                        double startCompactness = startCluster.GetCompactness();
                        double endCompactness = endCluster.GetCompactness();
                        double maxCompactness = std::max(startCompactness, endCompactness);
                        double minCompactness = std::min(startCompactness, endCompactness);

                        if(maxCompactness>minCompactness*fTPCLinesTrackFinderPset.CompactnessTol){
                            finalCompactLinearClustersV.push_back(startCluster);
                            finalCompactLinearClustersV.push_back(endCluster);
                        }
                        else{
                            finalCompactLinearClustersV.push_back(lCluster);
                        }
                    }
                    else{
                        finalCompactLinearClustersV.push_back(lCluster);
                    }
                    
                }
                else{
                    finalCompactLinearClustersV.push_back(lCluster);
                }
            }

            finalClustersV.clear();
            finalClustersV = finalCompactLinearClustersV;
        }

        //::DISPLAY
        if(fTPCLinesTrackFinderPset.Verbose>=2){
            fDisplay.Show(true, "DBSCAN clusters: compactness", hitList, pcaLine, hitPCATubeList, finalClustersV);
        }

    }

    // Vector of SLinearClusters to return
    std::vector<SLinearCluster> finalRecoTracks;
    if(finalClustersV.size()>0)
        finalRecoTracks = MakeTrack(finalClustersV);
    
    std::vector<SLinearCluster> recoTracks;
    recoTracks = finalRecoTracks;

    // Free hits for the next iteration
    std::vector<int> usedHitsIds;
    for(SLinearCluster &lCluster:recoTracks){
        std::vector<SHit> clusterHits = lCluster.GetHits();
        for(SHit &h:clusterHits){
            usedHitsIds.push_back(h.Id());
        }
    }

    std::vector<SHit> newFreeHits;
    for(SHit &h:hitList){
        if( std::find(usedHitsIds.begin(), usedHitsIds.end(), h.Id())==usedHitsIds.end()){
            newFreeHits.push_back(h);
        }
    }

    if(fTPCLinesTrackFinderPset.Verbose>=2)  std::cout<<" Remaining hits: "<<newFreeHits.size()<<std::endl;
    hitList.clear();
    for (SHit& h : newFreeHits) {
        hitList.push_back(h);
    }
    
    // return
    return {recoTracks};
}

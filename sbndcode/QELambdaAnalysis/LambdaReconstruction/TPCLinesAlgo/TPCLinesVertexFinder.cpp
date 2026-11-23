////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesAlgo.h
//
// \brief Definition of TPCLinesHough
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCLinesVertexFinder.h"

// ------- constructor
TPCLinesVertexFinder::TPCLinesVertexFinder(VertexFinderAlgorithmPsetType tpcLinesVertexFinderPset):
    fTPCLinesVertexFinderPset(tpcLinesVertexFinderPset)
{}


// ------- Get the angle of a vector in the 360 degrees space
double TPCLinesVertexFinder::GetAngle360(double x, double y) {
    if(x == 0){
        if(y>0) return 90;
        else return 270;
    }
    else{
        double a = std::atan(std::abs(y / x)) * 180.0 / M_PI;

        if (x > 0 && y < 0) { // quadrant 4
            a = 360 - a;
        } else if (x < 0 && y < 0) { // quadrant 3
            a = 180 + a;
        } else if (x < 0 && y > 0) { // quadrant 2
            a = 180 - a;
        }
        return a;
    }   
}


// ------- Check if the junction between a track and a triangle is contained within the triangle
bool TPCLinesVertexFinder::TrackTriangleJunctionContained(SOrigin ori, STriangle tri, double extraAngle){
    
    SPoint mainTrackVertex = ori.GetPoint();

    double juntionDirection[] = {
        tri.GetMainVertex().X() - mainTrackVertex.X(),
        tri.GetMainVertex().Y() - mainTrackVertex.Y()
    };

    double VDir1[] = {
        tri.GetVertexB().X() - tri.GetMainVertex().X(),
        tri.GetVertexB().Y() - tri.GetMainVertex().Y()
    };

    double VDir2[] = {
        tri.GetVertexC().X() - tri.GetMainVertex().X(),
        tri.GetVertexC().Y() - tri.GetMainVertex().Y()
    };

    double juntionDirectionAngle = GetAngle360(juntionDirection[0], juntionDirection[1]);
    double VDir1Angle = GetAngle360(VDir1[0], VDir1[1]);
    double VDir2Angle = GetAngle360(VDir2[0], VDir2[1]);
    
    if(fTPCLinesVertexFinderPset.Verbose>=1){
        std::cout << "    MainTrackVertex: " << mainTrackVertex;
        std::cout << "    V ANGLES --- Junction:  " << juntionDirectionAngle << ",1: " << VDir1Angle << ", 2: " << VDir2Angle << std::endl;
    }

    double minAngle = std::min(VDir1Angle, VDir2Angle);
    double maxAngle = std::max(VDir1Angle, VDir2Angle);

    bool junctionContained = false;
    
    if (maxAngle - minAngle < 180) {
        junctionContained = (minAngle-extraAngle<=juntionDirectionAngle) && (juntionDirectionAngle<=maxAngle+extraAngle);
    }
    else {
        junctionContained = (juntionDirectionAngle<=minAngle+extraAngle) || (juntionDirectionAngle>=maxAngle-extraAngle);
    }

    return junctionContained;
}


// ------- Get the number of hits between a track and a triangle
int TPCLinesVertexFinder::GetNHitsBetweenJunction(SLinearCluster mainTrack, STriangle triangle, std::vector<SLinearCluster> freeTrackList, double tol){


    SPoint p1(mainTrack.GetMinX(), mainTrack.GetYatMinX());
    double d1 = TPCLinesDistanceUtils::GetHitDistance(p1, triangle.GetMainVertex());
    SPoint p2(mainTrack.GetMaxX(), mainTrack.GetYatMaxX());
    double d2 = TPCLinesDistanceUtils::GetHitDistance(p2, triangle.GetMainVertex());
    SPoint mainTrackVertex = (d1 < d2) ? p1 : p2;


    // revisit this
    double juntionDirection[] = {
        triangle.GetMainVertex().X() - mainTrackVertex.X(),
        triangle.GetMainVertex().Y() - mainTrackVertex.Y()
    };

    double juntionEdges[2] = { std::min(mainTrackVertex.X(), triangle.GetMainVertex().X()), std::max(mainTrackVertex.X(), triangle.GetMainVertex().X()) };
    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "  Junction edges: " << juntionEdges[0] << ", " << juntionEdges[1] << std::endl;

    double juncSlope = juntionDirection[1] / juntionDirection[0];
    double juncIntercept = triangle.GetMainVertex().Y() - juncSlope * triangle.GetMainVertex().X();

    int nHitsInMiddle = 0;
    for (SLinearCluster & eTrack : freeTrackList) {

        if ((eTrack.GetMinX() > juntionEdges[1] && eTrack.GetMaxX() > juntionEdges[1]) ||
            (eTrack.GetMinX() < juntionEdges[0] && eTrack.GetMaxX() < juntionEdges[0])) {
            continue;
        }

        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "  Track" << eTrack.GetId() << " found in the junction\n";
        
        for (SHit& hit : eTrack.GetHits()) {
            double yHypo = juncSlope * hit.X() + juncIntercept;
            if (std::abs(yHypo - hit.Y()) < tol*hit.Width()) {
                nHitsInMiddle++;
                if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "   In middle hit: " << hit.X() << ", " << hit.Y() << std::endl;
            }
        }
    }

    return nHitsInMiddle;

}


// ------- Check if two vectors are equal
bool TPCLinesVertexFinder::areVectorsEqual(const std::vector<int>& vec1, const std::vector<int>& vec2) {
    if (vec1.size() != vec2.size()) {
        return false;
    }

    for (size_t i = 0; i < vec1.size(); ++i) {
        if (vec1[i] != vec2[i]) {
            return false;
        }
    }

    return true;
}


// ------- Get the intersection point between two tracks
SPoint TPCLinesVertexFinder::GetTracksIntersection(SLinearCluster track1, SLinearCluster track2, double dMax, bool useEdgeSlopes, bool useFit){
    SPoint intP(-1, -1);

    if (useFit) {
        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" Using fit not implemented\n";
        // Implement the useFit logic here.
        // The code related to evaluating track fits and finding intersects goes here.
    }
    else {
        LineEquation lineEq1 = track1.GetTrackEquation();
        LineEquation lineEq2 = track2.GetTrackEquation();
        double xInt = (lineEq1.Intercept() - lineEq2.Intercept()) / (lineEq2.Slope() - lineEq1.Slope());
        double yInt = lineEq1.Slope() * xInt + lineEq1.Intercept();

        intP = SPoint(xInt, yInt);

        if (useEdgeSlopes == true) {
            lineEq1 = track1.GetLineEquationAtX(intP.X());
            lineEq2 = track2.GetLineEquationAtX(intP.X());

            xInt = (lineEq1.Intercept() - lineEq2.Intercept()) / (lineEq2.Slope() - lineEq1.Slope());
            yInt = lineEq1.Slope() * xInt + lineEq1.Intercept();
            intP = SPoint(xInt, yInt);
        }
    }

    SHit intHit = SHit(intP.X(), intP.Y());
    double d1 = track1.GetHitCluster().GetMinDistanceToCluster(intHit);
    double d2 = track2.GetHitCluster().GetMinDistanceToCluster(intHit);

    if (d1 < dMax || d2 < dMax) {
        return intP;
    } else {
        return SPoint(-1, -1); // Return an appropriate "no intersection" value.
    }
}


// ------- Get hits contained in a given line equation  
int TPCLinesVertexFinder::GetHitsContainedInLineEquation(LineEquation trackEq, std::vector<SHit> hitList, float tol) {
    
    int nhits = 0;

    for (SHit& hit : hitList){
        float yHypo = trackEq.Slope() * hit.X() + trackEq.Intercept();
        if (std::abs(hit.Y() - yHypo) < tol * hit.Width()) {
            nhits++;
        }
    }

    return nhits;
}


// ----- Compare hits by min X
bool CompareHitsByMinX(SHit h1, SHit h2) {
    return h1.X() < h2.X();
}

// GetHitsContainedInHypo
// Input two tracks, get the NHits closer to a vertex
// Check how many hits of each track are contained in the line equation of the other track
std::vector<SHit> TPCLinesVertexFinder::GetMutualHitsContainedInHypo(SLinearCluster track1, SLinearCluster track2,  SPoint intP, int nHits, float tol) {

    // Get hits closest to the itneraction point
    int maxHits = std::min(nHits, track1.NHits());
    std::vector<SHit> track1Hits_ = track1.GetHits();
    //std::sort(track1Hits_.begin(), track1Hits_.end(), CompareHitsByMinX);
    std::vector<SHit> track1Hits;
    if (std::abs(intP.X() - track1.GetMinX()) < std::abs(intP.X() - track1.GetMaxX())) {
        track1Hits = std::vector<SHit>(track1Hits_.begin(), std::next(track1Hits_.begin(), maxHits) ); 
    }
    else {
        track1Hits = std::vector<SHit>(std::prev(track1Hits_.end(), maxHits), track1Hits_.end());
    }
    SLinearCluster trk1(track1Hits); 

    maxHits = std::min(nHits, track2.NHits());
    std::vector<SHit> track2Hits;
    std::vector<SHit> track2Hits_ = track2.GetHits();
    //std::sort(track2Hits_.begin(), track2Hits_.end(), CompareHitsByMinX);
    if (std::abs(intP.X() - track2.GetMinX()) < std::abs(intP.X() - track2.GetMaxX())) {
        track2Hits = std::vector<SHit>(track2Hits_.begin(), std::next(track2Hits_.begin(), maxHits) ); 
    }
    else {
        track2Hits = std::vector<SHit>(std::prev(track2Hits_.end(), maxHits), track2Hits_.end());
    }
    SLinearCluster trk2(track2Hits); 
   
    LineEquation track1Eq = trk1.GetTrackEquation();
    
    LineEquation track2Eq = trk2.GetTrackEquation();

    std::vector<SHit> containedHits;

    for (SHit& hit : track2.GetHits()) {
        float yHypo = track1Eq.Slope() * hit.X() + track1Eq.Intercept();
        if (std::abs(hit.Y() - yHypo) < tol * hit.Width()) {
            containedHits.push_back(hit);
        }
    }

    for (SHit& hit : track1.GetHits()) {
        float yHypo = track2Eq.Slope() * hit.X() + track2Eq.Intercept();
        if (std::abs(hit.Y() - yHypo) < tol * hit.Width()) {
            containedHits.push_back(hit);
        }
    }

    return containedHits;
}


// ----- Check arrow line intersection
SPoint TPCLinesVertexFinder::check_arrow_line_intersection(float Ax, float Ay, float Dx, float Dy, float line_slope, float line_intercept) {
    float m = line_slope, b = line_intercept;

    // Calculate t from the intersection equation
    float t = (m * Ax - Ay + b) / (Dy - m * Dx);

    float x = Ax + t * Dx;
    float y = Ay + t * Dy;

    SPoint intersection_point(x, y);

    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "ArrowDir intP " << intersection_point.X() << " " << intersection_point.Y() << std::endl;

    if ((x - Ax >= 0) == (Dx >= 0)) {
        return intersection_point;
    } else {
        return SPoint(-1, -1);
    }
}


// ------- Get the opposite point of a track
SPoint TPCLinesVertexFinder::GetTracksEquationOppositePoint(SLinearCluster track, std::vector<SLinearCluster> trackList, SPoint p){
    float minX1 = 1e3;
    float maxX1 = 0;
    for (SLinearCluster trk : trackList) {
        minX1 = std::min(minX1, trk.GetMinX());
        maxX1 = std::max(maxX1, trk.GetMaxX());
    }

    std::vector<float> Xpoints;
    Xpoints.push_back( std::min(p.X(), minX1) );
    Xpoints.push_back( std::max(p.X(), maxX1) );

    double m1 = track.GetTrackEquation().Slope();
    double n1 = track.GetTrackEquation().Intercept();

    std::vector<float> Ypoints;
    for (auto x : Xpoints) {
        Ypoints.push_back(m1 * x + n1);
    }

    SPoint pA(Xpoints[0], Ypoints[0]);
    SPoint pB(Xpoints[1], Ypoints[1]);

    double dA = TPCLinesDistanceUtils::GetHitDistance(pA, p);
    double dB = TPCLinesDistanceUtils::GetHitDistance(pB, p);


    return  (dA > dB) ? pA : pB;
}


// ------- Get the collinear tracks
std::vector<SLinearCluster> TPCLinesVertexFinder::GetCollinearTracks(SLinearCluster mainTrack, std::vector<SLinearCluster> trackList){

    std::vector<SLinearCluster> pointingTracksV;


    // Get directions pointing to the main track vertex
    for(auto &trk:trackList){
        // intersection point with the main track
        SPoint intP = GetTracksIntersection(trk, mainTrack, 10000);

        // check the intersection point is close enough to the mian track vertex
        double d = TPCLinesDistanceUtils::GetHitDistance(mainTrack.GetStartPoint(), intP);
        if(d<10){
            pointingTracksV.push_back(trk);
        }
        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"   Track "<<trk.GetId()<<" D="<<d<<" IntP: "<<intP;
    }

    // sort by minX
    std::sort(pointingTracksV.begin(), pointingTracksV.end(), [&](SLinearCluster& l1, SLinearCluster& l2) {return l1.GetMinX() < l2.GetMinX();} );    



    std::vector<std::vector<SLinearCluster>> groupedTracks;
    if(pointingTracksV.size()>0){
        // get consecutive segments
        groupedTracks.push_back({pointingTracksV[0]});
        SLinearCluster lastCluster = pointingTracksV[0];
        if(pointingTracksV.size()>=2){
            for (size_t i = 1; i < pointingTracksV.size(); ++i) {
                SLinearCluster currentCluster = pointingTracksV[i];
                
                // Consider the track if thje start and end are within 1 wire
                bool collinear = false;
                if ( std::abs(currentCluster.GetMinX() - lastCluster.GetMaxX()) <= 2){
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"   Contiguous "<<currentCluster.GetId()<<" to "<<lastCluster.GetId()<<std::endl;
                    // check is the edge hit for both tracks is compatible with the connectedness
                    SPoint endHit = lastCluster.GetEndPoint();
                    SPoint startHit = currentCluster.GetStartPoint();

                    double d12 = TPCLinesDistanceUtils::GetHitDistance(endHit, startHit);

                    double averageComp = (currentCluster.GetCompactness()+lastCluster.GetCompactness())/2.;
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"    AvComp: "<<averageComp<<" D12: "<<d12<<std::endl;
                    if(d12<5.*averageComp){
                        collinear=true;
                    }
                }

                if(collinear){
                    groupedTracks.back().push_back(currentCluster);
                }
                else{
                    groupedTracks.push_back({currentCluster});
                }

                lastCluster = currentCluster;

            }
         }
    }

    std::vector<SLinearCluster> finalTrackList;
    for (size_t i = 0; i < groupedTracks.size(); ++i) {
        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" Grouped tracks "<<i<<std::endl;
        std::vector<SHit> hits;
        for (size_t j = 0; j < groupedTracks[i].size(); ++j) {
            if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"  Trk: "<<groupedTracks[i][j].GetId()<<std::endl;
            std::vector<SHit> trk_hits = groupedTracks[i][j].GetHits();
            hits.insert(hits.end(), trk_hits.begin(), trk_hits.end());
        }

        int nHits = std::min(hits.size(), (size_t)5);

        std::sort(hits.begin(), hits.end(), [](SHit& a, SHit& b) {return a.X() < b.X();});
        std::vector<SHit> hitsStart = hits;
        hitsStart.resize(nHits);

        std::sort(hits.begin(), hits.end(), [](SHit& a, SHit& b) {return a.X() > b.X();});
        std::vector<SHit> hitsEnd = hits;
        //hitsEnd.resize(nHits);

        double avIntStart=0;
        double avIntEnd=0;
        for(int k=0; k<nHits; k++){
            if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"S:"<<hitsStart[k].Integral()<<hitsStart[k]<<" E:"<<hitsEnd[k].Integral()<<hitsEnd[k]<<std::endl;
            avIntStart+=hitsStart[k].Integral();
        }
        avIntStart/=nHits;

        for(size_t k=0; k<hitsEnd.size(); k++){
            avIntEnd+=hitsEnd[k].Integral();
        }
        
        avIntEnd/=hitsEnd.size();

        SLinearCluster finalTrack(hits);
        

        // get average hit signal at the beginning

        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<< "  Final track "<<finalTrack.GetMinX()<<" "<<finalTrack.GetMaxX()<<" Integral start/end "<<avIntStart<<" "<<avIntEnd<<std::endl;

        // add ig the energy deposition is larger at the beginning and the main track is larger
        if(avIntStart>1.25*avIntEnd && hits.size()){//<mainTrack.NHits()){
            finalTrackList.push_back(finalTrack);
        }
    }
  
    
                    
    return finalTrackList;
} 


// ------- Get the main direction
SLinearCluster TPCLinesVertexFinder::GetMainDirection(std::vector<SLinearCluster> & mainDirTrackList, 
                                std::vector<SLinearCluster>& freeTracksList,
                                std::vector<std::vector<SLinearCluster>> parallelTracks,
                                bool decideMainTrack,
                                int verbose)
{
    
    // longest
    std::vector<SLinearCluster> LongDirTrackList, LongFreeTracksList;
    SLinearCluster LongDirection = TPCLinesDirectionUtils::GetMainDirectionLongest(parallelTracks, LongDirTrackList, LongFreeTracksList);
    // downstream
    std::vector<SLinearCluster> DownDirTrackList, DownFreeTracksList;
    SLinearCluster DownDirection = TPCLinesDirectionUtils::GetMainDirectionDownstream(parallelTracks, DownDirTrackList, DownFreeTracksList);

    std::vector<int> LongDirectionIndexes;
    if(verbose>=1) std::cout << "Longest track: ";
    for (SLinearCluster &track : LongDirTrackList) {
        if(verbose>=1) std::cout << track.GetId() << " ";
        LongDirectionIndexes.push_back(track.GetId());
    }
    if(verbose>=1) std::cout << std::endl;

    if(verbose>=1) std::cout << "Downstream track: ";
    std::vector<int> DownDirectionIndexes;
    for (SLinearCluster &track : DownDirTrackList) {
        if(verbose>=1) std::cout << track.GetId() << " ";
        DownDirectionIndexes.push_back(track.GetId());
    }
    if(verbose>=1) std::cout << std::endl;

    bool downIsLong = areVectorsEqual(DownDirectionIndexes, LongDirectionIndexes);
    bool useLargest = true;

    if(decideMainTrack==true){
        if (!downIsLong) {
            double connTol = 3 * (LongDirection.GetConnectedness() + DownDirection.GetConnectedness()) / 2;
            double conn = TPCLinesDistanceUtils::GetClusterConnectedness(LongDirection.GetHitCluster(), DownDirection.GetHitCluster());
            SPoint intP = GetTracksIntersection(LongDirection, DownDirection, 20, true);

            if(verbose>=1) std::cout << "Long/main direction intersection: " << intP << " Connectedness: " << conn << " ConnTol: " << connTol << std::endl;

            if (conn > connTol) {
                useLargest = false;
                if(verbose>=1) std::cout << "  longest and most downstream are not connected, using the most downstream" << std::endl;
            }
            else {
                // if they are consecutive, still use the most downstream
                if(DownDirection.GetMaxX()<=LongDirection.GetMinX()){
                    useLargest = false;
                    if(verbose>=1) std::cout << "  longest is consecutive to the downstream, using the most downstream" << std::endl;
                }
                if(verbose>=1) std::cout << " connected, using the longest" << std::endl;
            }
        }
        else {
            if(verbose>=1) std::cout << "Longest/most downstream directions are the same" << std::endl;
        }
    }

    SLinearCluster MainDirection;
    mainDirTrackList.clear();
    freeTracksList.clear();
    if (!useLargest) {
        if(verbose>=1) std::cout << "Using the most downstream track" << std::endl;
        MainDirection = DownDirection;
        mainDirTrackList = DownDirTrackList;
        freeTracksList = DownFreeTracksList;
    }
    else{
        MainDirection = LongDirection;
        mainDirTrackList = LongDirTrackList;    
        freeTracksList = LongFreeTracksList;
    }

    return MainDirection;
}


// ------- Perform kinematic checks for the lambda decay
bool TPCLinesVertexFinder::LambdaDecayKinematicCheck(SOrigin mainOrigin, STriangle Triangle, SLinearCluster MainDirection, SLinearCluster track1, SLinearCluster track2, std::vector<SLinearCluster> FreeTracksList) {
    if (fTPCLinesVertexFinderPset.Verbose >= 1) {
        std::cout << " - - - Performing kinematic checks - - - " << track1.GetId() << " " << track2.GetId() << std::endl;
    }

    // -------- CHECK 1: Check if the triangle direction intersects the main direction
    /*SPoint start_point(Triangle.GetMainVertex().X(), Triangle.GetMainVertex().Y());
    SPoint direction_vector = Triangle.GetDirectorVector();
    double line_slope = MainDirection.GetTrackEquation().Slope();
    double line_intercept = MainDirection.GetTrackEquation().Intercept();
    SPoint intersection_point = check_arrow_line_intersection(start_point.X(), start_point.Y(),
                                                            direction_vector.X(), direction_vector.Y(),
                                                            line_slope, line_intercept);
    bool triangleIntersects = (intersection_point.X() != -1 && intersection_point.Y() != -1);

    if (fTPCLinesVertexFinderPset.Verbose >= 1) {
        std::cout << " - - - Check 1: Triangle direction intersects main direction - - - \n";
        std::cout << "Arrow line intersects the line equation at: (" << intersection_point.X() << ", " << intersection_point.Y() << ")" << std::endl;
        std::cout << "  ** Pass: " << triangleIntersects << std::endl;
    }*/

    // -------- CHECK 2: Junction between the main direction and the triangle vertex is contained within the triangle
    bool junctionContained = TrackTriangleJunctionContained(mainOrigin, Triangle, fTPCLinesVertexFinderPset.AngleTolerance);
    
    if (fTPCLinesVertexFinderPset.Verbose >= 1) {
        std::cout << " - - - Check 2: Junction between main direction origin and triangle vertex is contained - - - \n";
        std::cout << "  ** Pass: " << junctionContained << std::endl;
    }

    // -------- CHECK 3: Check that the junction doesn't cross other tracks and the vertex is not within the main track
    int nHitsInMiddle = GetNHitsBetweenJunction(MainDirection, Triangle, FreeTracksList, 1.5);
   
    bool passJunctionIsFree = (nHitsInMiddle <= 1);

    // set double precision limit
    double dVertexTrack = std::numeric_limits<double>::max();
    double distW_MainHit_MainTrack = std::numeric_limits<double>::max();

    double mainTrackWidthTol = 1. * MainDirection.GetAverageWidth();
    // if the triangle vertex is within the main track start/end, check the drift
    if(Triangle.GetMainVertex().X() > MainDirection.GetMinX() && Triangle.GetMainVertex().X() < MainDirection.GetMaxX()){
        // Get distance to the vertex
        LineEquation trackEq = MainDirection.GetTrackEquation();
        dVertexTrack = trackEq.GetDistance(Triangle.GetMainVertex());
        
        // OldMay passJunctionIsFree = passJunctionIsFree && (dVertexTrack > 1.5 * MainDirection.GetAverageWidth());

        SHit mainVertexHit = Triangle.GetMainVertexHit();
        distW_MainHit_MainTrack = MainDirection.GetHitCluster().GetMinDistanceToClusterW(mainVertexHit);
        passJunctionIsFree = passJunctionIsFree && (distW_MainHit_MainTrack > mainTrackWidthTol);
    }

    if (fTPCLinesVertexFinderPset.Verbose >= 1) {
        std::cout << " - - - Check 3: Junction does not intersect with other tracks - - - \n";
        std::cout << "NHits in the middle: " << nHitsInMiddle << " Pass? " << passJunctionIsFree << std::endl;
        std::cout << "Main track average width: " << MainDirection.GetAverageWidth() << " Distance to the vertex: " << dVertexTrack << std::endl;
        std::cout<<"   -- Vertex hits distanceW to main track: "<<distW_MainHit_MainTrack<<std::endl;
        std::cout << "  ** Pass: " << passJunctionIsFree << std::endl;
    }

    // -------- CHECK 4: EndPoint of the triangle
    
    SHit triangleEndHitA = Triangle.GetTrack1().GetOppositeEdgeHit(Triangle.GetMainVertex());
    SHit triangleEndHitB = Triangle.GetTrack2().GetOppositeEdgeHit(Triangle.GetMainVertex());
    double distW_EndHitA_MainTrack = MainDirection.GetHitCluster().GetMinDistanceToClusterW(triangleEndHitA);
    double distW_EndHitB_MainTrack = MainDirection.GetHitCluster().GetMinDistanceToClusterW(triangleEndHitB);

    bool passTriangleEndHits = (distW_EndHitA_MainTrack > mainTrackWidthTol) && (distW_EndHitB_MainTrack > mainTrackWidthTol);

    if (fTPCLinesVertexFinderPset.Verbose >= 1) {
        std::cout << " - - - Check 4: End hits of the triangle not contained in the main track - - - \n";
        std::cout << " End hit 1: " << triangleEndHitA.X() << " " << triangleEndHitA.Y() << " " << distW_EndHitA_MainTrack << std::endl;
        std::cout << " End hit 2: " << triangleEndHitB.X() << " " << triangleEndHitB.Y() << " " << distW_EndHitB_MainTrack << std::endl;
        std::cout << " Main track wdith tol: " << mainTrackWidthTol << std::endl;
        std::cout << " Pass? "<< passTriangleEndHits << std::endl;
    }

    // ---- CHECK 7: Sides of the triangle
    double lengthB = Triangle.GetSideLenghtB();
    double lengthC = Triangle.GetSideLenghtC();
    double lengthOppo = Triangle.GetOppositeSideLenght();
    double lenghtDeviation = 100 * ( (lengthB + lengthC) - lengthOppo )/ (lengthB + lengthC);
    bool passTriangleInequality = lenghtDeviation > fTPCLinesVertexFinderPset.TriangleInequalityTol;

    if (fTPCLinesVertexFinderPset.Verbose >= 1) {
        std::cout << " - - - Check 7: Sides of the triangle - - - \n";
        std::cout<<" Triangle lenghts B/C "<<lengthB<<" "<<lengthC<<" Oposite "<<lengthOppo<<std::endl;
        std::cout<< "B+C "<<lengthB+lengthC<<" Opposite "<<lengthOppo<<" Deviation: "<<lenghtDeviation<<std::endl;
        std::cout << "  Pass triangle inequality? " << passTriangleInequality << std::endl;
    }

    bool passKinematicChecks = junctionContained && passJunctionIsFree && passTriangleEndHits && passTriangleInequality;

    if (fTPCLinesVertexFinderPset.Verbose >= 1) {
        std::cout << " \n - - - PASS ALL CHECKS? - - - " << passKinematicChecks << std::endl;
    }

    return passKinematicChecks;
}


// ------- Main function: get the vertex candidates
std::vector<SOrigin> TPCLinesVertexFinder::GetAngleVertices(std::vector<SLinearCluster> trackList, SPoint ballVertex, std::vector<STriangle>& vertexList,  std::vector<SOrigin>& associatedOrigins,  SLinearCluster &mainDirection){

    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" In Origin finder\n";
    

    // ------- Fill sliding window line equations for each track
    for(SLinearCluster &trk:trackList){
        trk.FillSlidingWindowLineEquations(fTPCLinesVertexFinderPset.SlidingWindowN);
    }

    // ------- First look for all possible intersections
    // Vector to store the intersections
    std::vector<SOrigin> allIntersections;
    // ------ Loop 1
    if (trackList.size() > 0) {
        for (size_t ix = 0; ix < trackList.size(); ++ix) {
            SLinearCluster track1 = trackList[ix];

            // ------ Loop 2
            for (size_t jx = ix + 1; jx < trackList.size(); ++jx) {
                SLinearCluster track2 = trackList[jx];

                // ------ First check if the tracks are connected
                
                // Calculate connection based on the tracks connectedness
                float connTol =  5 * (track1.GetConnectedness() + track2.GetConnectedness()) / 2;
                float conn = TPCLinesDistanceUtils::GetClusterConnectedness(track1.GetHitCluster(), track2.GetHitCluster());
                bool connected = (conn < connTol); 

                // Don't waste time if tracks are not connected (very loose connected cut here)
                if (!connected) {
                    continue;
                }

                if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "\n----- Potential intersection " << track1.GetId() << " " << track2.GetId() << " Conn: " << conn << " Conn1: "<< track1.GetConnectedness() << " Conn2: " << track2.GetConnectedness()<<" Tol:" << connTol << std::endl;


                // ------ Look for the intersection points
                SPoint intP = GetTracksIntersection(track1, track2, 50, fTPCLinesVertexFinderPset.RefineVertexIntersection);

                if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"  Intersection Point: "<<intP;
                
                if(intP.X()==-1 and intP.Y()==-1) continue;

                // Get closest hit of each track to the intersection point
                std::pair<SHit, double> cloHit1Pair = track1.GetHitCluster().GetClosestHitToPoint(intP);
                SHit cloHit1 = cloHit1Pair.first;
                double dHit1 = cloHit1Pair.second;
                std::pair<SHit, double> cloHit2Pair = track2.GetHitCluster().GetClosestHitToPoint(intP);
                SHit cloHit2 = cloHit2Pair.first;
                double dHit2 = cloHit2Pair.second;
                SHit cloHit = (dHit1 < dHit2) ? cloHit1 : cloHit2;

                double intPYerror = cloHit.Width();

                // Get the edge hits of each track to the intersection point
                //SPoint edgeHit1 = ( std::abs(track1.GetStartPoint().X()-intP.X()) < std::abs(track1.GetEndPoint().X()-intP.X()) )? track1.GetStartPoint():track1.GetEndPoint();
                //SPoint edgeHit2 = ( std::abs(track2.GetStartPoint().X()-intP.X()) < std::abs(track2.GetEndPoint().X()-intP.X()) )? track2.GetStartPoint():track2.GetEndPoint();
                SPoint edgeHit1 = track1.GetEdgeHit(intP);
                SPoint edgeHit2 = track2.GetEdgeHit(intP);
                
                if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"Edge hits: \n"<<edgeHit1<<edgeHit2;
                

                // Get candidate vertex hits, i.e., hits from one track contained in the line hypothesis of the other track
                std::vector<SHit> vertexHits = GetMutualHitsContainedInHypo(track1, track2, intP, fTPCLinesVertexFinderPset.VertexHitsMinHits, fTPCLinesVertexFinderPset.VertexHitsTol);
                if(fTPCLinesVertexFinderPset.Verbose>=1){
                    std::cout << "    NVERTEX HITS " << vertexHits.size() << std::endl;
                    for(SHit & h:vertexHits) std::cout<<h;
                }

                // Check if the intersection point compactness is compatible with the two tracks
                double intPCompactness1 = track1.GetHitCluster().GetMinDistanceToCluster(intP);
                double intPCompactness2 = track2.GetHitCluster().GetMinDistanceToCluster(intP);
                bool intersectionPointCompactnessCompatible = intPCompactness1<fTPCLinesVertexFinderPset.VertexCompactnessTol * track1.GetCompactness() 
                                                            && intPCompactness2<fTPCLinesVertexFinderPset.VertexCompactnessTol * track2.GetCompactness();
                
                if(fTPCLinesVertexFinderPset.Verbose>=1){
                    std::cout<<" Compactness 1/2 "<<track1.GetCompactness()<<" "<<track2.GetCompactness()<<std::endl;
                    std::cout<<" IntPCompactness 1/2 "<<intPCompactness1<<" "<<intPCompactness2<<std::endl;
                    std::cout<<" Intersection point is compactness compatible: "<<intersectionPointCompactnessCompatible<<std::endl;
                }

                // Keep the intersection if
                // there are common vertex hits or
                // the intersection point is compatible with both track compactness
                if(vertexHits.size()==0 && intersectionPointCompactnessCompatible==false){
                    std::cout<<"SKIP... no vertex hit founds or intersection point not compactness compatible\n";
                    continue;
                }

                // Next check the intersection is near rhe edge hits of each track
                float dEdge1 = std::min(std::abs(cloHit1.X() - track1.GetMinX()), std::abs(cloHit1.X() - track1.GetMaxX()));
                float dEdge2 = std::min(std::abs(cloHit2.X() - track2.GetMinX()), std::abs(cloHit2.X() - track2.GetMaxX()));

                for (SHit& hit : vertexHits) {
                    //float min1 = std::min(std::abs(hit.X() - track1.GetMinX()), std::abs(hit.X() - track1.GetMaxX()));
                    //float min2 = std::min(std::abs(hit.X() - track2.GetMinX()), std::abs(hit.X() - track2.GetMaxX()));
                    float min1 = std::abs(hit.X() - edgeHit1.X());
                    if (min1 < dEdge1) {
                        dEdge1 = min1;
                    }

                    float min2 = std::abs(hit.X() - edgeHit2.X());
                    if (min2 < dEdge2) {
                        dEdge2 = min2;
                    }
                }

                if(fTPCLinesVertexFinderPset.Verbose>=1){
                    std::cout << " Closest hits to the intersection 1/2:\n";
                    std::cout << "   " << cloHit1 << "   " << cloHit2;
                    std::cout << " Closest hit " << cloHit;
                    std::cout << " DEdges: " << dEdge1 << " " << dEdge2 << std::endl;
                    std::cout << " Edge hits 1/2:\n";
                    std::cout << "  " << edgeHit1 << "  " << edgeHit2;
                }

                bool inEdgeOrigin;
                bool oneEdgeOrigin;
                
                float maxDEdge = fTPCLinesVertexFinderPset.MaxDistToEdge;

                if( (dEdge2 >= maxDEdge || dEdge1 >= maxDEdge)) {
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "   *** Intersection is not with edge hits" << std::endl;
                    inEdgeOrigin = false;
                    if(dEdge1 >= maxDEdge && dEdge2 >= maxDEdge){
                        oneEdgeOrigin = false;
                    }
                    else{
                        oneEdgeOrigin = true;
                    }
                }
                else{
                    inEdgeOrigin = true;
                }

                if(intersectionPointCompactnessCompatible){
                    maxDEdge = fTPCLinesVertexFinderPset.VertexCompactnessTol * std::max(track1.GetCompactness(), track2.GetCompactness());
                    float dEdge1Comp = std::hypot(edgeHit1.X() - intP.X(), edgeHit1.Y() - intP.Y());
                    float dEdge2Comp = std::hypot(edgeHit2.X() - intP.X(), edgeHit2.Y() - intP.Y());
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"  Max D Edge: "<<maxDEdge<<" Dist 1/2: "<<dEdge1Comp<<" "<<dEdge2Comp<<std::endl;
                    if(dEdge1Comp>maxDEdge || dEdge2Comp>maxDEdge){
                        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout << "   *** Intersection is not with edge hits" << std::endl;
                        inEdgeOrigin = false;
                        if(dEdge1Comp >= maxDEdge && dEdge2Comp >= maxDEdge){
                            oneEdgeOrigin = false;
                        }
                        else{
                            oneEdgeOrigin = true;
                        }
                    }
                }
                if(vertexHits.size()>0){
                    SLinearCluster vertexHitCluster(vertexHits);
            
                    float connVertexHitsCluster1 = TPCLinesDistanceUtils::GetClusterConnectedness(track1.GetHitCluster(), vertexHits);
                    float connVertexHitsCluster2 = TPCLinesDistanceUtils::GetClusterConnectedness(track2.GetHitCluster(), vertexHits);
                    float connOvVertexHitsCluster1 = TPCLinesDistanceUtils::GetClusterConnectednessOverlap(track1.GetHitCluster(), vertexHits);
                    float connOvVertexHitsCluster2 = TPCLinesDistanceUtils::GetClusterConnectednessOverlap(track2.GetHitCluster(), vertexHits);

                    float connXVertexHitsCluster1 = TPCLinesDistanceUtils::GetClusterMinDistanceX(track1.GetHitCluster(), vertexHits);
                    float connXVertexHitsCluster2 = TPCLinesDistanceUtils::GetClusterMinDistanceX(track2.GetHitCluster(), vertexHits);

                    if(fTPCLinesVertexFinderPset.Verbose>=1){
                        std::cout<<"   Vertex hit clister connectedneess: "<<vertexHitCluster.GetConnectedness()<<std::endl;
                        std::cout<<"   Connectedness HitsCluster 1/2: "<<connVertexHitsCluster1<<" "<<connVertexHitsCluster2<<std::endl;
                        std::cout<<"   Connectedness Overlap HitsCluster 1/2: "<<connOvVertexHitsCluster1<<" "<<connOvVertexHitsCluster2<<std::endl;
                        std::cout<<"   Connectedness X HitsCluster 1/2: "<<connXVertexHitsCluster1<<" "<<connXVertexHitsCluster2<<std::endl;
                    }
                    
                    // vertex hits must be connected to the tracks
                    if( connXVertexHitsCluster1>2 || connXVertexHitsCluster2>2){
                        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" VERTEX HITS NOT CONNECTED\n";
                        continue;
                    }   
                }

                // Finally define the vertex
                if(intersectionPointCompactnessCompatible==false){
                    float minD = 1e3;
                    SHit vertexHit(intP.X(), intP.Y());
                    for (SHit& hit : vertexHits) {
                        float d = std::hypot(hit.X() - intP.X(), hit.Y() - intP.Y());
                        if (d < minD) {
                            minD = d;
                            vertexHit = hit;
                        }
                    }
                    intP = SPoint(vertexHit.X(), vertexHit.Y());
                    intPYerror = vertexHit.Width();
                }
    
                if(fTPCLinesVertexFinderPset.Verbose>=1){
                    std::cout << "      Vertex set to: " << intP << std::endl;
                    std::cout << "      V TRACKs IDs " << track1.GetId() << " " << track2.GetId() << std::endl;
                }
                
                // check distance to the PANDORA vertex
                float dIntPBallVertex = std::hypot( 0.3*(intP.X() - ballVertex.X()), 0.075*(intP.Y() - ballVertex.Y()) );
                if(fTPCLinesVertexFinderPset.Verbose>=1){
                    std::cout << "     Distance to PANDORA vertex: " << dIntPBallVertex << std::endl;
                    std::cout << " Edge origin: " << inEdgeOrigin << " One edge origin: " << oneEdgeOrigin << std::endl;
                }
                
                // keep the vertex if is within the PANDORA ROI
                bool isIntersection= (intP.X()!=-1 && intP.Y()!=-1 && dIntPBallVertex<fTPCLinesVertexFinderPset.VertexDistanceROI);

                if(isIntersection && (inEdgeOrigin || oneEdgeOrigin)){
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"    Adding potential intersection..."<<std::endl;
                    SOrigin newOr = SOrigin(intP, {track1, track2}, inEdgeOrigin, intPYerror);
                    allIntersections.push_back(newOr);
                }

            }
        }
    }

    
    // ------ Second study the origin hierarchy
    // Reset and set variables
    std::vector<SOrigin> originList;
    vertexList.clear();
    // Variables for the origin assignment
    std::map<int, bool> usedTrack;
    for(SLinearCluster &trk:trackList) usedTrack[trk.GetId()]=false;
    std::vector<std::pair<int, int>> kinkTracks;

    // Longest track ID
    int longestTrackId = -1;
    int longestTrackNHits = -1;
    for(SLinearCluster &trk:trackList){
        if(trk.NHits()>longestTrackNHits){
            longestTrackNHits = trk.NHits();
            longestTrackId = trk.GetId();
        }
    }

    // First create origins for the intersection within the edges
    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"\n\n MERGING ORIGINS\n";
    for(SOrigin & ori:allIntersections){

        // Get the intersection attributes
        SLinearCluster track1 = ori.GetTrackEntry(0);
        SLinearCluster track2 = ori.GetTrackEntry(1);
        SPoint intP = ori.GetPoint();
        double intPYerror = ori.GetYError();
        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"  Intersection "<<track1.GetId()<<" "<<track2.GetId()<<" X/Y="<<ori.GetPoint().X()<<" "<<ori.GetPoint().Y()<<std::endl;

        // Edge intersection
        bool edgeIntersection = ori.IsEdgeOrigin();
        bool deltaRayIntersection = false;

        // Skip the intersection if it's not within the edge
        // Delta ray check
        if(!edgeIntersection){

            if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"  Intersection not within edges\n";
            
            // delta ray check
            int longTrackId = (track1.NHits()>track2.NHits())? track1.GetId():track2.GetId();
            int shortTrackNHits = (track1.NHits()<track2.NHits())? track1.NHits():track2.NHits();
            if(longTrackId==longestTrackId && shortTrackNHits<=5){
                deltaRayIntersection=true;
            }
            else{
                if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"   Not a potential delta ray...nHits short track:"<<shortTrackNHits<<"\n";
                continue;
            }
        }
        

        if(deltaRayIntersection){
            
            SLinearCluster deltaTrack = (track1.NHits()<track2.NHits())? track1:track2;
            SOrigin newOr(intP, {deltaTrack}, edgeIntersection, intPYerror);
            originList.push_back( newOr );
            if(fTPCLinesVertexFinderPset.Verbose>=1){
                std::cout<<"Adding delta ray origin\n";
                std::cout<<"   END = Adding new origin..."<<newOr;
            }
            
        }
        else{
            // ------ Origin assigment
            // Completeley new origin
            if(usedTrack[track1.GetId()]==false && usedTrack[track2.GetId()]==false){
                SOrigin newOr = SOrigin(intP, {track1, track2}, edgeIntersection, intPYerror);
                originList.push_back( newOr );
                // mark as used tracks
                usedTrack[track1.GetId()]=true;
                usedTrack[track2.GetId()]=true;
                if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"   END = Adding completely new origin..."<<newOr;
            }
            // If some of the tracks already used or both used
            else if(usedTrack[track1.GetId()]==false || usedTrack[track2.GetId()]==false || (usedTrack[track1.GetId()]==true && usedTrack[track2.GetId()]==true)){

                SLinearCluster newTrack = (usedTrack[track1.GetId()]==false) ? track1 : track2;
                SLinearCluster previousTrack = (usedTrack[track2.GetId()]==true) ? track1 : track2;

                bool merged=false;
                for(SOrigin & oriAux:originList){
                    
                    bool intersectionCompatible = std::abs(oriAux.GetPoint().X() - intP.X())<=2;
                    intersectionCompatible = intersectionCompatible && std::abs(oriAux.GetPoint().Y() - intP.Y())<2*oriAux.GetYError();
                    if(intersectionCompatible){
                        if(!oriAux.HasTrackIndex(newTrack.GetId())){
                            oriAux.AddTrack(newTrack, intP, intPYerror);
                            if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"   END = Adding new track to origin \n";
                        }
                        usedTrack[track1.GetId()]=true;
                        usedTrack[track2.GetId()]=true;
                        merged=true;
                    }
                    
                }

                if(!merged){
                    SOrigin newOr(intP, {track1, track2}, edgeIntersection, intPYerror);
                    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"   END = Adding new origin..."<<newOr;
                    originList.push_back( newOr );
                    usedTrack[track1.GetId()]=true;
                    usedTrack[track2.GetId()]=true;
                }

            }
        }
        
    }

    // Second create origins for the not within edges intersections
    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" MERGING KINKS\n";
    for(SOrigin & ori:allIntersections){
        // Skip the intersection if it's within the edges
        if(ori.IsEdgeOrigin()== true) continue;

        // get the intersection attributes
        SLinearCluster track1 = ori.GetTrackEntry(0);
        SLinearCluster track2 = ori.GetTrackEntry(1);

        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"  Potnetial Intersection Kink "<<ori.GetPoint().X()<<" "<<ori.GetPoint().Y()<<" Tracks: "<<ori.GetTrackEntry(0).GetId()<<" "<<ori.GetTrackEntry(1).GetId()<<std::endl;
        // check if the kink origin is comaptible with other origins
        bool merged=false;

        // if not merged, add the kink track to the shorter one
        if(!merged){
            // set the kink track to the shorter one
            int kinkTrackIx;
            int parentKinkTrackIx;
            SLinearCluster kinkTrack;
            SLinearCluster parentKinkTrack;
            if( track1.NHits()>track2.NHits() ){
                kinkTrackIx = track2.GetId();
                kinkTrack = track2;
                parentKinkTrack = track1;
                parentKinkTrackIx = track1.GetId();
            }
            else{
                kinkTrackIx = track1.GetId();
                kinkTrack = track1; 
                parentKinkTrack = track2;
                parentKinkTrackIx = track2.GetId();
            }

            SOrigin newOr(ori.GetPoint(), {kinkTrack, parentKinkTrack}, false, ori.GetYError(), parentKinkTrackIx);
            originList.push_back( newOr );
            usedTrack[kinkTrackIx]=true;
            if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"  Adding kink track "<<kinkTrack.GetId()<<" at "<<ori.GetPoint();
        }
    }
    
    // Third create origins for the unmatched tracks
    if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" MERGING UNMATCHED TRACKS\n";
    for(SLinearCluster & track:trackList){

        float d1 = std::hypot( 0.3*(track.GetStartPoint().X() - ballVertex.X()), 0.075*(track.GetStartPoint().Y() - ballVertex.Y()) );
        float d2 = std::hypot( 0.3*(track.GetEndPoint().X() - ballVertex.X()), 0.075*(track.GetEndPoint().Y() - ballVertex.Y()) );
        float d = std::min(d1, d2);
        SPoint edgePointClosest = (d1<d2)? track.GetStartPoint() : track.GetEndPoint();
        SPoint edgePointFarthest = (d1>d2)? track.GetStartPoint() : track.GetEndPoint();
        if(fTPCLinesVertexFinderPset.Verbose>=1){
            std::cout<<" Edge point closest: "<<edgePointClosest;
            std::cout<<" Edge point farthest: "<<edgePointFarthest;
        }

        double trackLength = std::hypot( 0.3*(track.GetStartPoint().X() - track.GetEndPoint().X()), 0.075*(track.GetStartPoint().Y()-track.GetEndPoint().Y()) );
    
        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<" Track "<<track.GetId()<<" "<<trackLength<<" "<<d<<"\n";
        if(d > fTPCLinesVertexFinderPset.VertexDistanceROI) continue;
        if(trackLength<2) continue;
        
        // add single origin if completely unmatched
        if(usedTrack[track.GetId()]==false){
            bool merged = false;
            for(SOrigin & ori:originList){
                bool intersectionCompatible = std::abs(ori.GetPoint().X() - edgePointClosest.X())<=2;
                intersectionCompatible = intersectionCompatible && std::abs(ori.GetPoint().Y() - edgePointClosest.Y())<2*ori.GetYError();
                if(intersectionCompatible){
                    ori.AddTrack(track, edgePointClosest, 0);
                    if(fTPCLinesVertexFinderPset.Verbose>=1)std::cout<<"   END = Adding to the existing one \n";
                    merged=true;
                }
            }
            if(!merged){
                if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"  Adding new origin for unmatched track "<<track.GetId()<<" at "<<edgePointClosest;
                originList.push_back( SOrigin(edgePointClosest, {track}, true, track.GetAverageWidth()) );
                
            }
            usedTrack[track.GetId()]=true;
            
        }
               
    }


    // ------ Sort origins by NHits
    std::sort(originList.begin(), originList.end(), [](SOrigin& obj1, SOrigin& obj2) {
        return obj1.TotalCharge() > obj2.TotalCharge();
    });

    // ------ Only check kinematics if the intersection is not within the MainDirection
    for(SOrigin &ori:originList){
        // only for multiplicity == 2 interactions 
        if(ori.Multiplicity()!=2) continue;
        if(ori.IsEdgeOrigin()==false) continue;        
        
        std::cout<<ori.Multiplicity()<<std::endl;
        SLinearCluster track1 = ori.GetTrackEntry(0);
        SLinearCluster track2 = ori.GetTrackEntry(1);

        if(fTPCLinesVertexFinderPset.Verbose>=1){
            std::cout<< "Looking for a triangle for tracks "<<track1.GetId()<<" and "<<track2.GetId()<<std::endl;
            std::cout<<"  Goodness: "<<track1.GetTrackEquation().Goodness()<<" "<<track2.GetTrackEquation().Goodness()<<std::endl;
        }
        
        // Only intersetctions between "good" tracks
        if(std::abs(track1.GetTrackEquation().Goodness())>fTPCLinesVertexFinderPset.MinTrackGoodness) continue;
        if(std::abs(track2.GetTrackEquation().Goodness())>fTPCLinesVertexFinderPset.MinTrackGoodness) continue;


        double occ1 = track1.GetOccupancy1D();
        double occ2 = track2.GetOccupancy1D();
        if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"    Occupancy: "<<occ1<<" "<<occ2<<std::endl;
        if(occ1<fTPCLinesVertexFinderPset.MinTrackOccupancy || occ2<fTPCLinesVertexFinderPset.MinTrackOccupancy){
            if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"    SKIP... low occupancy\n";
            continue;
        }

        for(SOrigin &ori1:originList){
            if(ori1.Multiplicity()!=1) continue;
            if(ori1.IsEdgeOrigin()==false) continue;
            
            SLinearCluster mainTrack = ori1.GetTrackEntry(0);

            // Only intersetctions between "good" tracks
            if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"  ** Testing main track "<<mainTrack.GetId()<<" goodness: "<<mainTrack.GetTrackEquation().Goodness()<<std::endl;
            if(std::abs(mainTrack.GetTrackEquation().Goodness())>fTPCLinesVertexFinderPset.MinTrackGoodness) continue;

            double occ3 = mainTrack.GetOccupancy1D();
            if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"    Occupancy 3: "<<occ3<<std::endl;
            if(occ3<fTPCLinesVertexFinderPset.MinTrackOccupancy){
                if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"    SKIP... low occupancy for long track\n";
                continue;
            }

            // check the track length
            double trackLength = std::hypot( 0.3*(mainTrack.GetStartPoint().X() - mainTrack.GetEndPoint().X()), 0.075*(mainTrack.GetStartPoint().Y()-mainTrack.GetEndPoint().Y()) );
            if(trackLength<2) continue;
            
            // Origins associated to main track and the V tracks
            std::vector<SOrigin> mainTrackOrigins;
            std::vector<SOrigin> VTracksOrigins;
            for(SOrigin &ori2:originList){
                for(int k=0; k<ori2.Multiplicity(); k++){
                    if(mainTrack.GetId()==ori2.GetTrackEntry(k).GetId()){
                        mainTrackOrigins.push_back(ori2);
                    }
                    else if( track1.GetId()==ori2.GetTrackEntry(k).GetId() || track2.GetId()==ori2.GetTrackEntry(k).GetId() ){
                        VTracksOrigins.push_back(ori2);
                    }
                }
            }

            std::unordered_set<int> mainConnectedTrackIDs;
            for(SOrigin &ori2:mainTrackOrigins){
                for(int k=0; k<ori2.Multiplicity(); k++){
                    mainConnectedTrackIDs.insert(ori2.GetTrackEntry(k).GetId());
                }
            }            
            std::unordered_set<int> VConnectedTrackIDs;
            for(SOrigin &ori2:VTracksOrigins){
                for(int k=0; k<ori2.Multiplicity(); k++){
                    VConnectedTrackIDs.insert(ori2.GetTrackEntry(k).GetId());
                }
            }


            if(fTPCLinesVertexFinderPset.Verbose>=1){
                std::cout<<"\n    Main track origins: ";
                for( auto &connTrkID:mainConnectedTrackIDs){
                    std::cout<<connTrkID<<" ";
                }
                std::cout<<"\n    V tracks delta ray origins: ";
                for( auto &connTrkID:VConnectedTrackIDs){
                    std::cout<<connTrkID<<" ";
                }
            }

            //main track and the V tracks cannot be connected
            std::unordered_set<int> intersectionOrigins;
            for(int connTrkID:mainConnectedTrackIDs){
                if(VConnectedTrackIDs.find(connTrkID)!=VConnectedTrackIDs.end()){
                    intersectionOrigins.insert(connTrkID);
                }
            }
            
            if(fTPCLinesVertexFinderPset.Verbose>=1){
                std::cout<<"\n    Intersection tracks: ";
                for( auto &connTrkID:intersectionOrigins){
                    std::cout<<connTrkID<<" ";
                }
            }

            bool connectedThroughOthers = !intersectionOrigins.empty();
            if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"    Connected through others "<<connectedThroughOthers<<std::endl;
            if(connectedThroughOthers==true)
                continue;
            

            // ------ Create the triangle object
            SPoint vertex1 = GetTracksEquationOppositePoint( track1, {track1}, ori.GetPoint() );
            SPoint vertex2 = GetTracksEquationOppositePoint( track2, {track2}, ori.GetPoint() );
            STriangle triangle = STriangle(ori.GetPoint(), vertex1, vertex2, SHit(ori.GetPoint().X(), ori.GetPoint().Y() ), track1, track2, mainTrack, track1.GetIntegral(), track2.GetIntegral());

            if(triangle.GetNWires()<=fTPCLinesVertexFinderPset.MinWires){
                if(fTPCLinesVertexFinderPset.Verbose>1) std::cout<<" Below the minimum number of wires\n";
                continue;
            }

            std::vector<SLinearCluster> otherTracks;
            for(SLinearCluster &trk:trackList){
                if(trk.GetId()==ori1.GetTrackEntry(0).GetId()) continue;
                if(trk.GetId()==track1.GetId()) continue;
                if(trk.GetId()==track2.GetId()) continue;
                otherTracks.push_back(trk);
            }
            bool passKinemanicChecks = LambdaDecayKinematicCheck(ori1, triangle, ori1.GetTrackEntry(0), track1, track2, otherTracks);
            
            if(passKinemanicChecks){
                if(fTPCLinesVertexFinderPset.Verbose>=1) std::cout<<"FOUND ORIGIN!\n";    
                vertexList.push_back(triangle);
                associatedOrigins.push_back(ori1);
            }
        
        }

    }


    // ------ Get main track for event display
    std::vector<std::vector<SLinearCluster>> parallelTracks;
    for(SLinearCluster &trk:trackList){
        parallelTracks.push_back({trk});
    }
    std::vector<SLinearCluster> MainDirTrackList;
    std::vector<SLinearCluster> FreeTracksList;
    mainDirection = GetMainDirection(MainDirTrackList, FreeTracksList, parallelTracks, fTPCLinesVertexFinderPset.DecideMainTrack, fTPCLinesVertexFinderPset.Verbose);

    // return origins vector
    return originList;

}











// CHECK 4: Check how many hits of triangle tracks are contained in the main direction hypothesis
    /*bool passAngleTracksNotInMain = true;
    int nhits_track1_inMainDir = GetHitsContainedInLineEquation(MainDirection.GetTrackEquation(), track1.GetHits());
    int nhits_track2_inMainDir = GetHitsContainedInLineEquation(MainDirection.GetTrackEquation(), track2.GetHits());
    passAngleTracksNotInMain = (1.0 * nhits_track1_inMainDir / track1.NHits() < fTPCLinesVertexFinderPset.MaxTrackFractionInMain &&
                                1.0 * nhits_track2_inMainDir / track2.NHits() < fTPCLinesVertexFinderPset.MaxTrackFractionInMain);

    if (fTPCLinesVertexFinderPset.Verbose >= 1) {
        std::cout << " - - - Check 4: Hits contained in the main direction - - - \n";
        std::cout << "NHits of track " << track1.GetId() << " in main direction: " << nhits_track1_inMainDir << std::endl;
        std::cout << "NHits of track " << track2.GetId() << " in main direction: " << nhits_track2_inMainDir << std::endl;
        std::cout << "Pass AngleTracksNotInMain: " << passAngleTracksNotInMain << std::endl;
        std::cout << "  ** Pass: " << passAngleTracksNotInMain << std::endl;
    }*/


    // CHECK 5: Check that the triangle tracks' start/end are not next to the main track edge hits
    /*bool passTriangleEdgesNotInMain = true;
    double distanceBMainStart = TPCLinesDistanceUtils::GetHitDistance(Triangle.GetVertexB(), MainDirection.GetStartPoint());
    double distanceBMainEnd = TPCLinesDistanceUtils::GetHitDistance(Triangle.GetVertexB(), MainDirection.GetEndPoint());
    double distanceCMainStart = TPCLinesDistanceUtils::GetHitDistance(Triangle.GetVertexC(), MainDirection.GetStartPoint());
    double distanceCMainEnd = TPCLinesDistanceUtils::GetHitDistance(Triangle.GetVertexC(), MainDirection.GetEndPoint());
    double meanComp = MainDirection.GetCompactness();
    double _tol = 1.0;

    passTriangleEdgesNotInMain = distanceBMainStart > _tol * meanComp &&
                                distanceBMainEnd > _tol * meanComp &&
                                distanceCMainStart > _tol * meanComp &&
                                distanceCMainEnd > _tol * meanComp;

    if (fTPCLinesVertexFinderPset.Verbose >= 1) {
        std::cout << " - - - Check 5: Triangle tracks are not the main track edge hits - - - \n";
        std::cout << "Triangle vertex BC: " << Triangle.GetVertexB().X() << " " << Triangle.GetVertexB().Y() << " " << Triangle.GetVertexC().X() << " " << Triangle.GetVertexC().Y()
        std::cout << "MainDirectionEdgeHits: " << MainDirection.GetStartPoint().X() << " " << MainDirection.GetStartPoint().Y() << " " << MainDirection.GetEndPoint().X() << " " << MainDirection.GetEndPoint().Y() << std::endl;
        std::cout << "MainDir compactness: " << meanComp << std::endl;
        std::cout << "Distances B: " << distanceBMainStart << " " << distanceBMainEnd << std::endl;
        std::cout << "Distances C: " << distanceCMainStart << " " << distanceCMainEnd << std::endl;
        std::cout << "  ** Pass: " << passTriangleEdgesNotInMain << std::endl;
    } */

    /* // ------ CHECK 6: Check that the opening angle is not 0
    if (fTPCLinesVertexFinderPset.Verbose >= 1) {
        std::cout << " - - - Check 6: Opening angle - - - \n";
        std::cout << " Opening Angle " << Triangle.GetOpeningAngle() << "\n";
    }
    bool passOpeningAngle = Triangle.GetOpeningAngle() > 2 && Triangle.GetOpeningAngle() < 178;
    if (fTPCLinesVertexFinderPset.Verbose >= 1) {
        std::cout << "Pass? " << passOpeningAngle << std::endl;
    }*/


    /* //------- CHECK 7: Triangle goodness
    double triangleMAE = Triangle.GetTriangleMAE();
    std::cout<<"  TRIANGLE MAE "<<triangleMAE<<std::endl;
    bool passTriangleGoodness = triangleMAE<fTPCLinesVertexFinderPset.MinTrackGoodness;
    if(fTPCLinesVertexFinderPset.Verbose>=1){
        std::cout<<" Pass triangle MAE? "<<passTriangleGoodness<<std::endl;
    }*/
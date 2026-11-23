////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleTriangle.cpp
//
// \brief Definition of SimpleTriangle
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCSimpleTriangles.h"

double get_momentum(double ke, double m) {
    return ke * std::sqrt(1 + 2. * m / ke);
}


double GetAngle360(double x, double y) {
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

STriangle::STriangle(){}

STriangle::STriangle(SPoint main_vertex, SPoint vertex_b, SPoint vertex_c, SHit mainhit, SLinearCluster track1, SLinearCluster track2, SLinearCluster mainTrack, double weight_b, double weight_c)
{
    fTrack1 = track1;
    fTrack2 = track2;
    fMainTrack = mainTrack;
    
    fMainVertex = main_vertex;
    fMainVertexHit = mainhit;
    
    if (vertex_b.X() < vertex_c.X()) {
        fVertexB = vertex_b;
        fVertexC = vertex_c;
    } else {
        fVertexB = vertex_c;
        fVertexC = vertex_b;
    }
    
    double xM = (fVertexB.X() + fVertexC.X()) / 2;
    double yM = (fVertexB.Y() + fVertexC.Y()) / 2;
    fMidPoint = SPoint(xM, yM);
    //std::cout << "MidPoint " << xM << " " << yM << std::endl;
    fDirectorVector = SPoint( fMainVertex.X() - xM, fMainVertex.Y() - yM );
    double slope = (yM - fMainVertex.Y()) / (xM - fMainVertex.X());
    double intercept = fMainVertex.Y() - slope * fMainVertex.X();
    fDirection = LineEquation(slope, intercept);

    double slopeB = (fVertexB.Y() - fMainVertex.Y()) / (fVertexB.X() - fMainVertex.X());
    double slopeC = (fVertexC.Y() - fMainVertex.Y()) / (fVertexC.X() - fMainVertex.X());

    double MProton = 938;
    double MPion = 130;

    double fConFactor = (1 / 0.0201293) * 23.6e-6;

    // Triangle directions
    fDir1 = LineEquation(slopeB, fVertexB.Y() - slopeB * fVertexB.X());
    fDir2 = LineEquation(slopeC, fVertexC.Y() - slopeC * fVertexC.X());

    // Hypothesis 1
    double slope1 = (MProton * slopeB + MPion * slopeC) / (MProton + MPion);
    double intercept1 = fMainVertex.Y() - slope1 * fMainVertex.X();
    fMomentumHypo1 = LineEquation(slope1, intercept1);

    // Hypothesis 2
    double slope2 = (MPion * slopeB + MProton * slopeC) / (MProton + MPion);
    double intercept2 = fMainVertex.Y() - slope2 * fMainVertex.X();
    fMomentumHypo2 = LineEquation(slope2, intercept2);

    // Calculate the opening angle in radians
    double VDir1[] = {
        GetVertexB().X() - GetMainVertex().X(),
        GetVertexB().Y() - GetMainVertex().Y()
    };

    double VDir2[] = {
        GetVertexC().X() - GetMainVertex().X(),
        GetVertexC().Y() - GetMainVertex().Y()
    };

    double VDir1Angle = GetAngle360(VDir1[0], VDir1[1]);
    double VDir2Angle = GetAngle360(VDir2[0], VDir2[1]);
    fOpeningAngle = std::min(360 - std::abs(VDir1Angle - VDir2Angle), std::abs(VDir1Angle - VDir2Angle));

    // Calculate momenta
    weight_b *= fConFactor;
    weight_c *= fConFactor;

    double Pb = get_momentum(weight_b, MProton);
    double Pc = get_momentum(weight_c, MPion);
    //std::cout << "Momentum weights hypo 1: " << Pb << ", " << Pc << std::endl;

    slope1 = (Pb * slopeB + Pc * slopeC) / (Pb + Pc);
    intercept1 = fMainVertex.Y() - slope1 * fMainVertex.X();
    fMomentumHypo1 = LineEquation(slope1, intercept1);


    Pb = get_momentum(weight_b, MPion);
    Pc = get_momentum(weight_c, MProton);
    //std::cout << "Momentum weights hypo 2: " << Pb << ", " << Pc << std::endl;

    slope2 = (Pb * slopeB + Pc * slopeC) / (Pb + Pc);
    intercept2 = fMainVertex.Y() - slope2 * fMainVertex.X();
    fMomentumHypo2 = LineEquation(slope2, intercept2);

}

// Get NHits functions
int STriangle::GetNHitsTriangle(){
    return fTrack1.NHits() + fTrack2.NHits();
}

int STriangle::GetNHitsTrack1(){
    return fTrack1.NHits();
}

int STriangle::GetNHitsTrack2(){
    return fTrack2.NHits();
}

int STriangle::GetNHitsMainTrack(){
    return fMainTrack.NHits();
}

// Get all hits function
std::vector<SHit> STriangle::GetAllHits(){
    
    std::vector<SHit> hitList;
    
    for(SHit &h:GetTrack1().GetHits()){
        hitList.push_back(h);
    }

    for(SHit &h:GetTrack2().GetHits()){
        hitList.push_back(h);
    }

    return hitList; 
}


// Get all hits function for main track
std::vector<SHit> STriangle::GetMainTrackAllHits(){
    
    std::vector<SHit> hitList;
    
    for(SHit &h:GetMainTrack().GetHits()){
        hitList.push_back(h);
    }

    return hitList;
}





// Get Length functions
double STriangle::GetLengthTrack1(){
    return fTrack1.GetLength();
}

double STriangle::GetLengthTrack2(){
    return fTrack2.GetLength();
}

double STriangle::GetLengthMainTrack(){
    return fMainTrack.GetLength();
}





double STriangle::GetDecayAngleDifference(){
    
    SPoint p1(fMainTrack.GetMinX(), fMainTrack.GetYatMinX());
    double d1 = TPCLinesDistanceUtils::GetHitDistance(p1, fMainVertex);
    SPoint p2(fMainTrack.GetMaxX(), fMainTrack.GetYatMaxX());
    double d2 = TPCLinesDistanceUtils::GetHitDistance(p2, fMainVertex);
    SPoint mainTrackVertex = (d1 < d2) ? p1 : p2;

    double juntionDirection[] = {
        fMainVertex.X() - mainTrackVertex.X(),
        fMainVertex.Y() - mainTrackVertex.Y()
    };

    double VDir1[] = {
        fVertexB.X() - fMainVertex.X(),
        fVertexB.Y() - fMainVertex.Y()
    };

    double VDir2[] = {
        fVertexC.X() - fMainVertex.X(),
        fVertexC.Y() - fMainVertex.Y()
    };

    double juntionDirectionAngle = GetAngle360(juntionDirection[0], juntionDirection[1]);
    double VDir1Angle = GetAngle360(VDir1[0], VDir1[1]);
    double VDir2Angle = GetAngle360(VDir2[0], VDir2[1]);

    double minAngle = std::min(VDir1Angle, VDir2Angle);
    double maxAngle = std::max(VDir1Angle, VDir2Angle);

    bool junctionContained = false;
    
    if (maxAngle - minAngle < 180) {
        junctionContained = (minAngle<=juntionDirectionAngle) && (juntionDirectionAngle<=maxAngle);
    }
    else {
        junctionContained = (juntionDirectionAngle<=minAngle) || (juntionDirectionAngle>=maxAngle);
    }

    if(junctionContained) return 0;
    else{
        double angle1 = std::min(360 - std::abs(juntionDirectionAngle - VDir1Angle), std::abs(juntionDirectionAngle - VDir1Angle));
        double angle2 = std::min(360 - std::abs(juntionDirectionAngle - VDir2Angle), std::abs(juntionDirectionAngle - VDir2Angle));
        return std::min(angle1, angle2);
    }

}


double STriangle::GetGap(){

    SPoint p1(fMainTrack.GetMinX(), fMainTrack.GetYatMinX());
    double d1 = TPCLinesDistanceUtils::GetHitDistance(p1, fMainVertex);
    SPoint p2(fMainTrack.GetMaxX(), fMainTrack.GetYatMaxX());
    double d2 = TPCLinesDistanceUtils::GetHitDistance(p2, fMainVertex);
    SPoint mainTrackVertex = (d1 < d2) ? p1 : p2;

    return std::hypot( 0.3*(mainTrackVertex.X()-fMainVertex.X()), 0.075*(mainTrackVertex.Y()-fMainVertex.Y()) );

}


int STriangle::GetNWires(){
    int minX = std::min(fTrack1.GetMinX(), fTrack2.GetMinX());
    int maxX = std::max(fTrack1.GetMaxX(), fTrack2.GetMaxX());
    return maxX - minX +1;
}


double distanceToLine(SHit hit, double slope, double intercept) {
    // Calculate the distance using the formula |Ax + By + C| / sqrt(A^2 + B^2)
    double A = slope;
    double B = -1; // For a line in the form y = mx + c, B is -1
    double C = intercept;
    
    double numerator = std::abs(A * hit.X() + B * hit.Y() + C);
    double denominator = std::sqrt(A * A + B * B);
    
    return numerator / denominator;
}

// Function to calculate Mean Absolute Error (MAE) for a linear regression model
double calculateSideMAE(const std::vector<SHit>& data, double slope, double intercept) {
    double sum = 0.0;
    for (const SHit& point : data) {
        sum += distanceToLine( point, slope, intercept );
    }
    return sum;
}

double STriangle::GetTriangleMAE(){

    std::vector<SHit> hits1 = fTrack1.GetHits();
    std::vector<SHit> hits2 = fTrack2.GetHits();

    double mae1 = calculateSideMAE(hits1, fDir1.Slope(), fDir1.Intercept());
    double mae2 = calculateSideMAE(hits2, fDir2.Slope(), fDir2.Intercept());
    double mae = (fTrack1.GetIntegral()*mae1+fTrack2.GetIntegral()*mae2) / (fTrack1.GetIntegral() + fTrack2.GetIntegral());

    double meanY = (fTrack1.GetIntegral()*fTrack1.GetMeanY()+fTrack2.GetIntegral()*fTrack2.GetMeanY()) / (fTrack1.GetIntegral() + fTrack2.GetIntegral());
    std::cout<<" MeanY: "<<meanY<<std::endl;
    double mae1Cte = calculateSideMAE(hits1, 0, meanY);
    double mae2Cte = calculateSideMAE(hits2, 0, meanY);
    double maeCte = (fTrack1.GetIntegral()*mae1Cte+fTrack2.GetIntegral()*mae2Cte) / (fTrack1.GetIntegral() + fTrack2.GetIntegral());

    std::cout<<" TriangleCTE "<<mae<<" cte="<<maeCte<<std::endl;

    return 0.5*(mae1+mae2);
}


// Bresenham's line drawing algorithm for steps in X
void drawLine(int x0, double y0, int x1, double y1, std::vector<SPoint>& pointList) {

    // slope and intercept
    double slope = (y1 - y0) / (x1 - x0);
    double intercept = y0 - slope * x0;

    // min and max x
    int minX = std::min(x0, x1);
    int maxX = std::max(x0, x1);

    // loop over x
    for (int x = minX; x <= maxX; x++) {
        double y = slope * x + intercept;
        pointList.push_back(SPoint((double)x, y));
    }
}

std::vector<SPoint> getPointsInTriangle(int x0, int y0, int x1, int y1, int x2, int y2) {
    std::vector<SPoint> trianglePoints;

    // Sort vertices in ascending order based on y-coordinate
    if (y0 > y1) std::swap(x0, x1), std::swap(y0, y1);
    if (y0 > y2) std::swap(x0, x2), std::swap(y0, y2);
    if (y1 > y2) std::swap(x1, x2), std::swap(y1, y2);

    // Draw lines for the three edges of the triangle
    drawLine(x0, y0, x1, y1, trianglePoints);
    drawLine(x1, y1, x2, y2, trianglePoints);
    drawLine(x2, y2, x0, y0, trianglePoints);

    return trianglePoints;
}

double getSegmentOverlap(double a1, double a2, double b1, double b2) {
    double overlapEnd=0;
    double overlapStart=0;
    
    bool overlap = false;
    //if ((a1 <= b1 && b1 <= a2) || (a1 <= b2 && b2 <= a2)) {
    if( (a1 <= b1 && b1 <= a2) || (b1 <= a1 && a1 <= b2)) {
        double overlap_start = std::max(a1, b1);
        double overlap_end = std::min(a2, b2);
        
        if (overlap_start < overlap_end) {
            overlapStart = overlap_start;
            overlapEnd = overlap_end;
            overlap = true;
        }
    }
    
    double overlapLength = (overlap)? std::abs(overlapEnd - overlapStart):0;

    //std::cout << "        Checking overlap: "<<a1<<" "<<a2<<" - "<<b1<<" "<<b2<<" LengthOverelap: "<<overlapLength<<std::endl; 

    return overlapLength;
}

double STriangle::ComputeCoveredArea() {
   
    // get all the points that define the triangle
    std::vector<SPoint> pointsInTriangle = getPointsInTriangle(fMainVertex.X(), fMainVertex.Y(), fVertexB.X(), fVertexB.Y(), fVertexC.X(), fVertexC.Y());

    // get min and max X
    int minX = std::min(fMainVertex.X(), std::min(fVertexB.X(), fVertexC.X()));
    int maxX = std::max(fMainVertex.X(), std::max(fVertexB.X(), fVertexC.X()));

    // create map with min/max Y value for each X step
    std::map<int, double> minY;
    std::map<int, double> maxY;

    // initialize min/max Y values for each X step
    for (int i = minX; i <= maxX; i++) {
        minY[i] = 1e9;
        maxY[i] = -1e9;
    }

    // loop over the points in the triangle
    for (size_t i = 0; i < pointsInTriangle.size(); i++) {
        int x = pointsInTriangle[i].X();
        double y = pointsInTriangle[i].Y();

        // update min/max Y values for each X step
        if (y < minY[x]) minY[x] = y;
        if (y > maxY[x]) maxY[x] = y;
    }

    // area of the triangle
    double totalArea = 0;
    for (int i = minX; i <= maxX; i++) {
        //std::cout << "X = " << i << ", minY = " << minY[i] << ", maxY = " << maxY[i] << std::endl;
        totalArea+= maxY[i] - minY[i];
    }


    // Get hits from the triangle tracks
    std::vector<SHit> hitList;
    std::vector<SHit> hitListAux = fTrack1.GetHits();
    hitList.insert(hitList.end(), hitListAux.begin(), hitListAux.end());
    hitListAux.clear();
    hitListAux = fTrack2.GetHits();
    hitList.insert(hitList.end(), hitListAux.begin(), hitListAux.end());

    // loop over the hits
    double totalOverlap = 0;
    for(SHit& h:hitList){
        double overlap = getSegmentOverlap(h.Y()-h.Width(), h.Y()+h.Width(), minY[h.X()], maxY[h.X()]);
        totalOverlap+=overlap;
    }


    return totalOverlap/totalArea;
}

double STriangle::GetMinX(){
    return std::min( fTrack1.GetMinX(), fTrack2.GetMinX() );
}

double STriangle::GetMaxX(){
    return std::max( fTrack1.GetMaxX(), fTrack2.GetMaxX() );
}

double STriangle::GetMinY(){
    return std::min( fTrack1.GetMinY(), fTrack2.GetMinY() );
}

double STriangle::GetMaxY(){
    return std::max( fTrack1.GetMaxY(), fTrack2.GetMaxY() );
}

double STriangle::GetOverlapWithMainTrack(){
    // min/max x of the maint track
    int minX = fMainTrack.GetMinX();
    int maxX = fMainTrack.GetMaxX();

    // min/max x of the triangle
    int minXTriangle = std::min(fMainVertex.X(), std::min(fVertexB.X(), fVertexC.X())  );
    int maxXTriangle = std::max(fMainVertex.X(), std::max(fVertexB.X(), fVertexC.X())  );

    // get overlap length in X
    return getSegmentOverlap(minX, maxX, minXTriangle, maxXTriangle);

}


std::vector<int> STriangle::GetClusterIDs(){

    std::vector<int> clusterIDs;
    std::vector<int> clusterIDsAux = fTrack1.GetHitCluster().GetClusterIDs();
    clusterIDs.insert(clusterIDs.end(), clusterIDsAux.begin(), clusterIDsAux.end());
    clusterIDsAux.clear();
    clusterIDsAux = fTrack2.GetHitCluster().GetClusterIDs();
    clusterIDs.insert(clusterIDs.end(), clusterIDsAux.begin(), clusterIDsAux.end());

    return clusterIDs;

}


std::vector<int> STriangle::GetMainTrackClusterIDs(){

    std::vector<int> clusterIDs;
    std::vector<int> clusterIDsAux = fMainTrack.GetHitCluster().GetClusterIDs();
    clusterIDs.insert(clusterIDs.end(), clusterIDsAux.begin(), clusterIDsAux.end());

    return clusterIDs;

}

std::vector<int> STriangle::GetTrack1ClusterIDs(){

    std::vector<int> clusterIDs;
    std::vector<int> clusterIDsAux = fTrack1.GetHitCluster().GetClusterIDs();
    clusterIDs.insert(clusterIDs.end(), clusterIDsAux.begin(), clusterIDsAux.end());

    return clusterIDs;

}

std::vector<int> STriangle::GetTrack2ClusterIDs(){

    std::vector<int> clusterIDs;
    std::vector<int> clusterIDsAux = fTrack2.GetHitCluster().GetClusterIDs();
    clusterIDs.insert(clusterIDs.end(), clusterIDsAux.begin(), clusterIDsAux.end());

    return clusterIDs;

}


void STriangle::GetVertexXYZ(double &x, double &y, double &z){

    x = 0;
    y = 0;
    z = 250;

    double minDToVertex = 1e9;
    for(SHit& h:GetAllHits()){
    
        if(h.SPZ()>0){
            double d = std::hypot( fMainVertex.X() - h.X(), fMainVertex.Y() - h.Y() );
        
            if(d<minDToVertex){
                minDToVertex = d;
                x = h.SPX();
                z = h.SPZ();
                y = h.SPY();
            }
        }

    }

    return;
}


void STriangle::GetPrimaryLeptonVertexXYZ(double &x, double &y, double &z){

    x = 0;
    y = 0;
    z = 250;

    SPoint p1(fMainTrack.GetMinX(), fMainTrack.GetYatMinX());
    double d1 = TPCLinesDistanceUtils::GetHitDistance(p1, fMainVertex);
    SPoint p2(fMainTrack.GetMaxX(), fMainTrack.GetYatMaxX());
    double d2 = TPCLinesDistanceUtils::GetHitDistance(p2, fMainVertex);
    SPoint mainTrackVertex = (d1 < d2) ? p1 : p2;

    double minDToVertex = 1e9;
    for(SHit& h:GetMainTrackAllHits()){
        
        if(h.SPZ()>0){
            double d = std::hypot( mainTrackVertex.X() - h.X(), mainTrackVertex.Y() - h.Y() );

            if(d<minDToVertex){
                minDToVertex = d;
                x = h.SPX();
                z = h.SPZ();
                y = h.SPY();
            }
        }

    }

    return;
}
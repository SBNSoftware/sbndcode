////////////////////////////////////////////////////////////////////////////
//
// \file DirectionRecoUtils.h
//
// \brief Definition of direction reconstruction
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef DIRECTIONRECO_UTILS_H
#define DIRECTIONRECO_UTILS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <map>


#if LAMBDAANA_LARSOFT == 1
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/SObjects/TPCSimpleHits.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/SObjects/TPCSimpleClusters.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/SObjects/TPCSimpleLines.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/SObjects/TPCLinesDistanceUtils.h"
#else
#include "TPCSimpleHits.h"
#include "TPCSimpleClusters.h"
#include "TPCSimpleLines.h"
#include "TPCLinesDistanceUtils.h"
#endif





namespace TPCLinesDirectionUtils{

    std::vector<SLinearCluster> SlopeTrackMerger(std::vector<SLinearCluster> trackList, double distTh, double slopeTh, int verbose=0);


    bool GetLineHypoDistance(SLinearCluster mergeTrack, SLinearCluster track, int tol=1, int verbose=0);


    bool FullTrackContained(SLinearCluster mergeTrack, SLinearCluster track, double tol = 1.0);


    std::pair<bool, int> GetNHitsInHypo(SLinearCluster track1, SLinearCluster track2, int fNHits=6, double fWidthTol=1.0);


    bool check_overlap_and_find_region(int a1, int a2, int b1, int b2, int & overlapStart, int & overlapEnd);


    std::vector<SLinearCluster> GetVertexTracks(
        std::vector<SLinearCluster> trackList,
        std::map<int, std::vector<int>>& shortToLongMap,
        std::map<int, int> & connectionsMap,
        int maxHitsShortTrack,
        float dTol,
        float connTolEps,
        int verbose=0
    );



    std::vector<std::vector<SLinearCluster>> GetParallelTracks(
        std::vector<SLinearCluster>& trackList, double dist1DTh, double fAngleTh, double fMaxDWires, int verbose=0);



    SLinearCluster GetMainDirectionMaxHits(std::vector<std::vector<SLinearCluster>> trackClusterList, std::vector<SLinearCluster> & selectedTracksList, std::vector<SLinearCluster> & freeTracksList );


    SLinearCluster GetMainDirectionLongest(std::vector<std::vector<SLinearCluster>> trackClusterList, std::vector<SLinearCluster> & selectedTracksList, std::vector<SLinearCluster> & freeTracksList, int minHits=6);

    SLinearCluster GetMainDirectionDownstream(std::vector<std::vector<SLinearCluster>> trackClusterList, std::vector<SLinearCluster> & selectedTracksList, std::vector<SLinearCluster> & freeTracksList, int minHits=6);


}





#endif // TPC_SIMPLE_CLUSTERS_H


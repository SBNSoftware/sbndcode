////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesAlgo.h
//
// \brief Definition of TPCLinesAlgo
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_LINES_ALGO_H
#define TPC_LINES_ALGO_H

#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <unordered_set>
#if LAMBDAANA_LARSOFT == 1
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/SObjects/TPCSimpleHits.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/SObjects/TPCSimpleClusters.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/SObjects/TPCSimpleTriangles.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/SObjects/TPCSimpleEvents.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/SObjects/LambdaAnaTree.h"

#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/TPCLinesAlgo/TPCLinesParameters.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/TPCLinesAlgo/TPCLinesHough.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/TPCLinesAlgo/TPCLinesTrackFinder.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/TPCLinesAlgo/TPCLinesVertexFinder.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/TPCLinesAlgo/TPCSimpleCalo.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/TPCLinesAlgo/TPCLinesDisplay.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/TPCLinesAlgo/TPCLinesDirectionRecoUtils.h"

#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/FRANS/ChargeDensity.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/src/FRANS/ChargeDensityPset.h"

#else
#include "TPCSimpleHits.h"
#include "TPCSimpleClusters.h"
#include "TPCSimpleTriangles.h"
#include "TPCSimpleEvents.h"
#include "LambdaAnaTree.h"

#include "TPCLinesParameters.h"
#include "TPCLinesHough.h"
#include "TPCLinesTrackFinder.h"
#include "TPCLinesVertexFinder.h"
#include "TPCSimpleCalo.h"
#include "TPCLinesDisplay.h"
#include "TPCLinesDirectionRecoUtils.h"

#include "ChargeDensity.h"
#include "ChargeDensityPset.h"


#endif





class TPCLinesAlgo {
    private:
        // Algorithm parameters
        TPCLinesAlgoPsetType fTPCLinesPset;

        // Input variables
        size_t fNTotalHits;
        std::vector<SHit> fHitList;
        std::vector<SHit> fHitListOutROI;
        std::vector<SHit> fUnmatchedHits;
        std::map<int, int> fClusterIdCounter;
        SVertex fVertex;
        SVertex fVertexTrue;
        int fMinX;
        double fMinY;
        int fMaxX;

        std::vector<SLinearCluster> fFinalTrackList;
        SLinearCluster fMainDirection;
        std::vector<STriangle> fAngleList;
        std::vector<SOrigin> fOrigins;

        // --- Algorithms ---
        // Hough algorithm
        TPCLinesHough fHoughAlgo;

        // Track finder and clustering algorithm
        TPCLinesTrackFinder fTrackFinder;

        // Vertex finder algorithm
        TPCLinesVertexFinder fVertexFinder;

        // FRANS
        ChargeDensity fFRANSAlgo;

        // Calorimetry
        STriangleCalo fTriangleCalo;

        bool fHasTriangle;


        // Display app directory
        TPCLinesDisplay fDisplay;
        std::string fDisplayAppPath;

        int GetWireDistance(int ch1, int ch2);

        std::vector<SHit> RemoveIsolatedHits(std::vector<SHit> hitListForHough, std::vector<SHit>& discardedHits, double maxD, int minNeighbours);

        std::vector<SLinearCluster> MergeIsolatedHits(std::vector<SLinearCluster> recoTrackList, std::vector<SHit> hitList, double dCleaning1D, std::vector<SHit> & discardedHits, double dTh = 3);

        std::vector<SLinearCluster> MergeOutOfROIHits(std::vector<SLinearCluster> recoTrackList, std::vector<SHit> hitList, std::vector<SHit> & remainingHitList, int xStepTol);

        std::vector<SLinearCluster> CreateOutOfROIClusters(std::vector<SLinearCluster> recoTrackList);

        std::vector<SHit> GetHitsInCluster(int clusterId);

        int GetNOverlappedChannelsWithAPABadChannels(int minCh, int maxCh);

        SEvent fRecoEvent;
        SPoint fMainVertex;

        // APA gap
        std::vector<int> fBadChannelsTPC0;
        std::vector<int> fBadChannelsTPC1;

    public:
        // --- Constructor ---
        TPCLinesAlgo(TPCLinesAlgoPsetType tpcLinesAlgoPset, FRANSPsetType fransPset);
        
        // --- Function to set the input variables ---
        bool SetHitList(int view,
                        std::vector<int>& vertex,
                        std::vector<int>& vertexTrue,
                        std::vector<SHit> hits);
        
        // --- Main analysis function ---
        void AnaView(std::string eventLabel);

        STriangle FillLambdaAnaTree(LambdaAnaTree &lambdaAnaTree);

        // --- Getters ---
        // Get average desnity of hits per wire
        double GetAverageHitDensity();
        
        // Get best view (best average chi2)
        int GetBestView(std::vector<int> *_Ch, std::vector<double> *_Chi2);
        
        // Get X-Y shifts
        int ShiftX(){ return fMinX; };
        double ShiftY(){ return fMinY; };
        
        // Get N input hits
        int GetNInputHits() {return fNTotalHits; };

        // Get reco event
        SEvent GetRecoEvent() { return fRecoEvent; };
        SPoint GetMainVertex() { return fMainVertex; };

        // --- Display function ---
        void Display(std::string name, TCanvas *canvas = nullptr, TCanvas *canvasCalo = nullptr, TCanvas *canvasFRANS = nullptr);

        
};

#endif // TPC_SIMPLE_LINES_H
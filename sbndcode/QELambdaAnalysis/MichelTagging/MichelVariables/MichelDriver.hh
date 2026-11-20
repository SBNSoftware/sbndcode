#include "sbndcode/QELambdaAnalysis/MichelTagging/MichelVariables/MichelObj.hh"
#include "sbndcode/QELambdaAnalysis/MichelTagging/MichelVariables/MichelGeo.hh"
#include "sbndcode/QELambdaAnalysis/MichelTagging/MichelVariables/MichelChi2.hh"
#include "sbndcode/QELambdaAnalysis/MichelTagging/MichelVariables/MichelAngles.hh"
#include "sbndcode/QELambdaAnalysis/MichelTagging/MichelVariables/MichelRadial.hh"
#include "sbndcode/QELambdaAnalysis/MichelTagging/MichelVariables/MichelBlips.hh"
#include "sbndcode/QELambdaAnalysis/MichelTagging/MichelVariables/MichelClassifier.hh"

// ROOT and C++ includes
#include <TTree.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TInterpreter.h"
#include "TTimeStamp.h"
#include "TMVA/RBDT.hxx"

#include <string.h>
#include <vector>
#include <limits>
#include <map>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <optional>


class MichelDriver {
private:
    MichelParams fParams;
    std::vector<MuTrack> fmuTracks;
    std::vector<MuHit> fmuHits;
    std::vector<MuCalo> fmuCalos;
    std::vector<MuBlip> fmuBlips;  

    MichelClassifier& fmichelReader;

    double frecoVtxX;
    double frecoVtxY;
    double frecoVtxZ; 

    // Reco features
    Bragg fbragg;
    PCAResult fpca;

    // Michel features
    MichelFeatures fmichelFeatures;
    float fMichelScore;

    // Angle features
    std::vector<double> fangles;

    // Radius features
    std::vector<double> fradii;
    std::vector<double> fcharge;
    std::vector<double> fdQdr;
    std::vector<double> fsmoothdQdr;
    int fdropIdx;
    double fmindQdr;

public:

    // Constructor
    MichelDriver(
        const MichelParams& params, 
        const std::vector<MuTrack>& muTracks, 
        const std::array<std::vector<MuHit>, 3>& muHits, 
        const std::array<std::vector<MuCalo>, 3>& muCalos, 
        const std::vector<MuBlip>& blips,
        const std::array<double, 3>& recoVtx,
        MichelClassifier& michelReader,
        int plane=2)
        : fParams(params)
        , fmuTracks(muTracks)
        , fmuHits(muHits[plane])
        , fmuCalos(muCalos[plane])
        , fmuBlips(blips) 
        , fmichelReader(michelReader)
        , frecoVtxX(recoVtx[0])
        , frecoVtxY(recoVtx[1])
        , frecoVtxZ(recoVtx[2]) 
    {}

    // --------- Getters --------------
    // std::vector<MuTrack> MuTracks() const { return fmuTracks; }
    // std::array<std::vector<MuHit>, 3> MuHits() const { return fmuHits; }
    // std::array<std::vector<MuCalo>, 3> MuCalos() const { return fmuCalo; }
    // std::vector<MuBlip> MuBlips() const { return fblips; }

    MichelFeatures Features() const {return fmichelFeatures; }
    Bragg getBragg() const { return fbragg; }
    PCAResult getPca() const { return fpca; }

    std::vector<double> getAngles() { return fangles; }

    std::vector<double> getRadii() { return fradii; }
    std::vector<double> getCharge() { return fcharge; }
    std::vector<double> getdQdr() { return fdQdr; }
    std::vector<double> getSmoothdQdr() { return fsmoothdQdr; }
    int getDropIdx() { return fdropIdx; }
    double getMindQdr() { return fmindQdr; }

    inline float getScore() { return fMichelScore; }

    // ---------- MAIN FUNCTION --------
    void Compute();    

};
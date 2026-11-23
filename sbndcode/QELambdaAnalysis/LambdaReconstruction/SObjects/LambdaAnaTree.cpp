////////////////////////////////////////////////////////////////////////////
//
// \file TPCLineParameters.h
//
// \brief Definition of LambdaAnaTree
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "LambdaAnaTree.h"

LambdaAnaTree::LambdaAnaTree()
: fTree(nullptr)
{}

LambdaAnaTree::LambdaAnaTree(TTree* tree, bool readMode)
: fTree(tree),
  fInputFileNameRead(new std::string)
{
    InitializeTree(readMode);
}

void LambdaAnaTree::SetTree(TTree* tree, bool readMode)
{
    fTree = tree;
    InitializeTree(readMode);
}

template <typename T>
void LambdaAnaTree::SetBranch(const char* name, T* variable, bool setAddress) {
    if (setAddress) {
        fTree->SetBranchAddress(name, variable);
        //std::cout << "Setting branch address for " << name << std::endl;
    } else {
        fTree->Branch(name, variable);
    }
}

// instance types
template void LambdaAnaTree::SetBranch(const char* name, int* variable, bool setAddress);
template void LambdaAnaTree::SetBranch(const char* name, float* variable, bool setAddress);
template void LambdaAnaTree::SetBranch(const char* name, double* variable, bool setAddress);
template void LambdaAnaTree::SetBranch(const char* name, bool* variable, bool setAddress);
template void LambdaAnaTree::SetBranch(const char* name, std::string* variable, bool setAddress);
template void LambdaAnaTree::SetBranch(const char* name, std::vector<float>* variable, bool setAddress);
template void LambdaAnaTree::SetBranch(const char* name, std::vector<double>* variable, bool setAddress);


void LambdaAnaTree::InitializeTree(bool readMode){
    // Set branch addresses for event information
    SetBranch("EventID", &fEventID, readMode);
    SetBranch("SubrunID", &fSubrunID, readMode);
    SetBranch("RunID", &fRunID, readMode);
    SetBranch("SliceID", &fSliceID, readMode);
    SetBranch("mctID", &fmctID, readMode);

    // Set branch addresses for true variables
    SetBranch("IntOrigin", &fIntOrigin, readMode);
    SetBranch("IntMode", &fIntMode, readMode);
    SetBranch("IntDirt", &fIntDirt, readMode);
    SetBranch("IntType", &fIntType, readMode);
    SetBranch("IntCCNC", &fIntCCNC, readMode);
    SetBranch("IntNuPDG", &fIntNuPDG, readMode);
    SetBranch("IntNuEnergy", &fIntNuEnergy, readMode);
    SetBranch("IntNProtons", &fIntNProtons, readMode);
    SetBranch("IntNNeutrons", &fIntNNeutrons, readMode);
    SetBranch("IntNPi0", &fIntNPi0, readMode);
    SetBranch("IntNPip", &fIntNPip, readMode);
    SetBranch("IntNPim", &fIntNPim, readMode);
    SetBranch("IntNMuonP", &fIntNMuonP, readMode);
    SetBranch("IntNMuonM", &fIntNMuonM, readMode);
    SetBranch("IntNElectronP", &fIntNElectronP, readMode);
    SetBranch("IntNElectronM", &fIntNElectronM, readMode);
    SetBranch("IntNLambda", &fIntNLambda, readMode);
    SetBranch("IntNSigma0", &fIntNSigma0, readMode);
    SetBranch("IntNSigmaPlus", &fIntNSigmaPlus, readMode);
    SetBranch("IntNSigmaMinus", &fIntNSigmaMinus, readMode);
    SetBranch("IntNHyperons", &fIntNHyperons, readMode);

    // Set branch addresses for true vertex information
    SetBranch("NuvE", &fNuvE, readMode);
    SetBranch("NuvT", &fNuvT, readMode);
    SetBranch("NuvX", &fNuvX, readMode);
    SetBranch("NuvY", &fNuvY, readMode);
    SetBranch("NuvZ", &fNuvZ, readMode);
    SetBranch("TruthIsFiducial", &fTruthIsFiducial, readMode);
    SetBranch("TruthIsAV", &fTruthIsAV, readMode);

    // Set branch addresses for slice information
    SetBranch("SliceCompleteness", &fSliceCompleteness, readMode);
    SetBranch("SlicePurity", &fSlicePurity, readMode);

    // Set branch addresses for main MC particle information
    SetBranch("MainMCParticlePDG", &fMainMCParticlePDG, readMode);
    SetBranch("MainMCParticleEnergy", &fMainMCParticleEnergy, readMode);
    SetBranch("MainMCParticleStartX", &fMainMCParticleStartX, readMode);
    SetBranch("MainMCParticleStartY", &fMainMCParticleStartY, readMode);
    SetBranch("MainMCParticleStartZ", &fMainMCParticleStartZ, readMode);
    SetBranch("MainMCParticleStartT", &fMainMCParticleStartT, readMode);
    SetBranch("MainMCParticleEndX", &fMainMCParticleEndX, readMode);
    SetBranch("MainMCParticleEndY", &fMainMCParticleEndY, readMode);
    SetBranch("MainMCParticleEndZ", &fMainMCParticleEndZ, readMode);
    SetBranch("MainMCParticleEndT", &fMainMCParticleEndT, readMode);
    SetBranch("MainMCParticleLength", &fMainMCParticleLength, readMode);
    SetBranch("MainMCParticleNDaughters", &fMainMCParticleNDaughters, readMode);
    SetBranch("MainMCParticleMichelDecay", &fMainMCParticleMichelDecay, readMode);

    // Set branch addresses for lambda true information
    SetBranch("Gap", &fGap, readMode);
    SetBranch("LambdaKE", &fLambdaKE, readMode);
    SetBranch("ProtonKE", &fProtonKE, readMode);
    SetBranch("PionKE", &fPionKE, readMode);
    SetBranch("LambdaBeta", &fLambdaBeta, readMode);
    SetBranch("LambdaGamma", &fLambdaGamma, readMode);
    SetBranch("LambdaDecayTransverseMomentum", &fLambdaDecayTransverseMomentum, readMode);
    SetBranch("LambdaDecayLongitudinalAsymmetry", &fLambdaDecayLongitudinalAsymmetry, readMode);
    SetBranch("LambdaX", &fLambdaX, readMode);
    SetBranch("LambdaY", &fLambdaY, readMode);
    SetBranch("LambdaZ", &fLambdaZ, readMode);
    SetBranch("HasLambda", &fHasLambda, readMode);
    SetBranch("HasLambdaVDecayed", &fHasLambdaVDecayed, readMode);

    // Set branch addresses for cosmic rejection
    SetBranch("IsClearCosmic", &fIsClearCosmic, readMode);
    SetBranch("NuScore", &fNuScore, readMode);
    SetBranch("bcfmScore", &fBcFMScore, readMode);

    // Set branch addresses for reco vertex information
    SetBranch("RecnuvX", &fRecnuvX, readMode);
    SetBranch("RecnuvY", &fRecnuvY, readMode);
    SetBranch("RecnuvZ", &fRecnuvZ, readMode);
    SetBranch("RecoIsFiducial", &fRecoIsFiducial, readMode);

    // Set branch addresses for FRANS PANDORA information
    SetBranch("FRANSScorePANDORA", &fFRANSScorePANDORA, readMode);
    SetBranch("NVertexTracks", &fNVertexTracks, readMode);
    SetBranch("ShowerEnergy", &fShowerEnergy, readMode);
    SetBranch("NShowers", &fNShowers, readMode);
    SetBranch("NShowerHits", &fNShowerHits, readMode);
    SetBranch("NShwTh75", &fNShwTh75, readMode);
    SetBranch("NShwTh100", &fNShwTh100, readMode);
    SetBranch("MainShowerEnergy", &fMainShowerEnergy, readMode);
    SetBranch("MainShowerScore", &fMainShowerScore, readMode);
    if(!readMode){
        SetBranch("ShowerEnergyVect", &fShowerEnergyVect, readMode);
        SetBranch("ShowerScoreVect", &fShowerScoreVect, readMode);
        SetBranch("ShowerNHitsVect", &fShowerNHitsVect, readMode);
        SetBranch("ShowerOutROINHitsVect", &fShowerOutROINHitsVect, readMode);
    }
    
    // Set branch addresses for # origins information
    SetBranch("NOrigins", &fNOrigins, readMode);
    SetBranch("NOriginsMult1", &fNOriginsMult1, readMode);
    SetBranch("NOriginsMult2", &fNOriginsMult2, readMode);
    SetBranch("NOriginsMultGT3", &fNOriginsMultGT3, readMode);
    SetBranch("NOriginsPairOneTwo", &fNOriginsPairOneTwo, readMode);

    // Set branch addresses for # unassociated origins information
    SetBranch("NUnOrigins", &fNUnOrigins, readMode);
    SetBranch("NUnOriginsMult1", &fNUnOriginsMult1, readMode);
    SetBranch("NUnOriginsMult2", &fNUnOriginsMult2, readMode);
    SetBranch("NUnOriginsMultGT3", &fNUnOriginsMultGT3, readMode);

    // Set branch addresses for angle information
    SetBranch("NAngles", &fNAngles, readMode);
    SetBranch("AngleFRANSScore", &fAngleFRANSScore, readMode);
    SetBranch("AngleGap", &fAngleGap, readMode);
    SetBranch("AngleNHits", &fAngleNHits, readMode);
    SetBranch("AngleNHitsTrack1", &fAngleNHitsTrack1, readMode);
    SetBranch("AngleNHitsTrack2", &fAngleNHitsTrack2, readMode);
    SetBranch("AngleMinNHits", &fAngleMinNHits, readMode);
    SetBranch("AngleNHitsMainTrack", &fAngleNHitsMainTrack, readMode);
    SetBranch("AngleLengthTrack1", &fAngleLengthTrack1, readMode);
    SetBranch("AngleLengthTrack2", &fAngleLengthTrack2, readMode);
    SetBranch("AngleLengthMainTrack", &fAngleLengthMainTrack, readMode);
    SetBranch("AngleLongestIsMain", &fAngleLongestIsMain, readMode);
    
    // Angle cleaness
    SetBranch("AngleCoveredArea", &fAngleCoveredArea, readMode);
    SetBranch("AngleDirtHits", &fAngleDirtHits, readMode);
    SetBranch("AngleDirtHitsRatio", &fAngleDirtHitsRatio, readMode);
    SetBranch("AngleDirtHitsWires", &fAngleDirtHitsWires, readMode);
    SetBranch("AngleDirtHitsWiresRatio", &fAngleDirtHitsWiresRatio, readMode);
    SetBranch("NFreeHits", &fNFreeHits, readMode);
    SetBranch("NUnassociatedHits", &fNUnassociatedHits, readMode);
    SetBranch("NMaxDirtUnassociatedHits", &fNMaxDirtUnassociatedHits, readMode);

    // Kinematics
    SetBranch("AngleDecayContainedDiff", &fAngleDecayContainedDiff, readMode);
    SetBranch("AngleOpeningAngle", &fAngleOpeningAngle, readMode);
    SetBranch("AngleMainTrackOverlap", &fAngleMainTrackOverlap, readMode);
    SetBranch("AnglePzSignMuon", &fAnglePzSignMuon, readMode);
    SetBranch("AnglePzSignLambda", &fAnglePzSignLambda, readMode);
    SetBranch("AnglePzSign", &fAnglePzSign, readMode);
    SetBranch("AngleGapOverlapWithAPAJuntion", &fAngleGapOverlapWithAPAJuntion, readMode);

    // Set branch addresses for charge ratio information
    // vertex
    SetBranch("AnglePassFit", &fAnglePassFit, readMode);
    SetBranch("AngleTwoLinesChi2", &fAngleTwoLinesChi2, readMode);
    SetBranch("AngleNVertexHits", &fAngleNVertexHits, readMode);
    SetBranch("AngleNBulkHits", &fAngleNBulkHits, readMode);
    SetBranch("AngleVertexHitIntegralRatio", &fAngleVertexHitIntegralRatio, readMode);
    SetBranch("AngleVertexHitIntegralDifference", &fAngleVertexHitIntegralDifference, readMode);
    SetBranch("AngleVertexHitIntegralRelativeDifference", &fAngleVertexHitIntegralRelativeDifference, readMode);
    SetBranch("AngleTrackLength1", &fAngleTrackLength1, readMode);
    SetBranch("AngleTrackLength2", &fAngleTrackLength2, readMode);
    SetBranch("AngleTrackLengthRatio", &fAngleTrackLengthRatio, readMode);
    // residual range
    SetBranch("AngleResidualRange1RMS", &fAngleResidualRange1RMS, readMode);
    SetBranch("AngleResidualRange2RMS", &fAngleResidualRange2RMS, readMode);
    SetBranch("AngleResidualRangeMinRMS", &fAngleResidualRangeMinRMS, readMode);
    SetBranch("AngleResidualRangeMaxRMS", &fAngleResidualRangeMaxRMS, readMode);
    SetBranch("AngleResidualRange1AngleRMS", &fAngleResidualRange1AngleRMS, readMode);
    SetBranch("AngleResidualRange2AngleRMS", &fAngleResidualRange2AngleRMS, readMode);
    SetBranch("AngleResidualRangeMinAngleRMS", &fAngleResidualRangeMinAngleRMS, readMode);
    SetBranch("AngleResidualRangeMaxAngleRMS", &fAngleResidualRangeMaxAngleRMS, readMode);
    // dQdx
    SetBranch("AngleChargeRatioAverage", &fAngleChargeRatioAverage, readMode);
    SetBranch("AngleChargeDifferenceAverage", &fAngleChargeDifferenceAverage, readMode);
    SetBranch("AngleChargeRelativeDifferenceAverage", &fAngleChargeRelativeDifferenceAverage, readMode);
    SetBranch("AngleChargeRatioIntegral", &fAngleChargeRatioIntegral, readMode);
    SetBranch("AngleChargeDifferenceIntegral", &fAngleChargeDifferenceIntegral, readMode);
    SetBranch("AngleChargeRelativeDifferenceIntegral", &fAngleChargeRelativeDifferenceIntegral, readMode);
    SetBranch("AnglePassChargeFit", &fAnglePassChargeFit, readMode);
    SetBranch("AngleBandOverlap", &fAngleBandOverlap, readMode);
    SetBranch("AngleBandCrossHits", &fAngleBandCrossHits, readMode);
    SetBranch("AngleChargeRatioFit", &fAngleChargeRatioFit, readMode);
    SetBranch("AngleChargeDifferenceFit", &fAngleChargeDifferenceFit, readMode);
    SetBranch("AngleChargeRelativeDifferenceFit", &fAngleChargeRelativeDifferenceFit, readMode);

    // Set branch addresses for calorimetry information
    SetBranch("NTracks", &fNTracks, readMode);
    SetBranch("NTracksLI", &fNTracksLI, readMode);
    SetBranch("NTracksHI", &fNTracksHI, readMode);
    SetBranch("KEHI", &fKEHI, readMode);
    SetBranch("KELI", &fKELI, readMode);
    SetBranch("OpeningAngle", &fOpeningAngle, readMode);
    SetBranch("InvariantMass", &fInvariantMass, readMode);
    SetBranch("KEDecayedMother", &fKEDecayedMother, readMode);
    SetBranch("EDecayedMother", &fEDecayedMother, readMode);
    SetBranch("RecoGap", &fRecoGap, readMode);
    SetBranch("RecoBeta", &fRecoBeta, readMode);
    SetBranch("RecoGamma", &fRecoGamma, readMode);
    SetBranch("RecoTauLab", &fRecoTauLab, readMode);
    SetBranch("RecoTau", &fRecoTau, readMode);
    SetBranch("RecoDecayTransverseMomentum", &fRecoDecayTransverseMomentum, readMode);
    SetBranch("RecoDecayLongitudinalAsymmetry", &fRecoDecayLongitudinalAsymmetry, readMode);
    SetBranch("RecoDecayTransverseMomentumFromGap", &fRecoDecayTransverseMomentumFromGap, readMode);
    SetBranch("RecoDecayLongitudinalAsymmetryFromGap", &fRecoDecayLongitudinalAsymmetryFromGap, readMode);
    SetBranch("RecoDecayGapAngleDifference", &fRecoDecayGapAngleDifference, readMode);

    // Final reco vertices
    SetBranch("RecoLeptonX", &fRecoLeptonX, readMode);
    SetBranch("RecoLeptonY", &fRecoLeptonY, readMode);
    SetBranch("RecoLeptonZ", &fRecoLeptonZ, readMode);
    SetBranch("RecoHyperonX", &fRecoHyperonX, readMode);
    SetBranch("RecoHyperonY", &fRecoHyperonY, readMode);
    SetBranch("RecoHyperonZ", &fRecoHyperonZ, readMode);

    if(readMode==false)
        SetBranch("InputFileName", &fInputFileName, readMode);
    else{
        SetBranch("InputFileName", &fInputFileNameRead, readMode);
    }

}

void LambdaAnaTree::ResetVars(){
    fEventID = -999;
    fSubrunID = -999;
    fRunID = -999;
    fSliceID = -999;
    fInputFileName = "";
    fmctID = -999;

    fIsClearCosmic = false;
    fNuScore = -999;
    fBcFMScore = -999;

    fIntOrigin = -999;
    fIntDirt = false;
    fIntMode = -999;
    fIntType = -999;
    fIntCCNC = -999;
    fIntNuPDG = -999;
    fIntNuEnergy = -999;
    fIntNProtons = -999;
    fIntNNeutrons = -999;
    fIntNPi0 = -999;
    fIntNPip = -999;
    fIntNPim = -999;
    fIntNMuonP = -999;
    fIntNMuonM = -999;
    fIntNElectronP = -999;
    fIntNElectronM = -999;
    fIntNLambda = -999;
    fIntNSigma0 = -999;
    fIntNSigmaPlus = -999;
    fIntNSigmaMinus = -999;
    fIntNHyperons = -999;

    fNuvE = -999;
    fNuvT = -999;
    fNuvX = -999;
    fNuvY = -999;
    fNuvZ = -999;
    fTruthIsFiducial = false;
    fTruthIsAV = false;

    fSliceCompleteness = -999;
    fSlicePurity = -999;

    fMainMCParticlePDG = -999;
    fMainMCParticleEnergy = -999;
    fMainMCParticleStartX = -999;
    fMainMCParticleStartY = -999;
    fMainMCParticleStartZ = -999;
    fMainMCParticleStartT = -999;
    fMainMCParticleEndX = -999;
    fMainMCParticleEndY = -999;
    fMainMCParticleEndZ = -999;
    fMainMCParticleEndT = -999;
    fMainMCParticleLength = -999;
    fMainMCParticleNDaughters = -999;
    fMainMCParticleMichelDecay = false;

    fGap = -999;
    fLambdaKE = -999;
    fProtonKE = -999;
    fPionKE = -999;
    fLambdaBeta = -999;
    fLambdaGamma = -999;
    fLambdaDecayTransverseMomentum = -999;
    fLambdaDecayLongitudinalAsymmetry = -999;
    fLambdaX = -999;
    fLambdaY = -999;
    fLambdaZ = -999;
    fHasLambda = false;
    fHasLambdaVDecayed = false;

    fRecnuvX = -999;
    fRecnuvY = -999;
    fRecnuvZ = -999;
    fRecoIsFiducial = false;

    fFRANSScorePANDORA = -999;
    fNVertexTracks = -999;
    fShowerEnergy = -999;
    fNShwTh75 = -999;
    fNShwTh100 = -999;
    fNShowers = -999;
    fNShowerHits = -999;
    fNShowersOutROI = -999;
    fMainShowerEnergy = -999;
    fMainShowerScore = -999;
    fShowerEnergyVect.clear();
    fShowerScoreVect.clear();
    fShowerNHitsVect.clear();
    fShowerOutROINHitsVect.clear();

    fNOrigins = -999;
    fNOriginsMult1 = -999;
    fNOriginsMult2 = -999;
    fNOriginsMultGT3 = -999;
    fNOriginsPairOneTwo = -999;

    fNUnOrigins = -999;
    fNUnOriginsMult1 = -999;
    fNUnOriginsMult2 = -999;
    fNUnOriginsMultGT3 = -999;

    // Angle information
    fNAngles = -999;

    fAngleFRANSScore = -999;
    fAngleFRANSDelta = -999;
    fAngleFRANSEta = -999;
    fAngleFRANSFitScore = -999;
    fAngleFRANSIota = -999;

    fAngleGap = -999;
    fAngleNHits = -999;
    fAngleNHitsTrack1 = -999;
    fAngleNHitsTrack2 = -999;
    fAngleMinNHits = -999;
    fAngleNHitsMainTrack = -999;
    fAngleLengthTrack1 = -999;
    fAngleLengthTrack2 = -999;
    fAngleLengthMainTrack = -999;
    fAngleLongestIsMain = false;
    
    // Angle cleaness
    fAngleCoveredArea = -999;
    fAngleDirtHits = -999;
    fAngleDirtHitsRatio = -999;
    fAngleDirtHitsWires = -999;
    fAngleDirtHitsWiresRatio = -999;
    fNFreeHits = -999;
    fNUnassociatedHits = -999;
    fNMaxDirtUnassociatedHits = -999;

    // Kinematics
    fAngleOpeningAngle = -999;
    fAngleDecayContainedDiff = -999;
    fAngleMainTrackOverlap = -999;
    fAnglePzSignMuon = -999;
    fAnglePzSignLambda = -999;
    fAnglePzSign = -999;
    fAngleGapOverlapWithAPAJuntion = -999;

    // Calorimetry
    fAnglePassFit = false;
    fAngleTwoLinesChi2 = -999;
    fAngleNVertexHits = -999;
    fAngleNBulkHits = -999;
    fAngleVertexHitIntegralRatio = -999;
    fAngleVertexHitIntegralDifference = -999;
    fAngleVertexHitIntegralRelativeDifference = -999;
    fAngleTrackLength1 = -999;
    fAngleTrackLength2 = -999;
    fAngleTrackLengthRatio = -999;
    fAngleResidualRange1RMS = -999;
    fAngleResidualRange2RMS = -999;
    fAngleResidualRangeMinRMS = -999;
    fAngleResidualRangeMaxRMS = -999;
    fAngleResidualRange1AngleRMS = -999;
    fAngleResidualRange2AngleRMS = -999;
    fAngleResidualRangeMinAngleRMS = -999;
    fAngleResidualRangeMaxAngleRMS = -999;
    fAngleChargeRatioAverage = -999;
    fAngleChargeDifferenceAverage = -999;
    fAngleChargeRelativeDifferenceAverage = -999;
    fAngleChargeRatioIntegral = -999;
    fAngleChargeDifferenceIntegral = -999;
    fAngleChargeRelativeDifferenceIntegral = -999;
    fAnglePassChargeFit = false;
    fAngleBandOverlap = -999;
    fAngleBandCrossHits = -999;
    fAngleChargeRatioFit = -999;
    fAngleChargeDifferenceFit = -999;
    fAngleChargeRelativeDifferenceFit = -999;

    // Calorimetry
    fNTracks = -999;
    fNTracksLI = -999;
    fNTracksHI = -999;
    fKEHI = -999;
    fKELI = -999;
    fOpeningAngle = -999;
    fInvariantMass = -999;
    fKEDecayedMother = -999;
    fEDecayedMother = -999;
    fRecoGap = -999;
    fRecoBeta = -999;
    fRecoGamma = -999;
    fRecoTauLab = -999;
    fRecoTau = -999;
    fRecoDecayTransverseMomentum = -999;
    fRecoDecayLongitudinalAsymmetry = -999;
    fRecoDecayTransverseMomentumFromGap = -999;
    fRecoDecayLongitudinalAsymmetryFromGap = -999;
    fRecoDecayGapAngleDifference = -999;

    // Final reco vertices
    fRecoLeptonX = -999;
    fRecoLeptonY = -999;
    fRecoLeptonZ = -999;
    fRecoHyperonX = -999;
    fRecoHyperonY = -999;
    fRecoHyperonZ = -999;

}


void LambdaAnaTree::GetEntry(int i){
    fTree->GetEntry(i);
}


// --- Print function ---
void LambdaAnaTree::PrintEventInfo() {
        // Event information
        std::cout << " ---- Lambda Ana Tree Summary ----\n";
        std::cout << "Run: " << fRunID << " Subrun: " << fSubrunID << " Event: " << fEventID << "\n";
        std::cout << " NuE: "<<fNuvE<<" NuT: "<<fNuvT<<" NuX: "<<fNuvX<<" NuY: "<<fNuvY<<" NuZ: "<<fNuvZ<<std::endl;
        std::cout << " Lambda gap: " << fGap << " LambdaKE: "<<fLambdaKE<<std::endl;
        std::cout << " ProtonKE: " << fProtonKE << " PionKE: "<<fPionKE<<std::endl;
        std::cout << " Has Lambda: " << fHasLambda << " Has Lambda V Decayed: "<<fHasLambdaVDecayed<<std::endl;

        // Slice information
        std::cout << " Slice Completeness: " << fSliceCompleteness << " Slice Purity: "<<fSlicePurity<<std::endl;

        // Reco vertex information
        std::cout << " Reco nuvX: " << fRecnuvX << " Reco nuvY: "<<fRecnuvY<<" Reco nuvZ: "<<fRecnuvZ<<std::endl;
        std::cout << " Reco Is Fiducial: " << fRecoIsFiducial << std::endl;

        // FRANS PANDORA information
        std::cout << " FRANS Score PANDORA: " << fFRANSScorePANDORA << std::endl;

        // # origins information
        std::cout << " # Origins: " << fNOrigins << " Mult 1: "<<fNOriginsMult1<<" Mult 2: "<<fNOriginsMult2<<" Mult >3: "<<fNOriginsMultGT3<<" # Origins Pair 1-2: "<<fNOriginsPairOneTwo<<std::endl;

        // # unassociated origins information
        std::cout << " # Unassociated Origins: " << fNUnOrigins << " Mult 1: "<<fNUnOriginsMult1<<" Mult 2: "<<fNUnOriginsMult2<<" Mult >3: "<<fNUnOriginsMultGT3<<std::endl;

        // Angle information
        std::cout << " --- Triangle Information --- " << std::endl;
        std::cout << " # Angles: " << fNAngles << std::endl;
        std::cout << " Angle FRANS Score: " << fAngleFRANSScore << std::endl;
        std::cout << " Angle Gap: " << fAngleGap << std::endl;
        std::cout << " Angle # Hits: " << fAngleNHits << std::endl;
        std::cout << " Angle # Hits Track 1: " << fAngleNHitsTrack1 << std::endl;
        std::cout << " Angle # Hits Track 2: " << fAngleNHitsTrack2 << std::endl;
        std::cout << " Angle Min # Hits: " << fAngleMinNHits << std::endl;
        std::cout << " Angle # Hits Main Track: " << fAngleNHitsMainTrack << std::endl;
        std::cout << " Angle Length Track 1: " << fAngleLengthTrack1 << std::endl;
        std::cout << " Angle Length Track 2: " << fAngleLengthTrack2 << std::endl;
        std::cout << " Angle Length Main Track: " << fAngleLengthMainTrack << std::endl;
        std::cout << " Angle Longest Is Main: " << fAngleLongestIsMain << std::endl;
    
        // Angle cleaness
        std::cout << " Angle Covered Area: " << fAngleCoveredArea << std::endl;
        std::cout << " Angle Dirt Hits: " << fAngleDirtHits << std::endl;
        std::cout << " Angle Dirt Hits Ratio: " << fAngleDirtHitsRatio << std::endl;
        std::cout << " Angle Dirt Hits Wires: " << fAngleDirtHitsWires << std::endl;
        std::cout << " Angle Dirt Hits Wires Ratio: " << fAngleDirtHitsWiresRatio << std::endl;
        std::cout << " Angle Unassociated Hits: " << fNUnassociatedHits << std::endl;
        std::cout << " Angle Free Hits: " << fNFreeHits << std::endl;
        std::cout << " Angle Max Dirt Unassociated Hits: " << fNMaxDirtUnassociatedHits << std::endl;

        // Kinematics
        std::cout << " --- Kinematics Information --- " << std::endl;
        std::cout << " Angle Opening Angle: " << fAngleOpeningAngle << std::endl;
        std::cout << " Angle Decay Contained Diff: " << fAngleDecayContainedDiff << std::endl;
        std::cout << " Angle Main Track Overlap: " << fAngleMainTrackOverlap << std::endl;
        std::cout << " Angle Pz Sign Muon: " << fAnglePzSignMuon << std::endl;
        std::cout << " Angle Pz Sign Lambda: " << fAnglePzSignLambda << std::endl;
        std::cout << " Angle Pz Sign: " << fAnglePzSign << std::endl;
        std::cout << " Angle Gap Overlap With APA Juntion: " << fAngleGapOverlapWithAPAJuntion << std::endl;

        // Calorimetry
        std::cout << " --- Calorimetry Information --- " << std::endl;
        std::cout << "    --- Vertex --- " << std::endl;
        std::cout << " Angle Pass Fit: " << fAnglePassFit << std::endl;
        std::cout << " Angle Two Lines Chi2: " << fAngleTwoLinesChi2 << std::endl;
        std::cout << " Angle # Vertex Hits: " << fAngleNVertexHits << std::endl;
        std::cout << " Angle # Bulk Hits: " << fAngleNBulkHits << std::endl;
        std::cout << " Angle Vertex Hit Integral Ratio: " << fAngleVertexHitIntegralRatio << std::endl;
        std::cout << " Angle Vertex Hit Integral Difference: " << fAngleVertexHitIntegralDifference << std::endl;
        std::cout << " Angle Vertex Hit Integral Relative Difference: " << fAngleVertexHitIntegralRelativeDifference << std::endl;
        std::cout << " Angle Track Length 1: " << fAngleTrackLength1 << std::endl;
        std::cout << " Angle Track Length 2: " << fAngleTrackLength2 << std::endl;
        std::cout << " Angle Track Length Ratio: " << fAngleTrackLengthRatio << std::endl;
        // residual range properties
        std::cout << "    --- Residual range --- " << std::endl;
        std::cout << " Angle Residual Range 1 RMS: " << fAngleResidualRange1RMS << std::endl;
        std::cout << " Angle Residual Range 2 RMS: " << fAngleResidualRange2RMS << std::endl;
        std::cout << " Angle Residual Range Min RMS: " << fAngleResidualRangeMinRMS << std::endl;
        std::cout << " Angle Residual Range Max RMS: " << fAngleResidualRangeMaxRMS << std::endl;
        std::cout << " Angle Residual Range 1 Angle RMS: " << fAngleResidualRange1AngleRMS << std::endl;
        std::cout << " Angle Residual Range 2 Angle RMS: " << fAngleResidualRange2AngleRMS << std::endl;
        std::cout << " Angle Residual Range Min Angle RMS: " << fAngleResidualRangeMinAngleRMS << std::endl;
        std::cout << " Angle Residual Range Max Angle RMS: " << fAngleResidualRangeMaxAngleRMS << std::endl;
        std::cout << "    --- dQdx --- " << std::endl;
        std::cout << " Angle Charge Ratio Average: " << fAngleChargeRatioAverage << std::endl;
        std::cout << " Angle Charge Difference Average: " << fAngleChargeDifferenceAverage << std::endl;
        std::cout << " Angle Charge Relative Difference Average: " << fAngleChargeRelativeDifferenceAverage << std::endl;
        std::cout << " Angle Charge Ratio Integral: " << fAngleChargeRatioIntegral << std::endl;
        std::cout << " Angle Charge Difference Integral: " << fAngleChargeDifferenceIntegral << std::endl;
        std::cout << " Angle Charge Relative Difference Integral: " << fAngleChargeRelativeDifferenceIntegral << std::endl;
        std::cout << " Angle Pass Charge Fit: " << fAnglePassChargeFit << std::endl;
        std::cout << " Angle Band Overlap: " << fAngleBandOverlap << std::endl;
        std::cout << " Angle Band Cross Hits: " << fAngleBandCrossHits << std::endl;
        std::cout << " Angle Charge Ratio Fit: " << fAngleChargeRatioFit << std::endl;
        std::cout << " Angle Charge Difference Fit: " << fAngleChargeDifferenceFit << std::endl;
        std::cout << " Angle Charge Relative Difference Fit: " << fAngleChargeRelativeDifferenceFit << std::endl;
        // PID
        std::cout << " --- PID Information --- " << std::endl;
        std::cout << " NTracks LI: " << fNTracksLI << std::endl;
        std::cout << " NTracks HI: " << fNTracksHI << std::endl;
        std::cout << " KE HI: " << fKEHI << std::endl;
        std::cout << " KE LI: " << fKELI << std::endl;
        std::cout << " Opening Angle: " << fOpeningAngle << std::endl;
        std::cout << " Invariant Mass: " << fInvariantMass << std::endl;
        std::cout << " KE Decayed Mother: " << fKEDecayedMother << std::endl;
        std::cout << " E Decayed Mother: " << fEDecayedMother << std::endl;
        std::cout << " Reco Gap: " << fRecoGap << std::endl;
        std::cout << " Reco Beta: " << fRecoBeta << std::endl;
        std::cout << " Reco Gamma: " << fRecoGamma << std::endl;
        std::cout << " Reco TauLab: " << fRecoTauLab << " Reco Tau: " << fRecoTau << std::endl;


    }
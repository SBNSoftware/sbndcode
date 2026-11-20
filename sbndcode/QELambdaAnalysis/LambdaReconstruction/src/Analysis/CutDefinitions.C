#include "CutEfficienciesDefinitions.C"

//--------- Cut definitions
double fCutCRUMBS = 0.;
double fCutMinNAngles = 1;
double fCutFRANS = 0.2;
double fCutNOrigins = 4;
double fCutNOriginsM3 = 0;
double fCutFRANSPANDORA = 0.2;
double fCutNShw = 1;
double fCutShwEnergy = 135;

//--------- Phase space definitions
std::vector<PlotDef> psLambdaKinematics = {
    {"NuvZ",  "0==0", CutType::kNone,  0, {-10,510,25}, false, "Z_{#nu} [cm]"}
    ,{"LambdaKE", "0==0", CutType::kNone,  0, {0, 0.4, 20}, false, "KE_{#Lambda} [GeV]"}
    ,{"NuvE", "0==0", CutType::kNone,  0, {0, 3, 10}, false, "E_{#nu} [GeV]"}
    ,{"Gap", "0==0", CutType::kNone, 0, {0, 20, 20}, false, "Gap [cm]"}
};

//---------Volume cuts
std::string fTruthInFV = "TruthIsFiducial==1 &&";
//std::string fTruthInAV = "abs(NuvX)<200 && abs(NuvY)<200 && NuvZ>0 && NuvZ<500 && ";
std::string fTruthInAV = "TruthIsAV && ";

//---------Counter cut
TCut fCounterCut = "SliceID==0";

std::vector<PlotDef> cutDefsTalk2 = {
    
    {"TruthIsFiducial || (IntOrigin==2 || IntDirt==1) ",  "TruthIsFiducial || (IntOrigin==2 || IntDirt==1)", CutType::kCenter, 1, {0,2,2}, true, "Truth in FV",  "Truth \\ in \\ FV"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}

    // Pre FRANS (all events)
    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kNone, fCutFRANSPANDORA, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score", true}
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score", true}
    
    // Secondary vertices
    ,{"NAngles", "NAngles", CutType::kRight, 1, {0, 10, 10}, true,  "# seconday vertices",  "\\ Secondary \\ vertex"}

    // FRANS
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, true, "FRANS score", "FRANS \\ score"}
    //,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kNone, fCutFRANSPANDORA, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score", true}

    // CLEANESS
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft,  1.5, {0, 20, 40}, true, "#Delta OpeningAngle [#circ]", "DeltaOpeningAngle"}
    ,{"NUnOrigins", "NUnOrigins", CutType::kLeftInt, 0, {0, 15, 15}, true,  "# extra vertices ",  "\\#  \\ extra \\ vertices"}    
    ,{"AnglePassChargeFit==1 && AngleResidualRangeMaxAngleRMS<=50", "(AnglePassChargeFit==1 && AngleResidualRangeMaxAngleRMS<=50)", CutType::kCenter, 1, {0, 2, 2}, true, "PreCalorimetryCut", "PreCalorimetryCut"}
    ,{"AnglePzSign", "AnglePzSign", CutType::kRight,  1, {-4, 4, 8}, true, "AnglePzSign", "AnglePzSign"}  
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 300, 60}, true, "Main track length [cm]", "Main \\ track \\ length \\ [cm]"}
    ,{"AngleDirtHits<=10 && NUnassociatedHits<=10", "(AngleDirtHits<=10 && NUnassociatedHits<=10)", CutType::kCenter, 1, {0, 2, 2}, true, "# unnassociated hits<10", "unnassociated hits<10"}
    //,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 10, {0, 400, 40}, true, "# unassociated hits", "\\# \\ unassociated \\ hits"}
   
    // After Cuts Distributions
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score"}
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 300, 60}, false, "Main track length [cm]", "Main \\ track \\ length \\ [cm]"}
    ,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 10, {0, 400, 40}, false, "# unassociated hits", "\\# \\ unassociated \\ hits"}
    ,{"AngleDirtHits", "AngleDirtHits", CutType::kLeft,  25, {0, 50, 25}, false, "#Dirt Hits", "Dirt \\ Hits"}
    ,{"AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS", CutType::kLeft, 50, {0, 360, 36}, false, "AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS"}
    ,{"AnglePassChargeFit", "AnglePassChargeFit", CutType::kRight, 1, {0, 2, 2}, false, "AnglePassChargeFit", "AnglePassChargeFit"}
    
}; 


std::vector<PlotDef> cutDefsTalk2Induction = {
    
    {"TruthIsFiducial || (IntOrigin==2 || IntDirt==1) ",  "TruthIsFiducial || (IntOrigin==2 || IntDirt==1)", CutType::kCenter, 1, {0,2,2}, true, "Truth in FV",  "Truth \\ in \\ FV"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}

    // FRANS (all events)
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score", true}
    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kNone, 0.1, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score", true}
    
    // Secondary vertices
    ,{"NAngles", "NAngles", CutType::kRight, 1, {0, 10, 10}, true,  "# seconday vertices",  "\\ Secondary \\ vertex"}

    // FRANS
    //,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, true, "FRANS score", "FRANS \\ score"}
    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, 0.1, {-.5,.5,40}, true, "FRANS score", "FRANS \\ score" }

    // CLEANESS
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft,  1.5, {0, 20, 40}, true, "#Delta OpeningAngle [#circ]", "DeltaOpeningAngle"}
    ,{"NUnOrigins", "NUnOrigins", CutType::kLeftInt, 0, {0, 15, 15}, true,  "# extra vertices ",  "\\#  \\ extra \\ vertices"}    
    ,{"AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS", CutType::kLeft, 50, {0, 300, 30}, true, "AngleRMS", "AngleRMS"}
    //,{"AnglePzSign", "AnglePzSign", CutType::kRight,  1, {-4, 4, 8}, true, "AnglePzSign", "AnglePzSign"}  
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 300, 60}, true, "Main track length [cm]", "Main \\ track \\ length \\ [cm]"}
    ,{"AngleDirtHits<=10 && NUnassociatedHits<=10", "(AngleDirtHits<=10 && NUnassociatedHits<=10)", CutType::kCenter, 1, {0, 2, 2}, true, "# unnassociated hits<10", "unnassociated hits<10"}
    //,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 10, {0, 400, 40}, true, "# unassociated hits", "\\# \\ unassociated \\ hits"}
   
    // After Cuts Distributions
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score"}
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 300, 60}, false, "Main track length [cm]", "Main \\ track \\ length \\ [cm]"}
    ,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 10, {0, 400, 40}, false, "# unassociated hits", "\\# \\ unassociated \\ hits"}
    ,{"AngleDirtHits", "AngleDirtHits", CutType::kLeft,  25, {0, 50, 25}, false, "#Dirt Hits", "Dirt \\ Hits"}
    ,{"AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS", CutType::kLeft, 50, {0, 360, 36}, false, "AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS"}
    ,{"AnglePassChargeFit", "AnglePassChargeFit", CutType::kRight, 1, {0, 2, 2}, false, "AnglePassChargeFit", "AnglePassChargeFit"}
    
}; 


std::vector<PlotDef> cutDefsPIDFull = {
    
    {"TruthIsFiducial || (IntOrigin==2 || IntDirt==1) ",  "TruthIsFiducial || (IntOrigin==2 || IntDirt==1)", CutType::kCenter, 1, {0,2,2}, true, "Truth in FV",  "Truth \\ in \\ FV"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}    

    // Secondary vertices
    ,{"NAngles", "NAngles", CutType::kRight, 1, {0, 10, 10}, true,  "# \\upsilonP{2}-\\upsilon^{1} pairs",  "\\upsilon^{2}-\\upsilon^{1} \\ pairs"}

    // FRANS
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, true, "FRANS score", "FRANS \\ score"}

    //PID
    ,{"NTracksLI>=1 && NTracksHI>=1", "(NTracksLI>=1 && NTracksHI>=1)", CutType::kCenter, 1, {0, 2, 2}, true, "NTracksLI>0 && NTracksHI>0", "PID"}
    
    // V characterization
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft,  2.5, {0, 20, 40}, true, "#Phi [#circ]", "Phi \\ [\\circ]"}
    ,{"AnglePzSign", "AnglePzSign", CutType::kRight,  1, {-4, 4, 8}, true, "Sign(Pz)", "Sign(Pz)"}
    ,{"NUnOrigins", "NUnOrigins", CutType::kLeftInt, 0, {0, 15, 15}, true,  "# extra vertices ",  "\\#  \\ extra \\ vertices"}    
    ,{"AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS", CutType::kLeft, 50, {0, 360, 36}, true, "VTracksRMS [#circ]", "VTracksRMS \\ [\\circ]"}
    ,{"NMaxDirtUnassociatedHits", "NMaxDirtUnassociatedHits", CutType::kLeft, 10, {0, 400, 40}, true, "# unassociated hits", "\\# \\ unassociated \\ hits"}
   
    //,{"NMaxDirtUnassociatedHits", "(NMaxDirtUnassociatedHits<10 && AngleLengthMainTrack>=10 && AngleGapOverlapWithAPAJuntion<=0.1)", CutType::kCenter, 1, {0, 400, 40}, true, "# unassociated hits", "\\# \\ unassociated \\ hits"}
   
    // Additional fiducial volume cuts
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 300, 60}, true, "Main track length [cm]", "Main \\ track \\ length \\ [cm]"}
    ,{"AngleGapOverlapWithAPAJuntion", "AngleGapOverlapWithAPAJuntion", CutType::kLeft,  0.1, {0, 1, 10}, true, "verlapWithAPAJuntion", "OverlapWithAPAJuntion"}

    ,{"InvariantMass", "InvariantMass", CutType::kCenter, 1115, {1000, 1500, 10}, false, "Invariant Mass [MeV]", "Invariant \\ Mass"}

    ,{"CRUMBSScore",  "CRUMBSScore",         CutType::kNone,   fCutCRUMBS, {-1,1,20}, true, "CRUMBSScore",  "CRUMBSScore"}
};


std::vector<PlotDef> cutDefsInd = {
    
    {"TruthIsFiducial || (IntOrigin==2 || IntDirt==1) ",  "TruthIsFiducial || (IntOrigin==2 || IntDirt==1)", CutType::kCenter, 1, {0,2,2}, true, "Truth in FV",  "Truth \\ in \\ FV"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}    

    // Secondary vertices
    ,{"NAngles", "NAngles", CutType::kRight, 1, {0, 10, 10}, true,  "# \\upsilon^{2}-\\upsilon^{1} pairs",  "\\upsilon^{2}-\\upsilon^{1} \\ pairs"}

    // FRANS
    //,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, 0.1, {-.5,.5,40}, true, "FRANS score", "FRANS \\ score"}
    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, 0.1, {-.5,.5,40}, true, "FRANS score", "FRANS \\ score" }

    // V characterization
    ,{"RecoDecayGapAngleDifference", "RecoDecayGapAngleDifference", CutType::kLeft, 0.5, {-1,4,20}, true, "#phi", "Phi"}
    ,{"RecoPzSign", "RecoPzSign", CutType::kRight,  1, {-4, 4, 8}, true, "Sign(Pz)", "Sign(Pz)"}
    ,{"NUnOrigins", "NUnOrigins", CutType::kLeftInt, 0, {0, 15, 15}, true,  "# extra vertices ",  "\\#  \\ extra \\ vertices"}    
    ,{"AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS", CutType::kLeft, 50, {0, 360, 36}, true, "AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS"}
    ,{"NMaxDirtUnassociatedHits", "NMaxDirtUnassociatedHits", CutType::kLeft, 10, {0, 400, 40}, true, "# unassociated hits", "\\# \\ unassociated \\ hits"}
   
    ,{"NTracksLI>=1 && NTracksHI>=1", "(NTracksLI>=1 && NTracksHI>=1)", CutType::kCenter, 1, {0, 2, 2}, true, "PID", "PID"}
    ,{"RecoLeptonPID==13", "(RecoLeptonPID==13)", CutType::kRight,  1, {0, 2, 2}, true, "MuonPID", "MuonPID"}

    // Additional fiducial volume cuts
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 300, 60}, true, "Main track length [cm]", "Main \\ track \\ length \\ [cm]"}
    //,{"AngleGapOverlapWithAPAJuntion", "AngleGapOverlapWithAPAJuntion", CutType::kLeft,  0.1, {0, 1, 10}, true, "verlapWithAPAJuntion", "OverlapWithAPAJuntion"}
};


std::vector<PlotDef> cutDefsCol = {
    
    {"TruthIsFiducial || (IntOrigin==2 || IntDirt==1) ",  "TruthIsFiducial || (IntOrigin==2 || IntDirt==1)", CutType::kCenter, 1, {0,2,2}, true, "Truth in FV",  "Truth \\ in \\ FV"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}    

    // Secondary vertices
    ,{"NAngles", "NAngles", CutType::kRight, 1, {0, 10, 10}, true,  "# \\upsilon^{2}-\\upsilon^{1} pairs",  "\\upsilon^{2}-\\upsilon^{1} \\ pairs"}

    ,{"RecoDecayGapAngleDifference", "RecoDecayGapAngleDifference", CutType::kLeft, 0.5, {-1,4,20}, false, "#phi", "Phi"}
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft, 0.5, {-10,100,20}, false, "#alpha", "Alpha"}

    ,{"RecoPzSign", "RecoPzSign", CutType::kRight,  1, {-4, 4, 8}, false, "SignReco(Pz)", "SignReco(Pz)"}
    ,{"AnglePzSign", "AnglePzSign", CutType::kRight,  1, {-4, 4, 8}, false, "Sign(Pz)", "Sign(Pz)"}

    // FRANS
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, true, "FRANS score", "FRANS \\ score"}
    //,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, 0.2, {-.5,.5,40}, true, "FRANS score", "FRANS \\ score" }

    ,{"NTracksLI>=1 && NTracksHI>=1", "(NTracksLI>=1 && NTracksHI>=1)", CutType::kCenter, 1, {0, 2, 2}, true, "PID", "PID"}
    ,{"RecoLeptonPID==13", "(RecoLeptonPID==13)", CutType::kRight,  1, {0, 2, 2}, true, "MuonPID", "MuonPID"}

    // V characterization
    ,{"RecoDecayGapAngleDifference", "RecoDecayGapAngleDifference", CutType::kLeft, 0.75, {-1,4,20}, true, "#phi", "Phi"}
    ,{"RecoPzSign", "RecoPzSign", CutType::kRight,  1, {-4, 4, 8}, true, "Sign(Pz)", "Sign(Pz)"}
    ,{"NUnOrigins", "NUnOrigins", CutType::kLeftInt, 0, {0, 15, 15}, true,  "# extra vertices ",  "\\#  \\ extra \\ vertices"}    
    ,{"AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS", CutType::kLeft, 50, {0, 360, 36}, true, "AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS"}
    ,{"NMaxDirtUnassociatedHits", "NMaxDirtUnassociatedHits", CutType::kLeft, 10, {0, 400, 40}, true, "# unassociated hits", "\\# \\ unassociated \\ hits"}
       

    // Additional fiducial volume cuts
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 300, 60}, true, "Main track length [cm]", "Main \\ track \\ length \\ [cm]"}
    //,{"AngleGapOverlapWithAPAJuntion", "AngleGapOverlapWithAPAJuntion", CutType::kLeft,  0.1, {0, 1, 10}, true, "verlapWithAPAJuntion", "OverlapWithAPAJuntion"}
};


std::vector<PlotDef> cutDefsVcharacterisation = {
    
    {"TruthIsFiducial || (IntOrigin==2 || IntDirt==1) ",  "TruthIsFiducial || (IntOrigin==2 || IntDirt==1)", CutType::kCenter, 1, {0,2,2}, true, "Truth in FV",  "Truth \\ in \\ FV"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}    

    // Secondary vertices
    ,{"NAngles", "NAngles", CutType::kRight, 1, {0, 10, 10}, true,  "# \\upsilon^{2}-\\upsilon^{1} pairs",  "\\upsilon^{2}-\\upsilon^{1} \\ pairs"}
    //PID
    ,{"NTracksLI>=1 && NTracksHI>=1", "(NTracksLI>=1 && NTracksHI>=1)", CutType::kCenter, 1, {0, 2, 2}, true, "NTracksLI>0 && NTracksHI>0", "PID"}
    // FRANS
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, true, "FRANS score", "FRANS \\ score", true}

    // V characterization
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft,  2.5, {0, 20, 40}, false, "#Phi [#circ]", "Phi \\ [\\circ]", true}
    ,{"RecoDecayGapAngleDifference", "RecoDecayGapAngleDifference", CutType::kLeft, 0.5, {-1,4,30}, true, "#phi_{3D}", "Phi_{3D}"}
    ,{"AnglePzSign", "AnglePzSign", CutType::kRight,  1, {-4, 4, 8}, false, "Sign(Pz)", "Sign(Pz)", true}
    ,{"NUnOrigins", "NUnOrigins", CutType::kLeftInt, 0, {0, 15, 15}, false,  "# extra vertices ",  "\\#  \\ extra \\ vertices", true}
    ,{"AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS", CutType::kLeft, 50, {0, 360, 36}, false, "VTracksRMS [#circ]", "VTracksRMS \\ [\\circ]", true}
    ,{"NMaxDirtUnassociatedHits", "NMaxDirtUnassociatedHits", CutType::kLeft, 10, {0, 400, 40}, false, "# unassociated hits", "\\# \\ unassociated \\ hits", true}
 
};


std::vector<PlotDef> cutDefsFinalSelection = {
    
    //{"TruthIsFiducial || (IntOrigin==2 || IntDirt==1) ",  "TruthIsFiducial || (IntOrigin==2 || IntDirt==1)", CutType::kCenter, 1, {0,2,2}, true, "Truth in FV",  "Truth \\ in \\ FV"}
   // {"IntOrigin==2 || IntDirt==1 || (IntNLambda>0 && IntMode==0 && abs(IntNuPDG)==14 && TruthIsFV) || (TruthIsAV && !(IntNLambda>0 && IntMode==0 && abs(IntNuPDG)==14))",  "IntOrigin==2 || IntDirt==1 || (IntNLambda>0 && IntMode==0 && abs(IntNuPDG)==14 && TruthIsFV) || (TruthIsAV && !(IntNLambda>0 && IntMode==0 && abs(IntNuPDG)==14))", CutType::kCenter, 1, {0,2,2}, true, "Truth in FV",  "Truth \\ in \\ FV"}

   {"IntOrigin==2 || IntDirt==1 || TruthIsAV",  "IntOrigin==2 || IntDirt==1 || TruthIsAV", CutType::kCenter, 1, {0,2,2}, true, "Truth in FV",  "Truth \\ in \\ FV"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}
    
    // Secondary vertices
    ,{"NAngles", "NAngles", CutType::kRight, 1, {0, 10, 10}, true,  "# \\upsilon^{2}-\\upsilon^{1} pairs",  "\\upsilon^{2}-\\upsilon^{1} \\ pairs"}

    // FRANS
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, true, "FRANS score", "FRANS \\ score"}

    //PID
    ,{"NTracksLI>=1 && NTracksHI>=1", "(NTracksLI>=1 && NTracksHI>=1)", CutType::kCenter, 1, {0, 2, 2}, true, "PID", "PID"}

    // Joint topological cut
    ,{"AnglePzSign>=1 && RecoDecayGapAngleDifference<=1. && NUnOrigins<=0 && AngleResidualRangeMaxAngleRMS<=50 && NMaxDirtUnassociatedHits<=10 && AngleLengthMainTrack>=10 && AngleGapOverlapWithAPAJuntion<=0.1 && abs(RecoLeptonX)<190 && abs(RecoLeptonY)<190 && RecoLeptonZ>10 && RecoLeptonZ<490", 
    "(AnglePzSign>=1 && RecoDecayGapAngleDifference<=1. && NUnOrigins<=0 && AngleResidualRangeMaxAngleRMS<=50 && NMaxDirtUnassociatedHits<=10 && AngleLengthMainTrack>=10 && AngleGapOverlapWithAPAJuntion<=0.1 && abs(RecoLeptonX)<190 && abs(RecoLeptonY)<190 && RecoLeptonZ>10 && RecoLeptonZ<490)", CutType::kCenter, 1, {0, 2, 2}, true, "V characterisation", "V characterisation"}

    //,{"CRUMBSScore",  "CRUMBSScore",         CutType::kRight,   fCutCRUMBS, {-1,1,20}, true, "CRUMBSScore",  "CRUMBSScore"}    
    ,{"InvariantMass", "InvariantMass", CutType::kCenter, 1115, {1000, 1500, 10}, false, "Invariant Mass [MeV]", "Invariant \\ Mass"}
}; 


std::vector<PlotDef> cutDefsPID = {
    
    {"TruthIsFiducial || (IntOrigin==2 || IntDirt==1) ",  "TruthIsFiducial || (IntOrigin==2 || IntDirt==1)", CutType::kCenter, 1, {0,2,2}, true, "Truth in FV",  "Truth \\ in \\ FV"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}
    
    //,{"CRUMBSScore",  "CRUMBSScore",         CutType::kNone,   fCutCRUMBS, {-1,1,20}, true, "CRUMBSScore",  "CRUMBSScore"}

    // Pre FRANS (all events)
    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kNone, fCutFRANSPANDORA, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score", true}
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score", true}
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score"}
    
    // Secondary vertices
    ,{"NAngles", "NAngles", CutType::kRight, 1, {0, 10, 10}, true,  "# \\upsilon^{2}-\\upsilon^{1} pairs",  "\\upsilon^{2}-\\upsilon^{1} \\ pairs"}

    // FRANS
    //,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, fCutFRANSPANDORA, {-.5,.5,40}, true, "FRANS score PANDORA", "FRANS \\ score \\ PANDORA"}
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, true, "FRANS score", "FRANS \\ score"}

    //PID
    ,{"NTracksLI>=1 && NTracksHI>=1", "(NTracksLI>=1 && NTracksHI>=1)", CutType::kCenter, 1, {0, 2, 2}, true, "PID", "PID"}

    // Joint topological cut
    ,{"AnglePzSign>=1 && AngleDecayContainedDiff<=2.5 && NUnOrigins<=0 && AngleResidualRangeMaxAngleRMS<=50 && NMaxDirtUnassociatedHits<=10 && AngleLengthMainTrack>=10 && AngleGapOverlapWithAPAJuntion<=0.1 && abs(RecoLeptonX)<190 && abs(RecoLeptonY)<190 && RecoLeptonZ>10 && RecoLeptonZ<490", 
    "(AnglePzSign>=1 && AngleDecayContainedDiff<=2.5 && NUnOrigins<=0 && AngleResidualRangeMaxAngleRMS<=50 && NMaxDirtUnassociatedHits<=10 && AngleLengthMainTrack>=10 && AngleGapOverlapWithAPAJuntion<=0.1 && abs(RecoLeptonX)<190 && abs(RecoLeptonY)<190 && RecoLeptonZ>10 && RecoLeptonZ<490)", CutType::kCenter, 1, {0, 2, 2}, true, "V characterisation", "V characterisation"}

    //,{"CRUMBSScore",  "CRUMBSScore",         CutType::kRight,   fCutCRUMBS, {-1,1,20}, true, "CRUMBSScore",  "CRUMBSScore"}    
    ,{"InvariantMass", "InvariantMass", CutType::kCenter, 1115, {1000, 1500, 10}, false, "Invariant Mass [MeV]", "Invariant \\ Mass"}
}; 


// --------- CUTS TO PLOT ORIGINS DISTRIBUTIONS ------------
std::vector<PlotDef> cutDefsOriginsDistributions = {
    {"TruthIsFiducial",   "0==0",            CutType::kNone,   0, {0, 2, 2}, true, "Truth in FV",  "No \\ cut"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0, 2, 2}, true, "Reco in FV",   "Reco \\ in \\ FV"}
    ,{"NOrigins",         "NOrigins",        CutType::kNone,   0, {0, 6, 6}, false, "# \\upsilon",         "\\#  \\ vertices", true}
    ,{"NOriginsMult1",    "NOriginsMult1",   CutType::kNone,   0, {0, 6, 6}, false, "# \\upsilon^{1}",     "\\#  \\ vertices \\ mult \\ 1", true}
    ,{"NOriginsMult2",    "NOriginsMult2",   CutType::kNone,   0, {0, 6, 6}, false, "# \\upsilon^{2}",     "\\#  \\ vertices \\ mult \\ 2", true}  
    ,{"NOriginsMultGT3",  "NOriginsMultGT3", CutType::kNone,   0, {0, 6, 6}, false, "# \\upsilon^{>=3}",   "\\# \\ vertices \\ mult \\ \\ > \\ 3", true}
};


// --------- CUTS TO PLOT RECO VERTEX XYZ ------------
std::vector<PlotDef> cutDefsRecoVertex = {
    {"RecnuvX",  "0==0", CutType::kNone,   0, {-200, 200, 20}, true, "Slice vertex X",  "No \\ cut", true}
    ,{"RecnuvY", "RecnuvY", CutType::kNone, 0, {-200, 200, 20}, true, "Slice vertex Y", "No \\ cut", true}
    ,{"RecnuvZ", "RecnuvZ", CutType::kNone, 0, {0, 500, 25}, true, "Slice vertex Z", "No \\ cut", true}
};


// --------- CUTS FOR BRAZIL CM ------------
std::vector<PlotDef> cutDefsBrazilCM = {
    
    {"TruthIsFiducial || (IntOrigin==2 || IntDirt==1) ",  "TruthIsFiducial || (IntOrigin==2 || IntDirt==1)", CutType::kCenter, 1, {0,2,2}, true, "Truth in FV",  "Truth \\ in \\ FV"}
    ,{"RecoIsFiducial",   "RecoIsFiducial",  CutType::kCenter, 1, {0,2,2}, true, "Reco in FV",   "Reco \\ in \\ FV"}

    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kNone, fCutFRANSPANDORA, {-.5,.5,40}, false, "FRANS score", "FRANS \\ score", true}
    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, fCutFRANSPANDORA, {-.5,.5,40}, true, "FRANS score", "FRANS \\ score"}
    
    ,{"NAngles",               "NAngles",             CutType::kRight, fCutMinNAngles , {0,5,5}, false, "# V", "\\# \\ V" }
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft, 1, {0, 15, 30}, false, "#Delta#alpha [#circ]", "\\Delta\\alpha\\ [\\circ]"}
    ,{"NUnOriginas",        "NUnOrigins",        CutType::kLeftInt,  0,   {0, 15, 15}, false,  "# extra vertices ",  "\\#  \\ extra \\ vertices"}    
    ,{"NShowers", "NShowers", CutType::kLeft, 1, {0, 10, 10}, false, "# showers", "\\# \\ showers"}
    ,{"NShowerHits", "NShowerHits", CutType::kLeft, 20, {0, 200, 20}, false, "# shower hits", "\\# \\ shower \\ hits"}

    ,{"NAngles>=1 && AngleDecayContainedDiff<=1 && NUnOrigins<=0", "(NAngles>=1 && AngleDecayContainedDiff<=1 && NUnOrigins<=0)", CutType::kCenter,  1,   {0, 2, 2}, true,  "topological cuts",  "\\ Topological \\ cuts"}

    ,{"NShowers", "NShowers", CutType::kLeft, 1, {0, 10, 10}, false, "# showers", "\\# \\ showers"}
    ,{"NShowerHits", "NShowerHits", CutType::kLeft, 20, {0, 200, 20}, false, "# shower hits", "\\# \\ shower \\ hits"}
}; 


// --------- CUTS FOR LOOP MACRO ------------
PlotDef startCutForLoop = {"TruthIsFiducial", "TruthIsFiducial", CutType::kCenter, 1,   {0, 2, 2}, true,  "fiducial volume",  "TruthFiducialVolume"};

PlotDef minimalCutForLoop = {
"RecoIsFiducial==1 && NAngles>=1 && AngleFRANSScore>0.2 && NTracksLI>=1 && NTracksHI>=1 && abs(RecoLeptonX)<190 && abs(RecoLeptonY)<190 && RecoLeptonZ>10 && RecoLeptonZ<490 && AngleLengthMainTrack>10",
"(RecoIsFiducial==1 && NAngles>=1 && AngleFRANSScore>0.2 && NTracksLI>=1 && NTracksHI>=1 && abs(RecoLeptonX)<190 && abs(RecoLeptonY)<190 && RecoLeptonZ>10 && RecoLeptonZ<490 && AngleLengthMainTrack>10)",
CutType::kCenter,  1,   {0, 2, 2}, true,  "minimal cut",  "MinimalCut"};

std::vector<PlotDef> cutDefsLoop{
    {"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft,  2.5, {0, 20, 40}, true, "#Phi [#circ]", "Phi \\ [\\circ]"}
    ,{"AnglePzSign", "AnglePzSign", CutType::kRight,  1, {-4, 4, 8}, true, "Sign(Pz)", "Sign(Pz)"}
    ,{"AngleResidualRangeMaxAngleRMS", "AngleResidualRangeMaxAngleRMS", CutType::kLeft, 50, {0, 360, 36}, true, "VTracksRMS [#circ]", "VTracksRMS \\ [\\circ]"}
    ,{"NMaxDirtUnassociatedHits<10 && NUnOrigins==0", "(NMaxDirtUnassociatedHits<10 && NUnOrigins==0)", CutType::kCenter, 1, {0, 2, 2}, true, "Additional activity", "Additional activity"}
    //,{"NUnOrigins", "NUnOrigins", CutType::kLeftInt, 0, {0, 15, 15}, true,  "# extra vertices ",  "\\#  \\ extra \\ vertices"}    
    //,{"NMaxDirtUnassociatedHits", "NMaxDirtUnassociatedHits", CutType::kLeft, 10, {0, 400, 40}, true, "# unassociated hits", "\\# \\ unassociated \\ hits"}
};


//--------- CUT REPORITORY DEFINITIONS ----------------
std::vector<PlotDef> cutDefsRepository = {

    // Generic
    {"TruthIsAV",  "TruthIsAV", CutType::kCenter, 1, {0,2,2}, true, "Truth in AV",  "Truth \\ in \\ AV"}
    ,{"CRUMBSScore",  "CRUMBSScore",         CutType::kRight,   fCutCRUMBS, {-1,1,20}, true, "CRUMBSScore",  "CRUMBSScore"}
    
    // # origins
    ,{"NOriginsPairOneTwo>0", "NOriginsPairOneTwo>0", CutType::kCenter, 1, {0,2,2}, true, "# origins mult 1>0, # origins mult2>0", "\\#\\ origins\\ mult \\ 1>0, \\# \\ origins \\ mult \\ 2>0"}
    ,{"NAngles",               "NAngles",             CutType::kRight, fCutMinNAngles , {0,5,5}, true, "# V", "\\# \\ V" }
    ,{"NUnOriginsMultGT3", "NUnOriginsMultGT3", CutType::kLeft,  0, {0, 5, 5}, false, "# origins mult 3", "\\# \\ origins \\ mult \\ 3"}
    ,{"NUnOrigins",        "NUnOrigins",        CutType::kLeft,  0,   {0, 15, 15}, true,  "# origins",  "\\#  \\ origins"}
    ,{"NVertexTracks", "NVertexTracks", CutType::kLeft, 3, {0, 10, 10}, true, "# vertex tracks", "NVertexTracks"}

    // FRANS
    ,{"FRANSScorePANDORA", "FRANSScorePANDORA", CutType::kRight, fCutFRANSPANDORA, {-.5,.5,40}, true, "FRANS score PANDORA", "FRANS \\ score \\ PANDORA"}
    ,{"AngleFRANSScore", "AngleFRANSScore", CutType::kRight, fCutFRANS, {-.5,.5,40}, true, "V FRANS score",   "V \\ FRANS \\ score"}
    
    // CLEANESS
    ,{"AngleCoveredArea", "AngleCoveredArea", CutType::kLeft, 0.8, {0, 2, 40}, false, "Covered Area", "Covered \\ Area"}
    ,{"AngleDirtHits", "AngleDirtHits", CutType::kLeft,  10, {0, 20, 40}, false, "#Dirt Hits", "Dirt \\ Hits"}
    ,{"AngleDirtHitsWires", "AngleDirtHitsWires", CutType::kLeft,  10, {0, 20, 40}, false, "#Dirt Hits Wires", "Dirt \\ Hits \\ Wires"}
    ,{"AngleNHitsTrack1", "AngleNHitsTrack1", CutType::kRight, 4, {0, 100, 50}, false, "track 1 # hits ", "Track1 \\ \\# \\ hits"}
    ,{"AngleNHitsTrack2", "AngleNHitsTrack2", CutType::kRight, 4, {0, 100, 50}, false, "track 2 # hits ", "Track2 \\ \\# \\ hits"}
    ,{"AngleLongestIsMain", "AngleLongestIsMain", CutType::kCenter,  1, {0, 2, 2}, false, "Longest is main", "Longest \\ is \\ main"}
    ,{"AngleNHits", "AngleNHits",  CutType::kRight, 10, {0, 200, 40}, false, "Angle # hits", "Angle \\ \\# \\ hits"}
    ,{"NFreeHits",         "NFreeHits",         CutType::kLeft,  30, {0, 200, 40}, false,  "# free hits",  "\\#  \\ free \\ hits"}
    ,{"NUnassociatedHits", "NUnassociatedHits", CutType::kLeft, 15, {0, 200, 40}, false, "# unassociated hits", "\\# \\ unassociated \\ hits"}
    
    // MAIN TRACK
    ,{"AngleLengthMainTrack", "AngleLengthMainTrack", CutType::kRight,  10, {0, 200, 40}, false, "Main track # hits", "Main \\ track \\ length \\ [cm]"}
    ,{"AngleMainTrackOverlap", "AngleMainTrackOverlap", CutType::kRight,  2, {0, 50, 50}, false, "OverlapWithMainTrack", "OverlapWithMainTrack"}
    ,{"AngleGapOverlapWithAPAJuntion", "AngleGapOverlapWithAPAJuntion", CutType::kLeft,  0.1, {0, 1, 10}, true, "AngleGapOverlapWithAPAJuntion", "AngleGapOverlapWithAPAJuntion"}

    // 2D CALORIMETRY VARIABLES
    ,{"AngleTwoLinesChi2", "AngleTwoLinesChi2", CutType::kNone, 1, {0, 1000, 100}, false, "Chi2", "Chi2"}
    ,{"AngleVertexHitIntegralRatio", "AngleVertexHitIntegralRatio", CutType::kNone, 1, {0, 3, 20}, false, "AngleVertexHitIntegralRatio", "AngleVertexHitIntegralRatio"}
    ,{"AngleVertexHitIntegralDifference", "AngleVertexHitIntegralDifference", CutType::kNone, 1, {-10, 10, 100}, false, "AngleVertexHitIntegralDifference", "AngleVertexHitIntegralDifference"}
    ,{"AngleTrackLengthRatio", "AngleTrackLengthRatio", CutType::kNone, 1, {0, 3, 20}, false, "AngleTrackLengthRatio", "AngleTrackLengthRatio"}
    ,{"AngleResidualRangeMinRMS", "AngleResidualRangeMinRMS", CutType::kNone, 1, {0, 2, 40}, false, "AngleResidualRangeMinRMS", "AngleResidualRangeMinRMS"}
    ,{"AngleChargeRatioFit", "AngleChargeRatioFit", CutType::kNone, 1, {0, 3, 20}, false, "AngleChargeRatioFit", "AngleChargeRatioFit"}
    ,{"AngleChargeDifferenceFit", "AngleChargeDifferenceFit", CutType::kNone, 1, {-10, 10, 100}, false, "AngleChargeDifferenceFit", "AngleChargeDifferenceFit"}
    ,{"AngleChargeRatioIntegral", "AngleChargeRatioIntegral", CutType::kNone, 1, {0, 3, 20}, false, "AngleChargeRatioIntegral", "AngleChargeRatioIntegral"}
    ,{"AngleChargeDifferenceIntegral", "AngleChargeDifferenceIntegral", CutType::kNone, 1, {-10, 10, 100}, false, "AngleChargeDifferenceIntegral", "AngleChargeDifferenceIntegral"}
    ,{"AngleNVertexHits", "AngleNVertexHits", CutType::kNone, 1, {0, 100, 20}, false, "AngleNVertexHits", "AngleNVertexHits"}
    ,{"AngleNBulkHits", "AngleNBulkHits", CutType::kNone, 1, {0, 100, 20}, false, "AngleNBulkHits", "AngleNBulkHits"}
    ,{"AngleBandOverlap", "AngleBandOverlap", CutType::kNone, 1, {0, 2, 50}, false, "AngleBandOverlap", "AngleBandOverlap"}
    ,{"AngleBandCrossHits", "AngleBandCrossHits", CutType::kNone, 1, {0, 2, 50}, false, "AngleBandCrossHits", "AngleBandCrossHits"}


    // SHOWERS
    ,{"ShowerEnergy", "ShowerEnergy", CutType::kLeft, fCutShwEnergy, {0, 300, 15}, false, "Shower Energy [MeV]", "ShowerEnergy [MeV]"}
    
    // 3D CALORIMETRY VARIABLES

    // TOPOLIGICAL CUTS
    ,{"AngleDecayContainedDiff", "AngleDecayContainedDiff", CutType::kLeft,  2.5, {0, 20, 40}, false, "#Delta OpeningAngle [#circ]", "DeltaOpeningAngle"}
    ,{"NAngles>=1 && AngleDecayContainedDiff<=1 && NUnOrigins<=0", "(NAngles>=1 && AngleDecayContainedDiff<=1 && NUnOrigins<=0)", CutType::kCenter,  1,   {0, 2, 2}, true,  "topological cuts",  "\\ Topological \\ cuts"}
    ,{"AnglePzSign>=1 && AngleDecayContainedDiff<=2.5 && NUnOrigins<=0 && AngleResidualRangeMaxAngleRMS<=50 && NMaxDirtUnassociatedHits<=10 && AngleLengthMainTrack>=10 && AngleGapOverlapWithAPAJuntion<=0.1", 
    "(AnglePzSign>=1 && AngleDecayContainedDiff<=2.5 && NUnOrigins<=0 && AngleResidualRangeMaxAngleRMS<=50 && NMaxDirtUnassociatedHits<=10 && AngleLengthMainTrack>=10 && AngleGapOverlapWithAPAJuntion<=0.1)", CutType::kCenter, 1, {0, 2, 2}, true, "V characterization", "V characterization"}
    ,{"NVertexTracks", "NVertexTracks", CutType::kLeft, 3, {0, 8, 8}, false, "# vertex tracks", "NVertexTracks"}
};

//--------- SIGNAL AND BG DEFINITIONS REPOSITORY ------------
std::string fLambdaQE = "(IntNLambda>0 && IntMode==0)";
std::string fLambdaQENuMu = "(IntNLambda>0 && IntMode==0 && abs(IntNuPDG)==14)";
std::string fLambdaQENuE = "(IntNLambda>0 && IntMode==0 && abs(IntNuPDG)==12)";
std::string fHyperonNoSignal = "( IntNHyperons>0 && !( IntMode==0 && IntNLambda>0 && abs(IntNuPDG)==14 ) )";
// -- Lambda QE
SampleDef fSaLambdaQENuMu(fTruthInAV+"IntOrigin==1 && IntDirt==0 && "+fLambdaQENuMu, "Signal", true, "Signal", 1);
SampleDef fSaLambdaQENuE(fTruthInAV+"IntOrigin==1 && IntDirt==0 && "+fLambdaQENuE, "#Lambda QE #bar{#nu}_{e}", false, "Lambda QE NuE");
SampleDef fSaLambdaNoSignal(fTruthInAV+"IntOrigin==1 && IntDirt==0 && IntNLambda>0 && !"+fLambdaQENuMu, "#Lambda Inclusive", false, "Lambda Inclusive");
SampleDef fSaLambdaResNuMu(fTruthInAV+"IntOrigin==1 && IntDirt==0 && IntNLambda>0 && IntMode==1 && abs(IntNuPDG)==14", "#Lambda Res #bar{#nu}_{#mu}", false, "Lambda Res NuMu");
SampleDef fSaLambdaDisNuMu(fTruthInAV+"IntOrigin==1 && IntDirt==0 && IntNLambda>0 && IntMode==2 && abs(IntNuPDG)==14", "#Lambda Dis #bar{#nu}_{#mu}", false, "Lambda Dis NuMu");
SampleDef fSaLambdaCohMecNuMu(fTruthInAV+"IntOrigin==1 && IntDirt==0 && IntNLambda>0 && (IntMode==3 || IntMode==10) && abs(IntNuPDG)==14", "#Lambda CohMec #bar{#nu}_{#mu}", false, "Lambda CohMec NuMu");
SampleDef fSaHyperonNoSignal(fTruthInAV+"IntOrigin==1 && IntDirt==0 && "+fHyperonNoSignal, "#splitline{Hyperon}{(other)}", false, "Hyperon (other)", 5);
// -- BNB inclusive
SampleDef fSaBNBInclusive(fTruthInAV+"IntOrigin==1 &&  IntDirt==0 && !"+fLambdaQENuMu, "BG  #nu", false, "BG BNB", 2);
SampleDef fSaBNBInclusiveNoHyperon(fTruthInAV+"IntOrigin==1 &&  IntDirt==0 && !"+fLambdaQENuMu+"&& !"+fHyperonNoSignal, "BG  #nu", false, "BG BNB", 2);
SampleDef fSaBNBInclusiveNoLambda(fTruthInAV+"IntOrigin==1 && IntDirt==0 && IntNLambda==0", "BG  #nu", false, "BG BNB", 2);
// -- Dirt
SampleDef fSaDirt("IntOrigin==1 && IntDirt==1", "Dirt", false, "Dirt", 3);
// -- Cosmic
SampleDef fSaCosmic("IntOrigin==2", "Cosmic", false, "Cosmic", 4);
// -- Dirt and Cosmics
SampleDef fSaCosmicDirt("IntOrigin==2 || (IntOrigin==1 && IntDirt==1)", "#splitline{Dirt+}{cosmic}", false, "Dirt/Cosmic", 3);
// --- BNB (exclusive)
// QE
SampleDef fSaBNBQENuMu(fTruthInAV+"IntOrigin==1 && IntDirt==0 && abs(IntNuPDG)==14 && IntDirt==0 && IntMode==0 && !"+fLambdaQENuMu, "QE", false, "BNB NuMu QE", 2);
SampleDef fSaBNBResNuMu(fTruthInAV+"IntOrigin==1 && IntDirt==0 && abs(IntNuPDG)==14 && IntDirt==0 && IntMode==1", "Res", false, "BNB NuMu Res", 3);
SampleDef fSaBNBDisNuMu(fTruthInAV+"IntOrigin==1 && IntDirt==0 && abs(IntNuPDG)==14 && IntDirt==0 && IntMode==2", "Dis", false, "BNB NuMu Dis", 4);
SampleDef fSaBNBMecNuMu(fTruthInAV+"IntOrigin==1 && IntDirt==0 && abs(IntNuPDG)==14 && IntDirt==0 && IntMode==3", "Mec", false, "BNB NuMu Mec", 5);
SampleDef fSaBNBCohNuMu(fTruthInAV+"IntOrigin==1 && IntDirt==0 && abs(IntNuPDG)==14 && IntDirt==0 && IntMode==10", "Coh", false, "BNB NuMu Coh", 7);
SampleDef fSaBNBNuE(fTruthInAV+"IntOrigin==1 && IntDirt==0 && abs(IntNuPDG)==12", "BNB #nu_{e}", false, "BNB NuE", 7);





std::vector<SampleDef> fSampleDefsByIntMode = {
    fSaLambdaQENuMu
    ,fSaBNBQENuMu
    ,fSaBNBResNuMu
    ,fSaBNBDisNuMu
    //,fSaBNBMecNuMu
    //,fSaBNBCohNuMu
    ,fSaBNBNuE
};


std::vector<SampleDef> fSampleDefsStd = {
    fSaLambdaQENuMu
    ,fSaBNBInclusive
    ,fSaCosmicDirt
    ,fSaHyperonNoSignal
};

std::vector<SampleDef> fSampleDefsCosmicDIrt = {
    fSaLambdaQENuMu
    ,fSaBNBInclusive
    ,fSaCosmic
    ,fSaDirt
};



std::vector<SampleDef> fSampleDefsOrigins = {
    fSaLambdaQENuMu
    ,fSaBNBQENuMu
    ,fSaBNBResNuMu
    ,fSaBNBDisNuMu
};

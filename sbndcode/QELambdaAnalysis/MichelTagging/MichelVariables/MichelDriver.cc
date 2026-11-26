#include "MichelDriver.hh"

void MichelDriver::Compute()
{
    // Set default Michel Vars

    fmichelFeatures.nMuonTracks = -1;
    fmichelFeatures.lowerChi2Ratio = -1;
    fmichelFeatures.angleQ80 = -1;
    fmichelFeatures.dropRadius = -1;
    fmichelFeatures.nBlips = -1;
    fmichelFeatures.score = -99;
    fmichelFeatures.featuresComputed = false;

    // 0. Get number of tracks
    fmichelFeatures.nMuonTracks = fmuTracks.size();

    // 1. Find bragg peak and perform PCA
    MichelGeo michelGeo(fParams.GeoPset, fmuCalos);
    michelGeo.Compute();
    fbragg = michelGeo.getBragg();
    fpca = michelGeo.getPCAResult();

    // std::cout << "Bragg found: " << fbragg.Found << std::endl;
    // std::cout << "PCA valid: " << fpca.Valid << std::endl;

    if (!fbragg.Found || !fpca.Valid){
        fmichelFeatures.featuresComputed = false;
        fmichelFeatures.score = -2;
        return;
    } 

    // 2. Chi2 var
    MichelChi2 michelChi2(fParams.GeoPset, fmuCalos, fbragg);
    fmichelFeatures.lowerChi2Ratio = michelChi2.LowerChi2Ratio();
    fmichelFeatures.allChi2Ratio = michelChi2.AllChi2Ratio();
    // std::cout << "lower chi2 ratio: " << fmichelFeatures.lowerChi2Ratio << std::endl;
    // std::cout << "all chi2 ratio: " << fmichelFeatures.allChi2Ratio << std::endl;

    // 3. Angle var
    MichelAngles michelAngles(fParams.AngleMaxRadius, fmuHits, fbragg, fpca);
    fmichelFeatures.angleQ80 = michelAngles.AnglesQ80();
    // std::cout << "angle q80: " << fmichelFeatures.angleQ80 << std::endl;

    fangles = michelAngles.getAngles();

    // 4. Radial var
    MichelRadial michelRadial(fParams.RadialPset, fmuHits, fbragg, fpca);
    michelRadial.Compute();
    fmichelFeatures.dropRadius = michelRadial.getDropRadius();
    // std::cout << "dropRadius: " << fmichelFeatures.dropRadius << std::endl;

    fradii = michelRadial.getRadii();
    fcharge = michelRadial.getCharge();
    fdQdr = michelRadial.getdQdr();
    fsmoothdQdr = michelRadial.getSmoothdQdr();
    fdropIdx = michelRadial.getDropIdx();
    fmindQdr = michelRadial.getMindQdr();
    
    
    // 5. Blips var
    MichelBlips michelBlips(
        fParams.BlipMaxWireSpan, 
        fParams.BlipMaxDist, 
        fmuBlips, fmuTracks, fbragg, 
        frecoVtxX, frecoVtxY, frecoVtxZ);
    fmichelFeatures.nBlips = michelBlips.NBlips();
    // std::cout << "nBlips: " << fmichelFeatures.nBlips << std::endl;
    
    fmichelFeatures.featuresComputed = true;


    if (fmichelFeatures.featuresComputed) {
        // std::cout << "Computing Michel Score" << std::endl;

        if (!fmichelReader.Loaded()) {
            // std::cout << "ERROR: Michel BDT not loaded" << std::endl;
            fmichelFeatures.score = -3;
            return;
        }

        fmichelFeatures.score = fmichelReader.Evaluate(fmichelFeatures);

    } else {
        fmichelFeatures.score = -2;
    }


}
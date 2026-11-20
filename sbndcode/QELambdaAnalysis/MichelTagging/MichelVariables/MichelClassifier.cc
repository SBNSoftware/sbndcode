#include "MichelClassifier.hh"


void MichelClassifier::Load(std::string const& weightsFile)
    {
        fReader.AddVariable("n_mu_tracks",      &fNMuTracks);
        fReader.AddVariable("lower_chi2_ratio", &fLowerChi2Ratio);
        fReader.AddVariable("angle_q80",        &fAngleQ80);
        fReader.AddVariable("drop_radius",      &fDropRadius);
        fReader.AddVariable("n_blips",          &fNBlips);

        fReader.BookMVA("BDT", weightsFile);
        fLoaded = true;
    }


float MichelClassifier::Evaluate(const MichelFeatures& michelFeatures)
{
    // CORRECT ORDER 
    //  n_mu_tracks 
    //  lower_chi2_ratio
    //  angle_q80
    //  drop_radius
    //  n_blips

    fNMuTracks      = michelFeatures.nMuonTracks;
    fLowerChi2Ratio = michelFeatures.lowerChi2Ratio;
    fAngleQ80       = michelFeatures.angleQ80;
    fDropRadius     = michelFeatures.dropRadius;
    fNBlips         = michelFeatures.nBlips;


    return fReader.EvaluateMVA("BDT");
}
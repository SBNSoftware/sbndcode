#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TString.h>
#include <vector>
#include <string>
#include <iostream>

void generate_sigma_u_tree(const char* outFileName = "/exp/sbnd/data/users/coackley/sigma_u_universes.root", int nUniverses = 100){

    std::vector<std::string> knobNames = {
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi",
        "NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi",

        "MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu",
        "MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar",
        "MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression",
        "MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio",
        "MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio",

        "MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement",

        "MINERvAE2p2h_SBN_v1_E2p2h_A_nu",
        "MINERvAE2p2h_SBN_v1_E2p2h_A_nubar",
        "MINERvAE2p2h_SBN_v1_E2p2h_B_nu",
        "MINERvAE2p2h_SBN_v1_E2p2h_B_nubar",

        "GENIEReWeight_SBN_v1_multisigma_AhtBY",
        "GENIEReWeight_SBN_v1_multisigma_BhtBY",
        "GENIEReWeight_SBN_v1_multisigma_CV1uBY",
        "GENIEReWeight_SBN_v1_multisigma_CV2uBY",
        "GENIEReWeight_SBN_v1_multisigma_CoulombCCQE",
        "GENIEReWeight_SBN_v1_multisigma_DecayAngMEC",
        "GENIEReWeight_SBN_v1_multisigma_EtaNCEL",
        "GENIEReWeight_SBN_v1_multisigma_FrAbs_N",
        "GENIEReWeight_SBN_v1_multisigma_FrAbs_pi",
        "GENIEReWeight_SBN_v1_multisigma_FrCEx_N",
        "GENIEReWeight_SBN_v1_multisigma_FrCEx_pi",
        "GENIEReWeight_SBN_v1_multisigma_FrInel_N",
        "GENIEReWeight_SBN_v1_multisigma_FrInel_pi",
        "GENIEReWeight_SBN_v1_multisigma_FrPiProd_N",
        "GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi",
        "GENIEReWeight_SBN_v1_multisigma_MFP_N",
        "GENIEReWeight_SBN_v1_multisigma_MFP_pi",
        "GENIEReWeight_SBN_v1_multisigma_MaCCRES",
        "GENIEReWeight_SBN_v1_multisigma_MaNCEL",
        "GENIEReWeight_SBN_v1_multisigma_MaNCRES",
        "GENIEReWeight_SBN_v1_multisigma_MvCCRES",
        "GENIEReWeight_SBN_v1_multisigma_MvNCRES",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi",
        "GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi",
        "GENIEReWeight_SBN_v1_multisigma_NormCCCOH",
        "GENIEReWeight_SBN_v1_multisigma_NormCCMEC",
        "GENIEReWeight_SBN_v1_multisigma_NormNCCOH",
        "GENIEReWeight_SBN_v1_multisigma_NormNCMEC",
        "GENIEReWeight_SBN_v1_multisigma_RDecBR1eta",
        "GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma",
        "GENIEReWeight_SBN_v1_multisigma_RPA_CCQE",
        "GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad",
        "GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi",
        "GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape",
        "GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE",
        "GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE",
        "GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE",
        "GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE"
    };

    int nKnobs = (int)knobNames.size();
    std::cout << "Generating sigma_u for " << nKnobs << " knobs, " << nUniverses << " universes each." << std::endl;

    TFile* outFile = new TFile(outFileName, "RECREATE");
    TTree* tree = new TTree("sigma_u_tree", "Deterministic sigma_u draws per knob per universe");

    int universeIndex;
    tree->Branch("universe_index", &universeIndex, "universe_index/I");

    std::vector<float> sigmaVals(nKnobs, 0.0f);
    for (int k = 0; k < nKnobs; ++k){
        TString branchName = knobNames[k] + "_sigmau";
        tree->Branch(branchName, &sigmaVals[k], branchName + "/F");
    }

    for (int u = 0; u < nUniverses; ++u){
        universeIndex = u;
        for (int k = 0; k < nKnobs; ++k){
            TString seedString = TString::Format("%s_%d", knobNames[k].c_str(), u);
            UInt_t seed = seedString.Hash();

            TRandom3 rng(seed);
            sigmaVals[k] = rng.Gaus(0.0, 1.0);
        }
        tree->Fill();
    }

    tree->Write();
    outFile->Close();

    std::cout << "Wrote sigma_u tree to " << outFileName << std::endl;
}

////////////////////////////////////////////////////////////////////////
// Class:       NuEWithSysts
// Plugin Type: analyzer (Unknown Unknown)
// File:        NuEWithSysts_module.cc
//
// Generated at Fri Oct 31 08:23:53 2025 by Rachel Coackley using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

// Sim Base
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

// Reco Base
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "sbnobj/Common/Reco/MVAPID.h"

// Tools
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

// LArSoft
#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
// #include "larcorealg/Geometry/ChannelMapStandardAlg.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcore/Geometry/WireReadout.h"

#include "sbnobj/Common/Reco/CRUMBSResult.h"

#include "larcoreobj/SummaryData/POTSummary.h"
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"

// C++
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <array>

// ROOT
#include "art_root_io/TFileService.h"
#include <TTree.h>
#include <TFile.h>

constexpr double def_double = -std::numeric_limits<double>::max();

namespace sbnd {
  class NuEWithSysts;
}


class sbnd::NuEWithSysts : public art::EDAnalyzer {
public:
  explicit NuEWithSysts(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuEWithSysts(NuEWithSysts const&) = delete;
  NuEWithSysts(NuEWithSysts&&) = delete;
  NuEWithSysts& operator=(NuEWithSysts const&) = delete;
  NuEWithSysts& operator=(NuEWithSysts&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void Slices(art::Event const& e);
  void trueSignal(art::Event const& e);
  void clearVectors();

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  unsigned int eventID; // Event num
  unsigned int runID; // Run num  
  unsigned int subRunID; // Subrun num

  int nuEScatter; // 0 = no nu+e scatter, 1 = nu+e scatter in event
  double nuEScatterTrueVX;
  double nuEScatterTrueVY;
  double nuEScatterTrueVZ;
  std::vector<float> nuEScatter_MCTruthFlux_weight_horncurrent;
  std::vector<float> nuEScatter_MCTruthFlux_weight_expskin;
  std::vector<float> nuEScatter_MCTruthFlux_weight_pioninexsec;
  std::vector<float> nuEScatter_MCTruthFlux_weight_pionqexsec;
  std::vector<float> nuEScatter_MCTruthFlux_weight_piontotxsec;
  std::vector<float> nuEScatter_MCTruthFlux_weight_nucleoninexsec;
  std::vector<float> nuEScatter_MCTruthFlux_weight_nucleonqexsec;
  std::vector<float> nuEScatter_MCTruthFlux_weight_nucleontotxsec;
  std::vector<float> nuEScatter_MCTruthFlux_weight_kplus;
  std::vector<float> nuEScatter_MCTruthFlux_weight_kmin;
  std::vector<float> nuEScatter_MCTruthFlux_weight_kzero;
  std::vector<float> nuEScatter_MCTruthFlux_weight_piplus;
  std::vector<float> nuEScatter_MCTruthFlux_weight_piminus;
  
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE;
  std::vector<float> nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE;
  

  TTree* NuEWeightsTree;
  std::vector<double>   reco_sliceID;
  std::vector<double>   reco_sliceInteraction;
  std::vector<double>   reco_sliceTrueVX;
  std::vector<double>   reco_sliceTrueVY;
  std::vector<double>   reco_sliceTrueVZ;
  std::vector<double>   reco_sliceOrigin;
  std::vector<double>   reco_sliceTrueCCNC;
  std::vector<double>   reco_sliceTrueNeutrinoType;

  std::vector<std::vector<float>>   reco_sliceMCTruthFlux_weight_horncurrent;
  std::vector<std::vector<float>>   reco_sliceMCTruthFlux_weight_expskin;
  std::vector<std::vector<float>>   reco_sliceMCTruthFlux_weight_pioninexsec;
  std::vector<std::vector<float>>   reco_sliceMCTruthFlux_weight_pionqexsec;
  std::vector<std::vector<float>>   reco_sliceMCTruthFlux_weight_piontotxsec;
  std::vector<std::vector<float>>   reco_sliceMCTruthFlux_weight_nucleoninexsec;
  std::vector<std::vector<float>>   reco_sliceMCTruthFlux_weight_nucleonqexsec;
  std::vector<std::vector<float>>   reco_sliceMCTruthFlux_weight_nucleontotxsec;
  std::vector<std::vector<float>>   reco_sliceMCTruthFlux_weight_kplus;
  std::vector<std::vector<float>>   reco_sliceMCTruthFlux_weight_kmin;
  std::vector<std::vector<float>>   reco_sliceMCTruthFlux_weight_kzero;
  std::vector<std::vector<float>>   reco_sliceMCTruthFlux_weight_piplus;
  std::vector<std::vector<float>>   reco_sliceMCTruthFlux_weight_piminus;
  
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE;
  std::vector<std::vector<float>> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE;

  art::ServiceHandle<cheat::ParticleInventoryService> particleInv;
  art::ServiceHandle<geo::Geometry> theGeometry;
  detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();

  const std::string PFParticleLabel;
  const std::string sliceLabel;
  const std::string sliceSCELabel;
  const std::string vertexLabel;
  const std::string nuGenModuleLabel;
  const std::string hitLabel;
  const std::string crumbsLabel;
  const std::string razzledLabel;
  const std::string trackLabel;
  const std::string showerLabel;
  const std::string MCTruthLabel;
  const std::string TruthLabel;
  const std::string spacePointLabel;
  double DLCurrent;
  const std::string clusterLabel;
  const std::string POTModuleLabel;
  double signal;
  const std::string fluxWeightLabel;
  const std::string genieWeightLabel;
  const std::string mcTruthFluxLabel;

  TFile *outputFile = TFile::Open("NuEWeightsAnalyserOutput.root","RECREATE");

};


sbnd::NuEWithSysts::NuEWithSysts(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  PFParticleLabel(p.get<std::string>("PFParticleLabel")),
  sliceLabel(p.get<std::string>("SliceLabel")),
  sliceSCELabel(p.get<std::string>("SliceSCELabel")),
  vertexLabel(p.get<std::string>("VertexLabel")),
  nuGenModuleLabel(p.get<std::string>("NuGenModuleLabel")),
  hitLabel(p.get<std::string>("HitLabel")),
  crumbsLabel(p.get<std::string>("CRUMBSLabel")),
  razzledLabel(p.get<std::string>("RazzledLabel")),
  trackLabel(p.get<std::string>("TrackLabel")),
  showerLabel(p.get<std::string>("ShowerLabel")),
  MCTruthLabel(p.get<std::string>("MCTruthLabel")),
  TruthLabel(p.get<std::string>("TruthLabel")),
  spacePointLabel(p.get<std::string>("SpacePointLabel")),
  DLCurrent(p.get<double>("DLCurrent")),
  clusterLabel(p.get<std::string>("ClusterLabel")),
  POTModuleLabel(p.get<std::string>("POTLabel")),
  signal(p.get<double>("Signal")),
  fluxWeightLabel(p.get<std::string>("FluxWeightLabel")),
  genieWeightLabel(p.get<std::string>("GenieWeightLabel")),
  mcTruthFluxLabel(p.get<std::string>("MCTruthFluxLabel"))
{

  art::ServiceHandle<art::TFileService> fs;  

  NuEWeightsTree = fs->make<TTree>("NuEWeights","");
  NuEWeightsTree->Branch("eventID", &eventID);
  NuEWeightsTree->Branch("runID", &runID);
  NuEWeightsTree->Branch("subRunID", &subRunID);
  NuEWeightsTree->Branch("DLCurrent", &DLCurrent);
  NuEWeightsTree->Branch("signal", &signal);
  
  NuEWeightsTree->Branch("nuEScatter", &nuEScatter);
  NuEWeightsTree->Branch("nuEScatterTrueVX", &nuEScatterTrueVX);
  NuEWeightsTree->Branch("nuEScatterTrueVY", &nuEScatterTrueVY);
  NuEWeightsTree->Branch("nuEScatterTrueVZ", &nuEScatterTrueVZ);
  
  NuEWeightsTree->Branch("nuEScatter_MCTruthFlux_weight_horncurrent", &nuEScatter_MCTruthFlux_weight_horncurrent);
  NuEWeightsTree->Branch("nuEScatter_MCTruthFlux_weight_expskin", &nuEScatter_MCTruthFlux_weight_expskin);
  NuEWeightsTree->Branch("nuEScatter_MCTruthFlux_weight_pioninexsec", &nuEScatter_MCTruthFlux_weight_pioninexsec);
  NuEWeightsTree->Branch("nuEScatter_MCTruthFlux_weight_pionqexsec", &nuEScatter_MCTruthFlux_weight_pionqexsec);
  NuEWeightsTree->Branch("nuEScatter_MCTruthFlux_weight_piontotxsec", &nuEScatter_MCTruthFlux_weight_piontotxsec);
  NuEWeightsTree->Branch("nuEScatter_MCTruthFlux_weight_nucleoninexsec", &nuEScatter_MCTruthFlux_weight_nucleoninexsec);
  NuEWeightsTree->Branch("nuEScatter_MCTruthFlux_weight_nucleonqexsec", &nuEScatter_MCTruthFlux_weight_nucleonqexsec);
  NuEWeightsTree->Branch("nuEScatter_MCTruthFlux_weight_nucleontotxsec", &nuEScatter_MCTruthFlux_weight_nucleontotxsec);
  NuEWeightsTree->Branch("nuEScatter_MCTruthFlux_weight_kplus", &nuEScatter_MCTruthFlux_weight_kplus);
  NuEWeightsTree->Branch("nuEScatter_MCTruthFlux_weight_kmin", &nuEScatter_MCTruthFlux_weight_kmin);
  NuEWeightsTree->Branch("nuEScatter_MCTruthFlux_weight_kzero", &nuEScatter_MCTruthFlux_weight_kzero);
  NuEWeightsTree->Branch("nuEScatter_MCTruthFlux_weight_piplus", &nuEScatter_MCTruthFlux_weight_piplus);
  NuEWeightsTree->Branch("nuEScatter_MCTruthFlux_weight_piminus", &nuEScatter_MCTruthFlux_weight_piminus);
 
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu", &nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar", &nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression", &nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio", &nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio", &nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement", &nuEScatter_MCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu", &nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar", &nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu", &nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar", &nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE);
  NuEWeightsTree->Branch("nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE);

  NuEWeightsTree->Branch("reco_sliceID", &reco_sliceID);
  NuEWeightsTree->Branch("reco_sliceInteraction", &reco_sliceInteraction);
  NuEWeightsTree->Branch("reco_sliceTrueVX", &reco_sliceTrueVX);
  NuEWeightsTree->Branch("reco_sliceTrueVY", &reco_sliceTrueVY);
  NuEWeightsTree->Branch("reco_sliceTrueVZ", &reco_sliceTrueVZ);
  NuEWeightsTree->Branch("reco_sliceOrigin", &reco_sliceOrigin);
  NuEWeightsTree->Branch("reco_sliceTrueCCNC", &reco_sliceTrueCCNC);
  NuEWeightsTree->Branch("reco_sliceTrueNeutrinoType", &reco_sliceTrueNeutrinoType);
  
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_horncurrent", &reco_sliceMCTruthFlux_weight_horncurrent);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_expskin", &reco_sliceMCTruthFlux_weight_expskin);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_pioninexsec", &reco_sliceMCTruthFlux_weight_pioninexsec);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_pionqexsec", &reco_sliceMCTruthFlux_weight_pionqexsec);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_piontotxsec", &reco_sliceMCTruthFlux_weight_piontotxsec);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_nucleoninexsec", &reco_sliceMCTruthFlux_weight_nucleoninexsec);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_nucleonqexsec", &reco_sliceMCTruthFlux_weight_nucleonqexsec);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_nucleontotxsec", &reco_sliceMCTruthFlux_weight_nucleontotxsec);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_kplus", &reco_sliceMCTruthFlux_weight_kplus);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_kmin", &reco_sliceMCTruthFlux_weight_kmin);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_kzero", &reco_sliceMCTruthFlux_weight_kzero);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_piplus", &reco_sliceMCTruthFlux_weight_piplus);
  NuEWeightsTree->Branch("reco_sliceMCTruthFlux_weight_piminus", &reco_sliceMCTruthFlux_weight_piminus);
 
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement", &reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu", &reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar", &reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu", &reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar", &reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE);
  NuEWeightsTree->Branch("reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE);
  
}

void sbnd::NuEWithSysts::analyze(art::Event const& e)
{
    clearVectors();

    eventID = e.id().event();
    runID = e.id().run();
    subRunID = e.id().subRun();

    //std::cout << "" << std::endl;
    //std::cout << "========================================================================================================" << std::endl; 
    //std::cout << "Run: " << runID << ", Subrun: " << subRunID << ", Event: " << eventID << ", DL/Current: " << DLCurrent << std::endl;
    Slices(e);
    trueSignal(e);

    NuEWeightsTree->Fill();
}

void sbnd::NuEWithSysts::trueSignal(art::Event const& e){
    //std::cout << "============= True NuE Scatters =============" << std::endl;

    art::Handle<std::vector<simb::MCTruth>> MCTruthHandle;
    std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
    if(e.getByLabel(TruthLabel, MCTruthHandle))
        art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

    art::FindManyP<std::map<std::string,std::vector<float>>> mcTruthFluxWeightMapAssns(MCTruthVec, e, fluxWeightLabel);
    art::FindManyP<std::map<std::string,std::vector<float>>> mcTruthGENIEWeightMapAssns(MCTruthVec, e, genieWeightLabel);

    int numNuEScatters = 0;

    double nuEScatterVX = -999999;
    double nuEScatterVY = -999999;
    double nuEScatterVZ = -999999;

    if(!MCTruthVec.empty()){
        for(auto &MCTruth : MCTruthVec){
            if(MCTruth->Origin() == simb::kBeamNeutrino){
                simb::MCNeutrino neutrino = MCTruth->GetNeutrino();
                simb::MCParticle neutrinoParticle = neutrino.Nu();
                if(neutrino.InteractionType() == 1098){
                    // This is a nu+e elastic scattering event
                    numNuEScatters++;

                    nuEScatterVX = neutrinoParticle.Vx();
                    nuEScatterVY = neutrinoParticle.Vy();
                    nuEScatterVZ = neutrinoParticle.Vz();

                    if(numNuEScatters == 1){
                        const std::vector<art::Ptr<std::map<std::string,std::vector<float>>>> mcTruthFluxWeightMaps(mcTruthFluxWeightMapAssns.at(MCTruth.key()));

                        //std::cout << "Size of mcTruthFluxWeightMaps vector = " << mcTruthFluxWeightMaps.size() << std::endl;
                        const art::Ptr<std::map<std::string, std::vector<float>>> &weightFluxMapPtr = mcTruthFluxWeightMaps.at(0);
                        const std::map<std::string, std::vector<float>> &weightFluxMap = *weightFluxMapPtr;

                        //std::cout << "Number of map entries = " << weightFluxMap.size() << std::endl;
                        
                        std::unordered_map<std::string, std::vector<float>*> targetMap = {
                            {"horncurrent_Flux",      &nuEScatter_MCTruthFlux_weight_horncurrent},
                            {"expskin_Flux",          &nuEScatter_MCTruthFlux_weight_expskin},
                            {"pioninexsec_Flux",      &nuEScatter_MCTruthFlux_weight_pioninexsec},
                            {"pionqexsec_Flux",       &nuEScatter_MCTruthFlux_weight_pionqexsec},
                            {"piontotxsec_Flux",      &nuEScatter_MCTruthFlux_weight_piontotxsec},
                            {"nucleoninexsec_Flux",   &nuEScatter_MCTruthFlux_weight_nucleoninexsec},
                            {"nucleonqexsec_Flux",    &nuEScatter_MCTruthFlux_weight_nucleonqexsec},
                            {"nucleontotxsec_Flux",   &nuEScatter_MCTruthFlux_weight_nucleontotxsec},
                            {"kplus_Flux",            &nuEScatter_MCTruthFlux_weight_kplus},
                            {"kminus_Flux",           &nuEScatter_MCTruthFlux_weight_kmin},
                            {"kzero_Flux",            &nuEScatter_MCTruthFlux_weight_kzero},
                            {"piplus_Flux",           &nuEScatter_MCTruthFlux_weight_piplus},
                            {"piminus_Flux",          &nuEScatter_MCTruthFlux_weight_piminus} 
                        };

                        for(const auto &[name, values] : weightFluxMap){
                            
                            //std::cout << "Map key = " << name << std::endl;
                            //std::cout << "Number of universes? = " << values.size() << std::endl;

                            //for(float v : values){
                                //std::cout << "  Weight = " << v << std::endl;
                            //}
                            
                            auto mapFill = targetMap.find(name);
                            if(mapFill == targetMap.end()) continue;

                            auto &target = *(mapFill->second);

                            target.assign(values.begin(), values.end());
                        }

                        const std::vector<art::Ptr<std::map<std::string,std::vector<float>>>> mcTruthGENIEWeightMaps(mcTruthGENIEWeightMapAssns.at(MCTruth.key()));

                        std::unordered_map<std::string, std::vector<float>*> targetGENIEMap = {
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi},
                            {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi", &nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi},
                            {"MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu", &nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu},
                            {"MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar", &nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar},
                            {"MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression", &nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression},
                            {"MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio", &nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio},
                            {"MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio", &nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio},
                            {"MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement", &nuEScatter_MCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement},
                            {"MINERvAE2p2h_SBN_v1_E2p2h_A_nu", &nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu},
                            {"MINERvAE2p2h_SBN_v1_E2p2h_A_nubar", &nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar},
                            {"MINERvAE2p2h_SBN_v1_E2p2h_B_nu", &nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu},
                            {"MINERvAE2p2h_SBN_v1_E2p2h_B_nubar", &nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar},
                            {"GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse},
                            {"GENIEReWeight_SBN_v1_multisim_COHVariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse},
                            {"GENIEReWeight_SBN_v1_multisim_CoulombCCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE},
                            {"GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse},
                            {"GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse},
                            {"GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse},
                            {"GENIEReWeight_SBN_v1_multisim_NCELVariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse},
                            {"GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse},
                            {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi},
                            {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi},
                            {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi},
                            {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi},
                            {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi},
                            {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi},
                            {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi},
                            {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi},
                            {"GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi},
                            {"GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi},
                            {"GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi},
                            {"GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi},
                            {"GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi},
                            {"GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi},
                            {"GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi},
                            {"GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi},
                            {"GENIEReWeight_SBN_v1_multisim_NormCCMEC", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC},
                            {"GENIEReWeight_SBN_v1_multisim_NormNCMEC", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC},
                            {"GENIEReWeight_SBN_v1_multisim_RDecBR1eta", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta},
                            {"GENIEReWeight_SBN_v1_multisim_RDecBR1gamma", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma},
                            {"GENIEReWeight_SBN_v1_multisim_RPA_CCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE},
                            {"GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse},
                            {"GENIEReWeight_SBN_v1_multisigma_AhtBY", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY},
                            {"GENIEReWeight_SBN_v1_multisigma_BhtBY", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY},
                            {"GENIEReWeight_SBN_v1_multisigma_CV1uBY", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY},
                            {"GENIEReWeight_SBN_v1_multisigma_CV2uBY", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY},
                            {"GENIEReWeight_SBN_v1_multisigma_CoulombCCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE},
                            {"GENIEReWeight_SBN_v1_multisigma_DecayAngMEC", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC},
                            {"GENIEReWeight_SBN_v1_multisigma_EtaNCEL", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL},
                            {"GENIEReWeight_SBN_v1_multisigma_FrAbs_N", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N},
                            {"GENIEReWeight_SBN_v1_multisigma_FrAbs_pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi},
                            {"GENIEReWeight_SBN_v1_multisigma_FrCEx_N", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N},
                            {"GENIEReWeight_SBN_v1_multisigma_FrCEx_pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi},
                            {"GENIEReWeight_SBN_v1_multisigma_FrInel_N", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N},
                            {"GENIEReWeight_SBN_v1_multisigma_FrInel_pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi},
                            {"GENIEReWeight_SBN_v1_multisigma_FrPiProd_N", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N},
                            {"GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi},
                            {"GENIEReWeight_SBN_v1_multisigma_MFP_N", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N},
                            {"GENIEReWeight_SBN_v1_multisigma_MFP_pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi},
                            {"GENIEReWeight_SBN_v1_multisigma_MaCCRES", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES},
                            {"GENIEReWeight_SBN_v1_multisigma_MaNCEL", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL},
                            {"GENIEReWeight_SBN_v1_multisigma_MaNCRES", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES},
                            {"GENIEReWeight_SBN_v1_multisigma_MvCCRES", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES},
                            {"GENIEReWeight_SBN_v1_multisigma_MvNCRES", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES},
                            {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi},
                            {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi},
                            {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi},
                            {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi},
                            {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi},
                            {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi},
                            {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi},
                            {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi},
                            {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi},
                            {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi},
                            {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi},
                            {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi},
                            {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi},
                            {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi},
                            {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi},
                            {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi},
                            {"GENIEReWeight_SBN_v1_multisigma_NormCCCOH", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH},
                            {"GENIEReWeight_SBN_v1_multisigma_NormCCMEC", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC},
                            {"GENIEReWeight_SBN_v1_multisigma_NormNCCOH", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH},
                            {"GENIEReWeight_SBN_v1_multisigma_NormNCMEC", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC},
                            {"GENIEReWeight_SBN_v1_multisigma_RDecBR1eta", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta},
                            {"GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma},
                            {"GENIEReWeight_SBN_v1_multisigma_RPA_CCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE},
                            {"GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad},
                            {"GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi},
                            {"GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape},
                            {"GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE},
                            {"GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE},
                            {"GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE},
                            {"GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE", &nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE},
                        };

                        for(const auto &weightGENIEMapPtr : mcTruthGENIEWeightMaps){
                            const std::map<std::string, std::vector<float>> &weightGENIEMap = *weightGENIEMapPtr;
                            for(const auto &[name, values] : weightGENIEMap){
                                auto mapFill = targetGENIEMap.find(name);
                                if(mapFill == targetGENIEMap.end()) continue;
                                auto &target = *(mapFill->second);
                                target.assign(values.begin(), values.end());
                            }
                        }

                        /*
                        for(size_t i = 0; i < mcTruthGENIEWeightMaps.size(); ++i){
                            const art::Ptr<std::map<std::string, std::vector<float>>> &weightGENIEMapPtr_3Jul = mcTruthGENIEWeightMaps.at(i);
                            const std::map<std::string, std::vector<float>> &weightGENIEMap_3Jul = *weightGENIEMapPtr_3Jul;
                            //std::cout << "  Number of map entries = " << weightGENIEMap_3Jul.size() << std::endl;
                            for(const auto &[name_3Jul, values_3Jul] : weightGENIEMap_3Jul){
                                //std::cout << "  Map key = " << name_3Jul << std::endl;
                                //std::cout << "  Number of universes = " << values_3Jul.size() << std::endl;
                                for(size_t j = 0; j < values_3Jul.size(); j++){
                                    //std::cout << "      Weight " << j+1 << ": " << values_3Jul.at(j) << std::endl;
                                }
                            }
                        }
                        */

                    }

                }
            }
        }
    }

    //if(numNuEScatters > 1) std::cout << "MORE THAN 1 NU+E ELASTIC SCATTER!!!!!!!!!" << std::endl;
    //std::cout << "Number of nu+e elastic scatters in the event = " << numNuEScatters << std::endl;
    //std::cout << "True neutrino vertex of nu+e elastic scatter = (" << nuEScatterVX << ", " << nuEScatterVY << ", " << nuEScatterVZ << ")" << std::endl;

    if(numNuEScatters != 0){
        nuEScatter = 1;
        nuEScatterTrueVX = nuEScatterVX; 
        nuEScatterTrueVY = nuEScatterVY; 
        nuEScatterTrueVZ = nuEScatterVZ;
    } else{
        nuEScatter = 0;
        nuEScatterTrueVX = nuEScatterVX; 
        nuEScatterTrueVY = nuEScatterVY; 
        nuEScatterTrueVZ = nuEScatterVZ;

        nuEScatter_MCTruthFlux_weight_horncurrent.push_back(-999999); 
        nuEScatter_MCTruthFlux_weight_expskin.push_back(-999999); 
        nuEScatter_MCTruthFlux_weight_pioninexsec.push_back(-999999); 
        nuEScatter_MCTruthFlux_weight_pionqexsec.push_back(-999999); 
        nuEScatter_MCTruthFlux_weight_piontotxsec.push_back(-999999); 
        nuEScatter_MCTruthFlux_weight_nucleoninexsec.push_back(-999999); 
        nuEScatter_MCTruthFlux_weight_nucleonqexsec.push_back(-999999); 
        nuEScatter_MCTruthFlux_weight_nucleontotxsec.push_back(-999999); 
        nuEScatter_MCTruthFlux_weight_kplus.push_back(-999999); 
        nuEScatter_MCTruthFlux_weight_kmin.push_back(-999999); 
        nuEScatter_MCTruthFlux_weight_kzero.push_back(-999999); 
        nuEScatter_MCTruthFlux_weight_piplus.push_back(-999999); 
        nuEScatter_MCTruthFlux_weight_piminus.push_back(-999999); 
    
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE.push_back(-999999);
        nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE.push_back(-999999);
    }
    //std::cout << "=============================================" << std::endl;
}

void sbnd::NuEWithSysts::Slices(art::Event const& e){
    const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
    //std::cout << "_________ Slices _________" << std::endl;

    int counter = 0;

    art::Handle<std::vector<recob::Slice>>  sliceHandle;
    std::vector<art::Ptr<recob::Slice>>     sliceVec;
    if(e.getByLabel(sliceLabel, sliceHandle))
        art::fill_ptr_vector(sliceVec, sliceHandle);

    art::Handle<std::vector<simb::MCTruth>> mcTruthHandle;
    std::vector<art::Ptr<simb::MCTruth>>    mcTruthVec;
    if(e.getByLabel(mcTruthFluxLabel, mcTruthHandle))
        art::fill_ptr_vector(mcTruthVec, mcTruthHandle);

    //std::cout << "------------ Number of slices = " << sliceVec.size() << std::endl;

    if(!sliceVec.empty()){
        art::Handle<std::vector<recob::Hit>>    hitHandle;
        std::vector<art::Ptr<recob::Hit>>       hitVec;

        if(e.getByLabel(hitLabel, hitHandle))
            art::fill_ptr_vector(hitVec, hitHandle);

        if(!hitVec.empty()){
            int sliceID(std::numeric_limits<int>::max());

            for(const art::Ptr<recob::Slice> &slice : sliceVec){
                sliceID = slice->ID();
                //std::cout << "Slice ID = " << sliceID << std::endl;
                if(sliceID == std::numeric_limits<int>::max()) continue;

                art::FindManyP<recob::Hit> sliceHitAssns(sliceVec, e, sliceLabel);
                const std::vector<art::Ptr<recob::Hit>> sliceHits(sliceHitAssns.at(slice.key()));


                // Gets the true particle ID of the truth particle who owns the most hits in the slice. True is for rollup.
                const int sliceTrueParticleID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, sliceHits, true);
                //std::cout << "The slice has most hits coming from true particle with ID = " << sliceTrueParticleID << std::endl;
                //std::cout << "Number of hits in slice = " << sliceHits.size() << std::endl;

                if(sliceTrueParticleID == std::numeric_limits<int>::min()) continue;

                // Get the MCParticle associated with the majority of the hits in the slice
                //const simb::MCParticle* sliceMCParticle = particleInv->TrackIdToParticle_P(sliceTrueParticleID);
                
                // Get the MCTruth associated with the majority of the hits in the slice
                const art::Ptr<simb::MCTruth> sliceMCTruth = particleInv->TrackIdToMCTruth_P(sliceTrueParticleID);

                // Get the MCNeutrino associated with the majority of the hits in the slice
                const simb::MCNeutrino sliceMCNeutrino = sliceMCTruth->GetNeutrino();
                const simb::MCParticle sliceMCNeutrinoParticle = sliceMCNeutrino.Nu();

                // Looking at flux systematics associated with MCTruth of slice
                art::FindManyP<std::map<std::string,std::vector<float>>> mcTruthFluxWeightMapAssns(mcTruthVec, e, fluxWeightLabel);
                const std::vector<art::Ptr<std::map<std::string,std::vector<float>>>> mcTruthFluxWeightMaps(mcTruthFluxWeightMapAssns.at(sliceMCTruth.key()));
        
                std::vector<float> reco_sliceMCTruthFlux_weight_horncurrent_vector;
                std::vector<float> reco_sliceMCTruthFlux_weight_expskin_vector;
                std::vector<float> reco_sliceMCTruthFlux_weight_pioninexsec_vector;
                std::vector<float> reco_sliceMCTruthFlux_weight_pionqexsec_vector;
                std::vector<float> reco_sliceMCTruthFlux_weight_piontotxsec_vector;
                std::vector<float> reco_sliceMCTruthFlux_weight_nucleoninexsec_vector;
                std::vector<float> reco_sliceMCTruthFlux_weight_nucleonqexsec_vector;
                std::vector<float> reco_sliceMCTruthFlux_weight_nucleontotxsec_vector;
                std::vector<float> reco_sliceMCTruthFlux_weight_kplus_vector;
                std::vector<float> reco_sliceMCTruthFlux_weight_kmin_vector;
                std::vector<float> reco_sliceMCTruthFlux_weight_kzero_vector;
                std::vector<float> reco_sliceMCTruthFlux_weight_piplus_vector;
                std::vector<float> reco_sliceMCTruthFlux_weight_piminus_vector;

                //std::cout << "============================" << std::endl;
                //std::cout << "Size of mcTruthFluxWeightMaps vector = " << mcTruthFluxWeightMaps.size() << std::endl;
                const art::Ptr<std::map<std::string, std::vector<float>>> &weightFluxMapPtr = mcTruthFluxWeightMaps.at(0);
                const std::map<std::string, std::vector<float>> &weightFluxMap = *weightFluxMapPtr;

                //std::cout << "Number of map entries = " << weightFluxMap.size() << std::endl;

                std::unordered_map<std::string, std::vector<float>*> targetMap = {
                    {"horncurrent_Flux",      &reco_sliceMCTruthFlux_weight_horncurrent_vector},
                    {"expskin_Flux",          &reco_sliceMCTruthFlux_weight_expskin_vector},
                    {"pioninexsec_Flux",      &reco_sliceMCTruthFlux_weight_pioninexsec_vector},
                    {"pionqexsec_Flux",       &reco_sliceMCTruthFlux_weight_pionqexsec_vector},
                    {"piontotxsec_Flux",      &reco_sliceMCTruthFlux_weight_piontotxsec_vector},
                    {"nucleoninexsec_Flux",   &reco_sliceMCTruthFlux_weight_nucleoninexsec_vector},
                    {"nucleonqexsec_Flux",    &reco_sliceMCTruthFlux_weight_nucleonqexsec_vector},
                    {"nucleontotxsec_Flux",   &reco_sliceMCTruthFlux_weight_nucleontotxsec_vector},
                    {"kplus_Flux",            &reco_sliceMCTruthFlux_weight_kplus_vector},
                    {"kminus_Flux",           &reco_sliceMCTruthFlux_weight_kmin_vector},
                    {"kzero_Flux",            &reco_sliceMCTruthFlux_weight_kzero_vector},
                    {"piplus_Flux",           &reco_sliceMCTruthFlux_weight_piplus_vector},
                    {"piminus_Flux",          &reco_sliceMCTruthFlux_weight_piminus_vector} 
                };

                for(const auto &[name, values] : weightFluxMap){
                    //std::cout << "Map key = " << name << std::endl;
                    //std::cout << "Number of universes? = " << values.size() << std::endl;

                    //for(float v : values){
                        //std::cout << "  Weight = " << v << std::endl;
                    //}

                    auto mapFill = targetMap.find(name);
                    if(mapFill == targetMap.end()) continue;

                    auto &target = *(mapFill->second);

                    target.assign(values.begin(), values.end());
                
                }
                
                //std::cout << "============================" << std::endl;

                // Looking at GENIE systematics associated with MCTruth of slice
                
                art::FindManyP<std::map<std::string, std::vector<float>>> mcTruthGENIEWeightMapAssns(mcTruthVec, e, genieWeightLabel);
                const std::vector<art::Ptr<std::map<std::string, std::vector<float>>>> mcTruthGENIEWeightMaps(mcTruthGENIEWeightMapAssns.at(sliceMCTruth.key()));
  
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE_vector;
                std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE_vector;
                
                std::unordered_map<std::string, std::vector<float>*> targetGENIEMap = {
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi_vector},
                    {"NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi", &reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi_vector},
                    {"MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu_vector},
                    {"MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar_vector},
                    {"MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression_vector},
                    {"MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio_vector},
                    {"MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio", &reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio_vector},
                    {"MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement", &reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement_vector},
                    {"MINERvAE2p2h_SBN_v1_E2p2h_A_nu", &reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu_vector},
                    {"MINERvAE2p2h_SBN_v1_E2p2h_A_nubar", &reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar_vector},
                    {"MINERvAE2p2h_SBN_v1_E2p2h_B_nu", &reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu_vector},
                    {"MINERvAE2p2h_SBN_v1_E2p2h_B_nubar", &reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar_vector},
                    {"GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse_vector},
                    {"GENIEReWeight_SBN_v1_multisim_COHVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse_vector},
                    {"GENIEReWeight_SBN_v1_multisim_CoulombCCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE_vector},
                    {"GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse_vector},
                    {"GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse_vector},
                    {"GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NCELVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NormCCMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC_vector},
                    {"GENIEReWeight_SBN_v1_multisim_NormNCMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC_vector},
                    {"GENIEReWeight_SBN_v1_multisim_RDecBR1eta", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta_vector},
                    {"GENIEReWeight_SBN_v1_multisim_RDecBR1gamma", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma_vector},
                    {"GENIEReWeight_SBN_v1_multisim_RPA_CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE_vector},
                    {"GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_AhtBY", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_BhtBY", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_CV1uBY", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_CV2uBY", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_CoulombCCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_DecayAngMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_EtaNCEL", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_FrAbs_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_FrAbs_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_FrCEx_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_FrCEx_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_FrInel_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_FrInel_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_FrPiProd_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_MFP_N", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_MFP_pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_MaCCRES", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_MaNCEL", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_MaNCRES", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_MvCCRES", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_MvNCRES", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NormCCCOH", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NormCCMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NormNCCOH", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_NormNCMEC", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_RDecBR1eta", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_RPA_CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE_vector},
                    {"GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE", &reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE_vector},
                };

                for(const auto &weightGENIEMapPtr : mcTruthGENIEWeightMaps){
                    const std::map<std::string, std::vector<float>> &weightGENIEMap = *weightGENIEMapPtr;
                    for(const auto &[name, values] : weightGENIEMap){
                        auto mapFill = targetGENIEMap.find(name);
                        if(mapFill == targetGENIEMap.end()) continue;
                        auto &target = *(mapFill->second);
                        target.assign(values.begin(), values.end());
                    }
                }


                /*
                //std::cout << "Size of mcTruthGENIEWeightMaps vector = " << mcTruthGENIEWeightMaps.size() << std::endl;
                for(size_t i = 0; i < mcTruthGENIEWeightMaps.size(); ++i){
                    //std::cout << "Vector " << i << ":" << std::endl;
                    const art::Ptr<std::map<std::string, std::vector<float>>> &weightGENIEMapPtr_3Jul = mcTruthGENIEWeightMaps.at(i);
                    const std::map<std::string, std::vector<float>> &weightGENIEMap_3Jul = *weightGENIEMapPtr_3Jul;
                    //std::cout << "  Number of map entries = " << weightGENIEMap_3Jul.size() << std::endl;
                    for(const auto &[name_3Jul, values_3Jul] : weightGENIEMap_3Jul){
                        //std::cout << "  Map key = " << name_3Jul << std::endl;
                        //std::cout << "  Number of universes = " << values_3Jul.size() << std::endl;
                        for(size_t j = 0; j < values_3Jul.size(); j++){
                            //std::cout << "      Weight " << j+1 << ": " << values_3Jul.at(j) << std::endl;
                        }
                    }
                }
                */

                counter++;
                reco_sliceID.push_back(sliceID);
                reco_sliceMCTruthFlux_weight_horncurrent.push_back(reco_sliceMCTruthFlux_weight_horncurrent_vector);
                reco_sliceMCTruthFlux_weight_expskin.push_back(reco_sliceMCTruthFlux_weight_expskin_vector);
                reco_sliceMCTruthFlux_weight_pioninexsec.push_back(reco_sliceMCTruthFlux_weight_pioninexsec_vector);
                reco_sliceMCTruthFlux_weight_pionqexsec.push_back(reco_sliceMCTruthFlux_weight_pionqexsec_vector);
                reco_sliceMCTruthFlux_weight_piontotxsec.push_back(reco_sliceMCTruthFlux_weight_piontotxsec_vector);
                reco_sliceMCTruthFlux_weight_nucleoninexsec.push_back(reco_sliceMCTruthFlux_weight_nucleoninexsec_vector);
                reco_sliceMCTruthFlux_weight_nucleonqexsec.push_back(reco_sliceMCTruthFlux_weight_nucleonqexsec_vector);
                reco_sliceMCTruthFlux_weight_nucleontotxsec.push_back(reco_sliceMCTruthFlux_weight_nucleontotxsec_vector);
                reco_sliceMCTruthFlux_weight_kplus.push_back(reco_sliceMCTruthFlux_weight_kplus_vector);
                reco_sliceMCTruthFlux_weight_kmin.push_back(reco_sliceMCTruthFlux_weight_kmin_vector);
                reco_sliceMCTruthFlux_weight_kzero.push_back(reco_sliceMCTruthFlux_weight_kzero_vector);
                reco_sliceMCTruthFlux_weight_piplus.push_back(reco_sliceMCTruthFlux_weight_piplus_vector);
                reco_sliceMCTruthFlux_weight_piminus.push_back(reco_sliceMCTruthFlux_weight_piminus_vector);

                reco_sliceMCTruthFlux_weight_horncurrent_vector.clear();
                reco_sliceMCTruthFlux_weight_expskin_vector.clear();
                reco_sliceMCTruthFlux_weight_pioninexsec_vector.clear();
                reco_sliceMCTruthFlux_weight_pionqexsec_vector.clear();
                reco_sliceMCTruthFlux_weight_piontotxsec_vector.clear();
                reco_sliceMCTruthFlux_weight_nucleoninexsec_vector.clear();
                reco_sliceMCTruthFlux_weight_nucleonqexsec_vector.clear();
                reco_sliceMCTruthFlux_weight_nucleontotxsec_vector.clear();
                reco_sliceMCTruthFlux_weight_kplus_vector.clear();
                reco_sliceMCTruthFlux_weight_kmin_vector.clear();
                reco_sliceMCTruthFlux_weight_kzero_vector.clear();
                reco_sliceMCTruthFlux_weight_piplus_vector.clear();
                reco_sliceMCTruthFlux_weight_piminus_vector.clear();

                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi_vector);
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi_vector);
                reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu.push_back(reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu_vector);
                reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar.push_back(reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar_vector);
                reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression.push_back(reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression_vector);
                reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio.push_back(reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio_vector);
                reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio.push_back(reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio_vector);
                reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement.push_back(reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement_vector);
                reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu.push_back(reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu_vector);
                reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar.push_back(reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar_vector);
                reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu.push_back(reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu_vector);
                reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar.push_back(reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE_vector);
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE_vector);

                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu_vector.clear();
                reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar_vector.clear();
                reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression_vector.clear();
                reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio_vector.clear();
                reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio_vector.clear();
                reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement_vector.clear();
                reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu_vector.clear();
                reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar_vector.clear();
                reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu_vector.clear();
                reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE_vector.clear();
                reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE_vector.clear();


                if(sliceMCTruth->Origin() == simb::kBeamNeutrino){
                    reco_sliceInteraction.push_back(sliceMCNeutrino.InteractionType());
                    reco_sliceTrueVX.push_back(sliceMCNeutrinoParticle.Vx());
                    reco_sliceTrueVY.push_back(sliceMCNeutrinoParticle.Vy());
                    reco_sliceTrueVZ.push_back(sliceMCNeutrinoParticle.Vz());
                    reco_sliceTrueCCNC.push_back(sliceMCNeutrino.CCNC());
                    reco_sliceTrueNeutrinoType.push_back(sliceMCNeutrinoParticle.PdgCode());
                    if(sliceMCNeutrino.InteractionType() == 1098){
                        reco_sliceOrigin.push_back(1);
                    } else{
                        reco_sliceOrigin.push_back(3);
                    }
                } else if(sliceMCTruth->Origin() == simb::kCosmicRay || sliceMCTruth->Origin() == 0){
                    reco_sliceInteraction.push_back(-100);        
                    reco_sliceTrueVX.push_back(-999999);
                    reco_sliceTrueVY.push_back(-999999);
                    reco_sliceTrueVZ.push_back(-999999);
                    reco_sliceOrigin.push_back(0);
                    reco_sliceTrueCCNC.push_back(-999999);
                    reco_sliceTrueNeutrinoType.push_back(-999999);
                } else {
                    reco_sliceInteraction.push_back(-999999);        
                    reco_sliceTrueVX.push_back(-999999);
                    reco_sliceTrueVY.push_back(-999999);
                    reco_sliceTrueVZ.push_back(-999999);
                    reco_sliceOrigin.push_back(-999999);
                    reco_sliceTrueCCNC.push_back(-999999);
                    reco_sliceTrueNeutrinoType.push_back(-999999);
                }
            }
        }
    }

    if(counter == 0){
        reco_sliceID.push_back(-999999);        
        reco_sliceInteraction.push_back(-999999);       
        reco_sliceTrueVX.push_back(-999999); 
        reco_sliceTrueVY.push_back(-999999); 
        reco_sliceTrueVZ.push_back(-999999); 
        reco_sliceOrigin.push_back(-999999);
        reco_sliceTrueCCNC.push_back(-999999);
        reco_sliceTrueNeutrinoType.push_back(-999999);
       
        std::vector<float> reco_sliceMCTruthFlux_weight_horncurrent_vector;
        std::vector<float> reco_sliceMCTruthFlux_weight_expskin_vector;
        std::vector<float> reco_sliceMCTruthFlux_weight_pioninexsec_vector;
        std::vector<float> reco_sliceMCTruthFlux_weight_pionqexsec_vector;
        std::vector<float> reco_sliceMCTruthFlux_weight_piontotxsec_vector;
        std::vector<float> reco_sliceMCTruthFlux_weight_nucleoninexsec_vector;
        std::vector<float> reco_sliceMCTruthFlux_weight_nucleonqexsec_vector;
        std::vector<float> reco_sliceMCTruthFlux_weight_nucleontotxsec_vector;
        std::vector<float> reco_sliceMCTruthFlux_weight_kplus_vector;
        std::vector<float> reco_sliceMCTruthFlux_weight_kmin_vector;
        std::vector<float> reco_sliceMCTruthFlux_weight_kzero_vector;
        std::vector<float> reco_sliceMCTruthFlux_weight_piplus_vector;
        std::vector<float> reco_sliceMCTruthFlux_weight_piminus_vector;

        reco_sliceMCTruthFlux_weight_horncurrent_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_expskin_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_pioninexsec_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_pionqexsec_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_piontotxsec_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_nucleoninexsec_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_nucleonqexsec_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_nucleontotxsec_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_kplus_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_kmin_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_kzero_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_piplus_vector.push_back(-999999);
        reco_sliceMCTruthFlux_weight_piminus_vector.push_back(-999999);

        reco_sliceMCTruthFlux_weight_horncurrent.push_back(reco_sliceMCTruthFlux_weight_horncurrent_vector);
        reco_sliceMCTruthFlux_weight_expskin.push_back(reco_sliceMCTruthFlux_weight_expskin_vector);
        reco_sliceMCTruthFlux_weight_pioninexsec.push_back(reco_sliceMCTruthFlux_weight_pioninexsec_vector);
        reco_sliceMCTruthFlux_weight_pionqexsec.push_back(reco_sliceMCTruthFlux_weight_pionqexsec_vector);
        reco_sliceMCTruthFlux_weight_piontotxsec.push_back(reco_sliceMCTruthFlux_weight_piontotxsec_vector);
        reco_sliceMCTruthFlux_weight_nucleoninexsec.push_back(reco_sliceMCTruthFlux_weight_nucleoninexsec_vector);
        reco_sliceMCTruthFlux_weight_nucleonqexsec.push_back(reco_sliceMCTruthFlux_weight_nucleonqexsec_vector);
        reco_sliceMCTruthFlux_weight_nucleontotxsec.push_back(reco_sliceMCTruthFlux_weight_nucleontotxsec_vector);
        reco_sliceMCTruthFlux_weight_kplus.push_back(reco_sliceMCTruthFlux_weight_kplus_vector);
        reco_sliceMCTruthFlux_weight_kmin.push_back(reco_sliceMCTruthFlux_weight_kmin_vector);
        reco_sliceMCTruthFlux_weight_kzero.push_back(reco_sliceMCTruthFlux_weight_kzero_vector);
        reco_sliceMCTruthFlux_weight_piplus.push_back(reco_sliceMCTruthFlux_weight_piplus_vector);
        reco_sliceMCTruthFlux_weight_piminus.push_back(reco_sliceMCTruthFlux_weight_piminus_vector);
  
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE_vector;
        std::vector<float> reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE_vector;
  
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE_vector.push_back(-999999);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE_vector.push_back(-999999);

        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi_vector);
        reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi.push_back(reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi_vector);
        reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu.push_back(reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu_vector);
        reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar.push_back(reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar_vector);
        reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression.push_back(reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression_vector);
        reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio.push_back(reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio_vector);
        reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio.push_back(reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio_vector);
        reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement.push_back(reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement_vector);
        reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu.push_back(reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu_vector);
        reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar.push_back(reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar_vector);
        reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu.push_back(reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu_vector);
        reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar.push_back(reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE_vector);
        reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE.push_back(reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE_vector);

    }

    //std::cout << "_________________________" << std::endl;
}

void sbnd::NuEWithSysts::clearVectors(){

    nuEScatter_MCTruthFlux_weight_horncurrent.clear();
    nuEScatter_MCTruthFlux_weight_expskin.clear();
    nuEScatter_MCTruthFlux_weight_pioninexsec.clear();
    nuEScatter_MCTruthFlux_weight_pionqexsec.clear();
    nuEScatter_MCTruthFlux_weight_piontotxsec.clear();
    nuEScatter_MCTruthFlux_weight_nucleoninexsec.clear();
    nuEScatter_MCTruthFlux_weight_nucleonqexsec.clear();
    nuEScatter_MCTruthFlux_weight_nucleontotxsec.clear();
    nuEScatter_MCTruthFlux_weight_kplus.clear();
    nuEScatter_MCTruthFlux_weight_kmin.clear();
    nuEScatter_MCTruthFlux_weight_kzero.clear();
    nuEScatter_MCTruthFlux_weight_piplus.clear();
    nuEScatter_MCTruthFlux_weight_piminus.clear();
  
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi.clear();
    nuEScatter_MCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi.clear();
    nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu.clear();
    nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar.clear();
    nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression.clear();
    nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio.clear();
    nuEScatter_MCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio.clear();
    nuEScatter_MCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement.clear();
    nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu.clear();
    nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar.clear();
    nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu.clear();
    nuEScatter_MCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE.clear();
    nuEScatter_MCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE.clear();

    reco_sliceID.clear();
    reco_sliceInteraction.clear();
    reco_sliceTrueVX.clear();
    reco_sliceTrueVY.clear();
    reco_sliceTrueVZ.clear();
    reco_sliceOrigin.clear();
    reco_sliceTrueCCNC.clear();
    reco_sliceTrueNeutrinoType.clear();
   
    reco_sliceMCTruthFlux_weight_horncurrent.clear();
    reco_sliceMCTruthFlux_weight_expskin.clear();
    reco_sliceMCTruthFlux_weight_pioninexsec.clear();
    reco_sliceMCTruthFlux_weight_pionqexsec.clear();
    reco_sliceMCTruthFlux_weight_piontotxsec.clear();
    reco_sliceMCTruthFlux_weight_nucleoninexsec.clear();
    reco_sliceMCTruthFlux_weight_nucleonqexsec.clear();
    reco_sliceMCTruthFlux_weight_nucleontotxsec.clear();
    reco_sliceMCTruthFlux_weight_kplus.clear();
    reco_sliceMCTruthFlux_weight_kmin.clear();
    reco_sliceMCTruthFlux_weight_kzero.clear();
    reco_sliceMCTruthFlux_weight_piplus.clear();
    reco_sliceMCTruthFlux_weight_piminus.clear();
  
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_2Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_CC_3Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_1Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_2Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_n_NC_3Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_np_CC_1Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_2Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_CC_3Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_1Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_2Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nu_p_NC_3Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_1Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_2Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_CC_3Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_1Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_2Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_n_NC_3Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_1Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_2Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_CC_3Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_1Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_2Pi.clear();
    reco_sliceMCTruthGENIE_weight_NOvAStyleNonResPionNorm_SBN_v1_NR_nubar_p_NC_3Pi.clear();
    reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nu.clear();
    reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_C12ToAr40_2p2hScaling_nubar.clear();
    reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_SPPLowQ2Suppression.clear();
    reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenuebar_xsec_ratio.clear();
    reco_sliceMCTruthGENIE_weight_MiscInteractionSysts_SBN_v1_nuenumu_xsec_ratio.clear();
    reco_sliceMCTruthGENIE_weight_MINERvAq0q3Weighting_SBN_v1_Mnv2p2hGaussEnhancement.clear();
    reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nu.clear();
    reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_A_nubar.clear();
    reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nu.clear();
    reco_sliceMCTruthGENIE_weight_MINERvAE2p2h_SBN_v1_E2p2h_B_nubar.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_COHVariationResponse.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_CoulombCCQE.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_DISBYVariationResponse.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_N_VariationResponse.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_FSI_pi_VariationResponse.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCELVariationResponse.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NCRESVariationResponse.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC1pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnCC2pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC1pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarnNC2pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC1pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpCC2pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC1pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvbarpNC2pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC1pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnCC2pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC1pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvnNC2pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC1pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpCC2pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC1pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NonRESBGvpNC2pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormCCMEC.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_NormNCMEC.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1eta.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RDecBR1gamma.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_RPA_CCQE.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisim_ZExpAVariationResponse.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_AhtBY.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_BhtBY.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV1uBY.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CV2uBY.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_CoulombCCQE.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_DecayAngMEC.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_EtaNCEL.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_N.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrAbs_pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_N.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrCEx_pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_N.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrInel_pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_N.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_FrPiProd_pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_N.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MFP_pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaCCRES.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCEL.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MaNCRES.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvCCRES.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_MvNCRES.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC1pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnCC2pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC1pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarnNC2pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC1pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpCC2pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC1pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvbarpNC2pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC1pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnCC2pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC1pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvnNC2pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC1pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpCC2pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC1pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NonRESBGvpNC2pi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCCOH.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormCCMEC.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCCOH.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_NormNCMEC.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1eta.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RDecBR1gamma.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_RPA_CCQE.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ThetaDelta2NRad.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_Theta_Delta2Npi.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_VecFFCCQEshape.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA1CCQE.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA2CCQE.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA3CCQE.clear();
    reco_sliceMCTruthGENIE_weight_GENIEReWeight_SBN_v1_multisigma_ZExpA4CCQE.clear();
}

void sbnd::NuEWithSysts::beginJob()
{
  // Implementation of optional member function here.
}

void sbnd::NuEWithSysts::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(sbnd::NuEWithSysts)

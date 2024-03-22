#ifndef SBND_NUANA_ALG_H
#define SBND_NUANA_ALG_H


#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/Geometry.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nugen/NuReweight/art/NuReweight.h" //GENIEReweight.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom.h"

#include <memory>

namespace sbnd{

  enum reweight {kNCEL, kQEMA, kQEVec, kResGanged, kCCRes, kNCRes,
        kCoh, kNonResRvp1pi, kNonResRvbarp1pi, kNonResRvp2pi,
        kNonResRvbarp2pi, kResDecay, kNC, kDIS, kDISnucl,
        kAGKY, kNReWeights};

  class NuAnaAlg
  {
  public:

    NuAnaAlg();
    // ~NuAnaAlg();

    void configureGeometry(art::ServiceHandle<geo::Geometry> );

    void configureReWeight(const std::vector<reweight> &,
                           const std::vector<std::vector<float>>&);
                           // std::vector<float>&,
                           // int);

    void calcWeight(        simb::MCTruth const&,
                            simb::GTruth const&,
                            std::vector<std::vector<float>>& );

    void packFluxWeight(    simb::MCFlux const& flux,
                            std::vector<std::vector<float>>&);

    unsigned int prepareSigmas(int, unsigned int,
                            std::vector<std::vector<float> > & );

    void parseWeights(const std::vector<std::string> &,
                            std::vector<reweight> &);

    // get the basic neutrino info:
    void packNeutrinoInfo(simb::MCNeutrino * neutrino,
                            int& nuchan,
                            int& inno,
                            double& enugen,
                            int& isCC,
                            int& mode,
                            double& thetaLep,
                            double& phiLep,
                            double& Elep,
                            std::vector<float>& neutMom,
                            std::vector<float>& vertex);

    // Pack up the flux info:
    void packFluxInfo(      simb::MCFlux const& flux,
                            int& ptype, int& tptype, int& ndecay,
                            std::vector<float>& neutVertexInWindow,
                            std::vector<float>& ParentVertex,
                            std::vector<float>& nuParentMomAtDecay,
                            std::vector<float>& nuParentMomAtProd,
                            std::vector<float>& nuParentMomTargetExit);

    // Pack up the genie info:
    void packGenieInfo(     simb::MCTruth const& truth,
                            std::vector<int> & GeniePDG,
                            std::vector<std::vector<float>>& GenieMomentum,
                            std::vector<std::string>& GenieProc,
                            int& NPi0FinalState,
                            int& NGamma,
                            int& NChargedPions);

    void packLarg4Info(std::vector<simb::MCParticle> const& mclarg4, int, int, int, int,
                            std::vector<std::vector<float> > & leptonPos,
                            std::vector<std::vector<float> > & leptonMom,
                            std::vector<std::vector<float> > & p1PhotonConversionPos,
                            std::vector<std::vector<float> > & p1PhotonConversionMom,
                            std::vector<std::vector<float> > & p2PhotonConversionPos,
                            std::vector<std::vector<float> > & p2PhotonConversionMom,
                            std::vector<std::vector<float> > & miscPhotonConversionPos,
                            std::vector<std::vector<float> > & miscPhotonConversionMom,
                            std::vector<std::vector<float> > & pionPos,
                            std::vector<std::vector<float> > & pionMom,
                            std::vector<std::vector<std::vector<float> > > & chargedPionPos,
                            std::vector<std::vector<std::vector<float> > > & chargedPionMom,
                            std::vector<int> & chargePionSign);

  private:

    simb::MCParticle const* getParticleByID(
            std::vector<simb::MCParticle> const& mclistLARG4, int ) const;
    simb::MCParticle const& getParticleByIDStrict(
            std::vector<simb::MCParticle> const& mclistLARG4, int ) const;


    void pack4Vector(const TLorentzVector& input, std::vector<float>& output) const;
    void pack3Vector(const TVector3& input, std::vector<float>& output) const;

    bool isInTPC(const TVector3 &) const;
    void GetPhotonConversionInfo( simb::MCParticle const& photon,
                                  TLorentzVector& ConversionPos,
                                  TLorentzVector& ConversionMom);
    // The reweighting utility class:
    // std::unique_ptr<rwgt::NuReweight> reweight;
    std::vector<std::vector<rwgt::NuReweight *> > reweightVector;

    // geometry boundaries:
    double xlow, xhigh, ylow, yhigh, zlow, zhigh;


    /* data */
  };
}

#endif

#ifndef LAR1ND_NUANA_ALG_H
#define LAR1ND_NUANA_ALG_H value

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "Geometry/Geometry.h"

#include "SimulationBase/MCParticle.h"
#include "NuReweight/art/NuReweight.h" //GENIEReweight.h"
#include "SimulationBase/MCNeutrino.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/GTruth.h"
#include "SimulationBase/MCFlux.h"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom.h"

#include <memory>

namespace lar1nd{

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

    void calcWeight(art::Ptr<simb::MCTruth>,
                            art::Ptr<simb::GTruth >,
                            std::vector<std::vector<float>>& );

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
    void packFluxInfo(art::Ptr<simb::MCFlux >  flux, 
                            int& ptype, int& tptype, int& ndecay,
                            std::vector<float>& neutVertexInWindow,
                            std::vector<float>& ParentVertex,
                            std::vector<float>& nuParentMomAtDecay,
                            std::vector<float>& nuParentMomAtProd,
                            std::vector<float>& nuParentMomTargetExit);

    // Pack up the genie info:
    void packGenieInfo(art::Ptr<simb::MCTruth>  truth,
                            std::vector<int> & GeniePDG,
                            std::vector<std::vector<float>>& GenieMomentum,
                            std::vector<std::string>& GenieProc,
                            int& NPi0FinalState,
                            int& NGamma,
                            int& NChargedPions);

    void packLarg4Info(art::Handle< std::vector<simb::MCParticle> > mclarg4, int, int, int, int,
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

    art::Ptr<simb::MCParticle> getParticleByID(
            art::Handle< std::vector<simb::MCParticle> > & mclistLARG4, int ) const;
    art::Ptr<simb::MCParticle> getParticleByPDG(
            art::Handle< std::vector<simb::MCParticle> > & mclistLARG4, int ) const;


    void pack4Vector(const TLorentzVector& input, std::vector<float>& output) const;
    void pack3Vector(const TVector3& input, std::vector<float>& output) const;

    bool isInTPC(TVector3 &) const;
    void GetPhotonConversionInfo( art::Ptr<simb::MCParticle> photon,
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

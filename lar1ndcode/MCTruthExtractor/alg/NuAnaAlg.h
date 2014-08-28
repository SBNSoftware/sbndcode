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

#include <memory>

namespace lar1nd{

  class NuAnaAlg
  {
  public:

    NuAnaAlg();
    // ~NuAnaAlg();
  
    void configureGeometry(art::ServiceHandle<geo::Geometry> );

    void configureReWeight();

    double calcWeight(art::Ptr<simb::MCTruth>  truth, art::Ptr<simb::GTruth >);

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
                          TLorentzVector& neutMom,
                          TVector3& vertex);

    // Pack up the flux info:
    void packFluxInfo(art::Ptr<simb::MCFlux >  flux, 
                            int& ptype, int& tptype, int& ndecay,
                            TVector3& ParentVertex,
                            TVector3& nuParentMomAtDecay,
                            TVector3& nuParentMomAtProd,
                            TVector3& nuParentMomTargetExit);

    // Pack up the genie info:
    void packGenieInfo(art::Ptr<simb::MCTruth>  truth,
                            std::vector<int> & GeniePDG,
                            std::vector<TLorentzVector>& GenieMomentum,
                            std::vector<std::string>& GenieProc,
                            int& NPi0FinalState,
                            int& NGamma,
                            int& NChargedPions);

    void packLarg4Info(art::Handle< std::vector<simb::MCParticle> > mclarg4, int, int, int,
                            std::vector<TLorentzVector> & leptonPos,
                            std::vector<TLorentzVector> & leptonMom,
                            std::vector<TLorentzVector> & p1PhotonConversionPos,
                            std::vector<TLorentzVector> & p1PhotonConversionMom,
                            std::vector<TLorentzVector> & p2PhotonConversionPos,
                            std::vector<TLorentzVector> & p2PhotonConversionMom,
                            std::vector<TLorentzVector> & miscPhotonConversionPos,
                            std::vector<TLorentzVector> & miscPhotonConversionMom,
                            std::vector<TLorentzVector> & pionPos,
                            std::vector<TLorentzVector> & pionMom,
                            std::vector<std::vector<TLorentzVector> > & chargedPionPos,
                            std::vector<std::vector<TLorentzVector> > & chargedPionMom,
                            std::vector<int> & chargePionSign);

  private:

    art::Ptr<simb::MCParticle> getParticleByID(
            art::Handle< std::vector<simb::MCParticle> > & mclistLARG4, int ) const;
    art::Ptr<simb::MCParticle> getParticleByPDG(
            art::Handle< std::vector<simb::MCParticle> > & mclistLARG4, int ) const;
    bool isInTPC(TVector3 &) const;
    void GetPhotonConversionInfo( art::Ptr<simb::MCParticle> photon,
                                  TLorentzVector& ConversionPos,
                                  TLorentzVector& ConversionMom);
    // The reweighting utility class:
    // std::unique_ptr<rwgt::NuReweight> reweight;
    rwgt::NuReweight * reweight;

    // geometry boundaries:
    double xlow, xhigh, ylow, yhigh, zlow, zhigh;


    /* data */
  };
}

#endif

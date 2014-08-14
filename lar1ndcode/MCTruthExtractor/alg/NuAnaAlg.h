#ifndef LAR1ND_NUANA_ALG_H
#define LAR1ND_NUANA_ALG_H value

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "Geometry/Geometry.h"

#include "SimulationBase/MCParticle.h"
#include "NuReweight/art/NuReweight.h" //GENIEReweight.h"
#include "SimulationBase/MCNeutrino.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCFlux.h"

#include "TVector3.h"
#include "TLorentzVector.h"


namespace lar1nd{

  class NuAnaAlg
  {
  public:

    NuAnaAlg();
    ~NuAnaAlg(){};
  
    void configureGeometry(art::ServiceHandle<geo::Geometry> );



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
                            int& NGamma);



  private:
    simb::MCParticle getParticleByID(
            art::Handle< std::vector<simb::MCParticle> > & mclistLARG4, int ) const;
    simb::MCParticle getParticleByPDG(
            art::Handle< std::vector<simb::MCParticle> > & mclistLARG4, int ) const;
    bool isInTPC(TVector3 &) const;
    void GetPhotonConversionInfo( simb::MCParticle& photon,
                                  TLorentzVector& ConversionPos,
                                  TLorentzVector& ConversionMom);
    // The reweighting utility class:
    // rwgt::NuReweight reweight;

    // geometry boundaries:
    double xlow, xhigh, ylow, yhigh, zlow, zhigh;


    /* data */
  };
}

#endif

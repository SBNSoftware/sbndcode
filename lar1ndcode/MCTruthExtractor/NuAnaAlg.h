#ifndef LAR1ND_NUANA_ALG_H
#define LAR1ND_NUANA_ALG_H value

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "Geometry/Geometry.h"


#include "SimulationBase/MCParticle.h"



namespace lar1nd{

  class NuAnaAlg
  {
  public:

    NuAnaAlg();
    ~NuAnaAlg();
  
    void configureGeometry(art::ServiceHandle<geo::Geometry> );
    void GetPhotonConversionInfo( simb::MCParticle& photon,
                                  TLorentzVector& ConversionPos,
                                  TLorentzVector& ConversionMom);

  private:
    simb::MCParticle getParticleByID(
            art::Handle< std::vector<simb::MCParticle> > & mclistLARG4, int );
    simb::MCParticle getParticleByPDG(
            art::Handle< std::vector<simb::MCParticle> > & mclistLARG4, int );
    bool isInTPC(TVector3) const;


    // geometry boundaries:
    double xlow, xhigh, ylow, yhigh, zlow, zhigh;


    /* data */
  };
}

#endif

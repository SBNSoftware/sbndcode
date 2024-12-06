#include <algorithm>

#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 

#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom()
#include "larcorealg/Geometry/GeometryCore.h"
#include "nusimdata/SimulationBase/MCParticle.h"

namespace filt{

  class LArG4FakeTriggerFilter : public art::EDFilter {
    public:

      explicit LArG4FakeTriggerFilter(fhicl::ParameterSet const & pset);
      virtual bool filter(art::Event& e) override;
      void reconfigure(fhicl::ParameterSet const& pset);
      virtual void beginJob() override;

    private:

      double fBeamTimeMin;  // Minimum time of beam window [us]
      double fBeamTimeMax;  // Maximum time of beam window [us]
      double fEnergyDeposit;// Minimum energy deposit in TPC for trigger [MeV]

      std::string fLArG4ModuleName;

      geo::GeometryCore const* fGeometryService;

      bool IsInterestingParticle(const simb::MCParticle& particle);
      double EnergyInTPC(const simb::MCParticle& particle);
  };


  LArG4FakeTriggerFilter::LArG4FakeTriggerFilter(fhicl::ParameterSet const & pset)
  : EDFilter(pset)
  {
    this->reconfigure(pset);

    fGeometryService = lar::providerFrom<geo::Geometry>();
  }


  void LArG4FakeTriggerFilter::reconfigure(fhicl::ParameterSet const& pset){

    fBeamTimeMin   = pset.get<double>("BeamTimeMin");
    fBeamTimeMax   = pset.get<double>("BeamTimeMax");
    fEnergyDeposit = pset.get<double>("EnergyDeposit");

    fLArG4ModuleName = pset.get<std::string>("LArG4ModuleName");

  }


  bool LArG4FakeTriggerFilter::filter(art::Event & e){

    auto const& particles = e.getProduct<std::vector<simb::MCParticle>>(fLArG4ModuleName);

    double e_dep = 0;

    for (simb::MCParticle const& particle : particles) {

      // Check particle time is within the beam time
      double time = particle.T() * 1e-3; // [us]
      if(time < fBeamTimeMin || time > fBeamTimeMax) continue;

      // Check particle is stable and charged
      if (!IsInterestingParticle(particle)) continue;
      
      // Add up the energy deposit inside the TPC
      e_dep += EnergyInTPC(particle); // [GeV]
    }

    // If the energy deposit within the beam time is greater than some limit then trigger the event
    return e_dep > fEnergyDeposit/1000.;

  }


  void LArG4FakeTriggerFilter::beginJob() {
  }


  // Check if particle is stable and charged
  bool LArG4FakeTriggerFilter::IsInterestingParticle(const simb::MCParticle& particle){

    // Stable final state
    if(particle.StatusCode() != 1) return false;

    // Charged: e, mu, pi, K, p
    int pdg = particle.PdgCode();
    if(!(pdg==11 || pdg==13 || pdg==211 || pdg==321 || pdg==2212)) return false;

    return true;

  }

  // Add up energy deposit in TPC
  double LArG4FakeTriggerFilter::EnergyInTPC(const simb::MCParticle& particle){

    geo::TPCGeo const& tpcGeo = fGeometryService->TPC({0, 0});
    double xmin = -2.0 * tpcGeo.HalfWidth();
    double xmax = 2.0 * tpcGeo.HalfWidth();
    double ymin = -tpcGeo.HalfHeight();
    double ymax = tpcGeo.HalfHeight();
    double zmin = 0.;
    double zmax = tpcGeo.Length();

    double e_dep = 0;

    int npts = particle.NumberTrajectoryPoints();
    for (int i = 1; i < npts; i++){
      TVector3 pos(particle.Vx(i), particle.Vy(i), particle.Vz(i));
      // Check if point is within the TPC
      if (pos[0] >= xmin && pos[0] <= xmax && pos[1] >= ymin && pos[1] <= ymax && pos[2] >= zmin && pos[2] <= zmax){
        e_dep += particle.E(i-1) - particle.E(i);
      }
    }

    return e_dep;
  }


  DEFINE_ART_MODULE(LArG4FakeTriggerFilter)

}

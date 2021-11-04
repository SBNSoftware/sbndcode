#include <algorithm>

#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 

#include "larcore/Geometry/Geometry.h"
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

      bool IsInterestingParticle(const art::Ptr<simb::MCParticle> particle);
      double EnergyInTPC(const art::Ptr<simb::MCParticle> particle);
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

    art::Handle<std::vector<simb::MCParticle> > particles;
    e.getByLabel(fLArG4ModuleName, particles);

    double e_dep = 0;

    for (unsigned int part_i = 0; part_i < particles->size(); part_i++){

      const art::Ptr<simb::MCParticle> particle(particles, part_i);

      // Check particle time is within the beam time
      double time = particle->T() * 1e-3; // [us]
      if(time < fBeamTimeMin || time > fBeamTimeMax) continue;

      // Check particle is stable and charged
      if (!IsInterestingParticle(particle)) continue;
      
      // Add up the energy deposit inside the TPC
      e_dep += EnergyInTPC(particle); // [GeV]
    }

    // If the energy deposit within the beam time is greater than some limit then trigger the event
    if(e_dep > fEnergyDeposit/1000.) return true;
    return false;

  }


  void LArG4FakeTriggerFilter::beginJob() {
  }


  // Check if particle is stable and charged
  bool LArG4FakeTriggerFilter::IsInterestingParticle(const art::Ptr<simb::MCParticle> particle){

    // Stable final state
    if(particle->StatusCode() != 1) return false;

    // Charged: e, mu, pi, K, p
    int pdg = particle->PdgCode();
    if(!(pdg==11 || pdg==13 || pdg==211 || pdg==321 || pdg==2212)) return false;

    return true;

  }

  // Add up energy deposit in TPC
  double LArG4FakeTriggerFilter::EnergyInTPC(const art::Ptr<simb::MCParticle> particle){

    double xmin = -2.0 * fGeometryService->DetHalfWidth();
    double xmax = 2.0 * fGeometryService->DetHalfWidth();
    double ymin = -fGeometryService->DetHalfHeight();
    double ymax = fGeometryService->DetHalfHeight();
    double zmin = 0.;
    double zmax = fGeometryService->DetLength();

    double e_dep = 0;

    int npts = particle->NumberTrajectoryPoints();
    for (int i = 1; i < npts; i++){
      TVector3 pos(particle->Vx(i), particle->Vy(i), particle->Vz(i));
      // Check if point is within the TPC
      if (pos[0] >= xmin && pos[0] <= xmax && pos[1] >= ymin && pos[1] <= ymax && pos[2] >= zmin && pos[2] <= zmax){
        e_dep += particle->E(i-1) - particle->E(i);
      }
    }

    return e_dep;
  }


  DEFINE_ART_MODULE(LArG4FakeTriggerFilter)

}

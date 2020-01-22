#include <algorithm>

#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 

#include "nusimdata/SimulationBase/MCParticle.h"

#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

namespace filt{

  class LArG4DirtFilter : public art::EDFilter {
    public:

      explicit LArG4DirtFilter(fhicl::ParameterSet const & pset);
      virtual bool filter(art::Event& e) override;
      void reconfigure(fhicl::ParameterSet const& pset);
      virtual void beginJob() override;

    private:

      double fEnergyDeposit;// Minimum energy deposit in TPC for trigger [MeV]
      bool fRemoveTpcOnly;
      double fTpcBuffer;

      std::string fLArG4ModuleName;

      sbnd::TPCGeoAlg fTpcGeo;

      bool IsInterestingParticle(const art::Ptr<simb::MCParticle> particle);
      double DaughterEDep(const art::Ptr<simb::MCParticle> particle, art::Handle<std::vector<simb::MCParticle> > particles);
      std::pair<bool, const art::Ptr<simb::MCParticle>> GetDaughter(int id, art::Handle<std::vector<simb::MCParticle> > particles);
  };


  LArG4DirtFilter::LArG4DirtFilter(fhicl::ParameterSet const & pset)
  : EDFilter(pset)
  {
    this->reconfigure(pset);

  }


  void LArG4DirtFilter::reconfigure(fhicl::ParameterSet const& pset){

    fEnergyDeposit = pset.get<double>("EnergyDeposit");
    fRemoveTpcOnly = pset.get<bool>("RemoveTpcOnly");
    fTpcBuffer = pset.get<double>("TpcBuffer");

    fLArG4ModuleName = pset.get<std::string>("LArG4ModuleName");

  }


  bool LArG4DirtFilter::filter(art::Event & event){

    art::Handle<std::vector<simb::MCParticle> > particles;
    event.getByLabel(fLArG4ModuleName, particles);

    double edep_tpc = 0;

    for (unsigned int part_i = 0; part_i < particles->size(); part_i++){
    
      const art::Ptr<simb::MCParticle> particle(particles, part_i);
      if(!IsInterestingParticle(particle)) continue;
      if(particle->Mother() != 0) continue;

      geo::Point_t vertex {particle->Vx(), 
                           particle->Vy(), 
                           particle->Vz()};

      if(fRemoveTpcOnly && fTpcGeo.InFiducial(vertex, -fTpcBuffer)) continue;

      edep_tpc += fTpcGeo.EDep(*particle);
      edep_tpc += DaughterEDep(particle, particles);

    }
    
    if(edep_tpc > fEnergyDeposit/1000.) return true;
    return false;

  }


  void LArG4DirtFilter::beginJob() {
  }


  // Check if particle is stable and charged
  bool LArG4DirtFilter::IsInterestingParticle(const art::Ptr<simb::MCParticle> particle){

    // Stable final state
    if(particle->StatusCode() != 1) return false;

    // Charged: e, mu, pi, K, p
    int pdg = particle->PdgCode();
    if(!(pdg==11 || pdg==13 || pdg==211 || pdg==321 || pdg==2212)) return false;

    return true;

  }
      
  double LArG4DirtFilter::DaughterEDep(const art::Ptr<simb::MCParticle> particle, art::Handle<std::vector<simb::MCParticle> > particles){

    if(particle->NumberDaughters() == 0) return 0;

    double edep = 0.;
    for(int j = 0; j < particle->NumberDaughters(); j++){
      int d_id = particle->Daughter(j);
      std::pair<bool, const art::Ptr<simb::MCParticle>> daughter = GetDaughter(d_id, particles);
      if(!daughter.first) continue;
      if(!IsInterestingParticle(daughter.second)) continue;
      edep += fTpcGeo.EDep(*daughter.second);
      edep += DaughterEDep(daughter.second, particles);
    }

    return edep;
  }

  std::pair<bool, const art::Ptr<simb::MCParticle>> LArG4DirtFilter::GetDaughter(int id, art::Handle<std::vector<simb::MCParticle> > particles){

    for (unsigned int part_i = 0; part_i < particles->size(); part_i++){
      const art::Ptr<simb::MCParticle> particle(particles, part_i);
      if(particle->TrackId() == id) return std::make_pair(true, particle);
    }
    const art::Ptr<simb::MCParticle> particle(particles, 0);
    return std::make_pair(false, particle);

  }

  DEFINE_ART_MODULE(LArG4DirtFilter)

}

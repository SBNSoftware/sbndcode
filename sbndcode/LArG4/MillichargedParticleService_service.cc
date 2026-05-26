#ifndef SBNDCODE_MILLICHARGEDPARTICLESERVICE_SERVICE_CC
#define SBNDCODE_MILLICHARGEDPARTICLESERVICE_SERVICE_CC

#include <fhiclcpp/ParameterSet.h>
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "TDatabasePDG.h"

struct MillichargedParticleService {
  double massFraction, charge;

  MillichargedParticleService(fhicl::ParameterSet const &p, art::ActivityRegistry &areg)
    : massFraction(p.get<double>("massFraction")),
      charge(p.get<double>("charge")) {
    std::cout << "CHARGE = " << charge << "!!!!!!" << std::endl;
    TDatabasePDG db;
    auto *muon = db.GetParticle(13);
    db.AddParticle(
		   "millicharged_name",
		   "millicharged_title",
		   muon->Mass() * massFraction,
		   muon->Stable(),
		   muon->Width(),
		   charge,
		   muon->ParticleClass(),
		   // 59 for "elaborate DM scenario", 52 for spin 0 DM                                                                    
		   5900052,
		   -1,
		   muon->TrackingCode()
		   );
  }
};

DECLARE_ART_SERVICE(MillichargedParticleService, LEGACY)
DEFINE_ART_SERVICE(MillichargedParticleService)


#endif //SBNDCODE_MILLICHARGEDPARTICLESERVICE_SERVICE_CC                                                                          




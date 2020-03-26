#ifndef SPACECHARGESERVICESBND_H
#define SPACECHARGESERVICESBND_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Run.h"
#include "sbndcode/SpaceCharge/SpaceChargeSBND.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"


namespace spacecharge
{
    class SpaceChargeServiceSBND : public SpaceChargeService
    {
    public:
	SpaceChargeServiceSBND(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

        void   reconfigure(fhicl::ParameterSet const& pset);
	void   preBeginRun(const art::Run& run);


	virtual const  provider_type* provider() const override
	{
	    return fProp.get();
	}

    private:

	std::unique_ptr<spacecharge::SpaceChargeSBND> fProp;

    }; // class SpaceChargeServiceSBND
} //namespace spacecharge
DECLARE_ART_SERVICE_INTERFACE_IMPL(spacecharge::SpaceChargeServiceSBND, spacecharge::SpaceChargeService, LEGACY)
#endif // SPACECHARGESERVICESBND_H

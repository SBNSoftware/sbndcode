////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SpaceChargeServiceSBND_service.cc; brief implementation of class for storing/accessing space charge distortions for SBND
// arbint@bnl.gov
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>

// LArSoft includes
#include "sbndcode/SpaceChargeServices/SpaceChargeServiceSBND.h"

// ROOT includes
#include "TMath.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

//-----------------------------------------------
spacecharge::SpaceChargeServiceSBND::SpaceChargeServiceSBND(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
{
    fProp.reset(new spacecharge::SpaceChargeSBND(pset));

    reg.sPreBeginRun.watch(this, &SpaceChargeServiceSBND::preBeginRun);
}

//----------------------------------------------
void spacecharge::SpaceChargeServiceSBND::preBeginRun(const art::Run& run)
{
    fProp->Update(run.id().run());
}

//------------------------------------------------
void spacecharge::SpaceChargeServiceSBND::reconfigure(fhicl::ParameterSet const& pset)
{
    fProp->Configure(pset);
    return;
}

//------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL(spacecharge::SpaceChargeServiceSBND, spacecharge::SpaceChargeService)

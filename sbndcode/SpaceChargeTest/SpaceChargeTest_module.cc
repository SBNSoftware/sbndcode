#ifndef SPACE_CHARGE_TEST
#define SPACE_CHARGE_TEST

/// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

// Larsoft includes
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

using namespace std;

namespace SpaceChargeTools
{
    class SpaceChargeTest;
}

class SpaceChargeTools::SpaceChargeTest : public art::EDAnalyzer
{
public:

    explicit SpaceChargeTest(fhicl::ParameterSet const & p);
    virtual ~SpaceChargeTest();

    virtual void beginJob();
    virtual void endJob();

    // Required function
    void analyze(art::Event const & e) override;
};

SpaceChargeTools::SpaceChargeTest::SpaceChargeTest(fhicl::ParameterSet const & p) : EDAnalyzer(p)
{
}

SpaceChargeTools::SpaceChargeTest::~SpaceChargeTest() {}
void SpaceChargeTools::SpaceChargeTest::beginJob() {}
void SpaceChargeTools::SpaceChargeTest::endJob() {}

void SpaceChargeTools::SpaceChargeTest::analyze(art::Event const & evt)
{
    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    cout << SCE->EnableSimSpatialSCE() << endl;
}

DEFINE_ART_MODULE(SpaceChargeTools::SpaceChargeTest)

#endif //SPACE_CHARGE_TEST

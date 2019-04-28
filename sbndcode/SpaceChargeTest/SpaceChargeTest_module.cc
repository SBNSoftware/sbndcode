#ifndef SPACE_CHARGE_TEST
#define SPACE_CHARGE_TEST

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

// Root includes
#include <TH1.h>

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

    virtual void beginJob() override;
    virtual void endJob() override;

    // Required function
    void analyze(art::Event const & e) override;

private:

    TH1D *hDx;
    TH1D *hDy;
    TH1D *hDz;
    TH1D *hEx;
    TH1D *hEy;
    TH1D *hEz;
};

SpaceChargeTools::SpaceChargeTest::SpaceChargeTest(fhicl::ParameterSet const & p) : EDAnalyzer(p)
{
}

SpaceChargeTools::SpaceChargeTest::~SpaceChargeTest() {}

void SpaceChargeTools::SpaceChargeTest::beginJob()
{
    // Access ART's TFileService, which will handle creating and writing histograms and n-tuples
    art::ServiceHandle<art::TFileService> fileServiceHandle;

    hDx = fileServiceHandle->make<TH1D>("hDx", "", 150, -1.0, 0.6);
    hDy = fileServiceHandle->make<TH1D>("hDy", "", 150, -6.0, 6.0);
    hDz = fileServiceHandle->make<TH1D>("hDz", "", 150, -6.0, 6.0);
    hEx = fileServiceHandle->make<TH1D>("hEx", "", 100, -0.06, 0.06);
    hEy = fileServiceHandle->make<TH1D>("hEy", "", 100, -0.04, 0.04);
    hEz = fileServiceHandle->make<TH1D>("hEz", "", 100, -0.04, 0.04);
}

void SpaceChargeTools::SpaceChargeTest::endJob() {}

void SpaceChargeTools::SpaceChargeTest::analyze(art::Event const & evt)
{
    int xMin = -206;
    int xMax = 206;
    int yMin = -210;
    int yMax = 210;
    int zMin = -10;
    int zMax = 510;

    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    cout << "Is Spatial SCE enabled? " << bool(SCE->EnableSimSpatialSCE()) << endl;
    cout << "Is E-field SCE enabled? " << bool(SCE->EnableSimEfieldSCE()) << endl;

    int nSkip = 5;
    for(int iX = xMin; iX <= xMax - nSkip; iX++)
        {
            iX = iX + nSkip;
            for(int iY = yMin; iY <= yMax - nSkip; iY++)
                {
                    iY = iY + nSkip;
                    for(int iZ = zMin; iZ <= zMax - nSkip; iZ++)
                        {
                            iZ = iZ + nSkip;
                            cout << iX << ", " << iY << ", " << iZ << endl;

                            geo::Point_t point = {double(iX), double(iY), double(iZ)};

                            geo::Vector_t spatialOffsets = SCE->GetPosOffsets(point);
                            if(!((spatialOffsets.X() == spatialOffsets.X()) &&
				 (spatialOffsets.Y() == spatialOffsets.Y()) &&
				 (spatialOffsets.Z() == 0.0)))
                                {
                                    hDx->Fill(spatialOffsets.X());
                                    hDy->Fill(spatialOffsets.Y());
                                    hDz->Fill(spatialOffsets.Z());
                                }

                            geo::Vector_t efieldOffsets = SCE->GetEfieldOffsets(point);
                            if(!((efieldOffsets.X() == efieldOffsets.X()) &&
				 (efieldOffsets.Y() == efieldOffsets.Y()) &&
				 (efieldOffsets.Z() == 0.0)))
                                {
                                    hEx->Fill(efieldOffsets.X());
                                    hEy->Fill(efieldOffsets.Y());
                                    hEz->Fill(efieldOffsets.Z());
                                }
                        }
                }
        }
}

DEFINE_ART_MODULE(SpaceChargeTools::SpaceChargeTest)

#endif //SPACE_CHARGE_TEST

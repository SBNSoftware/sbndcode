/**
 * \brief Unit tests for the CRT geometry
 *
 * \author Marco Del Tutto
 */

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"

#include "TGeoManager.h"
#include "TGeoNode.h"


class SBNDCRTGeometryTest;


class SBNDCRTGeometryTest : public art::EDAnalyzer {
public:
  explicit SBNDCRTGeometryTest(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SBNDCRTGeometryTest(SBNDCRTGeometryTest const&) = delete;
  SBNDCRTGeometryTest(SBNDCRTGeometryTest&&) = delete;
  SBNDCRTGeometryTest& operator=(SBNDCRTGeometryTest const&) = delete;
  SBNDCRTGeometryTest& operator=(SBNDCRTGeometryTest&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

};


SBNDCRTGeometryTest::SBNDCRTGeometryTest(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{}

void SBNDCRTGeometryTest::analyze(art::Event const& e)
{
  art::ServiceHandle<geo::Geometry const> geoService;

  std::cout << "Number of CRT AuxDets is " << geoService->NAuxDets() << std::endl;

  assert(geoService->NAuxDets() > 0);


  std::vector<int> used_copynumbers_auxdet;
  std::vector<int> used_copynumbers_sensitiveauxdet;

  for (size_t i = 0; i < geoService->NAuxDets(); i++) {
    const geo::AuxDetGeo& adGeo = geoService->AuxDet(i);

    std::set<std::string> volNames = { adGeo.TotalVolume()->GetName() };
    std::vector<std::vector<TGeoNode const*> > paths =
      geoService->FindAllVolumePaths(volNames);

    std::string path = "";
    for (size_t inode=0; inode<paths.at(0).size(); inode++) {
      path += paths.at(0).at(inode)->GetName();
      if (inode < paths.at(0).size() - 1) {
        path += "/";
      }
    }

    TGeoManager* manager = geoService->ROOTGeoManager();
    manager->cd(path.c_str());

    // We get the array of strips first, which is the AuxDet,
    // then from the AuxDet, we get the strip by picking the
    // daughter with the ID of the AuxDetSensitive, and finally
    // from the AuxDet, we go up and pick the module and tagger
    TGeoNode* nodeArray = manager->GetCurrentNode();
    TGeoNode* nodeModule = manager->GetMother(1);
    TGeoNode* nodeTagger = manager->GetMother(2);


    // Check every aux det has a different copynumber
    auto iter = std::find(used_copynumbers_auxdet.begin(), used_copynumbers_auxdet.end(),
                          nodeModule->GetNumber());
    assert(iter == used_copynumbers_auxdet.end());

    used_copynumbers_auxdet.push_back(nodeModule->GetNumber());


    // Ensure there are CRT strips
    assert(nodeArray->GetNdaughters() > 0);


    // Ensure that the number of AuxDetSensitive in the geometry service
    // is the same as the number of daughter nodes
    assert((int)geoService->NAuxDetSensitive(i) == nodeArray->GetNdaughters());

    std::cout << "Auxiliary detector ID " << i
              << " with copynumber " << nodeModule->GetNumber()
              << "\n\t strip array name: " << nodeArray->GetName()
              << "\n\t module name: " << nodeModule->GetName()
              << "\n\t tagger name: " << nodeTagger->GetName() << std::endl;

    // Loop over the AuxDetSensitive volumes
    for (int adsid = 0; adsid < nodeArray->GetNdaughters(); adsid++) {

      TGeoNode* nodeStrip = nodeArray->GetDaughter(adsid);

      std::cout << "\t\t has sensitive detector (strip) with ID " << adsid
                << ", name: " << nodeStrip->GetName()
                << " (copynumber " << nodeStrip->GetNumber() << ")" << std::endl;

      // Check every sensitive aux det has a different copynumber (this is used by LArG4)
      auto iter2 = std::find(used_copynumbers_sensitiveauxdet.begin(), used_copynumbers_sensitiveauxdet.end(),
                             nodeStrip->GetNumber());
      assert(iter2 == used_copynumbers_sensitiveauxdet.end());

      used_copynumbers_sensitiveauxdet.push_back(nodeStrip->GetNumber());
    }

  }


}

DEFINE_ART_MODULE(SBNDCRTGeometryTest)




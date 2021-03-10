////////////////////////////////////////////////////////////////////////
// Class:       Calibration
// Plugin Type: analyzer (art v3_05_01)
// File:        Calibration_module.cc
//
// Generated at Wed Dec 16 08:58:23 2020 by Vu Chi Lan Nguyen using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Addition framework
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"

//ROOT includes
#include <TTree.h>
#include <TH2D.h>
namespace sbnd {
  class Calibration;
}


class sbnd::Calibration : public art::EDAnalyzer {
public:
  explicit Calibration(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Calibration(Calibration const&) = delete;
  Calibration(Calibration&&) = delete;
  Calibration& operator=(Calibration const&) = delete;
  Calibration& operator=(Calibration&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  
  TTree *fTree;
  
  // variables to go in tree
  unsigned int 				fEventID;
  std::vector<geo::Point_t> 		*fPos;
  unsigned int 				fHitsCount;
  std::vector<float> 	 		fHits;

  //Histogram
  TH2D *fHitIntegralHist;

  //labels
  std::string fSpacePointLabel;
  std::string fHitLabel;
};


sbnd::Calibration::Calibration(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fSpacePointLabel = p.get<std::string>("SpacePointLabel");
  fHitLabel = p.get<std::string>("HitLabel");
}

void sbnd::Calibration::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  
  fEventID = e.id().event();
 
  // clear
  fPos->clear();
  fHits.clear();

  // initialise counter
  unsigned int nsp = 0;
 
  // Get SpacePoint information
  art::Handle< std::vector<recob::SpacePoint> > spacepointHandle;
  std::vector< art::Ptr<recob::SpacePoint> > spacepointList;

  if(e.getByLabel(fSpacePointLabel, spacepointHandle)){
    art::fill_ptr_vector(spacepointList, spacepointHandle);
  }

  //Find Association between SpacePoint and Hit
  art::FindManyP<recob::Hit> hitAssoc (spacepointList, e, fHitLabel);
   
  //Loop through space point, then loop through hits at each space point
  for(const art::Ptr<recob::SpacePoint> &sp: spacepointList){
     
    fPos->push_back(sp->position());
    nsp++;

    std::vector< art::Ptr<recob::Hit> > spHits = hitAssoc.at(sp.key());
    if(spHits.empty()) continue;
 
    //Get number of hits at this sp
    fHitsCount = spHits.size();     

    for(const art::Ptr<recob::Hit> &hit: spHits){

      if(!hit->WireID().isValid) continue;
   
      int planenum = hit->WireID().Plane;

      if(planenum != 2) continue;
     
      std::cout << "Event: " << e.id().event() << ", space point no.: " << nsp;
    
      fHits.push_back(hit->Integral());
   
      std::cout << ", fHits "  << hit->Integral()  << "\n";     

      fHitIntegralHist->Fill(sp->XYZ()[2], sp->XYZ()[1], hit->Integral());

     } // End of hits

  } //End of spacepoint


  fTree->Fill();

}

void sbnd::Calibration::beginJob()
{
  // Implementation of optional member function here.

  art::ServiceHandle<art::TFileService> tfs;

  fTree = tfs->make<TTree>("tree","Calibration Output Tree");

  fTree->Branch("eventID",&fEventID,"eventID/i");
  fTree->Branch("Position",&fPos);
  fTree->Branch("nHits", &fHitsCount,"nHits/i");
  fTree->Branch("HitIntegral", &fHits);

  fHitIntegralHist = tfs->make<TH2D>("HitIntegral","Hit Integral in Y-Z view", 200, 0, 500, 200, -200, 200);

}

void sbnd::Calibration::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(sbnd::Calibration)

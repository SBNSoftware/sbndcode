////////////////////////////////////////////////////////////////////////
// Class:       SBNDHitAnalyzerExample
// Plugin Type: analyzer (art v3_01_02)
// File:        SBNDHitAnalyzerExample_module.cc
//
// Generated at Sun Apr 21 19:27:27 2019 by Jonathan Asaadi using cetskelgen
// from cetlib version v3_05_01.
////////////////////////////////////////////////////////////////////////

// ##########################
// ### Framework Includes ###
// ##########################
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ############################################################
// ### Optional Includes to be able to write out histograms ###
// ############################################################
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

// ########################
// ### LArSoft Includes ###
// ########################
#include "lardataobj/RecoBase/Hit.h"

// #####################
// ### ROOT Includes ###
// #####################
#include "TTree.h"
#include <TH1F.h>

class SBNDHitAnalyzerExample;


class SBNDHitAnalyzerExample : public art::EDAnalyzer {
public:
  explicit SBNDHitAnalyzerExample(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SBNDHitAnalyzerExample(SBNDHitAnalyzerExample const&) = delete;
  SBNDHitAnalyzerExample(SBNDHitAnalyzerExample&&) = delete;
  SBNDHitAnalyzerExample& operator=(SBNDHitAnalyzerExample const&) = delete;
  SBNDHitAnalyzerExample& operator=(SBNDHitAnalyzerExample&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;


private:

  // ##############################################################
  // ### Declare the string for the module labels you will grab ###
  // ##############################################################
  std::string fHitModuleLabel;
  
  // #######################################################################
  // ### Declare the variables you will check and grab from the fcl file ###
  // #######################################################################
  double fminNumberHits;
  
  // ##############################################
  // ### Declare the histogram you want to fill ###
  // ##############################################
  TH1F* fnHits;

};

// ############################################################
// ### Where the parameters from the fcl file are picked up ###
// ############################################################
SBNDHitAnalyzerExample::SBNDHitAnalyzerExample(fhicl::ParameterSet const& p): 
   EDAnalyzer(p)  ,
   fHitModuleLabel          (p.get< std::string >("HitModuleLabel")),
   fminNumberHits           (p.get<   int       >("minNumberHits", 5))
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void SBNDHitAnalyzerExample::analyze(art::Event const& e)
{

// ###################################
// ### Getting the Hit Information ###
// ###################################
art::Handle< std::vector<recob::Hit> > hitListHandle; //<---Define hitListHandle as a vector of recob::Hit objects
std::vector<art::Ptr<recob::Hit> > hitlist; //<---Define hitlist as a pointer to recob::hits
   
// === Filling the hitlist from the hitlistHandle ===
if (e.getByLabel(fHitModuleLabel,hitListHandle))
   {art::fill_ptr_vector(hitlist, hitListHandle);}  


// === Declaring a useful variable ===   
int nhits = hitlist.size();

// === If there are too few hits in the event, skip it ===
if(nhits < fminNumberHits) {return;}

// === Filling the histogram ===
fnHits->Fill(nhits);


}

void SBNDHitAnalyzerExample::beginJob()
{

art::ServiceHandle<art::TFileService> tfs;
// ### Build your histogram ###
fnHits = tfs->make<TH1F>("nHits", "Number of hits per event", 8000, 0, 8000);

}


void SBNDHitAnalyzerExample::endJob()
{
  // Implementation of optional member function here.
}


DEFINE_ART_MODULE(SBNDHitAnalyzerExample)

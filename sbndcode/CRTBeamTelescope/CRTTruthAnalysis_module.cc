////////////////////////////////////////////////////////////////////////
// Class:       CRTTruthAnalysis
// Plugin Type: analyzer (Unknown Unknown)
// File:        CRTTruthAnalysis_module.cc
//
// Generated at Wed Oct 11 03:20:55 2023 by Jiaoyang Li using cetskelgen
// from  version .
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

#include "TTree.h"
#include "nusimdata/SimulationBase/MCTruth.h"

class CRTTruthAnalysis;


class CRTTruthAnalysis : public art::EDAnalyzer {
public:
  explicit CRTTruthAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTTruthAnalysis(CRTTruthAnalysis const&) = delete;
  CRTTruthAnalysis(CRTTruthAnalysis&&) = delete;
  CRTTruthAnalysis& operator=(CRTTruthAnalysis const&) = delete;
  CRTTruthAnalysis& operator=(CRTTruthAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  std::string _mctruth_label;

};


CRTTruthAnalysis::CRTTruthAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  _mctruth_label     = p.get<std::string>("MCTruthLabel", "generator");
  if(!_data_mode){
    _tree->Branch("mct_sp_pdg", &_mct_sp_pdg, "mct_sp_pdg/F");
  }
}

void CRTTruthAnalysis::analyze(art::Event const& e)
{
  // Implementation of required member function here.
}

void CRTTruthAnalysis::beginJob()
{
  // Implementation of optional member function here.
}

void CRTTruthAnalysis::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(CRTTruthAnalysis)

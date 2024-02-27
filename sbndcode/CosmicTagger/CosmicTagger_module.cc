////////////////////////////////////////////////////////////////////////
// Class:       CosmicTagger
// Plugin Type: producer (Unknown Unknown)
// File:        CosmicTagger_module.cc
//
// Generated at Fri Feb 17 15:50:27 2023 by Marco Del Tutto using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"

class CosmicTagger;


class CosmicTagger : public art::EDProducer {
public:
  explicit CosmicTagger(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CosmicTagger(CosmicTagger const&) = delete;
  CosmicTagger(CosmicTagger&&) = delete;
  CosmicTagger& operator=(CosmicTagger const&) = delete;
  CosmicTagger& operator=(CosmicTagger&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  std::string _mctruth_label;
  std::string _wire_label;

};


CosmicTagger::CosmicTagger(fhicl::ParameterSet const& p)
  : EDProducer{p}  // ,
  // More initializers here.
{
  _mctruth_label = p.get<std::string>("MCTruthLabel", "generator");
  _wire_label = p.get<std::string>("WireLaberl", "caldata");
}

void CosmicTagger::produce(art::Event& e)
{

  geo::GeometryCore const& geom = *(lar::providerFrom<geo::Geometry>());

  //
  // Get simb::MCTruth
  //
  art::Handle<std::vector<simb::MCTruth>> mct_h;
  e.getByLabel(_mctruth_label, mct_h);
  if(!mct_h.isValid()){
    std::cout << "MCTruth product " << _mctruth_label << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<simb::MCTruth>> mct_v;
  art::fill_ptr_vector(mct_v, mct_h);

  // Loop over neutrino interactions
  for (size_t i = 0; i < mct_v.size(); i++) {

    // Neutrino vertex
    double nu_x = mct_v[i]->GetNeutrino().Nu().Vx();
    double nu_y = mct_v[i]->GetNeutrino().Nu().Vy();
    double nu_z = mct_v[i]->GetNeutrino().Nu().Vz();
    std::cout << geom.NearestWireID({nu_x, nu_y, nu_z}, geo::PlaneID(0, 0, 0)) << std::endl;

  }

  //
  // Get recob::Wire
  //
  art::Handle<std::vector<recob::Wire>> wire_h;
  e.getByLabel(_wire_label, wire_h);
  if(!wire_h.isValid()){
    std::cout << "MCTruth product " << _wire_label << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<recob::Wire>> wire_v;
  art::fill_ptr_vector(wire_v, wire_h);

  for (auto const &wire : wire_v) {
    // Wires...
    std::cout << wire->Channel() << std::endl;
  }


}

DEFINE_ART_MODULE(CosmicTagger)

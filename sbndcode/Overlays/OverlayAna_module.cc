////////////////////////////////////////////////////////////////////////
// Class:       OverlayAna
// Plugin Type: analyzer (Unknown Unknown)
// File:        OverlayAna_module.cc
//
// Generated at Fri Aug 22 20:10:51 2025 by Marco Del Tutto using cetskelgen
// from cetlib version 3.18.02.
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
#include "art_root_io/TFileService.h"

#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "TTree.h"

class OverlayAna;


class OverlayAna : public art::EDAnalyzer {
public:
  explicit OverlayAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OverlayAna(OverlayAna const&) = delete;
  OverlayAna(OverlayAna&&) = delete;
  OverlayAna& operator=(OverlayAna const&) = delete;
  OverlayAna& operator=(OverlayAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  std::vector<std::string> _rawdigit_producers;

  TTree* _tree;

  int _run, _subrun, _event;
  std::vector<std::vector<float>> _tpc_waveforms_data;
  std::vector<std::vector<float>> _tpc_waveforms_mc;
  std::vector<std::vector<float>> _tpc_waveforms_overlay;

};


OverlayAna::OverlayAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  _rawdigit_producers = p.get<std::vector<std::string>>("RawDigits");

  if (_rawdigit_producers.size() != 3) {
    std::cout << "Need to provide three producer labels (Data, MC, Overlay)!" << std::endl;
    throw std::exception();
  }

  art::ServiceHandle<art::TFileService> fs;
  _tree = fs->make<TTree>("tree","");
  _tree->Branch("waveforms_data", "std::vector<std::vector<float>>", &_tpc_waveforms_data);
  _tree->Branch("waveforms_mc", "std::vector<std::vector<float>>", &_tpc_waveforms_mc);
  _tree->Branch("waveforms_overlay", "std::vector<std::vector<float>>", &_tpc_waveforms_overlay);

}

void OverlayAna::analyze(art::Event const& e)
{
  
  // art::ServiceHandle<geo::Geometry> geom;

  geo::WireReadoutGeom const& wr_geom = art::ServiceHandle<geo::WireReadout>()->Get();

  std::vector<std::vector<art::Ptr<raw::RawDigit>>> rawdigits;

  for(auto producer : _rawdigit_producers) {
    art::Handle<std::vector<raw::RawDigit>> rawdigit_h;
    e.getByLabel(producer, rawdigit_h);
    if(!rawdigit_h.isValid()){
      std::cout << "RawDigit product " << producer << " not found..." << std::endl;
      throw std::exception();
    }
    std::cout << "Found RawDigit product " << producer << std::endl;
    std::vector<art::Ptr<raw::RawDigit>> rawdigit_v;
    art::fill_ptr_vector(rawdigit_v, rawdigit_h);
    rawdigits.push_back(rawdigit_v);
  }


  // std::cout << "Number of available RawDigits " << rawdigit_v.size() << std::endl;
  
  _tpc_waveforms_data.clear();
  _tpc_waveforms_data.resize(wr_geom.Nchannels());
  _tpc_waveforms_mc.clear();
  _tpc_waveforms_mc.resize(wr_geom.Nchannels());
  _tpc_waveforms_overlay.clear();
  _tpc_waveforms_overlay.resize(wr_geom.Nchannels());

  for (int i = 0; i < 3; i++)
  {
    auto rawdigit_v = rawdigits[i];

    for (auto const &rawdigit : rawdigit_v)
    {

      unsigned int ch  = rawdigit->Channel();
      // float        ped = rawdigit->GetPedestal();

      auto adcs = rawdigit->ADCs();
      if(i == 0) _tpc_waveforms_data[ch].assign(adcs.begin(), adcs.end());
      if(i == 1) _tpc_waveforms_mc[ch].assign(adcs.begin(), adcs.end());
      if(i == 2) _tpc_waveforms_overlay[ch].assign(adcs.begin(), adcs.end());
    }
  }

  _tree->Fill();
}

DEFINE_ART_MODULE(OverlayAna)

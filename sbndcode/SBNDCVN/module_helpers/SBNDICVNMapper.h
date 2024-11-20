#ifndef LCVN_SBNDICVNMAPPER_H
#define LCVN_SBNDICVNMAPPER_H

#include  <iostream>
#include  <ostream>
#include  <list>
#include  <algorithm>
#include <numeric>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileDirectory.h"

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "larrecodnn/CVN/func/PixelMap.h"
#include "sbndcode/SBNDCVN/module_helpers/SBNDPixelMap.h"
#include "larrecodnn/CVN/interfaces/PixelMapProducer.h"

namespace lcvn {

  template <class T, class U>
  class SBNDICVNMapper : public art::EDProducer {
  public:
    explicit SBNDICVNMapper(fhicl::ParameterSet const& pset);

    void produce(art::Event& evt) override;

  protected:
    /// Module lablel for input clusters
    std::string fHitsModuleLabel;

    /// Instance lablel for cluster pixelmaps
    std::string fClusterPMLabel;

    /// Minimum number of hits for cluster to be converted to pixel map
    unsigned short fMinClusterHits;

    /// PixelMapProducer does the work for us
    T fProducer;
    
    bool fverbose;
    bool fUseSlice; 
    std::string fSliceLabel;
    std::string fPFParticleModuleLabel;
    std::string fT0Label;
    unsigned int fMapVecSize;
  };

}
#endif

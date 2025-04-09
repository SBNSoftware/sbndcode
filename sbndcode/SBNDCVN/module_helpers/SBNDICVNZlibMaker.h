#ifndef LCVN_SBNDICVNZLIBMAKER_H
#define LCVN_SBNDICVNZLIBMAKER_H

#include  <iostream>
#include  <ostream>
#include  <list>
#include  <algorithm>
#include <numeric>
#include <cstdlib>

// C/C++ includes
#include <cstdlib>
#include <iostream>

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileDirectory.h"
#include "boost/filesystem.hpp"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/Exception.h"

// Data products
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCTruth.h"
// #include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

// CVN includes
#include "larrecodnn/CVN/func/AssignLabels.h"
#include "larrecodnn/CVN/func/CVNImageUtils.h"
#include "larrecodnn/CVN/func/InteractionType.h"
#include "larrecodnn/CVN/func/LArTrainingData.h"
#include "larrecodnn/CVN/func/PixelMap.h"

// Compression
#include "math.h"
#include "zlib.h"

#include "TH1.h"

namespace lcvn {

  class SBNDICVNZlibMaker : public art::EDAnalyzer {
  public:
    explicit SBNDICVNZlibMaker(fhicl::ParameterSet const& pset);
    ~SBNDICVNZlibMaker();

    void beginJob() override;
    void analyze(const art::Event& evt) override {}
    void reconfigure(const fhicl::ParameterSet& pset);

  protected:
    std::string fOutputDir;
    std::string fPixelMapInput;
    bool fSetLog;
    std::vector<bool> fReverseViews;
    unsigned int fPlaneLimit;
    unsigned int fTDCLimit;
    bool fverbose;
    bool fUseSlice;
    std::string fSliceLabel;

    std::string out_dir;
    CVNImageUtils fImage;

    template <class T>
    void write_files(LArTrainingData<T> td, std::string evtid) = delete;
  };
}
#endif

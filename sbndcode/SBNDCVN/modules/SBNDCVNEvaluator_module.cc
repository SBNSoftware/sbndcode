////////////////////////////////////////////////////////////////////////
// \file    LArCVNEvaluator_module.cc
// \brief   Producer module creating CVN neural net results
// \author  Alexander Radovic - a.radovic@gmail.com
//          Saul Alonso Monsalve - saul.alonso.monsalve@cern.ch
////////////////////////////////////////////////////////////////////////

// C/C++ includes
#include <iostream>
#include <sstream>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larrecodnn/CVN/func/AssignLabels.h"
#include "larrecodnn/CVN/func/InteractionType.h"
#include "sbndcode/SBNDCVN/module_helpers/SBNDPixelMap.h"
#include "larrecodnn/CVN/func/Result.h"
#include "sbndcode/SBNDCVN/module_helpers/SBNDITFNetHandler.h"
#include "lardata/Utilities/AssociationUtil.h"

namespace lcvn {

  class SBNDCVNEvaluator : public art::EDProducer {
  public:
    explicit SBNDCVNEvaluator(fhicl::ParameterSet const& pset);

    void produce(art::Event& evt);

  private:
    /// Module label for input pixel maps
    std::string fPixelMapInput;
    std::string fPixelMapModuleLabel;
    std::string fResultLabel;

    /// Can use Caffe or Tensorflow
    std::string fCVNType;

    //lcvn::CaffeNetHandler fCaffeHandler;
    std::unique_ptr<lcvn::SBNDITFNetHandler> fTFHandler;

    /// Number of outputs fron neural net
    //unsigned int fNOutput;

    /// If there are multiple SBNDpixel maps per event can we use them?
    bool fMultiplePMs;
  };

  //.......................................................................
  SBNDCVNEvaluator::SBNDCVNEvaluator(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , fPixelMapInput(pset.get<std::string>("PixelMapInput"))
    , fPixelMapModuleLabel(pset.get<std::string>("PixelMapModuleLabel"))		    
    , fResultLabel(pset.get<std::string>("ResultLabel"))
    , fCVNType(pset.get<std::string>("CVNType"))
    ,
    //fCaffeHandler       (pset.get<fhicl::ParameterSet> ("CaffeNetHandler")),
    fTFHandler{art::make_tool<SBNDITFNetHandler>(pset.get<fhicl::ParameterSet>("SBNDTFHandler"))}
    ,
    //fNOutput       (fCaffeHandler.NOutput()),
    fMultiplePMs(pset.get<bool>("MultiplePMs"))
  {
    produces<std::vector<lcvn::Result>>(fResultLabel);
    produces<art::Assns<lcvn::SBNDPixelMap, lcvn::Result>>(fResultLabel);
  }
  //......................................................................
  void SBNDCVNEvaluator::produce(art::Event& evt)
  {

    /// Define containers for the things we're going to produce
    std::unique_ptr<std::vector<Result>> resultCol(new std::vector<Result>);
    auto assn = std::make_unique< art::Assns<lcvn::SBNDPixelMap, lcvn::Result> >();
    /// Load in the pixel maps
    std::vector<art::Ptr<lcvn::SBNDPixelMap>> pixelmaplist;
    art::InputTag itag1(fPixelMapModuleLabel, fPixelMapInput);
    auto pixelmapListHandle = evt.getHandle<std::vector<lcvn::SBNDPixelMap>>(itag1);
    if (pixelmapListHandle) art::fill_ptr_vector(pixelmaplist, pixelmapListHandle);

    if (fCVNType == "TF" || fCVNType == "Tensorflow" || fCVNType == "TensorFlow") {

      // If we have a pixel map then use the TF interface to give us a prediction
      if (pixelmaplist.size() > 0) {
        std::cout << "===================== size of the pixelmap list : " << pixelmaplist.size() << " ========================\n"; 
        std::vector<std::vector<float>> networkOutput = fTFHandler->Predict(*pixelmaplist[0]);
        // lcvn::Result can now take a vector of floats and works out the number of outputs
        resultCol->emplace_back(networkOutput);
        util::CreateAssn(*this, evt, *resultCol, pixelmaplist[0], *assn, fResultLabel);

        // Classify other pixel maps if they exist
        if (fMultiplePMs) {
          for (unsigned int p = 1; p < pixelmaplist.size(); ++p) {
            std::vector<std::vector<float>> output = fTFHandler->Predict(*pixelmaplist[p]);
            resultCol->emplace_back(output);
            util::CreateAssn(*this, evt, *resultCol, pixelmaplist[p], *assn, fResultLabel);
          }
        }
      }
    }
    else {
      mf::LogError("SBNDCVNEvaluator::produce")
        << "CVN Type not in the allowed list: Tensorflow" << std::endl;
      mf::LogError("SBNDCVNEvaluator::produce") << "Exiting without processing events" << std::endl;
      return;
    }

    evt.put(std::move(resultCol), fResultLabel);
    evt.put(std::move(assn), fResultLabel);
  }

  DEFINE_ART_MODULE(lcvn::SBNDCVNEvaluator)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////

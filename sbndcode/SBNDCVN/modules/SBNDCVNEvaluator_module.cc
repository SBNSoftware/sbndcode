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

#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larrecodnn/CVN/func/AssignLabels.h"
#include "larrecodnn/CVN/func/InteractionType.h"
#include "larrecodnn/CVN/func/Result.h"
#include "sbndcode/SBNDCVN/module_helpers/SBNDPixelMap.h"
#include "sbndcode/SBNDCVN/module_helpers/SBNDITFNetHandler.h"
#include "sbndcode/SBNDCVN/module_helpers/SBNDPixelMapProducer.h"

namespace lcvn {

  class SBNDCVNEvaluator : public art::EDProducer {
  public:
    explicit SBNDCVNEvaluator(fhicl::ParameterSet const& pset);

    void produce(art::Event& evt);

  private:
    /// Module label for input pixel maps
    std::string fPixelMapInput;
    std::string fPixelMapModuleLabel;
    art::InputTag fSliceLabel;
    art::InputTag fPFParticleModuleLabel;
    art::InputTag fT0Label;

    /// Can use Caffe or Tensorflow
    std::string fCVNType;

    //lcvn::CaffeNetHandler fCaffeHandler;
    std::unique_ptr<lcvn::SBNDITFNetHandler> fTFHandler;

    /// Number of outputs fron neural net
    //unsigned int fNOutput;

    SBNDPixelMapHitProducer fPMProducer;
  };

  //.......................................................................
  SBNDCVNEvaluator::SBNDCVNEvaluator(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , fPixelMapInput(pset.get<std::string>("PixelMapInput"))
    , fPixelMapModuleLabel(pset.get<std::string>("PixelMapModuleLabel"))
    , fSliceLabel(pset.get<art::InputTag>("SliceLabel"))
    , fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel"))
    , fT0Label(pset.get<art::InputTag>("T0Label"))
    , fCVNType(pset.get<std::string>("CVNType"))
    , fTFHandler{art::make_tool<SBNDITFNetHandler>(pset.get<fhicl::ParameterSet>("SBNDTFHandler"))}
    , fPMProducer(pset.get<fhicl::ParameterSet>("PixelMapProducer"))
  {
    produces<std::vector<lcvn::Result>>();
    if (fSliceLabel != ""){
      produces<art::Assns<recob::Slice, lcvn::Result>>();
    }
    else{
      produces<art::Assns<lcvn::SBNDPixelMap, lcvn::Result>>();
    }
  }
  //......................................................................
  void SBNDCVNEvaluator::produce(art::Event& evt)
  {

    /// Define containers for the things we're going to produce
    std::unique_ptr<std::vector<Result>> resultCol(new std::vector<Result>);
    auto assn1 = std::make_unique< art::Assns<recob::Slice, lcvn::Result> >();
    auto assn2 = std::make_unique< art::Assns<lcvn::SBNDPixelMap, lcvn::Result> >();

    if (fCVNType == "TF" || fCVNType == "Tensorflow" || fCVNType == "TensorFlow") {
      if (fSliceLabel != ""){ // by default use slice information
        std::vector < art::Ptr < recob::Slice > > slcList;
        auto slcHandle = evt.getHandle<std::vector<recob::Slice>>(fSliceLabel);
        if (!slcHandle.isValid()){
          throw cet::exception("SBNDCVNEvaluator") << "Unable to get slices using label "  << fSliceLabel;
        }
        else{
          art::fill_ptr_vector(slcList, slcHandle);
        }
        art::Handle< std::vector<recob::PFParticle> > PFPListHandle;
        std::vector< art::Ptr<recob::PFParticle> > PFPList;
        if( evt.getByLabel(fPFParticleModuleLabel,PFPListHandle))
          art::fill_ptr_vector(PFPList,PFPListHandle);

        art::FindManyP<recob::Hit> findManyHits(slcHandle, evt, fSliceLabel);
        art::FindManyP<recob::PFParticle> findManyPFPs(slcHandle, evt, fPFParticleModuleLabel);
        art::FindManyP<anab::T0> findManyT0s(PFPListHandle, evt, fT0Label);
        auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);
        for(auto const& slice : slcList){
	  std::cout << "********* " << evt.run() << "  " << evt.subRun() << "  " << evt.id().event() << "  " << slice->ID() << "  **************\n";
	  std::vector<float> pfp_T0_vec;
	  if(findManyPFPs.isValid()){
            std::vector<art::Ptr<recob::PFParticle>> slicePFPs = findManyPFPs.at(slice.key());
            if(slicePFPs.size()){
              for(auto const &pfp : slicePFPs){
                if(findManyT0s.isValid()){
                  std::vector<art::Ptr<anab::T0>> T0_vec = findManyT0s.at(pfp.key());
                  if(T0_vec.size()){
                    for(auto const& T0 : T0_vec){
                      pfp_T0_vec.push_back(T0->Time());
                    }
                  }
                }
              }
            }
	  }    
	  
	  float min_T0 = 0.;
	  if(pfp_T0_vec.size()){
            min_T0 = *min_element(pfp_T0_vec.begin(), pfp_T0_vec.end());
	  } 		      
	      
          if(findManyHits.isValid()){
            std::vector<art::Ptr<recob::Hit>> slicehits = findManyHits.at(slice.key());
            fPMProducer.Set_fT0_value(min_T0);
            SBNDPixelMap pm = fPMProducer.SBNDCreateMap(detProp, slicehits);
            auto nhits = fPMProducer.NROI();
            pm.SetTotHits(nhits);
            pm.fSliceID = slice->ID();
            std::vector<std::vector<float>> output = fTFHandler->Predict(pm);
            resultCol->emplace_back(output);
            util::CreateAssn(*this, evt, *resultCol, slice, *assn1);
          }
        }
      }
      else{ //Try to read pixel maps saved in the file
        /// Load in the pixel maps
        std::vector<art::Ptr<lcvn::SBNDPixelMap>> pixelmaplist;
        art::InputTag itag1(fPixelMapModuleLabel, fPixelMapInput);
        auto pixelmapListHandle = evt.getHandle<std::vector<lcvn::SBNDPixelMap>>(itag1);
        if (pixelmapListHandle) art::fill_ptr_vector(pixelmaplist, pixelmapListHandle);
        
        // If we have a pixel map then use the TF interface to give us a prediction
        if (pixelmaplist.size() > 0) {
          std::cout << "===================== size of the pixelmap list : " << pixelmaplist.size() << " ========================\n"; 
          std::vector<std::vector<float>> networkOutput = fTFHandler->Predict(*pixelmaplist[0]);
          // lcvn::Result can now take a vector of floats and works out the number of outputs
          for (unsigned int p = 0; p < pixelmaplist.size(); ++p) {
            std::vector<std::vector<float>> output = fTFHandler->Predict(*pixelmaplist[p]);
            resultCol->emplace_back(output);
            util::CreateAssn(*this, evt, *resultCol, pixelmaplist[p], *assn2);
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

    evt.put(std::move(resultCol));
    if (fSliceLabel != ""){
      evt.put(std::move(assn1));
    }
    else{
      evt.put(std::move(assn2));
    }
  }

  DEFINE_ART_MODULE(lcvn::SBNDCVNEvaluator)
} // end namespace cvn
////////////////////////////////////////////////////////////////////////

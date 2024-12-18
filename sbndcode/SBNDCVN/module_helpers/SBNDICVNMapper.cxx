#include "sbndcode/SBNDCVN/module_helpers/SBNDICVNMapper.h"
#include "lardataobj/RecoBase/Slice.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"

#include <algorithm> // std::min_element
#include <iterator>  // std::begin, std::end

namespace lcvn
{
  template <class T, class U> SBNDICVNMapper<T, U>::SBNDICVNMapper(fhicl::ParameterSet const& pset)
    : EDProducer{pset}
    , fHitsModuleLabel(pset.get<std::string>("HitsModuleLabel"))
    , fClusterPMLabel(pset.get<std::string>("ClusterPMLabel"))
    , fMinClusterHits(pset.get<unsigned short>("MinClusterHits"))
    , fProducer(pset.get<fhicl::ParameterSet>("PixelMapProducer"))
    , fverbose(pset.get<bool>("verbose"))
    , fUseSlice(pset.get<bool>("UseSlice"))
    , fSliceLabel(pset.get<std::string>("SliceLabel"))
    , fPFParticleModuleLabel(pset.get<std::string>("PFParticleModuleLabel"))
    , fT0Label(pset.get<std::string>("T0Label"))
    , fMapVecSize(pset.get<unsigned int>("MapVecSize"))
  {

    produces<std::vector<lcvn::PixelMap>>(fClusterPMLabel);
    produces<art::Assns<recob::Slice, lcvn::PixelMap>>(fClusterPMLabel);
  }
	
  //--------------------------------------------------------------------------------------------------------------------------	
	
template <class T, class U> void SBNDICVNMapper<T, U>::produce(art::Event& evt)
{
   if(fverbose) std::cout << "============ Calling the function SBNDICVNMapper::produce() ==============\n";
   
   if(fUseSlice){
      if(fverbose) std::cout << "============ Calling the function SBNDICVNMapper::produce() is using slices ==============\n";
      art::Handle< std::vector<recob::Slice> > SliceListHandle;
      std::vector< art::Ptr<recob::Slice> > SliceList;
      if( evt.getByLabel(fSliceLabel,SliceListHandle))
          art::fill_ptr_vector(SliceList,SliceListHandle);
      
      art::Handle< std::vector<recob::PFParticle> > PFPListHandle;
      std::vector< art::Ptr<recob::PFParticle> > PFPList;
      if( evt.getByLabel(fPFParticleModuleLabel,PFPListHandle))
          art::fill_ptr_vector(PFPList,PFPListHandle);
      
      art::FindManyP<U> findManyHits(SliceListHandle, evt, fSliceLabel);
      
      art::FindManyP<recob::PFParticle> findManyPFPs(SliceListHandle, evt, fPFParticleModuleLabel);
      art::FindManyP<anab::T0> findManyT0s(PFPListHandle, evt, fT0Label);
      
      std::unique_ptr< std::vector<PixelMap> > pmCol(new std::vector<PixelMap>);
      auto assn = std::make_unique< art::Assns<recob::Slice, lcvn::PixelMap> >();

      auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);
      
      for(auto const& slice : SliceList){
          if (fverbose) std::cout << "********* " << evt.run() << "  " << evt.subRun() << "  " << evt.id().event() << "  " << slice->ID() << "  **************\n";
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
	     std::vector<art::Ptr<U>> slicehits = findManyHits.at(slice.key());
	     fProducer.Set_fT0_value(min_T0);
	     PixelMap pm = fProducer.SBNDCreateMap(detProp, slicehits);
	     auto nhits = fProducer.NROI();
             pm.SetTotHits(nhits);
	     //pm.fSliceID = slice->ID();
             
	     if(nhits > fMinClusterHits && pmCol->size()<fMapVecSize){ 
	        pmCol->push_back(pm);
                util::CreateAssn(*this, evt, *pmCol, slice, *assn, fClusterPMLabel);
	     }
	  
          }
      }
      //std::cout<<pmCol->size()<<std::endl;
      evt.put(std::move(pmCol), fClusterPMLabel);
      evt.put(std::move(assn), fClusterPMLabel);
     
   }
   
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
   else{
       if(fverbose) std::cout << "============ Calling the function SBNDICVNMapper::produce() is using full event ==============\n";
       std::vector<art::Ptr<U>> hitlist;
       auto hitListHandle = evt.getHandle<std::vector<U>>(fHitsModuleLabel);
       if (hitListHandle) art::fill_ptr_vector(hitlist, hitListHandle);

       //Declaring containers for things to be stored in event
       std::unique_ptr<std::vector<PixelMap>> pmCol(new std::vector<PixelMap>);

       auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);
       PixelMap pm = fProducer.SBNDCreateMap(detProp, hitlist);
       auto nhits = fProducer.NROI();
       pm.SetTotHits(nhits);

       if (nhits > fMinClusterHits) pmCol->push_back(pm);

       evt.put(std::move(pmCol), fClusterPMLabel);
   }
  }
}

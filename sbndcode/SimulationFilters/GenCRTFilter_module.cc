#include <iostream>
#include <algorithm>

#include "TGeoManager.h"

#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace filt{

  class GenFilter : public art::EDFilter {
    public:
      explicit GenFilter(fhicl::ParameterSet const & pset);
      virtual ~GenFilter() {};
      virtual bool filter(art::Event& e);
      void reconfigure(fhicl::ParameterSet const& pset);
      void beginJob();

    private:

      std::vector<unsigned int> fTopHighCRTAuxDetIDs; 
      std::vector<unsigned int> fTopLowCRTAuxDetIDs; 
      std::vector<unsigned int> fBottomCRTAuxDetIDs; 
      std::vector<unsigned int> fFrontCRTAuxDetIDs; 
      std::vector<unsigned int> fBackCRTAuxDetIDs; 
      std::vector<unsigned int> fLeftCRTAuxDetIDs; 
      std::vector<unsigned int> fRightCRTAuxDetIDs; 

      bool fUseTopHighCRTs; 
      bool fUseTopLowCRTs; 
      bool fUseBottomCRTs; 
      bool fUseFrontCRTs; 
      bool fUseBackCRTs; 
      bool fUseLeftCRTs; 
      bool fUseRightCRTs; 
      std::vector<int> fPDGs;
      std::vector<double> fMinMomentums;
      std::vector<double> fMaxMomentums;
      double fCRTDimensionScaling;

      bool IsInterestingParticle(const simb::MCParticle &particle);
      void LoadCRTAuxDetIDs();
      bool UsesCRTAuxDets(const simb::MCParticle &particle, const std::vector<unsigned int> &crt_auxdet_vector);
      bool UsesCRTAuxDet(const simb::MCParticle &particle, geo::AuxDetGeo* const crt);
      bool RayIntersectsBox(TVector3 ray_origin, TVector3 ray_direction, TVector3 box_min_extent, TVector3 box_max_extent);
  };


  GenFilter::GenFilter::GenFilter(fhicl::ParameterSet const & pset)
  {
    this->reconfigure(pset);
  }


  void GenFilter::reconfigure(fhicl::ParameterSet const& pset){
    fUseTopHighCRTs = pset.get<bool>("UseTopHighCRTs");
    fUseTopLowCRTs = pset.get<bool>("UseTopLowCRTs");
    fUseBottomCRTs = pset.get<bool>("UseBottomCRTs");
    fUseFrontCRTs = pset.get<bool>("UseFrontCRTs");
    fUseBackCRTs = pset.get<bool>("UseBackCRTs");
    fUseLeftCRTs = pset.get<bool>("UseLeftCRTs");
    fUseRightCRTs = pset.get<bool>("UseRightCRTs");
    fPDGs = pset.get<std::vector<int> >("PDGs");
    fMinMomentums = pset.get<std::vector<double> >("MinMomentums");
    fMaxMomentums = pset.get<std::vector<double> >("MaxMomentums");
    fCRTDimensionScaling = pset.get<double>("CRTDimensionScaling");
  }


  bool GenFilter::filter(art::Event & e){
    std::vector< art::Handle< std::vector<simb::MCTruth> > > mclists;
    e.getManyByType(mclists);
    for (unsigned int i = 0; i < mclists.size() ; i++){
      for (unsigned int j = 0; j < mclists[i]->size(); j++){
        //Should have the truth record for the event now
        const art::Ptr<simb::MCTruth> mc_truth(mclists[i],j);
        for (int part = 0; part < mc_truth->NParticles(); part++){
          const simb::MCParticle particle = mc_truth->GetParticle(part);
          if (!IsInterestingParticle(particle)) continue;
          if (fUseTopHighCRTs){
            bool OK = UsesCRTAuxDets(particle,fTopHighCRTAuxDetIDs);
            if (!OK) return false;
            //std::cout<<"TopHighCRTs: " << OK << std::endl;
          }
          if (fUseTopLowCRTs){
            bool OK = UsesCRTAuxDets(particle,fTopLowCRTAuxDetIDs);
            if (!OK) return false;

            //std::cout<<"TopLowCRTs: " << OK << std::endl;
          }
          if (fUseBottomCRTs){
            bool OK = UsesCRTAuxDets(particle,fBottomCRTAuxDetIDs);
            if (!OK) return false;
            //std::cout<<"BottomCRTs: " << OK << std::endl;
          }
          if (fUseFrontCRTs){
            bool OK = UsesCRTAuxDets(particle,fFrontCRTAuxDetIDs);
            if (!OK) return false;
            //std::cout<<"FrontCRTs: " << OK << std::endl;
          }
          if (fUseBackCRTs){
            bool OK = UsesCRTAuxDets(particle,fBackCRTAuxDetIDs);
            if (!OK) return false;
            //std::cout<<"BackCRTs: " << OK << std::endl;
          }
          if (fUseLeftCRTs){
            bool OK = UsesCRTAuxDets(particle,fLeftCRTAuxDetIDs);
            if (!OK) return false;
            //std::cout<<"LeftCRTs: " << OK << std::endl;
          }
          if (fUseRightCRTs){
            bool OK = UsesCRTAuxDets(particle,fRightCRTAuxDetIDs);
            if (!OK) return false;
            //std::cout<<"RightCRTs: " << OK << std::endl;
          }
        }
      }
    }

    return true;
  }


  void GenFilter::beginJob() {
    LoadCRTAuxDetIDs();
  }


  bool GenFilter::IsInterestingParticle(const simb::MCParticle &particle){
    double mom = particle.Momentum().Vect().Mag();
    int pdg = particle.PdgCode();

    for (unsigned int particle_i = 0; particle_i < fPDGs.size(); particle_i++){
      //Check minimum momentum
      if (fMinMomentums[particle_i] > 0 && mom < fMinMomentums[particle_i]){
        return false;
      }
      //Check maximum momentum
      if (fMaxMomentums[particle_i] > 0 && mom > fMaxMomentums[particle_i]){
        return false;
      }
      //Check PDG
      if (fPDGs[particle_i]!=0 && pdg!=fPDGs[particle_i]){
        return false;
      }
    }
    return true;
  }


  void GenFilter::LoadCRTAuxDetIDs(){
    art::ServiceHandle<geo::Geometry> geom;

    for (unsigned int auxdet_i = 0; auxdet_i < geom->AuxDetGeoVec().size(); auxdet_i++){
      geo::AuxDetGeo* crt = geom->AuxDetGeoVec()[auxdet_i];
      const TGeoVolume* volModule = crt->TotalVolume();
      std::set<std::string> volNames = { volModule->GetName() };
      std::vector<std::vector<TGeoNode const*> > paths = geom->FindAllVolumePaths(volNames);

      std::string path = "";
      for (size_t inode=0; inode<paths.at(0).size(); inode++) {
        path += paths.at(0).at(inode)->GetName();
        if (inode < paths.at(0).size() - 1) {
            path += "/";
          }
      }

      TGeoManager* manager = geom->ROOTGeoManager();
      manager->cd(path.c_str());

      manager->GetCurrentNode();
      manager->GetMother(1);
      TGeoNode* nodeTagger = manager->GetMother(2);

      std::string taggerName = nodeTagger->GetName();

      //Now a bunch of if statements so that we can fill our CRT AuxDet ID vectors
      //BLEH
      if (taggerName.find("TopHigh")!=std::string::npos){
        fTopHighCRTAuxDetIDs.push_back(auxdet_i);
      }
      else if (taggerName.find("TopLow")!=std::string::npos){
        fTopLowCRTAuxDetIDs.push_back(auxdet_i);
      }
      else if (taggerName.find("Bot")!=std::string::npos){
        fBottomCRTAuxDetIDs.push_back(auxdet_i);
      }
      else if (taggerName.find("Front")!=std::string::npos){
        fFrontCRTAuxDetIDs.push_back(auxdet_i);
      }
      else if (taggerName.find("Back")!=std::string::npos){
        fBackCRTAuxDetIDs.push_back(auxdet_i);
      }
      else if (taggerName.find("Left")!=std::string::npos){
        fLeftCRTAuxDetIDs.push_back(auxdet_i);
      }
      else if (taggerName.find("Right")!=std::string::npos){
        fRightCRTAuxDetIDs.push_back(auxdet_i);
      }
      else {
        std::cout<<"Tagger with name: " << taggerName << " does not fit the logic.  This should not happen!!!!!!one"<<std::endl;
      }
    }

    std::cout<< "No. top high CRT AuxDets found: " << fTopHighCRTAuxDetIDs.size() << std::endl;
    std::cout<< "No. top low CRT AuxDets found: " << fTopLowCRTAuxDetIDs.size() << std::endl;
    std::cout<< "No. bottom CRT AuxDets found: " << fBottomCRTAuxDetIDs.size() << std::endl;
    std::cout<< "No. front CRT AuxDets found: " << fFrontCRTAuxDetIDs.size() << std::endl;
    std::cout<< "No. back CRT AuxDets found: " << fBackCRTAuxDetIDs.size() << std::endl;
    std::cout<< "No. left CRT AuxDets found: " << fLeftCRTAuxDetIDs.size() << std::endl;
    std::cout<< "No. right CRT AuxDets found: " << fRightCRTAuxDetIDs.size() << std::endl;
    return;
  }


  bool GenFilter::UsesCRTAuxDets(const simb::MCParticle &particle, const std::vector<unsigned int> &crt_auxdet_vector){
    //Loop over the aux dets, extract each one and then perform the test
    art::ServiceHandle<geo::Geometry> geom;

    for (unsigned int i = 0; i < crt_auxdet_vector.size(); i++){
      unsigned int auxdet_index = crt_auxdet_vector[i];
      geo::AuxDetGeo* crt = geom->AuxDetGeoVec()[auxdet_index];
      if (UsesCRTAuxDet(particle,crt)){
        return true;
      }
    }

    return false;
  }

  bool GenFilter::UsesCRTAuxDet(const simb::MCParticle &particle, geo::AuxDetGeo* const crt){
    //We need to prepare the particle's position and direction and construct a bounding box from the CRT for the ray-box collision algorithm
    //Start with the particle
    double particle_position_array[3], particle_direction_array[3], particle_local_position_array[3], particle_local_direction_array[3];
    particle_position_array[0] = particle.Position(0).X();
    particle_position_array[1] = particle.Position(0).Y();
    particle_position_array[2] = particle.Position(0).Z();
    particle_direction_array[0] = particle.Momentum(0).X()/particle.Momentum(0).Vect().Mag();
    particle_direction_array[1] = particle.Momentum(0).Y()/particle.Momentum(0).Vect().Mag();
    particle_direction_array[2] = particle.Momentum(0).Z()/particle.Momentum(0).Vect().Mag();
    crt->WorldToLocal(particle_position_array,particle_local_position_array);
    crt->WorldToLocalVect(particle_direction_array,particle_local_direction_array);
    TVector3 particle_local_position(particle_local_position_array[0],particle_local_position_array[1],particle_local_position_array[2]);
    TVector3 particle_local_direction(particle_local_direction_array[0],particle_local_direction_array[1],particle_local_direction_array[2]);
    double norm[3], lnorm[3];
    double center[3];
    crt->GetNormalVector(norm);
    crt->WorldToLocalVect(norm,lnorm);
    crt->GetCenter(center);

    //Now make the bounding box min and max extent from the CRT
    //In local coordinates, the normal of the CRT is parallel to the z-axis, the length is parallel to z, width parallel to x and height parallel to y
    //A gotchya: AuxDets pass half widths and half heights but FULL lengths.  So, divide length by 2
    TVector3 crt_min_extent(-1.*crt->HalfWidth1(),-1.*crt->HalfHeight(),-1.*crt->Length()/2.);
    //Scale the dimensions if the user wants them scaling
    crt_min_extent*=fCRTDimensionScaling;
    //Now make the max extent
    //Use symmetry to reduce change of a typo
    TVector3 crt_max_extent = crt_min_extent*-1.;

    //Now run the algorithm
    return RayIntersectsBox(particle_local_position,particle_local_direction,crt_min_extent,crt_max_extent);
  }


  bool GenFilter::RayIntersectsBox(TVector3 ray_origin, TVector3 ray_direction, TVector3 box_min_extent, TVector3 box_max_extent){
    double t1 = (box_min_extent.X() - ray_origin.X())/ray_direction.X();
    double t2 = (box_max_extent.X() - ray_origin.X())/ray_direction.X();
    double t3 = (box_min_extent.Y() - ray_origin.Y())/ray_direction.Y();
    double t4 = (box_max_extent.Y() - ray_origin.Y())/ray_direction.Y();
    double t5 = (box_min_extent.Z() - ray_origin.Z())/ray_direction.Z();
    double t6 = (box_max_extent.Z() - ray_origin.Z())/ray_direction.Z();

    double tmin = std::max(std::max(std::min(t1, t2), std::min(t3, t4)), std::min(t5, t6));
    double tmax = std::min(std::min(std::max(t1, t2), std::max(t3, t4)), std::max(t5, t6));

    if (tmin > tmax){
      //std::cout<<"Ray does not intersect the box"<<std::endl;
      return false;
    }
    if (tmax < 0){
      //std::cout<<"Ray path intersects but ray is aiming the wrong way"<<std::endl;
      return false;
    }


    //std::cout<<"Ray intersects the box"<<std::endl;
    return true;
  }


  DEFINE_ART_MODULE(GenFilter)

}

#include <iostream>
#include <algorithm>

#include "TGeoManager.h"

#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace filt{

  class GenFilter : public art::EDFilter {
    public:
      explicit GenFilter(fhicl::ParameterSet const & pset);
      virtual bool filter(art::Event& e) override;
      void reconfigure(fhicl::ParameterSet const& pset);
      virtual void beginJob() override;

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
      bool fUseReadoutWindow;

      geo::GeometryCore const* fGeometryService;
      detinfo::DetectorClocks const* fDetectorClocks;
      detinfo::DetectorProperties const* fDetectorProperties;
      double readoutWindow;
      double driftTime;

      bool IsInterestingParticle(const simb::MCParticle &particle);
      void LoadCRTAuxDetIDs();
      bool UsesCRTAuxDets(const simb::MCParticle &particle, const std::vector<unsigned int> &crt_auxdet_vector);
      bool UsesCRTAuxDet(const simb::MCParticle &particle, geo::AuxDetGeo const& crt);
      bool RayIntersectsBox(TVector3 ray_origin, TVector3 ray_direction, TVector3 box_min_extent, TVector3 box_max_extent);
      std::pair<double, double> XLimitsTPC(const simb::MCParticle &particle);
      std::pair<TVector3, TVector3> CubeIntersection(TVector3 min, TVector3 max, TVector3 start, TVector3 end);
  };


  GenFilter::GenFilter(fhicl::ParameterSet const & pset)
  : EDFilter(pset)
  {
    this->reconfigure(pset);
    
    fGeometryService = lar::providerFrom<geo::Geometry>();
    fDetectorClocks = lar::providerFrom<detinfo::DetectorClocksService>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    readoutWindow  = fDetectorClocks->TPCTick2Time((double)fDetectorProperties->ReadOutWindowSize()); // [us]
    driftTime = (2.*fGeometryService->DetHalfWidth())/fDetectorProperties->DriftVelocity(); // [us]
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
    fUseReadoutWindow = pset.get<bool>("UseReadoutWindow");
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

          double time = particle.T() * 1e-3; //[us]

          if (fUseReadoutWindow){
            if (time < -driftTime || time > readoutWindow) continue;
            // Get the minimum and maximum |x| position in the TPC
            std::pair<double, double> xLimits = XLimitsTPC(particle);
            // Calculate the expected time of arrival of those points
            double minTime = time + (2.0 * fGeometryService->DetHalfWidth() - xLimits.second)/fDetectorProperties->DriftVelocity();
            double maxTime = time + (2.0 * fGeometryService->DetHalfWidth() - xLimits.first)/fDetectorProperties->DriftVelocity();
            // If both times are below or above the readout window time then skip
            if((minTime < 0 && maxTime < 0) || (minTime > readoutWindow && maxTime > readoutWindow)) continue;
          }

          if (fUseTopHighCRTs){
            bool OK = UsesCRTAuxDets(particle, fTopHighCRTAuxDetIDs);
            if (!OK) continue;
            //std::cout<<"TopHighCRTs: " << OK << std::endl;
          }
          if (fUseTopLowCRTs){
            bool OK = UsesCRTAuxDets(particle, fTopLowCRTAuxDetIDs);
            if (!OK) continue;

            //std::cout<<"TopLowCRTs: " << OK << std::endl;
          }
          if (fUseBottomCRTs){
            bool OK = UsesCRTAuxDets(particle, fBottomCRTAuxDetIDs);
            if (!OK) continue;
            //std::cout<<"BottomCRTs: " << OK << std::endl;
          }
          if (fUseFrontCRTs){
            bool OK = UsesCRTAuxDets(particle, fFrontCRTAuxDetIDs);
            if (!OK) continue;
            //std::cout<<"FrontCRTs: " << OK << std::endl;
          }
          if (fUseBackCRTs){
            bool OK = UsesCRTAuxDets(particle, fBackCRTAuxDetIDs);
            if (!OK) continue;
            //std::cout<<"BackCRTs: " << OK << std::endl;
          }
          if (fUseLeftCRTs){
            bool OK = UsesCRTAuxDets(particle, fLeftCRTAuxDetIDs);
            if (!OK) continue;
            //std::cout<<"LeftCRTs: " << OK << std::endl;
          }
          if (fUseRightCRTs){
            bool OK = UsesCRTAuxDets(particle, fRightCRTAuxDetIDs);
            if (!OK) continue;
            //std::cout<<"RightCRTs: " << OK << std::endl;
          }

          return true;
        }
      }
    }

    return false;
  }


  void GenFilter::beginJob() {
    LoadCRTAuxDetIDs();
  }


  bool GenFilter::IsInterestingParticle(const simb::MCParticle &particle){
    double mom = particle.Momentum().Vect().Mag();
    int pdg = particle.PdgCode();

    for (unsigned int particle_i = 0; particle_i < fPDGs.size(); particle_i++){
      bool matched = true;
      //Check minimum momentum
      if (fMinMomentums[particle_i] > 0 && mom < fMinMomentums[particle_i]) matched = false;
      //Check maximum momentum
      if (fMaxMomentums[particle_i] > 0 && mom > fMaxMomentums[particle_i]) matched = false;
      //Check PDG
      if (fPDGs[particle_i]!=0 && pdg!=fPDGs[particle_i]) matched = false;
      if(matched) return true;
    }
    return false;
  }


  void GenFilter::LoadCRTAuxDetIDs(){
    art::ServiceHandle<geo::Geometry> geom;

    for (unsigned int auxdet_i = 0; auxdet_i < geom->NAuxDets(); auxdet_i++){
      geo::AuxDetGeo const& crt = geom->AuxDet(auxdet_i);
      const TGeoVolume* volModule = crt.TotalVolume();
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
      geo::AuxDetGeo const& crt = geom->AuxDet(auxdet_index);
      if (UsesCRTAuxDet(particle,crt)){
        return true;
      }
    }

    return false;
  }

  bool GenFilter::UsesCRTAuxDet(const simb::MCParticle &particle, geo::AuxDetGeo const& crt){
    //We need to prepare the particle's position and direction and construct a bounding box from the CRT for the ray-box collision algorithm
    //Start with the particle
    double particle_position_array[3], particle_direction_array[3], particle_local_position_array[3], particle_local_direction_array[3];
    particle_position_array[0] = particle.Position(0).X();
    particle_position_array[1] = particle.Position(0).Y();
    particle_position_array[2] = particle.Position(0).Z();
    particle_direction_array[0] = particle.Momentum(0).X()/particle.Momentum(0).Vect().Mag();
    particle_direction_array[1] = particle.Momentum(0).Y()/particle.Momentum(0).Vect().Mag();
    particle_direction_array[2] = particle.Momentum(0).Z()/particle.Momentum(0).Vect().Mag();
    crt.WorldToLocal(particle_position_array,particle_local_position_array);
    crt.WorldToLocalVect(particle_direction_array,particle_local_direction_array);
    TVector3 particle_local_position(particle_local_position_array[0],particle_local_position_array[1],particle_local_position_array[2]);
    TVector3 particle_local_direction(particle_local_direction_array[0],particle_local_direction_array[1],particle_local_direction_array[2]);
    double norm[3], lnorm[3];
    double center[3];
    crt.GetNormalVector(norm);
    crt.WorldToLocalVect(norm,lnorm);
    crt.GetCenter(center);

    //Now make the bounding box min and max extent from the CRT
    //In local coordinates, the normal of the CRT is parallel to the z-axis, the length is parallel to z, width parallel to x and height parallel to y
    //A gotchya: AuxDets pass half widths and half heights but FULL lengths.  So, divide length by 2
    TVector3 crt_min_extent(-1.*crt.HalfWidth1(),-1.*crt.HalfHeight(),-1.*crt.Length()/2.);
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

  std::pair<double, double> GenFilter::XLimitsTPC(const simb::MCParticle &particle){

    TVector3 start (particle.Position().X(), particle.Position().Y(), particle.Position().Z());
    TVector3 direction (particle.Momentum().X()/particle.Momentum().Vect().Mag(),
                        particle.Momentum().Y()/particle.Momentum().Vect().Mag(),
                        particle.Momentum().Z()/particle.Momentum().Vect().Mag());
    TVector3 end = start + direction;

    double minimum = 99999;
    double maximum = -99999;
    for(size_t c = 0; c < fGeometryService->Ncryostats(); c++){
      const geo::CryostatGeo& cryostat = fGeometryService->Cryostat(c);
      for(size_t t = 0; t < cryostat.NTPC(); t++){
        const geo::TPCGeo& tpcGeo = cryostat.TPC(t);
        // Find the intersection between the track and the TPC
        TVector3 min (tpcGeo.MinX(), tpcGeo.MinY(), tpcGeo.MinZ());
        TVector3 max (tpcGeo.MaxX(), tpcGeo.MaxY(), tpcGeo.MaxZ());
    
        std::pair<TVector3, TVector3> intersection = CubeIntersection(min, max, start, end);
        double x1 = std::abs(intersection.first.X());
        double x2 = std::abs(intersection.second.X());
        if(x1 != -99999 && x1 > maximum) maximum = x1;
        if(x1 != -99999 && x1 < minimum) minimum = x1;
        if(x2 != -99999 && x2 > maximum) maximum = x2;
        if(x2 != -99999 && x2 < minimum) minimum = x2;
      }
    }

    return std::make_pair(minimum, maximum);
  }

  std::pair<TVector3, TVector3> GenFilter::CubeIntersection(TVector3 min, TVector3 max, TVector3 start, TVector3 end){

    TVector3 dir = (end - start);
    TVector3 invDir (1./dir.X(), 1./dir.Y(), 1/dir.Z());
 
    double tmin, tmax, tymin, tymax, tzmin, tzmax;
 
    TVector3 enter (-99999, -99999, -99999);
    TVector3 exit (-99999, -99999, -99999);
 
    // Find the intersections with the X plane
    if(invDir.X() >= 0){
      tmin = (min.X() - start.X()) * invDir.X();
      tmax = (max.X() - start.X()) * invDir.X();
    }
    else{
      tmin = (max.X() - start.X()) * invDir.X();
      tmax = (min.X() - start.X()) * invDir.X();
    }
 
    // Find the intersections with the Y plane
    if(invDir.Y() >= 0){
      tymin = (min.Y() - start.Y()) * invDir.Y();
      tymax = (max.Y() - start.Y()) * invDir.Y();
    }
    else{
      tymin = (max.Y() - start.Y()) * invDir.Y();
      tymax = (min.Y() - start.Y()) * invDir.Y();
    }
 
    // Check that it actually intersects
    if((tmin > tymax) || (tymin > tmax)) return std::make_pair(enter, exit);
 
    // Max of the min points is the actual intersection
    if(tymin > tmin) tmin = tymin;
 
    // Min of the max points is the actual intersection
    if(tymax < tmax) tmax = tymax;
 
    // Find the intersection with the Z plane
    if(invDir.Z() >= 0){
      tzmin = (min.Z() - start.Z()) * invDir.Z();
      tzmax = (max.Z() - start.Z()) * invDir.Z();
    }
    else{
      tzmin = (max.Z() - start.Z()) * invDir.Z();
      tzmax = (min.Z() - start.Z()) * invDir.Z();
    }
 
    // Check for intersection
    if((tmin > tzmax) || (tzmin > tmax)) return std::make_pair(enter, exit);
 
    // Find final intersection points
    if(tzmin > tmin) tmin = tzmin;
 
    // Find final intersection points
    if(tzmax < tmax) tmax = tzmax;
 
    // Calculate the actual crossing points
    double xmin = start.X() + tmin * dir.X();
    double xmax = start.X() + tmax * dir.X();
    double ymin = start.Y() + tmin * dir.Y();
    double ymax = start.Y() + tmax * dir.Y();
    double zmin = start.Z() + tmin * dir.Z();
    double zmax = start.Z() + tmax * dir.Z();
 
    // Return pair of entry and exit points
    enter.SetXYZ(xmin, ymin, zmin);
    exit.SetXYZ(xmax, ymax, zmax);
    return std::make_pair(enter, exit);

  }

  DEFINE_ART_MODULE(GenFilter)

}

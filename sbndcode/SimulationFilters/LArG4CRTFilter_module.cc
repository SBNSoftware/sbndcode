#include <iostream>
#include <algorithm>

#include "TGeoManager.h"
#include "TVector3.h"

#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace filt{

  class LArG4CRTFilter : public art::EDFilter {
    public:
      explicit LArG4CRTFilter(fhicl::ParameterSet const & pset);
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
      std::string fLArG4ModuleName;
      bool fUseReadoutWindow;
      bool fUseTPC;

      geo::GeometryCore const* fGeometryService;
      double readoutWindow;
      double driftTime;
      
      bool IsInterestingParticle(const art::Ptr<simb::MCParticle> particle);
      void LoadCRTAuxDetIDs();
      bool UsesCRTAuxDets(const art::Ptr<simb::MCParticle> particle, const std::vector<unsigned int> &crt_auxdet_vector);
      bool EntersTPC(const art::Ptr<simb::MCParticle> particle);
      std::pair<double, double> XLimitsTPC(const art::Ptr<simb::MCParticle> particle);
  };


  LArG4CRTFilter::LArG4CRTFilter(fhicl::ParameterSet const & pset)
  : EDFilter(pset)
  {
    this->reconfigure(pset);

    fGeometryService = lar::providerFrom<geo::Geometry>();
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);
    readoutWindow  = clockData.TPCTick2Time((double)detProp.ReadOutWindowSize()); // [us]
    driftTime = (2.*fGeometryService->DetHalfWidth())/detProp.DriftVelocity(); // [us]
  }


  void LArG4CRTFilter::reconfigure(fhicl::ParameterSet const& pset){
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
    fLArG4ModuleName = pset.get<std::string>("LArG4ModuleName");
    fUseReadoutWindow = pset.get<bool>("UseReadoutWindow");
    fUseTPC = pset.get<bool>("UseTPC");
  }


  bool LArG4CRTFilter::filter(art::Event & e){
    art::Handle<std::vector<simb::MCParticle> > particles;
    e.getByLabel(fLArG4ModuleName,particles);

    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

    for (unsigned int part_i = 0; part_i < particles->size(); part_i++){
      const art::Ptr<simb::MCParticle> particle(particles,part_i);
      if (!IsInterestingParticle(particle)) continue;
      //std::cout<<"particlePDG: " << particle->PdgCode() << std::endl;
      //particle->Position().Vect().Print();
      //particle->Momentum().Vect().Print();
      double time = particle->T() * 1e-3; //[us]
      if (fUseReadoutWindow){
        if (time <= -driftTime || time >= readoutWindow) continue;
        // Get the minimum and maximum |x| position in the TPC
        std::pair<double, double> xLimits = XLimitsTPC(particle);
        // Calculate the expected time of arrival of those points
        double minTime = time + (2.0 * fGeometryService->DetHalfWidth() - xLimits.second)/detProp.DriftVelocity();
        double maxTime = time + (2.0 * fGeometryService->DetHalfWidth() - xLimits.first)/detProp.DriftVelocity();
        // If both times are below or above the readout window time then skip
        if((minTime < 0 && maxTime < 0) || (minTime > readoutWindow && maxTime > readoutWindow)) continue;
      }
      if (fUseTPC && !EntersTPC(particle)) continue;
      if (fUseTopHighCRTs){
        bool OK = UsesCRTAuxDets(particle,fTopHighCRTAuxDetIDs);
        if (!OK) continue;
        //std::cout<<"TopHighCRTs: " << OK << std::endl;
      }
      if (fUseTopLowCRTs){
        bool OK = UsesCRTAuxDets(particle,fTopLowCRTAuxDetIDs);
        if (!OK) continue;
        //std::cout<<"TopLowCRTs: " << OK << std::endl;
      }
      if (fUseBottomCRTs){
        bool OK = UsesCRTAuxDets(particle,fBottomCRTAuxDetIDs);
        if (!OK) continue;
        //std::cout<<"BottomCRTs: " << OK << std::endl;
      }
      if (fUseFrontCRTs){
        bool OK = UsesCRTAuxDets(particle,fFrontCRTAuxDetIDs);
        if (!OK) continue;
        //std::cout<<"FrontCRTs: " << OK << std::endl;
      }
      if (fUseBackCRTs){
        bool OK = UsesCRTAuxDets(particle,fBackCRTAuxDetIDs);
        if (!OK) continue;
        //std::cout<<"BackCRTs: " << OK << std::endl;
      }
      if (fUseLeftCRTs){
        bool OK = UsesCRTAuxDets(particle,fLeftCRTAuxDetIDs);
        if (!OK) continue;
        //std::cout<<"LeftCRTs: " << OK << std::endl;
      }
      if (fUseRightCRTs){
        bool OK = UsesCRTAuxDets(particle,fRightCRTAuxDetIDs);
        if (!OK) continue;
        //std::cout<<"RightCRTs: " << OK << std::endl;
      }
      //std::cout<<"Particle hit all CRTs"<<std::endl;
      return true;
    }

    return false;
  }


  void LArG4CRTFilter::beginJob() {
    LoadCRTAuxDetIDs();
  }


  bool LArG4CRTFilter::IsInterestingParticle(const art::Ptr<simb::MCParticle> particle){
    double mom = particle->Momentum().Vect().Mag();
    int pdg = particle->PdgCode();

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


  void LArG4CRTFilter::LoadCRTAuxDetIDs(){
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


  bool LArG4CRTFilter::UsesCRTAuxDets(const art::Ptr<simb::MCParticle> particle, const std::vector<unsigned int> &crt_auxdet_vector){
    //Loop over the aux dets, extract each one and then perform the test
    art::ServiceHandle<geo::Geometry> geom;
    //art:: ServiceHandle <geo:: AuxDetGeometry > adGeoService;
    //const  geo:: AuxDetGeometry* adG = &(* adGeoService);   // !!
    //const  geo:: AuxDetGeometryCore* adGeoCore = adG ->GetProviderPtr ();

    //Loop over the trajectory points made by this particle
    for (unsigned int pt_i = 0; pt_i < particle->NumberTrajectoryPoints(); pt_i++){
      //Get the position of the particle
      TLorentzVector position_lvector = particle->Position(pt_i);
      //We are going to use a function in the geometry service to see if there is a CRT at this particular position.  Said function requires an array, lets make one
      double position[3];
      position[0] = position_lvector.X();
      position[1] = position_lvector.Y();
      position[2] = position_lvector.Z();
      //The find the auxdet function throws a wobbler (an exception) if it can't find an auxdet.  Wrap what we want to do in a try catch statement pair
      try{
        unsigned int crt_id = geom->FindAuxDetAtPosition(position);
        //size_t adID , svID;
        //adGeoCore ->PositionToAuxDetChannel(position , adID , svID);
        //So we got an ID.  Lets compare it to all of the CRT IDs we are interested in
        for (unsigned int crt_i = 0; crt_i < crt_auxdet_vector.size(); crt_i++){
          if (crt_id == crt_auxdet_vector[crt_i]){
          //if (adID == crt_auxdet_vector[crt_i]){
            /*
            std::cout<<"Hit a CRT we wanted: DUMPING" <<std::endl;
            for (unsigned int i = 0; i < particle->NumberTrajectoryPoints(); i++){
              particle->Position(i).Print();
            }
            */
            //We found a CRT that we are interested in at this poisition!!!
            //We can leave now
            return true;
          }
        }
      }
      catch(...){}; //no CRT here, move along
    }
    return false;
  }

  bool LArG4CRTFilter::EntersTPC(const art::Ptr<simb::MCParticle> particle){
    bool enters = false;
    double xmin = -2.0 * fGeometryService->DetHalfWidth();
    double xmax = 2.0 * fGeometryService->DetHalfWidth();
    double ymin = -fGeometryService->DetHalfHeight();
    double ymax = fGeometryService->DetHalfHeight();
    double zmin = 0.;
    double zmax = fGeometryService->DetLength();

    int nTrajPoints = particle->NumberTrajectoryPoints();
    for (int traj_i = 0; traj_i < nTrajPoints; traj_i++){
      TVector3 trajPoint(particle->Vx(traj_i), particle->Vy(traj_i), particle->Vz(traj_i));
      // Check if point is within reconstructable volume
      if (trajPoint[0] >= xmin && trajPoint[0] <= xmax && trajPoint[1] >= ymin && trajPoint[1] <= ymax && trajPoint[2] >= zmin && trajPoint[2] <= zmax){
        enters = true;
      }
    }
    return enters;
  }

  std::pair<double, double> LArG4CRTFilter::XLimitsTPC(const art::Ptr<simb::MCParticle> particle){
    double xmin = -2.0 * fGeometryService->DetHalfWidth();
    double xmax = 2.0 * fGeometryService->DetHalfWidth();
    double ymin = -fGeometryService->DetHalfHeight();
    double ymax = fGeometryService->DetHalfHeight();
    double zmin = 0.;
    double zmax = fGeometryService->DetLength();

    double minimum = 99999;
    double maximum = -99999;

    int nTrajPoints = particle->NumberTrajectoryPoints();
    for (int traj_i = 0; traj_i < nTrajPoints; traj_i++){
      TVector3 trajPoint(particle->Vx(traj_i), particle->Vy(traj_i), particle->Vz(traj_i));
      // Check if point is within reconstructable volume
      if (trajPoint[0] >= xmin && trajPoint[0] <= xmax && trajPoint[1] >= ymin && trajPoint[1] <= ymax && trajPoint[2] >= zmin && trajPoint[2] <= zmax){
        if(std::abs(trajPoint[0]) < minimum) minimum = std::abs(trajPoint[0]);
        if(std::abs(trajPoint[0]) > maximum) maximum = std::abs(trajPoint[0]);
      }
    }
    return std::make_pair(minimum, maximum);
  }


  DEFINE_ART_MODULE(LArG4CRTFilter)

}

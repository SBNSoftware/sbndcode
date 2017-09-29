#include <iostream>
#include <algorithm>

#include "TGeoManager.h"

#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace filt{

  class LArG4CRTFilter : public art::EDFilter {
    public:
      explicit LArG4CRTFilter(fhicl::ParameterSet const & pset);
      virtual bool filter(art::Event& e) override;
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
      std::string fLArG4ModuleName;

      bool IsInterestingParticle(const art::Ptr<simb::MCParticle> particle);
      void LoadCRTAuxDetIDs();
      bool UsesCRTAuxDets(const art::Ptr<simb::MCParticle> particle, const std::vector<unsigned int> &crt_auxdet_vector);
  };


  LArG4CRTFilter::LArG4CRTFilter(fhicl::ParameterSet const & pset)
  {
    this->reconfigure(pset);
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
  }


  bool LArG4CRTFilter::filter(art::Event & e){
    art::Handle<std::vector<simb::MCParticle> > particles;
    e.getByLabel(fLArG4ModuleName,particles);

    for (unsigned int part_i = 0; part_i < particles->size(); part_i++){
      const art::Ptr<simb::MCParticle> particle(particles,part_i);
      if (!IsInterestingParticle(particle)) continue;
      //std::cout<<"particlePDG: " << particle->PdgCode() << std::endl;
      //particle->Position().Vect().Print();
      //particle->Momentum().Vect().Print();
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
        //So we got an ID.  Lets compare it to all of the CRT IDs we are interested in
        for (unsigned int crt_i = 0; crt_i < crt_auxdet_vector.size(); crt_i++){
          if (crt_id == crt_auxdet_vector[crt_i]){
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

  DEFINE_ART_MODULE(LArG4CRTFilter)

}

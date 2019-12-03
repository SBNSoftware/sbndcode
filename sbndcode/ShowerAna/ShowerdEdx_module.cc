#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//Larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

//Root Includes
#include "TMath.h"
#include "TTree.h"
#include "TRandom3.h"

//C++ Includes
#include <vector>
#include <iostream>


namespace ana {
  class ShowerdEdx;
}

class ana::ShowerdEdx: public art::EDAnalyzer {
  public:

    ShowerdEdx(const fhicl::ParameterSet& pset);

    void analyze(const art::Event& evt);
    void beginJob();

  private:
    //fcl Parameters
    std::string fGenieGenModuleLabel;
    std::string fLArGeantModuleLabel;
    std::string fMCShowerTag;

    float fTrackStubLength;
    int fGaussWidths;

    art::ServiceHandle<geo::Geometry> fGeom;
    art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
    art::ServiceHandle<art::TFileService> tfs;

    TTree* dEdxTree;
    float dEdx; // From MCShower
    float energy;
    float pdg;

    int trackID;
    float startX;
    float startY;
    float startZ;

    int numIDEs;
    int numShowers;
    int numMCShowers;
    int numParticles;

    int eventRun;
    int eventSubRun;
    int eventNumber;

    float thetaXZ;
    float thetaYZ;

    int width;
    float smear;

    std::string endProcess;

    float newdEdx; // Calculated in this module
    float dEdxMedian;
    std::vector<float> dEdxVector;
};


ana::ShowerdEdx::ShowerdEdx(const fhicl::ParameterSet& pset) : EDAnalyzer(pset){

  fGenieGenModuleLabel   = pset.get<std::string>("GenieGenModuleLabel");
  fLArGeantModuleLabel   = pset.get<std::string>("LArGeantModuleLabel");
  fMCShowerTag           = pset.get<std::string>("MCShowerTag");

  fTrackStubLength       = pset.get<float>("TrackStubLength");
  fGaussWidths           = pset.get<int>("GaussWidths");
}

void ana::ShowerdEdx::beginJob(){

  dEdxTree = tfs->make<TTree>("dEdxTree","dEdxTree");
  dEdxTree->Branch("dEdx",&dEdx,"dEdx/F");
  dEdxTree->Branch("energy",&energy,"energy/F");
  dEdxTree->Branch("pdg",&pdg,"pdg/F");

  dEdxTree->Branch("trackID",&trackID,"trackID/I");
  dEdxTree->Branch("startX",&startX,"startX/F");
  dEdxTree->Branch("startY",&startY,"startY/F");
  dEdxTree->Branch("startZ",&startZ,"startZ/F");

  dEdxTree->Branch("numIDEs",&numIDEs,"numIDEs/I");
  dEdxTree->Branch("numShowers",&numShowers,"numShowers/I");
  dEdxTree->Branch("numMCShowers",&numMCShowers,"numMCShowers/I");
  dEdxTree->Branch("numParticles",&numParticles,"numParticles/I");

  dEdxTree->Branch("eventRun",&eventRun,"eventRun/I");
  dEdxTree->Branch("eventSubRun",&eventSubRun,"eventSubRun/I");
  dEdxTree->Branch("eventNumber",&eventNumber,"eventNumber/I");

  dEdxTree->Branch("thetaXZ",&thetaXZ,"thetaXZ/F");
  dEdxTree->Branch("thetaYZ",&thetaYZ,"thetaYZ/F");

  dEdxTree->Branch("endProcess",&endProcess);
  dEdxTree->Branch("dEdxVector",&dEdxVector);
  dEdxTree->Branch("dEdxMedian",&dEdxMedian);

  dEdxTree->Branch("newdEdx",&newdEdx,"newdEdx/F");
  dEdxTree->Branch("width",&width,"width/I");
  dEdxTree->Branch("smear",&smear,"smear/F");

  gRandom = new TRandom3;

}

void ana::ShowerdEdx::analyze(const art::Event& evt){

  eventRun    = evt.run();
  eventSubRun = evt.subRun();
  eventNumber = evt.event();

  //Getting  MC truth information

  // Get all of the simChannels so we can access the IDEs
  art::Handle<std::vector<sim::SimChannel> > simChannelHandle;
  std::vector<art::Ptr<sim::SimChannel> > simchannels;
  if(evt.getByLabel(fLArGeantModuleLabel,simChannelHandle))
  {art::fill_ptr_vector(simchannels, simChannelHandle);}

  numShowers   = 0;
  numMCShowers = 0;
  numParticles = 0;

  // Create a map of MCShower track IDs to dEdx values, for comparison to MCParticle
  std::map<int,float> showerMap;
  auto const& mcshowers = *evt.getValidHandle<std::vector<sim::MCShower> >(fMCShowerTag);

  for (auto const &mcs: mcshowers) {
    if ((TMath::Abs(mcs.PdgCode())==11) || (TMath::Abs(mcs.PdgCode())==22)){
      int mcsTrackID = mcs.TrackID();
      float mcsdEdx    = mcs.dEdx();

      showerMap[mcsTrackID] = mcsdEdx;
      numMCShowers++;
    }
  }

  // Create a map of trackIDs to MCParticles
  std::map<int,const simb::MCParticle*> particleMap;
  const sim::ParticleList& particles = particleInventory->ParticleList();
  for (sim::ParticleList::const_iterator particleIt = particles.begin();
      particleIt != particles.end(); ++particleIt){
    const simb::MCParticle *particle = particleIt->second;
    particleMap[particle->TrackId()] = particle;
    numParticles++;
    std::cout << "TrackID: " << particle->TrackId() << " Pdg: " << particle->PdgCode() << " E: " << particle->E() << std::endl;
    if (TMath::Abs(particle->PdgCode()) == 11 || TMath::Abs(particle->PdgCode()) == 22) {
      numShowers++;
    }
  }

  // Look over the MCShowers
  for (std::map<int,float>::iterator showerIt = showerMap.begin();
      showerIt != showerMap.end(); showerIt++){

    // reinitialise tree values for each shower
    dEdx     = -99999;
    energy   = -99999;
    pdg      = -99999;

    trackID  = -99999;
    startX   = -99999;
    startY   = -99999;
    startZ   = -99999;

    thetaXZ  = -99999;
    thetaYZ  = -99999;

    numIDEs  = 0;
    newdEdx  = -99999;
    dEdxVector.clear();

    // Get the MCParticle that corresponds to the MCShower
    const simb::MCParticle *particle = particleMap.at((*showerIt).first);
    // Get MCParticle truth info for the TTree
    dEdx       = (*showerIt).second;
    energy     = particle->E();
    pdg        = particle->PdgCode();
    trackID    = particle->TrackId();
    endProcess = particle->EndProcess();

    std::vector<int> daughters;

    if(particle->Mother() != 0){continue;}

    // Find the start position of the shower for electrons
    TLorentzVector startPos;
    if (TMath::Abs(particle->PdgCode()) == 11) {
      startPos = particle->Position();
      daughters.push_back(trackID);
    } else if (TMath::Abs(particle->PdgCode()) == 22) {
      
      // Check for shower roll-up
      int numDaughters = particle->NumberDaughters();
      //      If shower is rolled up, take the photon ID
      if (numDaughters==0) {
        daughters.push_back(trackID);
      } else {
        // If the photon pair produces, take both daughters
        if (particle->EndProcess()=="conv"){
	  for (int daughterNum=0; daughterNum<numDaughters; daughterNum++){
            daughters.push_back(particle->Daughter(daughterNum));
	  
          }
        } else { // otherwise (compton) take first daughter only (big electron?)
          daughters.push_back(particle->Daughter(0));
	          }
      }

      std::cout << "particle->EndProcess(): " << particle->EndProcess() << std::endl;

      //For no only keep conv events where poth particles go past the lenght cut 
      if(daughters.size() < 2){std::cout << "daughters not enough" << std::endl;return;}

      //Only keep daughters which pass the length cut
      for(auto const& daughter_id: daughters){
	if(particleMap.find(daughter_id) == particleMap.end()){continue;}
	const simb::MCParticle *daughter = particleMap.at(daughter_id);
	int numTrajPoints = daughter->NumberTrajectoryPoints();
	TLorentzVector startPos_daughter = daughter->Position(0);
	TLorentzVector endPos_daughter = daughter->Position(numTrajPoints-1);
	for (int i=0; i< numTrajPoints; i++){
	  if ((daughter->E(i)< 0.01)){
	    endPos_daughter = daughter->Position(i);
	  }
	}
	if((startPos_daughter.Vect()-endPos_daughter.Vect()).Mag() < fTrackStubLength){std::cout << "track length not big enough" << std::endl; return;}
      }


      // For photons we want the end of the particle as this is
      // when the photon interacts and the showers start
      int numTrajPoints = particle->NumberTrajectoryPoints();
      bool firstPoint = true;
      for (int i=0; i< numTrajPoints; i++){
        if ((particle->E(i)<0.9*energy) && firstPoint){
          startPos = particle->Position(i);
          firstPoint = false;
        }
      }
    } else {
      std::cout<<"PDG: "<<particle->PdgCode()<<", Not a shower!!"<<std::endl;
      return;
    }

    startX = startPos.Px();
    startY = startPos.Py();
    startZ = startPos.Pz();

    // Calculate the angle of the track
    TVector3 pVector(particle->Px(),particle->Py(),particle->Pz());
    TVector3 pVectorNorm = pVector.Unit();
    float pX = pVectorNorm.X();
    float pY = pVectorNorm.Y();
    float pZ = pVectorNorm.Z();

    thetaXZ = TMath::ATan(pX/pZ);
    thetaYZ = TMath::ATan(pY/pZ);

    float dESum   = 0;

    // Loop over the simChannels to get the IDEs
    std::map<geo::PlaneID, std::vector<float> > planedEdxMap;
    
    for (const auto channelIt: simchannels){

      // std::cout<<"Channel: "<<channelIt->Channel()<<std::endl;

      std::map<int,int> daughterIDEs;
      std::map<int,float> daughterSum;
    
      for (const auto daughterID: daughters){ // loop over shower daughters
        daughterIDEs[daughterID] = 0;
      }

      int wireIDEs  = 0;
      float wireSum = 0;

      // Get a map of TDCs to IDEs, then loop over the IDEs
      auto tdc_ide_map = channelIt->TDCIDEMap();
      // std::cout<<"TDC IDE size: "<<tdc_ide_map.size()<<std::endl;
      for(auto const& tdc_ide_pair : tdc_ide_map) {
        // std::cout<<"Time: "<<tdc_ide_pair.first << std::endl;
        auto const& ide_v = tdc_ide_pair.second;
        // std::cout<<"IDE vector size: "<<ide_v.size()<<std::endl;
        for(auto const& ide : ide_v) {

	  if(particleMap.find(ide.trackID) == particleMap.end()){continue;}
	  if(TMath::Abs(particleMap[ide.trackID]->PdgCode()) != 11){continue;}
	  
	  //	  for (const auto daughterID: daughters){ // loop over shower daughters
	    

            // Calculate distances from ides to MCParticle starting position
            TVector3 ds((ide.x - startX), (ide.y - startY), (ide.z - startZ));

            // Only use IDEs within the distance cut of true vertex
            // and along the direction of the particle
            if (ds.Mag() < fTrackStubLength && ds.Dot(pVectorNorm) > 0){

              //if (TMath::Abs(ide.trackID) == TMath::Abs(daughterID)){
	      wireSum += ide.energy;
	      wireIDEs++;
	      daughterIDEs[ide.trackID]++;
	      //              }§
            }// daughters§
	} // ide_v
      } // tdc_ide_pair
    


     // tdc_ide_map
    
      //check all the daughters contributed to the wire
      bool allDaughtersHit = true;
      for (const auto daughterID: daughters){ // loop over shower daughters
	if (daughterIDEs[daughterID]==0){
	  allDaughtersHit = false;
        }
      }
      if(allDaughtersHit == false){continue;}

      if (wireIDEs>0 && allDaughtersHit){
        numIDEs += wireIDEs;
	dESum += wireSum;

        //Calculate the pitch
        std::vector<geo::WireID> wires = fGeom->ChannelToWire(channelIt->Channel());
        geo::PlaneID planeID = wires.front().planeID();
        float wirePitch = fGeom->WirePitch(planeID);
        float angleToVert = fGeom->WireIDToWireGeo(wires.front()).ThetaZ() - 0.5*TMath::Pi();
        float cosGamma = std::abs(sin(angleToVert)*pY+cos(angleToVert)*pZ);

        float pitch = wirePitch/cosGamma;
        float wiredEdx = wireSum/pitch;

        planedEdxMap[planeID].push_back(wiredEdx);
      }
    }

    if (numIDEs==0){
      std::cout<<"No IDEs, Returning"<<std::endl;
      return;
    }

    geo::PlaneID bestPlane;
    unsigned int maxHits = 0;
    for (const auto planedEdx: planedEdxMap){
      // std::cout<<planedEdx.first.toString()<<" and "<<planedEdx.second.size()<<std::endl;
      if (planedEdx.second.size()>maxHits){
        bestPlane = planedEdx.first;
        maxHits = planedEdx.second.size();
      }
    }
    dEdxVector = planedEdxMap[bestPlane];

    // std::cout<<"Best Plane: "<<bestPlane.toString()<<" with hits: "<<maxHits<<std::endl;
    // Check we have some IDEs within cut, else return -99999

    std::cout<<"Vector Size: "<<dEdxVector.size()<<" median "<<dEdxMedian<<std::endl;

    //Calculate the median without the first and last 2 points
    if (maxHits>4){
      dEdxVector.erase(dEdxVector.begin(), dEdxVector.begin()+2);
      dEdxVector.erase(dEdxVector.end()-2, dEdxVector.end());
      dEdxMedian = TMath::Median(dEdxVector.size(), &dEdxVector[0]);
    } else {
      dEdxVector.clear();
      dEdxMedian = -999;
    }


    // for (const auto dEdx:dEdxVector){
    //   std::cout<<dEdx<<std::endl;
    // }

    newdEdx = dESum/fTrackStubLength; // Calculate dE/dx
    newdEdx = newdEdx/3; // We are triple counting planes in the IDEs

    float dEdxNoSmear = newdEdx;
    for (int gaussWidth=0; gaussWidth<fGaussWidths; gaussWidth++){
      width = gaussWidth;
      smear = gRandom->Gaus(1,(float)width/100.);
      newdEdx = dEdxNoSmear*smear;

      std::cout<<"End process: "<<endProcess<<" Vector Size: "<<dEdxVector.size()<<" median "<<dEdxMedian<<std::endl;

      dEdxTree->Fill();
    }
  }// loop over mcshowers
}
DEFINE_ART_MODULE(ana::ShowerdEdx)

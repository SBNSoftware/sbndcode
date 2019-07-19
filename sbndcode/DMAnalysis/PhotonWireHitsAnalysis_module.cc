//

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/BeamInfo.h"

#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"

#include "sbndcode/RecoUtils/RecoUtils.h"

#include <cstring>
#include <vector>
#include <map>
#include <iterator>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional>
#include <typeinfo>

#include "TTree.h"
#include "TTimeStamp.h"


namespace sbndwire {
  class PhotonWireHitsAnalysis; 
}


class sbndwire::PhotonWireHitsAnalysis : public art::EDAnalyzer {
public:
  explicit PhotonWireHitsAnalysis(fhicl::ParameterSet const & p);

  PhotonWireHitsAnalysis(PhotonWireHitsAnalysis const &) = delete;
  PhotonWireHitsAnalysis(PhotonWireHitsAnalysis &&) = delete;
  PhotonWireHitsAnalysis & operator = (PhotonWireHitsAnalysis const &) = delete;
  PhotonWireHitsAnalysis & operator = (PhotonWireHitsAnalysis &&) = delete;

  void analyze(art::Event const & e) override;
  void beginJob();
  void endJob();

private:

  TTree *fOutputTree;
  TTree *fEventTree;

  unsigned int fEvent;
  std::string fTruthLabel;
  std::string fHitProducerLabel;
  float fhit_RMS;
  float fhit_PeakAmplitude;
  float fhit_Integral;
  int fdistancetotwohits;
  float ftheta;
  float fenergy;
  float favSigmaChange;

};

sbndwire::PhotonWireHitsAnalysis::PhotonWireHitsAnalysis(fhicl::ParameterSet const & p) : EDAnalyzer(p) {

  fTruthLabel = p.get<std::string>("TruthLabel");
  fHitProducerLabel = p.get<std::string>("HitLabel");

}

void sbndwire::PhotonWireHitsAnalysis::analyze(art::Event const & e) {

  fEvent = e.id().event();
  std::cout<< "NEW EVENT "<<fEvent<<std::endl;
  if (!e.isRealData()) {
    int hitcount(0);
    unsigned int wirenum(0);
    unsigned int firstwire(0);
    float peaktime(0);
    fdistancetotwohits = 0; 
    favSigmaChange = 0;
    //int tpc(0);

    art::ValidHandle<std::vector<simb::MCParticle>> mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel);
    art::ValidHandle<std::vector<recob::Hit>> hitsVec = e.getValidHandle<std::vector<recob::Hit>>(fHitProducerLabel);

    // Get Energy and Separation angle of electron pair
    fenergy = 0; ftheta = 0;
    

    if (mcParticles.isValid()) {
      for (unsigned int t=0; t<mcParticles->size(); ++t) { // Loop through truth particles
        const simb::MCParticle mcpart = mcParticles->at(t);

        if (mcpart.Process() == "primary" && mcpart.PdgCode() == 22) {
          fenergy = mcpart.E();
          ftheta = mcpart.Momentum().Theta();
        }
      }
    }
    

    if (fenergy>0.05) { // only do hit analysis if both electrons >50MeV
      float sigma1(0), sigma2(0);
      float splitwire(0);
      // Loop over hits and find which tpc has greatest number
      int hits0(0), hits1(0);
      unsigned int tpc;
      std::vector<recob::Hit> const& hits = *hitsVec;
      for (recob::Hit const& hit: hits) {
        if (hit.WireID().TPC==0) { hits0++; }
        if (hit.WireID().TPC==1) { hits1++; }
      }
      std::cout<<"hits tpc 0: "<<hits0<<" hits tpc 1: "<<hits1<<std::endl;
      if (hits0>hits1) { tpc=0; }
      else { tpc=1; }
      std::cout<<"Look in tpc "<<tpc<<std::endl<< std::endl;

      // Start analysis to find first wire with two hits
      auto const* geom = lar::providerFrom<geo::Geometry>();

      for (geo::PlaneID const& pid: geom->IteratePlaneIDs()) {
        std::cout<< "new plane "<< pid <<std::endl;
        //std::vector<recob::Hit> const& hits = *hitsVec; 
        bool twohits(false);
        for (recob::Hit const& hit: hits) { // Loop over hits
          
          if (hit.WireID().planeID() != pid) continue;

          if (hit.WireID().TPC == tpc && hit.WireID().Cryostat == 0  &&  hit.WireID().Plane == 2 && hit.PeakAmplitude()>35 && twohits==false) {
          // require amplitude of hits to be greater than 35 adc units
            std::cout << std::endl << "plane: " << hit.WireID().planeID() << " wire: " << hit.WireID().Wire << std::endl;
            std::cout << "rms: " << hit.RMS() << " integral: " << hit.Integral() << std::endl;

            if (tpc==0 && hit.WireID().TPC==1) {hitcount=0;} // first hit in second tpc
            tpc = hit.WireID().TPC;
  
            hitcount++;
            if (hitcount==1) { // define first wire number 
              std::cout<< "FIRST HIT "<<hit.WireID().Wire<<std::endl;
              firstwire = hit.WireID().Wire; 
              wirenum = hit.WireID().Wire;
              peaktime = hit.PeakTime();
              sigma1 = hit.SigmaPeakTime();
            }
  
            else { // not first wire
              if (hit.WireID().Wire == wirenum) { // if wire num is same as previous hit
                std::cout<<"TWO HITS ON SAME WIRE "<<wirenum<<std::endl;

                if ((hit.PeakTime()-peaktime)<20) {  // require this second hit to be within 20 ticks of first
                  fdistancetotwohits = wirenum-firstwire;
                  std::cout << "num wires "<< fdistancetotwohits<<std::endl;
                  twohits=true; 
                  splitwire = hit.WireID().Wire; 
                }
                else { std::cout<< "too far away"<<std::endl; }
              }
  
              else { // wire is not the same number as previous hit
                if (hit.WireID().Wire != wirenum+1) { 
                  std::cout<< "GAP IN HITS, no event"<<std::endl;
                  hitcount=1; firstwire = hit.WireID().Wire;
                  wirenum = hit.WireID().Wire;
                  peaktime = hit.PeakTime();
                  sigma1 = hit.SigmaPeakTime();
                  //break;
                }
                else { // wire must be wirenum+1 ie carry on to next one
                  wirenum = hit.WireID().Wire;
                  peaktime = hit.PeakTime(); 
                  std::cout<<"next hit, on "<<wirenum+1<<std::endl;
                }
              } // if only one hit on previous wire
            }
          }

          fhit_RMS = hit.RMS();
          fhit_PeakAmplitude = hit.PeakAmplitude();
          fhit_Integral = hit.Integral();
  
          fOutputTree->Fill(); // Fills once per hit
        }
      }

      for (recob::Hit const& hit: hits) {
        if (splitwire!=0) {
          if (hit.WireID().Wire == splitwire-1) {sigma2 = hit.SigmaPeakTime(); }
        }
      }
      if (sigma1!=0 && sigma2!=0) { favSigmaChange = (sigma2-sigma1)/fdistancetotwohits; }


    }
    if (fdistancetotwohits!=0) {
      fEventTree->Fill(); // Fills once per event
    }
  }
}

void sbndwire::PhotonWireHitsAnalysis::beginJob() {

  // use TFileService to create tree and then add branches
  art::ServiceHandle<art::TFileService> tfs;
  fOutputTree = tfs->make<TTree>("mytree","My Tree");
  fOutputTree->Branch("event",&fEvent,"event/i");
  fOutputTree->Branch("hit_RMS",&fhit_RMS,"hit_RMS");
  fOutputTree->Branch("hit_PeakAmplitude",&fhit_PeakAmplitude,"hit_PeakAmplitude");
  fOutputTree->Branch("hit_Integral",&fhit_Integral,"hit_Integral");


  fEventTree = tfs->make<TTree>("myeventtree","My Event Tree");
  fEventTree->Branch("event",&fEvent,"event/i");
  fEventTree->Branch("distancetotwohits",&fdistancetotwohits,"distancetotwohits/i");
  fEventTree->Branch("energy",&fenergy,"energy/f");
  fEventTree->Branch("theta",&ftheta,"theta/f");
  fEventTree->Branch("avSigmaChange",&favSigmaChange,"favSigmaChange/f");
}

void sbndwire::PhotonWireHitsAnalysis::endJob() {

}

DEFINE_ART_MODULE(sbndwire::PhotonWireHitsAnalysis)

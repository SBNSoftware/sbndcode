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
#include "TGraph.h"
#include "TTimeStamp.h"


namespace sbndwire {
  class WireHitsAnalysis; 
}


class sbndwire::WireHitsAnalysis : public art::EDAnalyzer {
public:
  explicit WireHitsAnalysis(fhicl::ParameterSet const & p);

  WireHitsAnalysis(WireHitsAnalysis const &) = delete;
  WireHitsAnalysis(WireHitsAnalysis &&) = delete;
  WireHitsAnalysis & operator = (WireHitsAnalysis const &) = delete;
  WireHitsAnalysis & operator = (WireHitsAnalysis &&) = delete;

  void analyze(art::Event const & e) override;
  void beginJob();
  void endJob();

private:

  TTree *fOutputTree;
  TTree *fEventTree;

  //TGraph *fRatioGraph1;
  //TGraph *fRatioGraph2;

  std::vector<double> fAmplitudeArea;
  std::vector<double> fWire;

  unsigned int fEvent;
  std::string fTruthLabel;
  std::string fHitProducerLabel;
  float fhit_RMS;
  float fhit_PeakAmplitude;
  float fhit_Integral;
  int fdistancetotwohits2;
  int fdistancetotwohits1;
  int fdistancetotwohits0;

  float fsepAngle;
  float fpairE;
  float ftheta_one;
  float ftheta_two;
  float favSigmaChange;
  int numpoints;

  float ftheta;
  float fenergy;
  double fAvChange;
  double fAvChange_firstfive;

  float ffirstfivechange_amparea;
  float ftotalchange_amparea;
  float fKE_high;
  float fKE_low;

};

sbndwire::WireHitsAnalysis::WireHitsAnalysis(fhicl::ParameterSet const & p) : EDAnalyzer(p) {

  fTruthLabel = p.get<std::string>("TruthLabel");
  fHitProducerLabel = p.get<std::string>("HitLabel");

}

void sbndwire::WireHitsAnalysis::analyze(art::Event const & e) {
  numpoints = 0;
  fEvent = e.id().event();
  std::cout<< "NEW EVENT "<<fEvent<<std::endl;
  if (!e.isRealData()) {
    int hitcount(0);
    unsigned int wirenum(0);
    unsigned int firstwire(0);
    float peaktime(0);
    fdistancetotwohits0 = -1; 
    fdistancetotwohits1 = -1;
    fdistancetotwohits2 = -1;
    fKE_high = -1; fKE_low = -1;
    ftheta_one = 0; ftheta_two = 0;
    favSigmaChange = 0;
    bool energyPass(false); int elecCounter(0);
    fenergy = 0; ftheta = 0;

    art::ValidHandle<std::vector<simb::MCParticle>> mcParticles = e.getValidHandle<std::vector<simb::MCParticle>>(fTruthLabel);
    art::ValidHandle<std::vector<recob::Hit>> hitsVec = e.getValidHandle<std::vector<recob::Hit>>(fHitProducerLabel);

    // Get Energy and Separation angle of electron pair
    float px_one(0), py_one(0), pz_one(0), mom_one(0), KE_one(0);
    float px_two(0), py_two(0), pz_two(0), mom_two(0), KE_two(0);
    fpairE = 0; fsepAngle = 0;
    

    if (mcParticles.isValid()) {
      for (unsigned int t=0; t<mcParticles->size(); ++t) { // Loop through truth particles
        const simb::MCParticle mcpart = mcParticles->at(t);

        // look for photons
        if (mcpart.Process() == "primary" && mcpart.PdgCode() == 22) {
          fenergy = mcpart.E();
          ftheta = mcpart.Momentum().Theta();
        }

        // look for electron pair
        if (mcpart.Process() == "primary" && mcpart.PdgCode() == 11) {
          elecCounter++;
          //KE.push_back(mcpart.E()-0.000511);
          //Vx.push_back(mcpart.Vx()); Vy.push_back(mcpart.Vy()); Vz.push_back(mcpart.Vz());

          if (elecCounter==1) {
            px_one = mcpart.Px(); py_one = mcpart.Py(); pz_one = mcpart.Pz();
            mom_one = mcpart.P();
            KE_one = mcpart.E()-0.000511;
            ftheta_one = mcpart.Momentum().Theta();
          }
          if (elecCounter==2) {
            px_two = mcpart.Px(); py_two = mcpart.Py(); pz_two = mcpart.Pz();
            mom_two = mcpart.P();
            KE_two = mcpart.E()-0.000511;
            ftheta_two = mcpart.Momentum().Theta();
          }
        }
      }

      if (elecCounter == 2) { //two electrons in the event
        double dotProd = px_one*px_two + py_one*py_two + pz_one*pz_two;
        fsepAngle = acos(dotProd/(mom_one*mom_two));
        fpairE = KE_one + KE_two;
        if (KE_one>0.05 && KE_two>0.05) { energyPass=true; }
        if (KE_one>KE_two) { fKE_high = KE_one; fKE_low = KE_two; }
        else { fKE_high = KE_two; fKE_low = KE_one; }
      }

    }

    if (energyPass || fenergy>0.1) { // only do hit analysis if both electrons >50MeV
       // or photon energy is greater than 100MeV

      // Loop over hits and find which tpc has greatest number
      int hits0(0), hits1(0);
      //float sigma1(0), sigma2(0);
      float splitwire(0);
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
        bool twohits0(false); bool twohits1(false); bool twohits2(false);
       
        // get vector of number of hits on each wire in plane 2
        // use this later to check that wires after the split position have at least 2 hits 
        unsigned int numHits[1666] = {};
        for (recob::Hit const& hit: hits) { // Loop over hits
          if (hit.WireID().TPC == tpc && hit.WireID().Cryostat == 0  && hit.PeakAmplitude()>35 && hit.WireID().Plane == 2) {
            unsigned int temp = hit.WireID().Wire;
            numHits[temp] = numHits[temp]+1; // increase counter for that wire by 1
          }
        }

        for (recob::Hit const& hit: hits) { // Loop over hits
          
          if (hit.WireID().planeID() != pid) continue;

          if (hit.WireID().TPC == tpc && hit.WireID().Cryostat == 0  && hit.PeakAmplitude()>35) {
            if (hit.WireID().Plane == 2 && twohits2==false) { // require amplitude of hits to be greater than 35 adc units
            std::cout << std::endl << "plane: " << hit.WireID().planeID() << " wire: " << hit.WireID().Wire << std::endl;
            std::cout << "rms: " << hit.RMS() << " integral: " << hit.Integral() << std::endl;

            hitcount++;
            if (hitcount==1) { // define first wire number 
              std::cout<< "FIRST HIT "<<hit.WireID().Wire<<std::endl;
              firstwire = hit.WireID().Wire; 
              wirenum = hit.WireID().Wire;
              peaktime = hit.PeakTime();
            }
  
            else { // not first wire
              if (hit.WireID().Wire == wirenum) { // if wire num is same as previous hit
                std::cout<<"TWO HITS ON SAME WIRE "<<wirenum<<std::endl;

                if ((hit.PeakTime()-peaktime)<20) {  // require this second hit to be within 20 ticks of first

                  if (numHits[wirenum+1]>=2 && numHits[wirenum+2]>=2 && numHits[wirenum+3]>=2 && numHits[wirenum+4]>=2) {

                    fdistancetotwohits2 = wirenum-firstwire;
                    std::cout << "num wires "<< fdistancetotwohits2<<std::endl;
                    twohits2=true;
                    hitcount=0;   
                    splitwire = hit.WireID().Wire;
                  }
                  else {std::cout<<"failed check for continuous second hit. noise?"<<std::endl;}
                }
                else { std::cout<< "too far away"<<std::endl; }
              }
  
              else { // wire is not the same number as previous hit
                if (hit.WireID().Wire != wirenum+1) { 
                  std::cout<< "GAP IN HITS, no event"<<std::endl;
                  hitcount=1; firstwire = hit.WireID().Wire; // this becomes the new 'first hit'
                  wirenum = hit.WireID().Wire;
                  peaktime = hit.PeakTime();
                }
                else { // wire must be wirenum+1 ie carry on to next one
                  wirenum = hit.WireID().Wire;
                  peaktime = hit.PeakTime(); 
                  std::cout<<"next hit, on "<<wirenum+1<<std::endl;
                }
              } // if only one hit on previous wire
            }
            }

            if (hit.WireID().Plane == 1 && twohits1==false) { // require amplitude of hits to be greater than 35 adc units
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
              }
  
              else { // not first wire
                if (hit.WireID().Wire == wirenum) { // if wire num is same as previous hit
                  std::cout<<"TWO HITS ON SAME WIRE "<<wirenum<<std::endl;
  
                  if ((hit.PeakTime()-peaktime)<20) {  // require this second hit to be within 20 ticks of first
                    fdistancetotwohits1 = wirenum-firstwire;
                    std::cout << "num wires "<< fdistancetotwohits1<<std::endl;
                    twohits1=true;
                    hitcount = 0;
                    splitwire = hit.WireID().Wire;
                  }
                  else { std::cout<< "too far away"<<std::endl; }
                }
  
                else { // wire is not the same number as previous hit
                  if (hit.WireID().Wire != wirenum+1) {
                    std::cout<< "GAP IN HITS, no event"<<std::endl;
                    hitcount=1; firstwire = hit.WireID().Wire; // this becomes the new 'first hit'
                    wirenum = hit.WireID().Wire;
                    peaktime = hit.PeakTime();
                  }
                  else { // wire must be wirenum+1 ie carry on to next one
                    wirenum = hit.WireID().Wire;
                    peaktime = hit.PeakTime();
                    std::cout<<"next hit, on "<<wirenum+1<<std::endl;
                  }
                } // if only one hit on previous wire
              }
            }



            if (hit.WireID().Plane == 0 && twohits0==false) { // require amplitude of hits to be greater than 35 adc units
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
              }

              else { // not first wire
                if (hit.WireID().Wire == wirenum) { // if wire num is same as previous hit
                  std::cout<<"TWO HITS ON SAME WIRE "<<wirenum<<std::endl;

                  if ((hit.PeakTime()-peaktime)<20) {  // require this second hit to be within 20 ticks of first
                    fdistancetotwohits0 = wirenum-firstwire;
                    std::cout << "num wires "<< fdistancetotwohits0<<std::endl;
                    twohits0=true;
                    hitcount = 0;
                    splitwire = hit.WireID().Wire;
                  }
                  else { std::cout<< "too far away"<<std::endl; }
                }

                else { // wire is not the same number as previous hit
                  if (hit.WireID().Wire != wirenum+1) {
                    std::cout<< "GAP IN HITS, no event"<<std::endl;
                    hitcount=1; firstwire = hit.WireID().Wire; // this becomes the new 'first hit'
                    wirenum = hit.WireID().Wire;
                    peaktime = hit.PeakTime();
                  }
                  else { // wire must be wirenum+1 ie carry on to next one
                    wirenum = hit.WireID().Wire;
                    peaktime = hit.PeakTime();
                    std::cout<<"next hit, on "<<wirenum+1<<std::endl;
                  }
                } // if only one hit on previous wire
              }
            }

       
          }

          fhit_RMS = hit.RMS();
          fhit_PeakAmplitude = hit.PeakAmplitude();
          fhit_Integral = hit.Integral();
  
          fOutputTree->Fill(); // Fills once per hit
        }
        
      }
      int n(0); int m(0);
      float AmpIntChange(0); float AmpIntChange_firstfive(0);
      std::cout<< "firstwire: "<<firstwire <<" splitwire: "<<splitwire<<std::endl;
      double lastwire(-1); 
      double lastpeaktime(-1);
      float firstamparea(-1); 
      float fifthamparea(-1);
      float lastamparea(-1);
      // loop through hits again (now knowing where split pos is)
      for (recob::Hit const& hit: hits) {
        // do this for plane 2 (collection)
        if (hit.WireID().TPC == tpc && hit.WireID().Cryostat == 0  &&  hit.WireID().Plane == 2 && hit.PeakAmplitude()>35) {
          if (hit.WireID().Wire>=firstwire && hit.WireID().Wire<=splitwire+5) {
            n++;
            fAmplitudeArea.push_back(hit.PeakAmplitude()/hit.Integral());
            fWire.push_back(hit.WireID().Wire);
            std::cout<<n<<" wire "<<hit.WireID().Wire<<" AmplitudeArea ratio "<<hit.PeakAmplitude()/hit.Integral()<<std::endl;

          }
          
          if (hit.WireID().Wire>=firstwire && hit.WireID().Wire<splitwire) {
            if (lastwire!=-1 && (hit.PeakTime()-lastpeaktime)<20) {
              m++;
              AmpIntChange+=(hit.PeakAmplitude()/hit.Integral())-lastwire; 
            }
            // get average of first five changes in amp/area
            if (lastwire!=-1 && (splitwire-firstwire)>6 && (hit.PeakTime()-lastpeaktime)<20) {
              if (hit.WireID().Wire-firstwire<6) {
                AmpIntChange_firstfive+=(hit.PeakAmplitude()/hit.Integral())-lastwire;
              }
            }
            
            // get the amp/area for first fifth and last wire before splitting
            if (hit.WireID().Wire==firstwire) { firstamparea = hit.PeakAmplitude()/hit.Integral(); }
            if (hit.WireID().Wire==firstwire+5 && (hit.PeakTime()-lastpeaktime)<20) { 
              fifthamparea = hit.PeakAmplitude()/hit.Integral();
             }
            if (hit.WireID().Wire==splitwire-1 && (hit.PeakTime()-lastpeaktime)<20) { 
              lastamparea = hit.PeakAmplitude()/hit.Integral(); 
            }
            // only update the amp/area to compare and the peaktime if its the first wire OR if it is less than 20 ticks from previous
             if (hit.WireID().Wire==firstwire || (hit.PeakTime()-lastpeaktime)<20) {
              lastwire = hit.PeakAmplitude()/hit.Integral();
              lastpeaktime = hit.PeakTime();
            }

            
          }


        }


      
      }
      double scale =static_cast<double>(m);
      fAvChange = AmpIntChange/scale;
      fAvChange_firstfive = AmpIntChange_firstfive/static_cast<double>(5);
      ftotalchange_amparea = lastamparea-firstamparea;
      ffirstfivechange_amparea = fifthamparea - firstamparea;
    }
    if (fdistancetotwohits2!=-1 && fdistancetotwohits1!=-1 && fdistancetotwohits0!=-1) {
      fEventTree->Fill(); // Fills once per event
    }
  }
}

void sbndwire::WireHitsAnalysis::beginJob() {

  // use TFileService to create tree and then add branches
  art::ServiceHandle<art::TFileService> tfs;
  fOutputTree = tfs->make<TTree>("mytree","My Tree");
  fOutputTree->Branch("event",&fEvent,"event/i");
  fOutputTree->Branch("hit_RMS",&fhit_RMS,"hit_RMS");
  fOutputTree->Branch("hit_PeakAmplitude",&fhit_PeakAmplitude,"hit_PeakAmplitude");
  fOutputTree->Branch("hit_Integral",&fhit_Integral,"hit_Integral");


  fEventTree = tfs->make<TTree>("myeventtree","My Event Tree");
  fEventTree->Branch("event",&fEvent,"event/i");
  fEventTree->Branch("distancetotwohits_plane2",&fdistancetotwohits2,"distancetotwohits2/i");
  fEventTree->Branch("distancetotwohits_plane1",&fdistancetotwohits1,"distancetotwohits1/i");
  fEventTree->Branch("distancetotwohits_plane0",&fdistancetotwohits0,"distancetotwohits0/i");
  fEventTree->Branch("AvChange",&fAvChange,"AvChange/d");
  fEventTree->Branch("AvChange_firstfive",&fAvChange_firstfive,"AvChange_firstfive/d");

  fEventTree->Branch("pairE",&fpairE,"pairE/f");
  fEventTree->Branch("separationAngle",&fsepAngle,"separationAngle/f");
  fEventTree->Branch("theta_one",&ftheta_one,"theta_one/f");
  fEventTree->Branch("theta_two",&ftheta_two,"theta_two/f");
  fEventTree->Branch("avSigmaChange",&favSigmaChange,"favSigmaChange/f");

  fEventTree->Branch("photonEnergy",&fenergy,"energy/f");
  fEventTree->Branch("photonTheta",&ftheta,"theta/f");

  fEventTree->Branch("KE_high",&fKE_high,"KE_high/f");
  fEventTree->Branch("KE_low",&fKE_low,"KE_low/f");
  fEventTree->Branch("totalchange_amparea",&ftotalchange_amparea,"totalchange_amparea/f");
  fEventTree->Branch("firstfivechange_amparea",&ffirstfivechange_amparea,"firstfivechange_amparea/f");


  //fRatioGraph1   = tfs->makeAndRegister<TGraph>("channelGraph", "Graph title;channel;# PE",numpoints,&fAmplitudeArea[0],fWire[0]);
}

void sbndwire::WireHitsAnalysis::endJob() {
}

DEFINE_ART_MODULE(sbndwire::WireHitsAnalysis)

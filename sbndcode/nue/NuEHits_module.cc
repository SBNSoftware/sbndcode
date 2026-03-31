////////////////////////////////////////////////////////////////////////
// Class:       NuEHits
// Plugin Type: analyzer (Unknown Unknown)
// File:        NuEHits_module.cc
//
// Generated at Mon Jan 19 02:34:51 2026 by Rachel Coackley using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

// Sim Base
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"

// Reco Base
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

// Tools
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

// LArSoft
#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
// #include "larcorealg/Geometry/ChannelMapStandardAlg.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcore/Geometry/WireReadout.h"

#include "sbnobj/Common/Reco/CRUMBSResult.h"

#include "larcoreobj/SummaryData/POTSummary.h"
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"

// C++
#include <string>
#include <vector>
#include <map>
#include <cmath>

// ROOT
#include "art_root_io/TFileService.h"
#include <TTree.h>
#include <TFile.h>

constexpr double def_double = -std::numeric_limits<double>::max();

namespace sbnd {
  class NuEHits;
}


class sbnd::NuEHits : public art::EDAnalyzer {
public:
  explicit NuEHits(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NuEHits(NuEHits const&) = delete;
  NuEHits(NuEHits&&) = delete;
  NuEHits& operator=(NuEHits const&) = delete;
  NuEHits& operator=(NuEHits&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;
  void PFPs(art::Event const& e); 
  void MCParticles(art::Event const& e);  
  void Hits(art::Event const& e);
  void SpacePoints(art::Event const& e);
  void clearVectors();
 
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  unsigned int eventID;
  unsigned int runID;
  unsigned int subRunID;

  TTree* NuEHitsTree;
  std::vector<double>   truth_neutrinoVX;
  std::vector<double>   truth_neutrinoVY;
  std::vector<double>   truth_neutrinoVZ;

  std::vector<double>   reco_neutrinoVX;
  std::vector<double>   reco_neutrinoVY;
  std::vector<double>   reco_neutrinoVZ;
  std::vector<double>   reco_neutrinoSliceID;
  
  std::vector<double>   hit_plane;
  std::vector<double>   hit_x;
  std::vector<double>   hit_uvz;
  std::vector<double>   hit_sliceID;
  std::vector<double>   hit_PFPID;

  std::vector<double> spacepoint_x;
  std::vector<double> spacepoint_y;
  std::vector<double> spacepoint_z;
  std::vector<double> spacepoint_sliceID;
  std::vector<double> spacepoint_PFPID;

  art::ServiceHandle<cheat::ParticleInventoryService> particleInv;
  art::ServiceHandle<geo::Geometry> theGeometry;
  detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();

  const std::string PFParticleLabel;
  const std::string sliceLabel;
  const std::string sliceSCELabel;
  const std::string vertexLabel;
  const std::string hitLabel;
  const std::string showerLabel;
  const std::string TruthLabel;
  const std::string spacePointLabel;
  double DLCurrent;
  double signal;

  TFile *outputFile = TFile::Open("NuEHitsAnalyserOutput.root","RECREATE");

};


sbnd::NuEHits::NuEHits(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  PFParticleLabel(p.get<std::string>("PFParticleLabel")),
  sliceLabel(p.get<std::string>("SliceLabel")),
  sliceSCELabel(p.get<std::string>("SliceSCELabel")),
  vertexLabel(p.get<std::string>("VertexLabel")),
  hitLabel(p.get<std::string>("HitLabel")),
  showerLabel(p.get<std::string>("ShowerLabel")),
  TruthLabel(p.get<std::string>("TruthLabel")),
  spacePointLabel(p.get<std::string>("SpacePointLabel")),
  DLCurrent(p.get<double>("DLCurrent")),
  signal(p.get<double>("Signal"))
{
  art::ServiceHandle<art::TFileService> fs;  
  NuEHitsTree = fs->make<TTree>("NuEHits","");
  NuEHitsTree->Branch("eventID", &eventID);
  NuEHitsTree->Branch("runID", &runID);
  NuEHitsTree->Branch("subRunID", &subRunID);
  NuEHitsTree->Branch("DLCurrent", &DLCurrent);
  NuEHitsTree->Branch("signal", &signal);
  NuEHitsTree->Branch("truth_neutrinoVX", &truth_neutrinoVX);
  NuEHitsTree->Branch("truth_neutrinoVY", &truth_neutrinoVY);
  NuEHitsTree->Branch("truth_neutrinoVZ", &truth_neutrinoVZ);
  NuEHitsTree->Branch("reco_neutrinoVX", &reco_neutrinoVX);
  NuEHitsTree->Branch("reco_neutrinoVY", &reco_neutrinoVY);
  NuEHitsTree->Branch("reco_neutrinoVZ", &reco_neutrinoVZ);
  NuEHitsTree->Branch("reco_neutrinoSliceID", &reco_neutrinoSliceID);
  NuEHitsTree->Branch("hit_plane", &hit_plane);
  NuEHitsTree->Branch("hit_x", &hit_x);
  NuEHitsTree->Branch("hit_uvz", &hit_uvz);
  NuEHitsTree->Branch("hit_sliceID", &hit_sliceID);
  NuEHitsTree->Branch("hit_PFPID", &hit_PFPID);
  NuEHitsTree->Branch("spacepoint_x", &spacepoint_x);
  NuEHitsTree->Branch("spacepoint_y", &spacepoint_y);
  NuEHitsTree->Branch("spacepoint_z", &spacepoint_z);
  NuEHitsTree->Branch("spacepoint_sliceID", &spacepoint_sliceID);
  NuEHitsTree->Branch("spacepoint_PFPID", &spacepoint_PFPID);
    
}

void sbnd::NuEHits::analyze(art::Event const& e)
{
    clearVectors();    

    eventID = e.id().event();
    runID = e.id().run();
    subRunID = e.id().subRun();

    MCParticles(e);
    PFPs(e);
    Hits(e);
    SpacePoints(e);

    NuEHitsTree->Fill();
}

void sbnd::NuEHits::SpacePoints(art::Event const& e){
    std::cout << "_________ SpacePoints _________" << std::endl;
    art::Handle<std::vector<recob::SpacePoint>> spacePointHandle;
    std::vector<art::Ptr<recob::SpacePoint>> spacePointVec;
    if(e.getByLabel(spacePointLabel, spacePointHandle))
        art::fill_ptr_vector(spacePointVec, spacePointHandle);

    art::Handle<std::vector<recob::PFParticle>> PFPHandle;
    std::vector<art::Ptr<recob::PFParticle>> PFPVec;
    if(e.getByLabel(PFParticleLabel, PFPHandle))
        art::fill_ptr_vector(PFPVec, PFPHandle);

    std::cout << "Number of spacepoints in event = " << spacePointVec.size() << std::endl;
    if(!spacePointVec.empty()){
        // Get association between spacepoint and PFP
        art::FindManyP<recob::PFParticle> spacePointPFPAssns(spacePointVec, e, spacePointLabel);
        
        // Get association between PFP and Slice
        art::FindManyP<recob::Slice> PFPSliceAssns(PFPVec, e, sliceSCELabel);

        // Loop through the spacepoints
        for(const art::Ptr<recob::SpacePoint> &spacePoint : spacePointVec){
            // Getting the PFP associated with the SpacePoint
            const std::vector<art::Ptr<recob::PFParticle>> spacePointPFPs(spacePointPFPAssns.at(spacePoint.key()));
            if(spacePointPFPs.size() != 0){
                // There is a PFP associated with the SpacePoint
                art::Ptr<recob::PFParticle> spacePointPFP = spacePointPFPs.at(0);
                spacepoint_PFPID.push_back(spacePointPFP->Self());

                const std::vector<art::Ptr<recob::Slice>> spacePointSlices(PFPSliceAssns.at(spacePointPFP.key()));
                if(spacePointSlices.size() != 0){
                    // There is a slice associated with the SpacePoint
                    art::Ptr<recob::Slice> spacePointSlice = spacePointSlices.at(0);
                    spacepoint_sliceID.push_back(spacePointSlice->ID());
                } else{
                    spacepoint_sliceID.push_back(-999999);
                }
            } else{
                spacepoint_PFPID.push_back(-999999);
                spacepoint_sliceID.push_back(-999999);
            }

            spacepoint_x.push_back(spacePoint->XYZ()[0]);
            spacepoint_y.push_back(spacePoint->XYZ()[1]);
            spacepoint_z.push_back(spacePoint->XYZ()[2]);

        }
    } else{
        std::cout << "No SpacePoints in the event" << std::endl;
        spacepoint_x.push_back(-999999);    
        spacepoint_y.push_back(-999999);    
        spacepoint_z.push_back(-999999);    
        spacepoint_sliceID.push_back(-999999);    
        spacepoint_PFPID.push_back(-999999);    
    }

    std::cout << "_______________________________" << std::endl;
}

void sbnd::NuEHits::PFPs(art::Event const& e){
    std::cout << "_________ PFParticles _________" << std::endl;
    art::Handle<std::vector<recob::PFParticle>> pfpHandle;
    std::vector<art::Ptr<recob::PFParticle>> pfpVec;
    if(e.getByLabel(PFParticleLabel, pfpHandle))
        art::fill_ptr_vector(pfpVec, pfpHandle);

    int neutrinoCounter = 0;

    if(!pfpVec.empty()){
        art::FindManyP<recob::Vertex> pfpVertexAssns(pfpVec, e, vertexLabel);
        art::FindOneP<recob::Slice> pfpSliceAssns(pfpVec, e, sliceSCELabel);    
                
        for(const art::Ptr<recob::PFParticle> &pfp : pfpVec){
            if(!(pfp->PdgCode() == std::numeric_limits<int>::max())){
                const art::Ptr<recob::Slice> pfpSlice(pfpSliceAssns.at(pfp.key()));
                const std::vector<art::Ptr<recob::Vertex>> pfpVertexs(pfpVertexAssns.at(pfp.key()));
                
                if(pfpVertexs.size() > 0){
                    if(!(pfp->IsPrimary() && (pfp->PdgCode() == 12 || pfp->PdgCode() == 14))){
                        // PFP is not the reco neutrino
                    } else{
                        // This is the reco neutrino
                        neutrinoCounter++;
                        const art::Ptr<recob::Vertex> &pfpVertex(pfpVertexs.front());
                    
                        reco_neutrinoVX.push_back(pfpVertex->position().X());
                        reco_neutrinoVY.push_back(pfpVertex->position().Y());
                        reco_neutrinoVZ.push_back(pfpVertex->position().Z());
                        reco_neutrinoSliceID.push_back(pfpSlice->ID());
                    }
                }
            }
        } 
    }

    if(neutrinoCounter == 0){
        reco_neutrinoVX.push_back(-999999);
        reco_neutrinoVY.push_back(-999999);
        reco_neutrinoVZ.push_back(-999999);
        reco_neutrinoSliceID.push_back(-999999);
    }
    std::cout << "_________________________" << std::endl;
}


void sbnd::NuEHits::MCParticles(art::Event const& e){
    art::Handle<std::vector<simb::MCTruth>> MCTruthHandle;
    std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
    if(e.getByLabel(TruthLabel, MCTruthHandle))
        art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

    std::cout << "_________ MCTruth _________" << std::endl;
    int counter = 0;
    if(!MCTruthVec.empty()){
        for(auto &MCTruth : MCTruthVec){
            std::cout << "MCTruth: Origin = " << MCTruth->Origin() << ", Number of Particles = " << MCTruth->NParticles() << std::endl;
            if(MCTruth->Origin() == simb::kBeamNeutrino){
                simb::MCNeutrino neutrino = MCTruth->GetNeutrino();
                simb::MCParticle neutrinoParticle = neutrino.Nu();
                std::cout << "Neutrino: Track ID = " << neutrinoParticle.TrackId() << ", CCNC = " << neutrino.CCNC() << ", Interaction Type = " << neutrino.InteractionType() << "Vertex = (" << neutrinoParticle.Vx() << ", " << neutrinoParticle.Vy() << ", " << neutrinoParticle.Vz() << "), Number of Daughters = " << neutrinoParticle.NumberDaughters() << std::endl; 
                int neutrinoTrackID = neutrinoParticle.TrackId();

                if(neutrino.InteractionType() == 1098){
                    // This is a Nu+E elastic scattering event
                    for(int i = 0; i < MCTruth->NParticles(); i++){
                        simb::MCParticle particle = MCTruth->GetParticle(i);
                        if(particle.TrackId() != neutrinoTrackID){
                            // This is not the neutrino
                            if(particle.StatusCode() == 1 && (particle.PdgCode() == 11 || particle.PdgCode() == -11)){
                                // StatusCode = 1: G4 is tracking this particle, it is outgoing
                                // Particle has to be an electron or positron
                                // This is the recoil electron!
                                
                                counter++;
                                truth_neutrinoVX.push_back(neutrinoParticle.Vx());
                                truth_neutrinoVY.push_back(neutrinoParticle.Vy());
                                truth_neutrinoVZ.push_back(neutrinoParticle.Vz());
                            }

                        }
                    }
                } else{
                    // This is a BNB or Cosmic event
                }
            }
        }
    }
    
    // There is no true recoil electron from a Nu+E elastic scattering interaction in the event. 
    if(counter == 0){
        truth_neutrinoVX.push_back(-999999);
        truth_neutrinoVY.push_back(-999999);
        truth_neutrinoVZ.push_back(-999999);
    }

    std::cout << "_________________________" << std::endl;
}


void sbnd::NuEHits::Hits(art::Event const& e){
    std::cout << "_________ Hits _________" << std::endl;
    
    detinfo::DetectorPropertiesData propD = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
    
    art::Handle<std::vector<recob::Hit>> hitHandle;
    std::vector<art::Ptr<recob::Hit>> hitVec;

    if(e.getByLabel(hitLabel, hitHandle))
        art::fill_ptr_vector(hitVec, hitHandle);

    std::cout << "Number of hits in event = " << hitVec.size() << std::endl;
    if(!hitVec.empty()){
        // There is a hit in the event
        
        // Get the slices
        art::Handle<std::vector<recob::Slice>> sliceHandle;
        std::vector<art::Ptr<recob::Slice>> sliceVec;
        if(e.getByLabel(sliceLabel, sliceHandle))
            art::fill_ptr_vector(sliceVec, sliceHandle);
       
        // Get the showers
        art::Handle<std::vector<recob::Shower>> showerHandle;
        std::vector<art::Ptr<recob::Shower>> showerVec;
        if(e.getByLabel(showerLabel, showerHandle))
            art::fill_ptr_vector(showerVec, showerHandle);

        // Get the PFPs
        art::Handle<std::vector<recob::PFParticle>> pfpHandle;
        std::vector<art::Ptr<recob::PFParticle>> pfpVec;
        if(e.getByLabel(PFParticleLabel, pfpHandle))
            art::fill_ptr_vector(pfpVec, pfpHandle);

        for(const art::Ptr<recob::Hit> &hit : hitVec){
            // Fill Hit branches:
            hit_plane.push_back(hit->View());

            if(hit->View() == 0 || hit->View() == 1 || hit->View() == 2){
                const geo::WireID& wireID(hit->WireID());
                const double xpos = propD.ConvertTicksToX(hit->PeakTime(), wireID.Plane, wireID.TPC, wireID.Cryostat);
                auto const& channelMapAlg =  art::ServiceHandle<geo::WireReadout const>()->Get();
                auto const wire = channelMapAlg.Plane(geo::PlaneID(wireID.Cryostat, wireID.TPC, wireID.Plane)).Wire(wireID);
                geo::Point_t start = wire.GetStart();
                geo::Point_t end = wire.GetEnd();

                const double ay(start.Y());
                const double az(start.Z());
                const double by(end.Y());
                const double bz(end.Z());

                const double ny(by - ay);
                const double nz(bz - az);
                const double n2(ny * ny + nz * nz);

                const double ry(ay - (ay * ny + az * nz) * ny / n2);
                const double rz(az - (ay * ny + az * nz) * nz / n2);
                const double sign((rz > 0.0) ? +1.0 : -1.0);
                const double uvz = sign * std::sqrt(ry * ry + rz * rz); // this is the U/V/Z coordinate
                
                hit_x.push_back(xpos);
                hit_uvz.push_back(uvz);
            } else{
                hit_x.push_back(-999999);
                hit_uvz.push_back(-999999);
            }
            
            // Getting the slice associated with the hit
            art::FindManyP<recob::Slice> hitSliceAssns(hitVec, e, sliceLabel);
            const std::vector<art::Ptr<recob::Slice>> hitSlices(hitSliceAssns.at(hit.key()));
            if(hitSlices.size() != 0){
                // There is a slice associated with the hit
                art::Ptr<recob::Slice> slice = hitSlices.at(0);
                hit_sliceID.push_back(slice->ID());
            } else{
                // There is no slice associated with the hit
                hit_sliceID.push_back(-999999);
            }

            // Getting the shower associated with the hit
            art::FindManyP<recob::Shower> hitShowerAssns(hitVec, e, showerLabel);
            const std::vector<art::Ptr<recob::Shower>> hitShowers(hitShowerAssns.at(hit.key()));

            if(hitShowers.size() != 0){
                // Getting the PFP associated with the shower
                art::FindManyP<recob::PFParticle> showerPFPAssns(hitShowers, e, showerLabel);
                const std::vector<art::Ptr<recob::PFParticle>> showerPFP(showerPFPAssns.at(0));
                art::Ptr<recob::PFParticle> pfp = showerPFP.at(0);
                hit_PFPID.push_back(pfp->Self());
            } else{
                // There is no PFP associated with the hit
                hit_PFPID.push_back(-999999);
            }
        }
    }
}

void sbnd::NuEHits::clearVectors(){
  truth_neutrinoVX.clear();
  truth_neutrinoVY.clear();
  truth_neutrinoVZ.clear();
  
  reco_neutrinoVX.clear();
  reco_neutrinoVY.clear();
  reco_neutrinoVZ.clear();
  reco_neutrinoSliceID.clear();

  hit_plane.clear();
  hit_x.clear();
  hit_uvz.clear();
  hit_sliceID.clear();
  hit_PFPID.clear();
  
  spacepoint_x.clear();
  spacepoint_y.clear();
  spacepoint_z.clear();
  spacepoint_sliceID.clear();
  spacepoint_PFPID.clear();
}

void sbnd::NuEHits::beginJob()
{
  // Implementation of optional member function here.
}

void sbnd::NuEHits::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(sbnd::NuEHits)

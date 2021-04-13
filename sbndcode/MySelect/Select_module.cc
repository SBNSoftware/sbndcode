////////////////////////////////////////////////////////////////////////
// Class:       Select
// Plugin Type: analyzer (art v3_05_01)
// File:        Select_module.cc
//
// Generated at Wed Jan 20 11:32:02 2021 by Vu Chi Lan Nguyen using cetskelgen
// from cetlib version v3_10_00.
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

// Additional framework includes
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"

// ROOT includes
#include <TTree.h>

namespace sbnd {
  class Select;
}


class sbnd::Select : public art::EDAnalyzer {

public:
  explicit Select(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Select(Select const&) = delete;
  Select(Select&&) = delete;
  Select& operator=(Select const&) = delete;
  Select& operator=(Select&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Service Handle
  art::ServiceHandle<art::TFileService> tfs;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

  //labels
  std::string fPFParticleLabel;
  std::string fT0Label;
  std::string fTrackLabel;
  std::string fSpacePointLabel;
  std::string fCalorimetryLabel;
  std::string fHitLabel;
  
  // Declare trees
  TTree *pfpTree; 	// tree for selected CR
  TTree *trueTree;	//tree for all true CR

  // Things to go in pfptree
  unsigned int		pfpEventID;
  size_t 		pfpParticleID;
  bool			pfpIsPrimary;
  int		        pfpParticlePDG;
  double		pfpTrackLength;
  double 		pfpYZAngle;
  double		pfpXZAngle;
  std::vector<float>	pfpdQdx_0;
  std::vector<float>	pfpX_0;
  std::vector<float>	pfpY_0;
  std::vector<float>	pfpZ_0;
  std::vector<float>	pfpdQdx_1;
  std::vector<float>	pfpX_1;
  std::vector<float>	pfpY_1;
  std::vector<float>	pfpZ_1;
  std::vector<float>	pfpdQdx_2;
  std::vector<float>	pfpX_2;
  std::vector<float>	pfpY_2;
  std::vector<float>	pfpZ_2;
  int			pfpTrueTrackID;
  int 		        pfpTruePDG;
  int			pfpTrueMother;
  bool 			pfpTrueCrossesApa;
  bool 			pfpTrueCrossesCpa;  

  // analysis class
  TPCGeoAlg fTPCGeo;

  //Things go into true tree
  unsigned int	mcEventID;
  int 		mcTrackID;
  int  		mcTruePDG;
  bool  	mcCrossesApa;
  bool 		mcCrossesCpa; 
  int 		mcTrueMother;
  int 		mcPfpParticleID;

  // Map to hold track truth id and particle id of a pfp
  std::map<int, bool> trueTrackIDpfpID;

  //Additional member functions
  bool TrackCrossesApa(const art::Ptr<recob::Track> trk) const;

  bool TrackCrossesCpa(const art::Ptr<recob::Track> trk) const;
  
  std::vector< art::Ptr<recob::PFParticle> > GetCrossingCR_T0(
	const std::vector< art::Ptr<recob::PFParticle> > &pfps, 
	const art::FindManyP<recob::Track> &trackAssoc,
 	const art::FindManyP<anab::T0> &t0Assoc) const;
  
  std::vector< art::Ptr<recob::PFParticle> > GetCrossingCR(
	const std::vector< art::Ptr<recob::PFParticle> > &pfps, 
	const art::FindManyP<recob::Track> &trackAssoc) const;
  
  bool MCParticleCrossesCpa(
	const simb::MCParticle& particle, 
	const sbnd::TPCGeoAlg TPCGeo) const; 

  double LengthTPC0(const simb::MCParticle& particle) const;
  
  double LengthTPC1(const simb::MCParticle& particle) const;
  
  bool MinLengthSquared(double lengthTPC0, double lengthTPC1) const;
  
  void clearPfpTree();

  void clearTrueTree();
};


sbnd::Select::Select(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.

  fPFParticleLabel 	= p.get<std::string>("PFParticleLabel");
  fT0Label 		= p.get<std::string>("T0Label");
  fTrackLabel 		= p.get<std::string>("TrackLabel");
  fSpacePointLabel 	= p.get<std::string>("SpacePointLabel");
  fCalorimetryLabel 	= p.get<std::string>("CalorimetryLabel");
  fHitLabel 		= p.get<std::string>("HitLabel");
}

void sbnd::Select::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  //std::cout << std::endl;
  //std::cout << "Reading event: " << e.id().event() << std::endl;
  
  // and clock data for event
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);

  // Define handle and fill vector
 
  art::Handle< std::vector<recob::PFParticle> > PFParticleHandle;
  std::vector< art::Ptr<recob::PFParticle> > PFParticleList;
  if(e.getByLabel(fPFParticleLabel, PFParticleHandle)){
    art::fill_ptr_vector(PFParticleList, PFParticleHandle);
  }

  art::Handle< std::vector<recob::Track> > TrackHandle;
  std::vector< art::Ptr<recob::Track> > TrackList;
  if(e.getByLabel(fTrackLabel, TrackHandle)){
    art::fill_ptr_vector(TrackList, TrackHandle);
  }

  // Define association
 
  art::FindManyP<recob::Track> trackAssoc (PFParticleList, e, fTrackLabel);
  art::FindManyP<anab::T0> t0Assoc (PFParticleList, e, fT0Label);

  art::FindManyP<anab::Calorimetry> calAssoc (TrackList, e, fCalorimetryLabel);
  art::FindManyP<recob::Hit> hitAssoc (TrackList, e , fHitLabel);

  //--------------------------Selection------------------------------//

  if(PFParticleList.empty()) return;

  std::vector< art::Ptr<recob::PFParticle> > CrossingCR = this->GetCrossingCR_T0(PFParticleList, trackAssoc, t0Assoc);

  //std::cout << "No. of cathode-anode crossing muon: " << CrossingCR.size() << std::endl;

 
  //-------------------------Fill pfptree----------------------------// 
  // Save event every pfparticle as a single row in the tree

  for(const art::Ptr<recob::PFParticle> &pfp: CrossingCR){

    clearPfpTree();

    //Event ID
    
    pfpEventID = e.id().event();
    
    // PFParticle    
    
    pfpParticleID = pfp->Self();
    pfpIsPrimary = pfp->IsPrimary();
    pfpParticlePDG = pfp->PdgCode();
    
    //std::cout << "Saving " << "Paritcle ID: " << pfpParticleID;
    //std::cout << std::endl;
    
    // Track
    
    std::vector< art::Ptr<recob::Track> > pfpTrack = trackAssoc.at(pfp.key());

    if(pfpTrack.empty()) continue;
     
    for(const art::Ptr<recob::Track> &trk: pfpTrack){

      pfpTrackLength = trk->Length();     
      pfpYZAngle = trk->ZenithAngle();
      pfpXZAngle = trk->AzimuthAngle();
      
      // Calorimetry
    
      std::vector< art::Ptr<anab::Calorimetry> > trkCal = calAssoc.at(trk.key());
     
      if(trkCal.empty()) continue;

      for(const art::Ptr<anab::Calorimetry> &cal: trkCal){
   
        if(!cal->PlaneID().isValid) continue;
       

        int planenum = cal->PlaneID().Plane;

        if( planenum == 0 ){
	  for(unsigned int i =0; i < cal->dQdx().size(); i++){ 
            pfpdQdx_0.push_back(cal->dQdx()[i]);
            pfpX_0.push_back(cal->XYZ()[i].X());
            pfpY_0.push_back(cal->XYZ()[i].Y());
            pfpZ_0.push_back(cal->XYZ()[i].Z());
          }
        }

        if( planenum == 1 ){
	  for(unsigned int i =0; i < cal->dQdx().size(); i++){ 
            pfpdQdx_1.push_back(cal->dQdx()[i]);
            pfpX_1.push_back(cal->XYZ()[i].X());
            pfpY_1.push_back(cal->XYZ()[i].Y());
            pfpZ_1.push_back(cal->XYZ()[i].Z());
          }
        }
      
        if( planenum == 2 ){
	  for(unsigned int i =0; i < cal->dQdx().size(); i++){ 
            pfpdQdx_2.push_back(cal->dQdx()[i]);
            pfpX_2.push_back(cal->XYZ()[i].X());
            pfpY_2.push_back(cal->XYZ()[i].Y());
            pfpZ_2.push_back(cal->XYZ()[i].Z());
          }
        }
     } 

      // MC Truth Info
     
      if(!e.isRealData()){
  
        std::vector< art::Ptr<recob::Hit> > trkHit = hitAssoc.at(trk.key());

        if(trkHit.empty()) continue;

	int trkidtruth = TruthMatchUtils::TrueParticleIDFromTotalTrueEnergy(clockData, trkHit, true);
        bool valid = TruthMatchUtils::Valid(trkidtruth);  

        if(!valid) continue;
    
        //std::cout << "Track ID: " << trkidtruth << " is valid: " << valid << std::endl;
       
        const simb::MCParticle *particle = particleInventory->TrackIdToParticle_P(trkidtruth);
       
	if(!particle){          
	  std::cout << "MCParticle is nullptr!" << std::endl;
 	  continue;
        }
 
        pfpTrueTrackID 	= particle->TrackId();
        pfpTruePDG     	= particle->PdgCode();
        pfpTrueMother 	= particle->Mother(); 
          
        //std::cout << "Getting MC data..\n";
        //std::cout << "True PDG code: " << pfpTruePDG;  
        //std::cout << ", True TrackID: " << pfpTrueTrackID;
        //std::cout << ", True Mother: " << pfpTrueMother;
        //std::cout << std::endl;
          
        pfpTrueCrossesApa = fTPCGeo.CrossesApa(*particle);
        pfpTrueCrossesCpa = this->MCParticleCrossesCpa(*particle, fTPCGeo);
        trueTrackIDpfpID[pfpTrueTrackID] = pfpParticleID;
 
      } // End of MC loop
     
    } // End of track loop  
 
    pfpTree->Fill(); 
  
  } // End of pfparticle loop

  //---------------------Get MC Truth to fill in true tree------------------------//
  const sim::ParticleList& MCParticleList = particleInventory->ParticleList();

  for(auto const &MCIter: MCParticleList){
    
    clearTrueTree();

    const simb::MCParticle* mcp = MCIter.second;

    mcEventID = e.id().event();
    mcTrackID = mcp->TrackId();
    mcCrossesApa = fTPCGeo.CrossesApa(*mcp);
    mcCrossesCpa = this->MCParticleCrossesCpa(*mcp, fTPCGeo);
    mcTruePDG = mcp->PdgCode();  
    mcTrueMother = mcp->Mother();  


    if(mcTrueMother != 0) continue;
    if(std::abs(mcTruePDG) != 13) continue;
    if(!mcCrossesCpa) continue;
    if(!mcCrossesApa) continue;
    if(!MinLengthSquared(this->LengthTPC0(*mcp), this->LengthTPC1(*mcp))) continue;

    for(auto &[k,v] :trueTrackIDpfpID){
     
      if(mcTrackID == k){
        mcPfpParticleID = v;
      }      
    }

    trueTree->Fill();
  } // End of MCParticle loop


}// End  

void sbnd::Select::beginJob()
{
  // Implementation of optional member function here.
  
  pfpTree = tfs->make<TTree>("pfptree", "Selected Reconstructed Crossing Muons");
  trueTree = tfs->make<TTree>("truetree", "True MC Crossing Muons");

  pfpTree->Branch("eventID", &pfpEventID, "eventID/i");
  pfpTree->Branch("ParticleID", &pfpParticleID, "ParticleID/i");
  pfpTree->Branch("IsPrimary", &pfpIsPrimary, "IsPrimary/O");
  pfpTree->Branch("ParticlePDG", &pfpParticlePDG, "ParticlePDG/I");
  pfpTree->Branch("TrackLength", &pfpTrackLength, "TrackLength/D");
  pfpTree->Branch("YZAngle", &pfpYZAngle, "YZAngle/D");
  pfpTree->Branch("XZAngle", &pfpXZAngle, "XZAngle/D");
  pfpTree->Branch("dQdx_0", &pfpdQdx_0);
  pfpTree->Branch("X_0", &pfpX_0);
  pfpTree->Branch("Y_0", &pfpY_0);
  pfpTree->Branch("Z_0", &pfpZ_0);
  pfpTree->Branch("dQdx_1", &pfpdQdx_1);
  pfpTree->Branch("X_1", &pfpX_1);
  pfpTree->Branch("Y_1", &pfpY_1);
  pfpTree->Branch("Z_1", &pfpZ_1);
  pfpTree->Branch("dQdx_2", &pfpdQdx_2);
  pfpTree->Branch("X_2", &pfpX_2);
  pfpTree->Branch("Y_2", &pfpY_2);
  pfpTree->Branch("Z_2", &pfpZ_2);
  pfpTree->Branch("TrueTrackID", &pfpTrueTrackID, "TrueTrackID/I");
  pfpTree->Branch("TruePDG", &pfpTruePDG, "TruePDG/I");
  pfpTree->Branch("TrueMother", &pfpTrueMother, "TrueMother/I"); 
  pfpTree->Branch("TrueCrossesApa", &pfpTrueCrossesApa, "TrueCrossesApa/O");
  pfpTree->Branch("TrueCrossesCpa", &pfpTrueCrossesCpa, "TrueCrossesCpa/O");

  trueTree->Branch("eventID", &mcEventID, "eventID/i");
  trueTree->Branch("TrackID", &mcTrackID, "TrackID/I");
  trueTree->Branch("TruePDG", &mcTruePDG, "TruePDG/I");
  trueTree->Branch("CrossesApa", &mcCrossesApa, "CrossesApa/O");
  trueTree->Branch("CrossesCpa", &mcCrossesCpa, "CrossesCpa/O");
  trueTree->Branch("TrueMother", &mcTrueMother, "TrueMother/I");
  trueTree->Branch("PfpParticleID", &mcPfpParticleID, "PfpParticleID/I");
}


  //----------------------Member Functions------------------------------//


//-------------------------------------------------------------------
// Get pfparticle that crosses cathode-anode

std::vector< art::Ptr<recob::PFParticle> > sbnd::Select::GetCrossingCR(const std::vector< art::Ptr<recob::PFParticle> > &pfps, const art::FindManyP<recob::Track> &trackAssoc) const{

  std::vector< art::Ptr<recob::PFParticle> > crosser;
  crosser.clear();

  for(const art::Ptr<recob::PFParticle> &pfp: pfps){
   
    if(!pfp->IsPrimary()) continue; //check if it's a cosmic
   
    if(!(std::abs(pfp->PdgCode()) == 13)) continue; //check if  it's muon
   
    std::vector< art::Ptr<recob::Track> > pfpTrack = trackAssoc.at(pfp.key());

    if(pfpTrack.empty()) continue;

    for(const art::Ptr<recob::Track> &trk: pfpTrack){

      if(!(this->TrackCrossesCpa(trk))) continue;  
      if(!(this->TrackCrossesApa(trk))) continue;
      crosser.push_back(pfp);
    }
  }
  return crosser;
}  

//-------------------------------------------------------------------
// Get pfparticle that crosses cathode-anode

std::vector< art::Ptr<recob::PFParticle> > sbnd::Select::GetCrossingCR_T0(const std::vector< art::Ptr<recob::PFParticle> > &pfps, const art::FindManyP<recob::Track> &trackAssoc, const art::FindManyP<anab::T0> &t0Assoc) const{

  std::vector< art::Ptr<recob::PFParticle> > crosser;
  crosser.clear();

  for(const art::Ptr<recob::PFParticle> &pfp: pfps){
   
    if(!pfp->IsPrimary()) continue; //check if it's a cosmic
   
    if(!(std::abs(pfp->PdgCode()) == 13)) continue; //check if  it's muon
   
    std::vector< art::Ptr<anab::T0> > pfpT0 = t0Assoc.at(pfp.key());

    if(pfpT0.empty()) continue; //T0 stiches cosmic track that crosses cathode -> look for cathode crosser 
 
    std::vector< art::Ptr<recob::Track> > pfpTrack = trackAssoc.at(pfp.key());

    if(pfpTrack.empty()) continue;

    for(const art::Ptr<recob::Track> &trk: pfpTrack){

      if(!(this->TrackCrossesApa(trk))) continue;
      crosser.push_back(pfp);
    }
  }
  return crosser;
}  
//-----------------------------------------------------------------
// Determine if a reconstructed track crosses the cpa

bool sbnd::Select::TrackCrossesCpa(const art::Ptr<recob::Track> trk) const{
  
  for(size_t i = 0; i < trk->NPoints()-1; i++){

    geo::Point_t pos = trk->LocationAtPoint(i);
    double x = pos.X();
 
    geo::Point_t pos1 = trk->LocationAtPoint(i+1);
    double x1 = pos1.X();

    if(x <= -0.5 && x1 >= 0.5) return true;
    if(x >= 0.5 && x1 <= -0.5) return true;
 
    }
 
  return false;
}
//------------------------------------------------------------------
// Determine if a reconstructed track crosses the apa 

bool sbnd::Select::TrackCrossesApa(const art::Ptr<recob::Track> trk) const{

  geo::Point_t Start = trk->Start();    
  geo::Point_t End = trk->End();
  
  float distStart = 200 - abs(Start.X()); 
  float distEnd = 200 - abs(End.X());
   
  if(distStart < 10 || distEnd < 10){
    return true;
  }
  return false;
} 

//------------------------------------------------------------------
// Determine if a true particle crosses the CPA within the detector volume

bool sbnd::Select::MCParticleCrossesCpa(const simb::MCParticle& particle, const sbnd::TPCGeoAlg TPCGeo) const{  
  for(size_t i = 0; i < particle.NumberTrajectoryPoints()-1; i++){
    double x = particle.Vx(i);
    double y = particle.Vy(i);
    double z = particle.Vz(i);
    double x1 = particle.Vx(i+1);
    double y1 = particle.Vy(i+1);
    double z1 = particle.Vz(i+1);

    if( y >= TPCGeo.MinY() && z >= TPCGeo.MinZ() && y <= TPCGeo.MaxY() && z <= TPCGeo.MaxZ() &&
        y1 >= TPCGeo.MinY() && z1 >= TPCGeo.MinZ() && y1 <= TPCGeo.MaxY() && z1 <= TPCGeo.MaxZ()){
      if(x <= 0 && x1 >= 0) return true;
      if(x >= 0 && x1 <= 0) return true;
    }
  }
 return false;
}

//------------------------------------------------------------------
// Determine if a true particle, that crosses both APA & CPA, has min length squared > 50cm in TPC0

double sbnd::Select::LengthTPC0(const simb::MCParticle& particle) const{

  TPCGeoAlg TPCGeo;
  bool first = true;
  double length = 0; 
  
  TVector3 start (-99999, -99999, -99999);
  TVector3 end (-99999, -99999, -99999); 

  for(size_t i = 0; i < particle.NumberTrajectoryPoints()-1; i++){
    double x = particle.Vx(i);
    double y = particle.Vy(i);
    double z = particle.Vz(i);

    if( x > TPCGeo.MinX() && y > TPCGeo.MinY() && z > TPCGeo.MinZ() && x < 0 && y < TPCGeo.MaxY() && z < TPCGeo.MaxZ() ){
      if(first){
        first = false;
        start.SetXYZ(x, y, z);
      }
      end.SetXYZ(x, y, z);
    }
  } 
  
  if(start.X() != -99999 && end.X() != -99999){
    length = (start-end).Mag();
  }
  
  return length;
} 

//------------------------------------------------------------------
// Determine if a true particle, that crosses both APA & CPA, has min length squared > 50cm in TPC1

double sbnd::Select::LengthTPC1(const simb::MCParticle& particle) const{

  TPCGeoAlg TPCGeo;
  bool first = true;
  double length = 0; 
  
  TVector3 start (-99999, -99999, -99999);
  TVector3 end (-99999, -99999, -99999); 

  for(size_t i = 0; i < particle.NumberTrajectoryPoints()-1; i++){
    double x = particle.Vx(i);
    double y = particle.Vy(i);
    double z = particle.Vz(i);

    if( x > 0 && y > TPCGeo.MinY() && z > TPCGeo.MinZ() && x < TPCGeo.MaxX() && y < TPCGeo.MaxY() && z < TPCGeo.MaxZ() ){
      if(first){
        first = false;
        start.SetXYZ(x, y, z);
      }
      end.SetXYZ(x, y, z);
    }
  } 
  
  if(start.X() != -99999 && end.X() != -99999){
    length = (start-end).Mag();
  }
  
  return length;
} 
//------------------------------------------------------------------
// Determine if a true particle, that crosses both APA & CPA, has min length squared > 50cm in both TPC
bool sbnd::Select::MinLengthSquared(double lengthTPC0, double lengthTPC1) const{

  //Stitching Cosmics rays requires min length squared 50cm on both TPC

  if(lengthTPC0*lengthTPC0 > 50 && lengthTPC1*lengthTPC1 > 50) return true;  
  return false;
}

//-------------------------------------------------------------------
void sbnd::Select::clearPfpTree()
{
  pfpEventID = 99999;
  pfpParticleID = 99999;
  pfpIsPrimary = 99999;
  pfpParticlePDG = 99999;
  pfpTrackLength = 99999;
  pfpYZAngle = 99999;
  pfpXZAngle = 99999;
  pfpdQdx_0.clear();
  pfpX_0.clear();
  pfpY_0.clear();
  pfpZ_0.clear();
  pfpdQdx_1.clear();
  pfpX_1.clear();
  pfpY_1.clear();
  pfpZ_1.clear();
  pfpdQdx_2.clear();
  pfpX_2.clear();
  pfpY_2.clear();
  pfpZ_2.clear();
  pfpTrueTrackID = 99999;
  pfpTruePDG = 99999;
  pfpTrueMother = 99999;
  pfpTrueCrossesApa = 99999;
  pfpTrueCrossesCpa = 99999;
}

void sbnd::Select::clearTrueTree()
{
  mcEventID = 99999;
  mcTrackID = 99999;
  mcTruePDG = 99999;
  mcCrossesApa = 99999;
  mcCrossesCpa = 99999;
  mcTrueMother = 99999;
  mcPfpParticleID = 99999;
}
DEFINE_ART_MODULE(sbnd::Select)

////////////////////////////////////////////////////////////////////////
// Class:       analysis
// Plugin Type: analyzer (art v3_05_01)
// File:        analysis_module.cc
//
// Generated at Tue Sep 29 05:57:52 2020 by Henry Lay using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

//Standard
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//Sim Base
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/Simulation/SimChannel.h"

//Reco Base
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

//Tools
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

//LArSoft
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/GeometryUtilities.h"
#include "larcorealg/Geometry/ChannelMapStandardAlg.h"
#include "larcore/Geometry/Geometry.h"

//Root
#include "art_root_io/TFileService.h"
#include "TTree.h"
#include "TVector3.h"
#include "TGraph2D.h"
#include "TCanvas.h"

class analysis;


class analysis : public art::EDAnalyzer {
public:
  explicit analysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  analysis(analysis const&) = delete;
  analysis(analysis&&) = delete;
  analysis& operator=(analysis const&) = delete;
  analysis& operator=(analysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;


private:

  // Declare member data here.
  void ClearData();
  void SetupMaps(art::Event const &e);
  void SimulationProcessor(art::Event const &e);
  void MCParticleProcessor(art::Event const &e, 
			   art::Ptr<simb::MCParticle> const &particle,
			   bool const &isCosmic);
  
  void ReconstructionProcessor(art::Event const &e);
  void TrackProcessor(art::Event const &e,
		      art::Ptr<recob::PFParticle> const &pfp,
		      art::Ptr<recob::Track> const &track);
  void ShowerProcessor(art::Event const &e,
		       art::Ptr<recob::PFParticle> const &pfp,
		       art::Ptr<recob::Shower> const &shower);
  void HitProcessor(art::Event const &e, int const &nuID);
  void SliceHits(art::Event const &e, std::vector<art::Ptr<recob::Hit> > const &sliceHits);
  void VertexProcessor(art::Event const &e);

  const art::Ptr<simb::MCParticle> GetMCParticle(int const &trackID);
  const art::Ptr<recob::PFParticle> GetPrimaryPFP();
  const art::Ptr<recob::PFParticle> GetPFP(unsigned int const &index);
  const std::vector<art::Ptr<recob::PFParticle> > GetPrimaryNeutrinoPFPs();
  const std::vector<art::Ptr<recob::PFParticle> > GetSlicePrimary(std::vector<art::Ptr<recob::PFParticle> > const &slicePFPs);
  const art::Ptr<simb::MCParticle> GetTrueParticle(std::vector<art::Ptr<recob::Hit> > const &trackHits);
  const art::Ptr<simb::MCParticle> GetTrueShowerParticle(std::vector<art::Ptr<recob::Hit> > const &showerHits);
  unsigned int HierarchyPrimary(art::Ptr<recob::PFParticle> const &pfp);

  float HitPurity(art::Event const &e, art::Ptr<recob::Shower> const &shower);
  float HitPurity(art::Event const &e, art::Ptr<recob::Track> const &track);
  float Completeness(art::Event const &e, art::Ptr<recob::Shower> const &shower);
  float Completeness(art::Event const &e, art::Ptr<recob::Track> const &track);

  float TruedEdx(art::Ptr<simb::MCParticle> const &particle);
  float TrueTrackLength(art::Ptr<simb::MCParticle> const &particle);
  float TrackLength(art::Ptr<recob::Track> const &track);
  float MeanScatter(art::Ptr<recob::Track> const &track);
  float StdDevScatter(art::Ptr<recob::Track> const &track, const float &meanScatter);

  void TestingPlace(art::Event const &e);


  art::Handle<std::vector<simb::MCTruth> > eHandleNeutrinos;
  art::Handle<std::vector<simb::MCTruth> > eHandleCosmics;
  art::Handle<std::vector<simb::MCParticle> > eHandleParticles;
  art::Handle<std::vector<recob::Hit> > eHandleHits;
  art::Handle<std::vector<recob::Track> > eHandleTracks;
  art::Handle<std::vector<recob::Shower> > eHandleShowers;
  art::Handle<std::vector<recob::PFParticle> > eHandlePFPs;
  art::Handle<std::vector<sim::SimChannel> > eHandleSimChannels;
  art::Handle<std::vector<recob::Vertex> > eHandleVertices;
  art::Handle<std::vector<recob::Slice> > eHandleSlices;

  fhicl::ParameterSet fhiclP;

  detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataForJob();
  detinfo::DetectorPropertiesData propData = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataForJob();
  art::ServiceHandle<geo::Geometry> GeoHandler;
  geo::GeometryCore const* geom = art::ServiceHandle<geo::Geometry>()->provider();

  std::string fNuGenModuleLabel, fCosmicGenModuleLabel, fLArGeantModuleLabel,
    fPFParticleModuleLabel, fTrackModuleLabel, fShowerModuleLabel,
    fVertexModuleLabel, fHitsModuleLabel, fClusterModuleLabel,
    fSpacePointModuleLabel, fParticleIDModuleLabel, fSliceModuleLabel;
  bool fProcessCosmics;

  TTree *fEventTree;

  int nNeutrinos, nMCParticles, nUsedSlices, nSlices, nTracks, nShowers, nVertices;

  unsigned int eRunID, eSubRunID, eEventID;

  std::vector<int> CCNC, mode, interactionType, target, hitNuc, hitQuark;
  std::vector<bool> nuMu, nuE, nuMuBar, nuEBar;
  std::vector<float> W, X, Y, qSquared, pt, theta, nuEn, leptonP,
    trueVT, trueVX, trueVY, trueVZ, nuSurvival, nuLoss, cosContam, xCorrection,
    nuPx, nuPy, nuPz;
  std::map<int,int> MCPDGMap, hitsMap;

  std::vector<int> mc_trackID, mc_statusCode, mc_PDG, mc_mother, mc_nDaughters,
    mc_nTrajectoryPoints, mc_nuID;
  std::vector<bool> mc_isPrimary, mc_isCosmic;
  std::vector<std::vector<int> > mc_daughters, mc_daughtersPDG;
  std::vector<float> mc_x0, mc_y0, mc_z0, mc_xEnd, mc_yEnd, mc_zEnd, mc_pX0,
    mc_pY0, mc_pZ0, mc_energy0, mc_momentum, mc_pXEnd, mc_pYEnd, mc_pZEnd,
    mc_energyEnd, mc_mass, mc_phi, mc_theta, mc_length;
  std::vector<std::string> mc_process, mc_endProcess;

  std::vector<int> sl_primaryIndex, sl_nTracks, sl_nShowers;
  std::vector<float> sl_VX, sl_VY, sl_VZ, sl_completeness, sl_purity;

  std::vector<int> tr_index, tr_truePDG, tr_trueTrackID, tr_parent, tr_nHits,
    tr_sliceID;
  std::vector<float> tr_purity, tr_completeness, tr_x0, tr_y0, tr_z0, tr_xEnd,
    tr_yEnd, tr_zEnd, tr_momentum, tr_phi, tr_theta, tr_length, tr_trackPFOScore,
    tr_chi2Proton, tr_chi2Pion, tr_chi2Muon, tr_meanScatter, tr_stdDevScatter;
  std::vector<std::vector<long unsigned int> > tr_daughters;
  std::vector<bool> tr_isPrimary;
    
  std::vector<int> sh_index, sh_truePDG, sh_trueTrackID, sh_trueMotherPDG, 
    sh_trueMotherTrackID, sh_parent, sh_nHits, sh_sliceID;
  std::vector<float> sh_purity, sh_completeness, sh_x0, sh_y0, sh_z0, sh_energy, 
    sh_phi, sh_theta, sh_length, sh_openAngle, sh_dEdx, sh_trackPFOScore,
    sh_showerVertexSep; 
  std::vector<std::vector<long unsigned int> > sh_daughters;
  std::vector<bool> sh_isPrimary;
  std::vector<TVector3> sh_direction;

  std::vector<float> v_X, v_Y, v_Z;


  std::vector<unsigned int> PrimaryPFPCodes;
  std::map<int,int> PFPPrimaryMap, mcNuMap;
  std::map<int,bool> mcCosmicMap;
  std::set<int> nuPFPs, cosPFPs;

  
};


analysis::analysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fNuGenModuleLabel (p.get<std::string>("NuGenModuleLabel")),
  fCosmicGenModuleLabel (p.get<std::string>("CosmicGenModuleLabel")),
  fLArGeantModuleLabel (p.get<std::string>("LArGeantModuleLabel")),
  fPFParticleModuleLabel (p.get<std::string>("PFParticleModuleLabel")),
  fTrackModuleLabel (p.get<std::string>("TrackModuleLabel")),
  fShowerModuleLabel (p.get<std::string>("ShowerModuleLabel")),
  fVertexModuleLabel (p.get<std::string>("VertexModuleLabel")),
  fHitsModuleLabel (p.get<std::string>("HitsModuleLabel")),
  fClusterModuleLabel (p.get<std::string>("ClusterModuleLabel")),
  fSpacePointModuleLabel (p.get<std::string>("SpacePointModuleLabel")),
  fParticleIDModuleLabel (p.get<std::string>("ParticleIDModuleLabel")),
  fSliceModuleLabel (p.get<std::string>("SliceModuleLabel")),
  fProcessCosmics (p.get<bool>("ProcessCosmics"))
  // More initializers here.
  {
    // Call appropriate consumes<>() for any products to be retrieved by this module.
    art::ServiceHandle<art::TFileService> tfs;
    fhiclP = p;
  

    fEventTree = tfs->make<TTree>("EventTree","Data tree");

  // Eventwide data
  fEventTree->Branch("RunID",&eRunID);
  fEventTree->Branch("SubRunID",&eSubRunID);
  fEventTree->Branch("EventID",&eEventID);

  
  // True event data
  fEventTree->Branch("nNeutrinos",&nNeutrinos);
  fEventTree->Branch("CCNC",&CCNC);
  fEventTree->Branch("mode",&mode);
  fEventTree->Branch("interactionType",&interactionType);
  fEventTree->Branch("target",&target);
  fEventTree->Branch("hitNuc",&hitNuc);
  fEventTree->Branch("hitQuark",&hitQuark);
  fEventTree->Branch("nuMu",&nuMu);
  fEventTree->Branch("nuE",&nuE);
  fEventTree->Branch("nuMuBar",&nuMuBar);
  fEventTree->Branch("nuEBar",&nuEBar);
  fEventTree->Branch("W",&W);
  fEventTree->Branch("X",&X);
  fEventTree->Branch("Y",&Y);
  fEventTree->Branch("qSquared",&qSquared);
  fEventTree->Branch("pt",&pt);
  fEventTree->Branch("theta",&theta);
  fEventTree->Branch("nuEn",&nuEn);
  fEventTree->Branch("leptonP",&leptonP);
  fEventTree->Branch("trueVT",&trueVT);
  fEventTree->Branch("trueVX",&trueVX);
  fEventTree->Branch("trueVY",&trueVY);
  fEventTree->Branch("trueVZ",&trueVZ);
  fEventTree->Branch("nuPx",&nuPx);
  fEventTree->Branch("nuPy",&nuPy);
  fEventTree->Branch("nuPz",&nuPz);
  fEventTree->Branch("nuSurvival",&nuSurvival);
  fEventTree->Branch("nuLoss",&nuLoss);
  fEventTree->Branch("cosContam",&cosContam);
  fEventTree->Branch("mcPDGMap",&MCPDGMap);
  fEventTree->Branch("hitsMap",&hitsMap);
  fEventTree->Branch("xCorrection",&xCorrection);

  // MCParticle data
  fEventTree->Branch("nMCParticles",&nMCParticles);
  fEventTree->Branch("mc_trackID",&mc_trackID);
  fEventTree->Branch("mc_statusCode",&mc_statusCode);
  fEventTree->Branch("mc_PDG",&mc_PDG);
  fEventTree->Branch("mc_isPrimary",&mc_isPrimary);
  fEventTree->Branch("mc_isCosmic",&mc_isCosmic);
  fEventTree->Branch("mc_mother",&mc_mother);
  fEventTree->Branch("mc_nDaughters",&mc_nDaughters);
  fEventTree->Branch("mc_daughters",&mc_daughters);
  fEventTree->Branch("mc_daughtersPDG",&mc_daughtersPDG);
  fEventTree->Branch("mc_nTrajectoryPoints",&mc_nTrajectoryPoints);
  fEventTree->Branch("mc_x0",&mc_x0);
  fEventTree->Branch("mc_y0",&mc_y0);
  fEventTree->Branch("mc_z0",&mc_z0);
  fEventTree->Branch("mc_xEnd",&mc_xEnd);
  fEventTree->Branch("mc_yEnd",&mc_yEnd);
  fEventTree->Branch("mc_zEnd",&mc_zEnd);
  fEventTree->Branch("mc_pX0",&mc_pX0);
  fEventTree->Branch("mc_pY0",&mc_pY0);
  fEventTree->Branch("mc_pZ0",&mc_pZ0);
  fEventTree->Branch("mc_energy0",&mc_energy0);
  fEventTree->Branch("mc_momentum",&mc_momentum);
  fEventTree->Branch("mc_pXEnd",&mc_pXEnd);
  fEventTree->Branch("mc_pYEnd",&mc_pYEnd);
  fEventTree->Branch("mc_pZEnd",&mc_pZEnd);
  fEventTree->Branch("mc_energyEnd",&mc_energyEnd);
  fEventTree->Branch("mc_mass",&mc_mass);
  fEventTree->Branch("mc_phi",&mc_phi);
  fEventTree->Branch("mc_theta",&mc_theta);
  fEventTree->Branch("mc_length",&mc_length);
  fEventTree->Branch("mc_process",&mc_process);
  fEventTree->Branch("mc_endProcess",&mc_endProcess);
  fEventTree->Branch("mc_nuID",&mc_nuID);

  // Slice data
  fEventTree->Branch("nSlices",&nSlices);
  fEventTree->Branch("nUsedSlices",&nUsedSlices);
  fEventTree->Branch("sl_primaryIndex",&sl_primaryIndex);
  fEventTree->Branch("sl_VX",&sl_VX);
  fEventTree->Branch("sl_VY",&sl_VY);
  fEventTree->Branch("sl_VZ",&sl_VZ);
  fEventTree->Branch("sl_completeness",&sl_completeness);
  fEventTree->Branch("sl_purity",&sl_purity);
  fEventTree->Branch("sl_nTracks",&sl_nTracks);
  fEventTree->Branch("sl_nShowers",&sl_nShowers);
    
  // Track data
  fEventTree->Branch("nTracks",&nTracks);
  fEventTree->Branch("tr_index",&tr_index);
  fEventTree->Branch("tr_truePDG",&tr_truePDG);
  fEventTree->Branch("tr_trueTrackID",&tr_trueTrackID);
  fEventTree->Branch("tr_purity",&tr_purity);
  fEventTree->Branch("tr_completeness",&tr_completeness);
  fEventTree->Branch("tr_parent",&tr_parent);
  fEventTree->Branch("tr_daughters",&tr_daughters);
  fEventTree->Branch("tr_isPrimary",&tr_isPrimary);
  fEventTree->Branch("tr_x0",&tr_x0);
  fEventTree->Branch("tr_y0",&tr_y0);
  fEventTree->Branch("tr_z0",&tr_z0);
  fEventTree->Branch("tr_xEnd",&tr_xEnd);
  fEventTree->Branch("tr_yEnd",&tr_yEnd);
  fEventTree->Branch("tr_zEnd",&tr_zEnd);
  fEventTree->Branch("tr_momentum",&tr_momentum);
  fEventTree->Branch("tr_phi",&tr_phi);
  fEventTree->Branch("tr_theta",&tr_theta);
  fEventTree->Branch("tr_length",&tr_length);
  fEventTree->Branch("tr_nHits",&tr_nHits);
  fEventTree->Branch("tr_trackPFOScore",&tr_trackPFOScore);
  fEventTree->Branch("tr_chi2Proton",&tr_chi2Proton);
  fEventTree->Branch("tr_chi2Pion",&tr_chi2Pion);
  fEventTree->Branch("tr_chi2Muon",&tr_chi2Muon);
  fEventTree->Branch("tr_meanScatter",&tr_meanScatter);
  fEventTree->Branch("tr_stdDevScatter",&tr_stdDevScatter);
  fEventTree->Branch("tr_sliceID",&tr_sliceID);

  // Shower data
  fEventTree->Branch("nShowers",&nShowers);
  fEventTree->Branch("sh_index",&sh_index);
  fEventTree->Branch("sh_truePDG",&sh_truePDG);
  fEventTree->Branch("sh_trueTrackID",&sh_trueTrackID);
  fEventTree->Branch("sh_trueMotherPDG",&sh_trueMotherPDG);
  fEventTree->Branch("sh_trueMotherTrackID",&sh_trueMotherTrackID);
  fEventTree->Branch("sh_purity",&sh_purity);
  fEventTree->Branch("sh_completeness",&sh_completeness);
  fEventTree->Branch("sh_parent",&sh_parent);
  fEventTree->Branch("sh_daughters",&sh_daughters);
  fEventTree->Branch("sh_isPrimary",&sh_isPrimary);
  fEventTree->Branch("sh_x0",&sh_x0);
  fEventTree->Branch("sh_y0",&sh_y0);
  fEventTree->Branch("sh_z0",&sh_z0);
  fEventTree->Branch("sh_energy",&sh_energy);
  fEventTree->Branch("sh_phi",&sh_phi);
  fEventTree->Branch("sh_theta",&sh_theta);
  fEventTree->Branch("sh_length",&sh_length);
  fEventTree->Branch("sh_direction",&sh_direction,32000,0);
  fEventTree->Branch("sh_openAngle",&sh_openAngle);
  fEventTree->Branch("sh_dEdx",&sh_dEdx);
  fEventTree->Branch("sh_nHits",&sh_nHits);
  fEventTree->Branch("sh_trackPFOScore",&sh_trackPFOScore);
  fEventTree->Branch("sh_showerVertexSep",&sh_showerVertexSep);
  fEventTree->Branch("sh_sliceID",&sh_sliceID);

  // Vertex data
  fEventTree->Branch("nVertices",&nVertices);
  fEventTree->Branch("v_X",&v_X);
  fEventTree->Branch("v_Y",&v_Y);
  fEventTree->Branch("v_Z",&v_Z);
}

void analysis::SimulationProcessor(art::Event const &e)
{
  eRunID = e.run(); 
  eSubRunID = e.subRun();
  eEventID = e.event();
  
  nNeutrinos = eHandleNeutrinos->size();

  bool isCosmic = false;
  
  for(int nu_i = 0; nu_i < nNeutrinos; ++nu_i){
    const art::Ptr<simb::MCTruth> truthNeutrino(eHandleNeutrinos,nu_i);
    if(truthNeutrino.isNull()) continue;

    if(truthNeutrino->Origin() == 1) {
      const simb::MCNeutrino neutrino = truthNeutrino->GetNeutrino();
      const simb::MCParticle neutrinoParticle = neutrino.Nu();
      const simb::MCParticle lepton = neutrino.Lepton();
      
      detinfo::DetectorPropertiesData propD = art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e);
      util::GeometryUtilities geomU = util::GeometryUtilities(*geom,clockData,propD);

      for(auto pID : geom->IteratePlaneIDs()){
	xCorrection.push_back( propD.ConvertTicksToX(clockData.TPCG4Time2Tick(neutrinoParticle.T())- clockData.Time2Tick(clockData.BeamGateTime()),pID) - propD.ConvertTicksToX(0,pID) );
	break;
      }
      CCNC.push_back(neutrino.CCNC());
      mode.push_back(neutrino.Mode());
      interactionType.push_back(neutrino.InteractionType());
      target.push_back(neutrino.Target());
      hitNuc.push_back(neutrino.HitNuc());
      hitQuark.push_back(neutrino.HitQuark());

      int neutrinoType = neutrinoParticle.PdgCode();
      bool NuMu = false, NuE = false, NuMuBar = false, NuEBar = false;
      if(neutrinoType == 14) NuMu = true;
      else if(neutrinoType == 12) NuE = true;
      else if(neutrinoType == -14) NuMuBar = true;
      else if(neutrinoType == -12) NuEBar = true;
      
      nuMu.push_back(NuMu);
      nuE.push_back(NuE);
      nuMuBar.push_back(NuMuBar);
      nuEBar.push_back(NuEBar);
      
      theta.push_back(neutrino.Theta());
      W.push_back(neutrino.W());
      X.push_back(neutrino.X());
      Y.push_back(neutrino.Y());
      qSquared.push_back(neutrino.QSqr());
      pt.push_back(neutrino.Pt());
      nuEn.push_back(neutrinoParticle.E());
      leptonP.push_back(lepton.E());
  
      trueVT.push_back(neutrinoParticle.T());
      trueVX.push_back(neutrinoParticle.Vx());
      trueVY.push_back(neutrinoParticle.Vy());
      trueVZ.push_back(neutrinoParticle.Vz());
      nuPx.push_back(neutrinoParticle.Px());
      nuPy.push_back(neutrinoParticle.Py());
      nuPz.push_back(neutrinoParticle.Pz());

    }
  
    art::FindManyP<simb::MCParticle> nuParticleAssn(eHandleNeutrinos,e,fLArGeantModuleLabel);
    std::vector<art::Ptr<simb::MCParticle> > particles = nuParticleAssn.at(truthNeutrino.key());
    
    for(unsigned int part_i = 0; part_i < particles.size(); ++part_i) {
      const art::Ptr<simb::MCParticle> particle = particles[part_i];
      mcNuMap[particle->TrackId()] = truthNeutrino->GetNeutrino().Nu().TrackId();
      MCParticleProcessor(e,particle,isCosmic);
    }
    HitProcessor(e,truthNeutrino->GetNeutrino().Nu().TrackId());
  }

  if(fProcessCosmics){
    isCosmic = true;
    
    art::FindManyP<simb::MCParticle> cosmicParticleAssn(eHandleCosmics,e,fLArGeantModuleLabel);
    
    for(unsigned int cosm_i = 0; cosm_i < eHandleCosmics->size(); ++cosm_i){
      const art::Ptr<simb::MCTruth> cosmic(eHandleCosmics,cosm_i);
      std::vector<art::Ptr<simb::MCParticle> > cosParticles = cosmicParticleAssn.at(cosmic.key());
      
      for(unsigned int part_i = 0; part_i < cosParticles.size(); ++part_i) {
	const art::Ptr<simb::MCParticle> particle = cosParticles[part_i];
	MCParticleProcessor(e,particle,isCosmic);
      }
    }
  }
}

void analysis::MCParticleProcessor(art::Event const &e, 
				   art::Ptr<simb::MCParticle> const &particle, 
				   bool const &isCosmic) 
{
  bool isPrimary = false;

  mc_trackID.push_back(particle->TrackId());
  mc_statusCode.push_back(particle->StatusCode());
  mc_PDG.push_back(particle->PdgCode());
  mc_mother.push_back(particle->Mother());

  mcCosmicMap[mc_trackID.back()] = isCosmic;
  mc_isCosmic.push_back(isCosmic);

  if(mc_mother.back() == 0) isPrimary = true;
  mc_isPrimary.push_back(isPrimary);

  mc_nDaughters.push_back(particle->NumberDaughters());
  std::vector<int> dts, dtsPDG;
  for(int i = 0; i < mc_nDaughters.back(); ++i) {
    int code = particle->Daughter(i);
    dts.push_back(code);
    dtsPDG.push_back(MCPDGMap[code]);
  }
  mc_daughters.push_back(dts);
  mc_daughtersPDG.push_back(dtsPDG);

  mc_nTrajectoryPoints.push_back(particle->NumberTrajectoryPoints());
  mc_x0.push_back(particle->Vx());
  mc_y0.push_back(particle->Vy());
  mc_z0.push_back(particle->Vz());
  mc_xEnd.push_back(particle->EndX());
  mc_yEnd.push_back(particle->EndY());
  mc_zEnd.push_back(particle->EndZ());
  mc_pX0.push_back(particle->Px());
  mc_pY0.push_back(particle->Py());
  mc_pZ0.push_back(particle->Pz());
  mc_energy0.push_back(particle->E());
  mc_momentum.push_back(particle->P());
  mc_pXEnd.push_back(particle->EndPx());
  mc_pYEnd.push_back(particle->EndPy());
  mc_pZEnd.push_back(particle->EndPz());
  mc_energyEnd.push_back(particle->EndE());
  mc_mass.push_back(particle->Mass());
  mc_phi.push_back(particle->Momentum().Phi());
  mc_theta.push_back(particle->Momentum().Theta());
  mc_length.push_back(TrueTrackLength(particle));
  //  mc_dEdx.push_back(TruedEdx(particle));

  mc_process.push_back(particle->Process());
  mc_endProcess.push_back(particle->EndProcess());

  mc_nuID.push_back(mcNuMap[particle->TrackId()]);

  ++nMCParticles;
}

void analysis::ReconstructionProcessor(art::Event const &e)
{
  art::FindManyP<recob::Vertex> pfpVertexAssoc(eHandlePFPs,e,fVertexModuleLabel);
  art::FindManyP<recob::Track> pfpTrackAssoc(eHandlePFPs,e,fTrackModuleLabel);
  art::FindManyP<recob::Shower> pfpShowerAssoc(eHandlePFPs,e,fShowerModuleLabel);
  art::FindManyP<recob::PFParticle> slicePFPAssoc(eHandleSlices,e,fSliceModuleLabel);
  art::FindManyP<recob::Hit> sliceHitAssoc(eHandleSlices,e,fSliceModuleLabel);
  art::FindManyP<larpandoraobj::PFParticleMetadata> fmpfpmd(eHandlePFPs,e,fPFParticleModuleLabel);
  
  nSlices = eHandleSlices->size();
  
  for(unsigned int slice_i = 0; slice_i < eHandleSlices->size(); ++slice_i){
    const art::Ptr<recob::Slice> slice(eHandleSlices,slice_i);
    std::vector<art::Ptr<recob::PFParticle> > slicePFPs = slicePFPAssoc.at(slice.key());
    std::vector<art::Ptr<recob::Hit> > sliceHits = sliceHitAssoc.at(slice.key());

    std::vector<art::Ptr<recob::PFParticle> > primaryPFPs = GetSlicePrimary(slicePFPs);
    if(primaryPFPs.size() == 0) continue;
    else if(primaryPFPs.size() > 1) {std::cout << "Multiple primaries!!!" << std::endl; continue;}

    SliceHits(e,sliceHits);
    
    art::Ptr<recob::PFParticle> nuPrimary = primaryPFPs[0];

    std::vector<art::Ptr<recob::Vertex> > pfpVertices = pfpVertexAssoc.at(nuPrimary.key());
    if(pfpVertices.size() != 1) {std::cout << "Multiple vertices!!!" << std::endl; return;}
    
    sl_primaryIndex.push_back(nuPrimary->Self());
    sl_VX.push_back(pfpVertices[0]->position().X());
    sl_VY.push_back(pfpVertices[0]->position().Y());
    sl_VZ.push_back(pfpVertices[0]->position().Z());
    
    std::vector<long unsigned int> nuDaughters = nuPrimary->Daughters();
    int sliceTracks = 0, sliceShowers = 0;
    
    for(auto id : nuDaughters){
      const art::Ptr<recob::PFParticle> daughter = GetPFP(id);

      float pfoScore = -999;

      const std::vector< art::Ptr<larpandoraobj::PFParticleMetadata> > pfpMetaVec = fmpfpmd.at(daughter->Self());
      if(pfpMetaVec.size() != 1) std::cout << "PFP meta data > 1" << std::endl;
      for (auto const& pfpMeta: pfpMetaVec)  {
	larpandoraobj::PFParticleMetadata::PropertiesMap propertiesMap = pfpMeta->GetPropertiesMap();
	auto const pfpTrackScoreIter = propertiesMap.find("TrackScore");
	pfoScore = pfpTrackScoreIter->second;
      }
      
      if(daughter->PdgCode() == 13) {
	std::vector<art::Ptr<recob::Track> > daughterTracks = pfpTrackAssoc.at(daughter.key());
	if(daughterTracks.size() != 1) continue;
	else{
	  tr_sliceID.push_back(nUsedSlices);
	  tr_trackPFOScore.push_back(pfoScore);
	  TrackProcessor(e,daughter,daughterTracks[0]);
	  ++sliceTracks;
	}
      }
      else if(daughter->PdgCode() == 11) {
	std::vector<art::Ptr<recob::Shower> > daughterShowers = pfpShowerAssoc.at(daughter.key());
	if(daughterShowers.size() != 1) continue;
	else {
	  sh_sliceID.push_back(nUsedSlices);
	  sh_trackPFOScore.push_back(pfoScore);
	  ShowerProcessor(e,daughter,daughterShowers[0]);
	  ++sliceShowers;
	}
      }
    }
    sl_nTracks.push_back(sliceTracks);
    sl_nShowers.push_back(sliceShowers);
    ++nUsedSlices;
  }
}

void analysis::TrackProcessor(art::Event const &e,
			      art::Ptr<recob::PFParticle> const &pfp,
			      art::Ptr<recob::Track> const &track)
{
  bool isPrimary = false;

  art::FindManyP<recob::Hit> trackHitsAssoc(eHandleTracks,e,fTrackModuleLabel);
  std::vector<art::Ptr<recob::Hit> > trackHits = trackHitsAssoc.at(track.key());
  
  art::FindManyP<anab::ParticleID>pidTrackAssoc(eHandleTracks,e,fParticleIDModuleLabel);
  std::vector<art::Ptr<anab::ParticleID> > trackPIDs = pidTrackAssoc.at(track.key());
  art::Ptr<anab::ParticleID> trackPID = trackPIDs[2];

  art::Ptr<simb::MCParticle> mc = GetTrueParticle(trackHits);
  if(mc.isNull()) return;

  tr_index.push_back(pfp->Self());
  tr_truePDG.push_back(mc->PdgCode());
  tr_trueTrackID.push_back(mc->TrackId());
  tr_purity.push_back(HitPurity(e,track));
  tr_completeness.push_back(Completeness(e,track));
  tr_parent.push_back(pfp->Parent());
  tr_daughters.push_back(pfp->Daughters());
  for(unsigned int i=0; i<PrimaryPFPCodes.size(); ++i) {
    if(pfp->Parent() == PrimaryPFPCodes[i]) isPrimary = true;
  }
  tr_isPrimary.push_back(isPrimary);
  tr_x0.push_back(track->Start().X());
  tr_y0.push_back(track->Start().Y());
  tr_z0.push_back(track->Start().Z());
  tr_xEnd.push_back(track->End().X());
  tr_yEnd.push_back(track->End().Y());
  tr_zEnd.push_back(track->End().Z());
  tr_momentum.push_back(track->StartMomentum());
  tr_phi.push_back(track->StartDirection().Phi());
  tr_theta.push_back(track->StartDirection().Theta());
  tr_length.push_back(TrackLength(track));
  tr_nHits.push_back(trackHits.size());
  tr_chi2Proton.push_back(trackPID->Chi2Proton());
  tr_chi2Pion.push_back(trackPID->Chi2Pion());
  tr_chi2Muon.push_back(trackPID->Chi2Muon());
  tr_meanScatter.push_back(MeanScatter(track));
  tr_stdDevScatter.push_back(StdDevScatter(track,tr_meanScatter.back()));

  ++nTracks;
}

void analysis::ShowerProcessor(art::Event const &e,
			       art::Ptr<recob::PFParticle> const &pfp,
			       art::Ptr<recob::Shower> const &shower)
{
  bool isPrimary = false;
  
  art::FindManyP<recob::Hit> showerHitsAssoc(eHandleShowers,e,fShowerModuleLabel);
  std::vector<art::Ptr<recob::Hit> > showerHits = showerHitsAssoc.at(shower.key());

  art::Ptr<simb::MCParticle> mc = GetTrueShowerParticle(showerHits);
  if(mc.isNull()) return;

  sh_index.push_back(pfp->Self());
  sh_truePDG.push_back(mc->PdgCode());
  sh_trueTrackID.push_back(mc->TrackId());

  art::Ptr<simb::MCParticle> mcMother = GetMCParticle(mc->Mother());
  if(mcMother.isNonnull()){
    sh_trueMotherPDG.push_back(mcMother->PdgCode());
    sh_trueMotherTrackID.push_back(mcMother->TrackId());
  }
  else{
    sh_trueMotherPDG.push_back(-999);
    sh_trueMotherTrackID.push_back(-999);
  }

  sh_purity.push_back(HitPurity(e,shower));
  sh_completeness.push_back(Completeness(e,shower));
  sh_parent.push_back(pfp->Parent());
  sh_daughters.push_back(pfp->Daughters());
  for(unsigned int i=0; i<PrimaryPFPCodes.size(); ++i) {
    if(pfp->Parent() == PrimaryPFPCodes[i]) isPrimary = true;
  }
  sh_isPrimary.push_back(isPrimary);
  sh_x0.push_back(shower->ShowerStart().X());
  sh_y0.push_back(shower->ShowerStart().Y());
  sh_z0.push_back(shower->ShowerStart().Z());
      
  int best_plane = shower->best_plane();
  if(best_plane < 0 || best_plane > 2) best_plane = 2;
  std::vector<double> showerEnergy = shower->Energy();
  std::vector<double> showerdEdx = shower->dEdx();
  if(showerEnergy.size() == 0) sh_energy.push_back(-999);
  else sh_energy.push_back(showerEnergy[best_plane]);
  if(showerdEdx.size() == 0) sh_dEdx.push_back(-999);
  else sh_dEdx.push_back(showerdEdx[best_plane]);

  sh_phi.push_back(shower->Direction().Phi());
  sh_theta.push_back(shower->Direction().Theta());
  sh_length.push_back(shower->Length());
  sh_direction.push_back(shower->Direction());
  sh_openAngle.push_back(shower->OpenAngle());
  sh_nHits.push_back(showerHits.size());
  
  TVector3 vertex(sl_VX.back(),sl_VY.back(),sl_VZ.back());
  sh_showerVertexSep.push_back((vertex - shower->ShowerStart()).Mag());

  ++nShowers;
}

void analysis::HitProcessor(art::Event const &e, int const &nuID)
{
  art::FindManyP<recob::Track> hitTrackAssoc(eHandleHits,e,fTrackModuleLabel);
  art::FindManyP<recob::Shower> hitShowerAssoc(eHandleHits,e,fShowerModuleLabel);
  art::FindManyP<recob::PFParticle> trackPFOAssoc(eHandleTracks,e,fTrackModuleLabel);
  art::FindManyP<recob::PFParticle> showerPFOAssoc(eHandleShowers,e,fShowerModuleLabel);
  

  int hitLocations[6] = {0,0,0,0,0,0};

  for(unsigned int hit_i = 0; hit_i < eHandleHits->size(); ++hit_i) {
    const art::Ptr<recob::Hit> hit(eHandleHits,hit_i);
    std::vector<art::Ptr<recob::Track> > track = hitTrackAssoc.at(hit.key());
    std::vector<art::Ptr<recob::Shower> > shower = hitShowerAssoc.at(hit.key());
    int trueTrackID = -999, objectIndex = -999;

    trueTrackID = TruthMatchUtils::TrueParticleID(clockData,hit,true);

    if(shower.size() == 0 && track.size() != 0){
      std::vector<art::Ptr<recob::PFParticle> > pfo = trackPFOAssoc.at(track[0].key());
      objectIndex = pfo[0]->Self();      
    }
    else if(track.size() == 0 && shower.size() != 0){
      std::vector<art::Ptr<recob::PFParticle> > pfo = showerPFOAssoc.at(shower[0].key());
      objectIndex = pfo[0]->Self();
    }
    else if(shower.size() == 0 && track.size() == 0) {}


    if((!mcCosmicMap[trueTrackID]) && mcNuMap[trueTrackID] == nuID
       && nuPFPs.count(objectIndex)==1) hitLocations[0]++;
    else if((!mcCosmicMap[trueTrackID]) && mcNuMap[trueTrackID] == nuID
	    && cosPFPs.count(objectIndex)==1) hitLocations[1]++;
    else if((!mcCosmicMap[trueTrackID]) && mcNuMap[trueTrackID] == nuID) 
      hitLocations[2]++;

    else if((mcCosmicMap[trueTrackID]) 
	    && nuPFPs.count(objectIndex)==1) hitLocations[3]++;
    else if((mcCosmicMap[trueTrackID]) 
	    && cosPFPs.count(objectIndex)==1) hitLocations[4]++;
    else if((mcCosmicMap[trueTrackID])) hitLocations[5]++;

  }

  nuSurvival.push_back(hitLocations[0]/static_cast<float>(hitLocations[0]+hitLocations[1]+hitLocations[2]));
  nuLoss.push_back(hitLocations[1]/static_cast<float>(hitLocations[0]+hitLocations[1]+hitLocations[2]));
  cosContam.push_back(hitLocations[3]/static_cast<float>(hitLocations[0]+hitLocations[3]));
  
}

void analysis::SliceHits(art::Event const &e, std::vector<art::Ptr<recob::Hit> > const &sliceHits)
{
  int nuTrueHits = 0, nuRecoHits = 0;

  for(unsigned int hit_i = 0; hit_i < eHandleHits->size(); ++hit_i) {
    const art::Ptr<recob::Hit> hit(eHandleHits,hit_i);
    int trueTrackID = -999;
    trueTrackID = TruthMatchUtils::TrueParticleID(clockData,hit,true);
    
    if(!mcCosmicMap[trueTrackID]) nuTrueHits++;
  }

  for(auto hit : sliceHits){
    int trueTrackID = -999;
    trueTrackID = TruthMatchUtils::TrueParticleID(clockData,hit,true);
    
    if(!mcCosmicMap[trueTrackID]) nuRecoHits++;
  }   

  sl_completeness.push_back(nuRecoHits/static_cast<float>(nuTrueHits));
  sl_purity.push_back(nuRecoHits/static_cast<float>(sliceHits.size()));  
}

void analysis::VertexProcessor(art::Event const &e)
{
  nVertices = eHandleVertices->size();
  
  for(int v_i = 0; v_i < nVertices; ++v_i) {
    const art::Ptr<recob::Vertex> vertex(eHandleVertices,v_i);
    
    v_X.push_back(vertex->position().X());
    v_Y.push_back(vertex->position().Y());
    v_Z.push_back(vertex->position().Z());
  }
}

void analysis::SetupMaps(art::Event const &e) 
{
  MCPDGMap.clear();
  hitsMap.clear();
  mcCosmicMap.clear();
  PFPPrimaryMap.clear();
  nuPFPs.clear();
  cosPFPs.clear();
  mcNuMap.clear();

  for(unsigned int part_i = 0; part_i < eHandleParticles->size(); ++part_i) {
    const art::Ptr<simb::MCParticle> particle(eHandleParticles,part_i);
    MCPDGMap[particle->TrackId()] = particle->PdgCode();
  }
  
  for(unsigned hit_i = 0; hit_i < eHandleHits->size(); ++hit_i) {
    const art::Ptr<recob::Hit> hit(eHandleHits,hit_i);
    hitsMap[TruthMatchUtils::TrueParticleID(clockData,hit,true)]++;
  }

  for(unsigned int pfp_i = 0; pfp_i < eHandlePFPs->size(); ++pfp_i) {
    const art::Ptr<recob::PFParticle> pfp(eHandlePFPs,pfp_i);
    if(pfp->IsPrimary()){
      if(pfp->PdgCode() == 12 || pfp->PdgCode() == -12 
	 || pfp->PdgCode() == 14 || pfp->PdgCode() == -14){
	PrimaryPFPCodes.push_back(pfp->Self());
	nuPFPs.insert(pfp->Self());
	PFPPrimaryMap[pfp->Self()] = pfp->Self();
      }
      else cosPFPs.insert(pfp->Self());
    }
  }
  
  bool regress = true;
  
  while(regress) {
    regress = false;
  
    for(unsigned int pfp_i = 0; pfp_i < eHandlePFPs->size(); ++pfp_i) {
      const art::Ptr<recob::PFParticle> pfp(eHandlePFPs,pfp_i);
      if(pfp->IsPrimary()) continue;
      if(nuPFPs.count(pfp->Parent()) == 1) {nuPFPs.insert(pfp->Self()); PFPPrimaryMap[pfp->Self()] = PFPPrimaryMap[pfp->Parent()];}
      else if(cosPFPs.count(pfp->Parent()) == 1) cosPFPs.insert(pfp->Self());
      else regress = true;
    }
  }
}

const art::Ptr<simb::MCParticle> analysis::GetMCParticle(int const &trackID)
{
  for(unsigned int part_i = 0; part_i < eHandleParticles->size(); ++part_i) {
    const art::Ptr<simb::MCParticle> particle(eHandleParticles,part_i);
    if(particle->TrackId() == trackID) return particle;
  }

  const art::Ptr<simb::MCParticle> nullReturn;
  return nullReturn;
}

const art::Ptr<recob::PFParticle> analysis::GetPrimaryPFP() 
{
  for(unsigned int pfp_i = 0; pfp_i < eHandlePFPs->size(); ++pfp_i) {
    const art::Ptr<recob::PFParticle> pfp(eHandlePFPs,pfp_i);
    if(pfp->IsPrimary()) return pfp;
  }
  
  art::Ptr<recob::PFParticle> nullReturn;
  return nullReturn;
}

const art::Ptr<recob::PFParticle> analysis::GetPFP(unsigned int const &index) 
{
  for(unsigned int pfp_i = 0; pfp_i < eHandlePFPs->size(); ++pfp_i) {
    const art::Ptr<recob::PFParticle> pfp(eHandlePFPs,pfp_i);
    if(pfp->Self()==index) return pfp;
  }
  
  art::Ptr<recob::PFParticle> nullReturn;
  return nullReturn;
}

const std::vector<art::Ptr<recob::PFParticle> > analysis::GetPrimaryNeutrinoPFPs() 
{
  std::vector<art::Ptr<recob::PFParticle> > primaries;
  for(unsigned int pfp_i = 0; pfp_i < eHandlePFPs->size(); ++pfp_i) {
    const art::Ptr<recob::PFParticle> pfp(eHandlePFPs,pfp_i);
    if(pfp->IsPrimary()){
      if(pfp->PdgCode() == 12 || pfp->PdgCode() == -12 
	 || pfp->PdgCode() == 14 || pfp->PdgCode() == -14){
	primaries.push_back(pfp);
      }
    }
  }
  return primaries;
}

const std::vector<art::Ptr<recob::PFParticle> > analysis::GetSlicePrimary(std::vector<art::Ptr<recob::PFParticle> > const &slicePFPs)
{
  std::vector<art::Ptr<recob::PFParticle> > primaries;
  for(auto pfp : slicePFPs) {
    if(pfp->IsPrimary()){
      primaries.push_back(pfp);
    }
  }
  return primaries;
}

const art::Ptr<simb::MCParticle> analysis::GetTrueParticle(std::vector<art::Ptr<recob::Hit> > const &trackHits)
{
  int particleID = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData,trackHits,true);

  for(unsigned int i = 0; i < eHandleParticles->size(); i++) {
      const::art::Ptr<simb::MCParticle> particle(eHandleParticles,i);
      if(particle->TrackId() == particleID) return particle;
    }

  art::Ptr<simb::MCParticle> nullReturn;
  return nullReturn;
}

const art::Ptr<simb::MCParticle> analysis::GetTrueShowerParticle(std::vector<art::Ptr<recob::Hit> > const &showerHits)
{
  art::Ptr<simb::MCParticle> original = GetTrueParticle(showerHits);
  if(original.isNull()) return original;

  int trackID = original->Mother();

  while(MCPDGMap[trackID] == 11 || 
	MCPDGMap[trackID] == -11 || 
	MCPDGMap[trackID] == 22) {
    original = GetMCParticle(trackID);
    trackID = original->Mother();
  }
  return original;
}

unsigned int analysis::HierarchyPrimary(art::Ptr<recob::PFParticle> const &pfp)
{
  bool regress = true;
  art::Ptr<recob::PFParticle> parent = pfp;
  
  while(regress){
    if(parent->IsPrimary()) return parent->Self();
    else parent = GetPFP(parent->Parent());
  }
  return 999999;  
}

float analysis::HitPurity(art::Event const &e, art::Ptr<recob::Shower> const &shower)
{
  art::FindManyP<recob::Hit> showerHitsAssn(eHandleShowers,e,fShowerModuleLabel);
  std::vector<art::Ptr<recob::Hit> > showerHits = showerHitsAssn.at(shower.key());
  art::Ptr<simb::MCParticle> mc = GetTrueShowerParticle(showerHits);
  if(mc.isNull()) return -999;

  int trackID = mc->TrackId();

  std::map<int,int> showerHitsMap;

  for(unsigned int i=0; i<showerHits.size(); ++i) {
    showerHitsMap[TruthMatchUtils::TrueParticleID(clockData,showerHits[i],true)]++;
  }
  return showerHitsMap[trackID]/static_cast<float>(showerHits.size());
}

float analysis::HitPurity(art::Event const &e, art::Ptr<recob::Track> const &track)
{
  art::FindManyP<recob::Hit> trackHitsAssn(eHandleTracks,e,fTrackModuleLabel);
  std::vector<art::Ptr<recob::Hit> > trackHits = trackHitsAssn.at(track.key());
  art::Ptr<simb::MCParticle> mc = GetTrueParticle(trackHits);
  if(mc.isNull()) return -999;
  int trackID = mc->TrackId();

  std::map<int,int> trackHitsMap;

  for(unsigned int i=0; i<trackHits.size(); ++i) {
    trackHitsMap[TruthMatchUtils::TrueParticleID(clockData,trackHits[i],true)]++;
  }
  return trackHitsMap[trackID]/static_cast<float>(trackHits.size());
}

float analysis::Completeness(art::Event const &e, art::Ptr<recob::Shower> const &shower)
{
  art::FindManyP<recob::Hit> showerHitsAssn(eHandleShowers,e,fShowerModuleLabel);
  std::vector<art::Ptr<recob::Hit> > showerHits = showerHitsAssn.at(shower.key());

  art::Ptr<simb::MCParticle> mc = GetTrueShowerParticle(showerHits);
  if(mc.isNull()) return -999;
  int trackID = mc->TrackId();

  std::map<int,int> showerHitsMap;

  for(unsigned int i=0; i<showerHits.size(); ++i) {
    showerHitsMap[TruthMatchUtils::TrueParticleID(clockData,showerHits[i],true)]++;
  }
  
  return showerHitsMap[trackID]/static_cast<float>(hitsMap[trackID]);
}

float analysis::Completeness(art::Event const &e, art::Ptr<recob::Track> const &track)
{
  art::FindManyP<recob::Hit> trackHitsAssn(eHandleTracks,e,fTrackModuleLabel);
  std::vector<art::Ptr<recob::Hit> > trackHits = trackHitsAssn.at(track.key());

  art::Ptr<simb::MCParticle> mc = GetTrueParticle(trackHits);
  if(mc.isNull()) return -999;
  int trackID = mc->TrackId();

  std::map<int,int> trackHitsMap;

  for(unsigned int i=0; i<trackHits.size(); ++i) {
    trackHitsMap[TruthMatchUtils::TrueParticleID(clockData,trackHits[i],true)]++;
  }
  
  return trackHitsMap[trackID]/static_cast<float>(hitsMap[trackID]);
}

float analysis::TruedEdx(art::Ptr<simb::MCParticle> const &particle)
{
  const int tID = particle->TrackId();
  std::vector<TVector3> stubPoints;

  float length = 0, energy = 0;
  unsigned int nTrajPoints = particle->NumberTrajectoryPoints();

  stubPoints.push_back(particle->Position(0).Vect());
  // particle->Position(0).Vect().Print();
  
  for(unsigned int point = 1; point < nTrajPoints; ++point) {
    TVector3 diff = particle->Position(point).Vect() - particle->Position(point-1).Vect();
    length += diff.Mag();
    // particle->Position(point).Vect().Print();
    if(length > 3) break;
    stubPoints.push_back(particle->Position(point).Vect());
  }

  for(unsigned int sim_i = 0; sim_i < eHandleSimChannels->size(); ++sim_i){
    const art::Ptr<sim::SimChannel> simChannel(eHandleSimChannels,sim_i);
    auto simIDEMap = simChannel->TDCIDEMap();    
    auto chID = simChannel->Channel();
    art::ServiceHandle<geo::Geometry> geom;
    auto wireIDVector = geom->ChannelToWire(chID);

    for(auto x : simIDEMap){
      auto ideVec = x.second;
      for(auto ide : ideVec){
	if(ide.trackID == tID) {
	  energy += ide.energy;
	  // std::cout <<  "IDE: ";
	  // std::cout << wireIDVector[0].toString();
	  // std::cout << " -- Energy: " << ide.energy
	  // 	    << " x: " << ide.x << " y: " << ide.y << " z: " << ide.z << std::endl;
	}
      }
    }  
  }
  // if(isPrimary && PDG == 13 && energy != 0) std::cout << "TRACKID: " << tID << "\tEnergy: " << energy << "\tTrueEnergy: " << particle->Momentum(0).Vect().Mag() << std::endl;
  return 0.0;
}

float analysis::TrueTrackLength(art::Ptr<simb::MCParticle> const &particle)
{
  float length = 0;
  unsigned int nTrajPoints = particle->NumberTrajectoryPoints();

  if(nTrajPoints < 2) return length;

  for(unsigned int point = 1; point < nTrajPoints; ++point) {
    TVector3 l = particle->Position(point).Vect();
    if(l.X() > 200 || l.X() < -200 || l.Y() > 200 || l.Y() < -200 ||
       l.Z() > 500 || l.Z() < 0) break;

    TVector3 diff = particle->Position(point).Vect() - particle->Position(point-1).Vect();
    length += TMath::Sqrt(diff.Mag2());
  }

  return length;
}

float analysis::TrackLength(art::Ptr<recob::Track> const &track)
{
  float length = 0;
  int nTrajPoints = track->NumberTrajectoryPoints();

  if(nTrajPoints < 2) return length;

  for(int point = 1; point < nTrajPoints; ++point) {
    TVector3 l = track->LocationAtPoint<TVector3>(point);
    if(l.X() == -999 || l.Y() == -999 || l.Z() == -999) break;

    TVector3 diff = track->LocationAtPoint<TVector3>(point) - track->LocationAtPoint<TVector3>(point-1);
    length += TMath::Sqrt(diff.Mag2());
  }
  return length;
}

float analysis::MeanScatter(art::Ptr<recob::Track> const &track)
{
  float sumScatAngle = 0;
  int N = 0;
  int nTrajPoints = track->NumberTrajectoryPoints();

  if(nTrajPoints < 2) return sumScatAngle;

  for(int point = 1; point < nTrajPoints; ++point) {
    TVector3 iDir = track->DirectionAtPoint<TVector3>(point-1);
    TVector3 fDir = track->DirectionAtPoint<TVector3>(point);
    if(fDir.X() == -999 || fDir.Y() == -999 || fDir.Z() == -999) break;

    sumScatAngle += TMath::RadToDeg() * fDir.Angle(iDir);
    N+=1;
  }
  float mean = -999;
  if(N>0) mean = sumScatAngle / static_cast<float>(N);

  return mean;
}

float analysis::StdDevScatter(art::Ptr<recob::Track> const &track, const float &meanScatter)
{
  float sumSqDiff = 0;
  int N = 0;
  int nTrajPoints = track->NumberTrajectoryPoints();

  if(nTrajPoints < 2) return sumSqDiff;

  for(int point = 1; point < nTrajPoints; ++point) {
    TVector3 iDir = track->DirectionAtPoint<TVector3>(point-1);
    TVector3 fDir = track->DirectionAtPoint<TVector3>(point);
    if(fDir.X() == -999 || fDir.Y() == -999 || fDir.Z() == -999) break;

    float scatAngle = TMath::RadToDeg() * fDir.Angle(iDir);
    
    sumSqDiff += TMath::Power((scatAngle-meanScatter),2);
    N+=1;
  }
  float stdDev = -999;
  if(N>1){
    float var = sumSqDiff / static_cast<float>(N-1);
    stdDev = TMath::Sqrt(var);
  }
  return stdDev;
}

void analysis::TestingPlace(art::Event const& e) 
{
}

void analysis::ClearData()
{
  eRunID = -999; eSubRunID = -999; eEventID = -999;
  nNeutrinos = 0; nMCParticles = 0; nUsedSlices = 0; nSlices = 0;
  nTracks = 0; nShowers = 0; nVertices = 0;

  CCNC.clear(); mode.clear(); interactionType.clear(); target.clear(); hitNuc.clear(); 
  hitQuark.clear();
  nuMu.clear(); nuE.clear(); nuMuBar.clear(); nuEBar.clear();
  W.clear(); X.clear(); Y.clear(); qSquared.clear(); pt.clear(); theta.clear(); 
  nuEn.clear(); leptonP.clear(); trueVT.clear(); trueVX.clear(); trueVY.clear(); trueVZ.clear(); 
  nuSurvival.clear(); nuLoss.clear(); cosContam.clear(); xCorrection.clear();
  nuPx.clear(); nuPy.clear(); nuPz.clear(); 

  mc_trackID.clear(); mc_statusCode.clear(); mc_PDG.clear(); mc_mother.clear(); 
  mc_nDaughters.clear(); mc_nTrajectoryPoints.clear(); mc_nuID.clear();
  mc_isPrimary.clear(); mc_isCosmic.clear();
  mc_daughters.clear(); mc_daughtersPDG.clear();
  mc_x0.clear(); mc_y0.clear(); mc_z0.clear(); mc_xEnd.clear(); mc_yEnd.clear(); 
  mc_zEnd.clear(); mc_pX0.clear(); mc_pY0.clear(); mc_pZ0.clear(); mc_energy0.clear(); 
  mc_momentum.clear(); mc_pXEnd.clear(); mc_pYEnd.clear(); mc_pZEnd.clear();
  mc_energyEnd.clear(); mc_mass.clear(); mc_phi.clear(); mc_theta.clear(); 
  mc_length.clear();
  mc_process.clear(); mc_endProcess.clear();

  sl_primaryIndex.clear(); sl_nTracks.clear(); sl_nShowers.clear();
  sl_VX.clear(); sl_VY.clear(); sl_VZ.clear(); sl_completeness.clear(); sl_purity.clear();

  tr_index.clear(); tr_truePDG.clear(); tr_trueTrackID.clear(); tr_parent.clear(); 
  tr_nHits.clear(); tr_sliceID.clear();
  tr_purity.clear(); tr_completeness.clear(); tr_x0.clear(); tr_y0.clear(); tr_z0.clear();
  tr_xEnd.clear(); tr_yEnd.clear(); tr_zEnd.clear(); tr_momentum.clear(); tr_phi.clear();
  tr_theta.clear(); tr_length.clear(); tr_trackPFOScore.clear(); tr_chi2Proton.clear(); 
  tr_chi2Pion.clear(); tr_chi2Muon.clear(); tr_meanScatter.clear(); 
  tr_stdDevScatter.clear();
  tr_daughters.clear();
  tr_isPrimary.clear();

  sh_index.clear(); sh_truePDG.clear(); sh_trueTrackID.clear(); sh_trueMotherPDG.clear();
  sh_trueMotherTrackID.clear(); sh_parent.clear(); sh_nHits.clear(); sh_sliceID.clear();
  sh_purity.clear(); sh_completeness.clear(); sh_x0.clear(); sh_y0.clear(); sh_z0.clear();
  sh_energy.clear(); sh_phi.clear(); sh_theta.clear(); sh_length.clear(); 
  sh_openAngle.clear(); sh_dEdx.clear(); sh_trackPFOScore.clear(); 
  sh_showerVertexSep.clear();
  sh_daughters.clear();
  sh_isPrimary.clear();
  sh_direction.clear();

  v_X.clear(); v_Y.clear(); v_Z.clear();
}


void analysis::analyze(art::Event const& e)
{
  e.getByLabel(fNuGenModuleLabel,eHandleNeutrinos);
  if(fProcessCosmics) e.getByLabel(fCosmicGenModuleLabel,eHandleCosmics);
  e.getByLabel(fLArGeantModuleLabel,eHandleParticles);
  e.getByLabel(fHitsModuleLabel,eHandleHits);
  e.getByLabel(fTrackModuleLabel,eHandleTracks);
  e.getByLabel(fShowerModuleLabel,eHandleShowers);
  e.getByLabel(fPFParticleModuleLabel,eHandlePFPs);
  e.getByLabel(fLArGeantModuleLabel,eHandleSimChannels);
  e.getByLabel(fVertexModuleLabel,eHandleVertices);
  e.getByLabel(fSliceModuleLabel,eHandleSlices);

  clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  ClearData();
  SetupMaps(e);
  SimulationProcessor(e);
  ReconstructionProcessor(e);
  VertexProcessor(e);

  fEventTree->Fill();
}

DEFINE_ART_MODULE(analysis)

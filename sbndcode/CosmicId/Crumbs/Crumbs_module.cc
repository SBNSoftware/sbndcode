////////////////////////////////////////////////////////////////////////
// Class:       Crumbs
// Plugin Type: analyzer
// File:        Crumbs_module.cc
//
// Generated at Wed Jan  5 08:25:29 2022 by Henry Lay using cetskelgen
// from  version .
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

#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "lardataobj/RecoBase/Slice.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "art_root_io/TFileService.h"
#include "TTree.h"
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larcoreobj/SummaryData/POTSummary.h"

class Crumbs;


class Crumbs : public art::EDAnalyzer {
public:
  explicit Crumbs(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  Crumbs(Crumbs const&) = delete;
  Crumbs(Crumbs&&) = delete;
  Crumbs& operator=(Crumbs const&) = delete;
  Crumbs& operator=(Crumbs&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginSubRun(const art::SubRun& sr) override;
  void endJob() override;

  void ClearMaps();
  void ResetVars();
  void SetupMaps(art::Event const& e);
  std::vector<art::Ptr<recob::Hit> > GetAllSliceHits(art::Event const& e, const art::Ptr<recob::PFParticle> pPrimary);
  std::map<int,float> SlicePurity(art::Event const& e, std::vector<art::Ptr<recob::Hit> > sliceHits);
  float SliceCompleteness(art::Event const& e, std::vector<art::Ptr<recob::Hit> > sliceHits, std::vector<art::Ptr<recob::Hit> > allHits, const int matchedGenID);
  int SliceTruthId(std::map<int, float> purities);
  std::vector<art::Ptr<anab::T0> > GetCRTTrackT0s(art::Event const& e, const art::Ptr<recob::PFParticle> pPrimary);
  std::vector<art::Ptr<anab::T0> > GetCRTHitT0s(art::Event const& e, const art::Ptr<recob::PFParticle> pPrimary);


private:

  // Declare member data here.

  bool fVerbose, fProcessUnambiguousSlices, fProcessNeutrinos, fProcessCosmics;

  std::string fMCParticlesModuleLabel, fPFParticlesModuleLabel, fHitsModuleLabel, fTrackModuleLabel, fShowerModuleLabel, fSliceModuleLabel,
    fGeneratorModuleLabel, fCosmicModuleLabel, fFlashMatchModuleLabel, fCRTTrackMatchModuleLabel, fCRTHitMatchModuleLabel;

  std::map<int, int> fTrackToGenMap;
  std::map<int, std::string> fGenTypeMap;

  sbnd::TPCGeoAlg fTpcGeo;

  TTree *fNuSliceTree, *fNotNuSliceTree, *fRunTree;

  double tpc_NuScore, tpc_CRFracHitsInLongestTrack, tpc_CRLongestTrackDeflection, tpc_CRLongestTrackDirY, tpc_CRNHitsMax,
    tpc_NuEigenRatioInSphere, tpc_NuNFinalStatePfos, tpc_NuNHitsTotal, tpc_NuNSpacePointsInSphere, tpc_NuVertexY, tpc_NuWeightedDirZ,
    pds_FMTotalScore, pds_FMYScore, pds_FMZScore, pds_FMRRScore, pds_FMRatioScore, pds_FMPE, crt_TrackScore, crt_HitScore;

  unsigned eventID, subRunID, runID, slicePDG, sliceIndex, matchedIndex;
  std::string matchedType;
  double matchedPurity, matchedCompleteness;

  float fTotalPOT;
};


Crumbs::Crumbs(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fVerbose                    (p.get<bool>("Verbose",false)),
  fProcessUnambiguousSlices   (p.get<bool>("ProcessUnambiguousSlices",false)),
  fProcessNeutrinos           (p.get<bool>("ProcessNeutrinos",true)),
  fProcessCosmics             (p.get<bool>("ProcessCosmics",true)),
  fMCParticlesModuleLabel     (p.get<std::string>("MCParticlesModuleLabel")),
  fPFParticlesModuleLabel     (p.get<std::string>("PFParticlesModuleLabel")),
  fHitsModuleLabel            (p.get<std::string>("HitsModuleLabel")),
  fTrackModuleLabel           (p.get<std::string>("TrackModuleLabel")),
  fShowerModuleLabel          (p.get<std::string>("ShowerModuleLabel")),
  fSliceModuleLabel           (p.get<std::string>("SliceModuleLabel")),
  fGeneratorModuleLabel       (p.get<std::string>("GeneratorModuleLabel")),
  fCosmicModuleLabel          (p.get<std::string>("CosmicModuleLabel")),
  fFlashMatchModuleLabel      (p.get<std::string>("FlashMatchModuleLabel")),
  fCRTTrackMatchModuleLabel   (p.get<std::string>("CRTTrackMatchModuleLabel")),
  fCRTHitMatchModuleLabel     (p.get<std::string>("CRTHitMatchModuleLabel"))

{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  art::ServiceHandle<art::TFileService> tfs;
  fNuSliceTree = tfs->make<TTree>("NuSliceTree","True Nu Slice data TTree");

  fNuSliceTree->Branch("tpc_NuScore",&tpc_NuScore);
  fNuSliceTree->Branch("tpc_CRFracHitsInLongestTrack",&tpc_CRFracHitsInLongestTrack);
  fNuSliceTree->Branch("tpc_CRLongestTrackDeflection",&tpc_CRLongestTrackDeflection);
  fNuSliceTree->Branch("tpc_CRLongestTrackDirY",&tpc_CRLongestTrackDirY);
  fNuSliceTree->Branch("tpc_CRNHitsMax",&tpc_CRNHitsMax);
  fNuSliceTree->Branch("tpc_NuEigenRatioInSphere",&tpc_NuEigenRatioInSphere);
  fNuSliceTree->Branch("tpc_NuNFinalStatePfos",&tpc_NuNFinalStatePfos);
  fNuSliceTree->Branch("tpc_NuNHitsTotal",&tpc_NuNHitsTotal);
  fNuSliceTree->Branch("tpc_NuNSpacePointsInSphere",&tpc_NuNSpacePointsInSphere);
  fNuSliceTree->Branch("tpc_NuVertexY",&tpc_NuVertexY);
  fNuSliceTree->Branch("tpc_NuWeightedDirZ",&tpc_NuWeightedDirZ);

  fNuSliceTree->Branch("pds_FMTotalScore",&pds_FMTotalScore);
  fNuSliceTree->Branch("pds_FMYScore",&pds_FMYScore);
  fNuSliceTree->Branch("pds_FMZScore",&pds_FMZScore);
  fNuSliceTree->Branch("pds_FMRRScore",&pds_FMRRScore);
  fNuSliceTree->Branch("pds_FMRatioScore",&pds_FMRatioScore);
  fNuSliceTree->Branch("pds_FMPE",&pds_FMPE);

  fNuSliceTree->Branch("crt_TrackScore",&crt_TrackScore);
  fNuSliceTree->Branch("crt_HitScore",&crt_HitScore);

  fNuSliceTree->Branch("eventID",&eventID);
  fNuSliceTree->Branch("subRunID",&subRunID);
  fNuSliceTree->Branch("runID",&runID);
  fNuSliceTree->Branch("slicePDG",&slicePDG);
  fNuSliceTree->Branch("sliceIndex",&sliceIndex);
  fNuSliceTree->Branch("matchedIndex",&matchedIndex);
  fNuSliceTree->Branch("matchedType",&matchedType);
  fNuSliceTree->Branch("matchedPurity",&matchedPurity);
  fNuSliceTree->Branch("matchedCompleteness",&matchedCompleteness);

  fNotNuSliceTree = tfs->make<TTree>("NotNuSliceTree","True Non-Nu Slice data TTree");

  fNotNuSliceTree->Branch("tpc_NuScore",&tpc_NuScore);
  fNotNuSliceTree->Branch("tpc_CRFracHitsInLongestTrack",&tpc_CRFracHitsInLongestTrack);
  fNotNuSliceTree->Branch("tpc_CRLongestTrackDeflection",&tpc_CRLongestTrackDeflection);
  fNotNuSliceTree->Branch("tpc_CRLongestTrackDirY",&tpc_CRLongestTrackDirY);
  fNotNuSliceTree->Branch("tpc_CRNHitsMax",&tpc_CRNHitsMax);
  fNotNuSliceTree->Branch("tpc_NuEigenRatioInSphere",&tpc_NuEigenRatioInSphere);
  fNotNuSliceTree->Branch("tpc_NuNFinalStatePfos",&tpc_NuNFinalStatePfos);
  fNotNuSliceTree->Branch("tpc_NuNHitsTotal",&tpc_NuNHitsTotal);
  fNotNuSliceTree->Branch("tpc_NuNSpacePointsInSphere",&tpc_NuNSpacePointsInSphere);
  fNotNuSliceTree->Branch("tpc_NuVertexY",&tpc_NuVertexY);
  fNotNuSliceTree->Branch("tpc_NuWeightedDirZ",&tpc_NuWeightedDirZ);

  fNotNuSliceTree->Branch("pds_FMTotalScore",&pds_FMTotalScore);
  fNotNuSliceTree->Branch("pds_FMYScore",&pds_FMYScore);
  fNotNuSliceTree->Branch("pds_FMZScore",&pds_FMZScore);
  fNotNuSliceTree->Branch("pds_FMRRScore",&pds_FMRRScore);
  fNotNuSliceTree->Branch("pds_FMRatioScore",&pds_FMRatioScore);
  fNotNuSliceTree->Branch("pds_FMPE",&pds_FMPE);

  fNotNuSliceTree->Branch("crt_TrackScore",&crt_TrackScore);
  fNotNuSliceTree->Branch("crt_HitScore",&crt_HitScore);

  fNotNuSliceTree->Branch("eventID",&eventID);
  fNotNuSliceTree->Branch("subRunID",&subRunID);
  fNotNuSliceTree->Branch("runID",&runID);
  fNotNuSliceTree->Branch("slicePDG",&slicePDG);
  fNotNuSliceTree->Branch("sliceIndex",&sliceIndex);
  fNotNuSliceTree->Branch("matchedIndex",&matchedIndex);
  fNotNuSliceTree->Branch("matchedType",&matchedType);
  fNotNuSliceTree->Branch("matchedPurity",&matchedPurity);
  fNotNuSliceTree->Branch("matchedCompleteness",&matchedCompleteness);

  fTotalPOT = 0;

  fRunTree = tfs->make<TTree>("RunTree","Run data TTree");
  fRunTree->Branch("totalPOT",&fTotalPOT);
}

void Crumbs::beginSubRun(const art::SubRun& sr)
{
  art::Handle< sumdata::POTSummary > potListHandle;

  if(sr.getByLabel(fGeneratorModuleLabel,potListHandle))
    fTotalPOT+=potListHandle->totpot;
  else
    fTotalPOT+=0.;
}

void Crumbs::endJob()
{
  fRunTree->Fill();
}

void Crumbs::ClearMaps() 
{
  fTrackToGenMap.clear();
  fGenTypeMap.clear();
}

void Crumbs::ResetVars()
{
  tpc_NuScore = -999999.; tpc_CRFracHitsInLongestTrack = -999999.; tpc_CRLongestTrackDeflection = -999999.; tpc_CRLongestTrackDirY = -999999.; tpc_CRNHitsMax = -999999.;
  tpc_NuEigenRatioInSphere = -999999.; tpc_NuNFinalStatePfos = -999999.; tpc_NuNHitsTotal = -999999.; tpc_NuNSpacePointsInSphere = -999999.; tpc_NuVertexY = -999999.;
  tpc_NuWeightedDirZ = -999999.;

  pds_FMTotalScore = -999999.; pds_FMYScore = -999999.; pds_FMZScore = -999999.; pds_FMRRScore = -999999.; pds_FMRatioScore = -999999.; pds_FMPE = -999999.;

  crt_TrackScore = -4.; crt_HitScore = -4.;

  slicePDG = 999999; sliceIndex = 999999; matchedIndex = 999999;
  matchedType = "";
  matchedPurity = -999999.; matchedCompleteness = -999999.;

}

void Crumbs::SetupMaps(art::Event const& e)
{

  unsigned nNu(0), nCos(0);

  if(fProcessNeutrinos)
    {
      art::Handle<std::vector<simb::MCTruth> > handleMCTruthNu;
      e.getByLabel(fGeneratorModuleLabel, handleMCTruthNu);
      art::FindManyP<simb::MCParticle> truthNuMCPAssn(handleMCTruthNu,e,fMCParticlesModuleLabel);

      for (unsigned int i = 0; i < handleMCTruthNu->size(); ++i){
	const art::Ptr<simb::MCTruth> mcTruth(handleMCTruthNu, i);
	const simb::MCParticle nu = mcTruth->GetNeutrino().Nu();

	if(!fTpcGeo.InVolume(nu))
	  fGenTypeMap[i] = "DirtNu";
	else
	  fGenTypeMap[i] = "Nu";
    
	const std::vector<art::Ptr<simb::MCParticle> > particles = truthNuMCPAssn.at(mcTruth.key());
    
	for (auto const& particle : particles)
	  {
	    fTrackToGenMap[particle->TrackId()] = i;
	  }
	++nNu;
      }
    }

  if(fProcessCosmics)
    {
      art::Handle<std::vector<simb::MCTruth> > handleMCTruthCosmic;
      e.getByLabel(fCosmicModuleLabel, handleMCTruthCosmic);

      art::FindManyP<simb::MCParticle> truthCosmicMCPAssn(handleMCTruthCosmic,e,fMCParticlesModuleLabel);

      for (unsigned int i = 0; i < handleMCTruthCosmic->size(); ++i){
	const art::Ptr<simb::MCTruth> mcTruth(handleMCTruthCosmic, i);

	fGenTypeMap[i + nNu] = "Cosmic";
    
	const std::vector<art::Ptr<simb::MCParticle> > particles = truthCosmicMCPAssn.at(mcTruth.key());
    
	for (auto const& particle : particles)
	  {
	    fTrackToGenMap[particle->TrackId()] = i + nNu;
	  }
	++nCos;
      }
    }
  
  for(auto const& [id, type] : fGenTypeMap){
    if(fVerbose)
      std::cout << "MCTruth " << id << " " << type << std::endl;
  }

  eventID = e.event();
  subRunID = e.subRun();
  runID = e.run();
}

void Crumbs::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  this->ClearMaps();
  this->SetupMaps(e);

  art::Handle<std::vector<recob::PFParticle> > handlePFPs;
  e.getByLabel(fPFParticlesModuleLabel, handlePFPs);

  art::Handle<std::vector<recob::Slice> > handleSlices;
  e.getByLabel(fSliceModuleLabel, handleSlices);

  art::Handle<std::vector<recob::Hit> > handleHits;
  e.getByLabel(fHitsModuleLabel, handleHits);
  std::vector<art::Ptr<recob::Hit> > allHits;
  art::fill_ptr_vector(allHits, handleHits);

  art::Handle<std::vector<sbn::SimpleFlashMatch> > handleFMs;
  e.getByLabel(fFlashMatchModuleLabel, handleFMs);

  std::vector<art::Ptr<recob::PFParticle> > unambigPrimariesVector, nuPrimariesVector;

  for (unsigned int i = 0; i < handlePFPs->size(); ++i){
    const art::Ptr<recob::PFParticle> pParticle(handlePFPs, i);
    if(pParticle->IsPrimary()) {
      int pdg = std::abs(pParticle->PdgCode());
      if(pdg == 13) unambigPrimariesVector.push_back(pParticle);
      else if(pdg == 12 || pdg == 14) nuPrimariesVector.push_back(pParticle);
    }
  }

  art::FindManyP<larpandoraobj::PFParticleMetadata> pfpMetadataAssoc(handlePFPs, e, fPFParticlesModuleLabel);
  art::FindManyP<sbn::SimpleFlashMatch> pfpFMAssoc(handlePFPs, e, fFlashMatchModuleLabel);

  unsigned nSlices(0);
  
  for (auto const& primary : nuPrimariesVector)
    {
      this->ResetVars();

      if(fVerbose)
	std::cout << "\nNeutrino Slice" << std::endl;

      std::vector<art::Ptr<recob::Hit> > sliceHits = this->GetAllSliceHits(e, primary);
      std::map<int, float> puritiesMap = this->SlicePurity(e, sliceHits);
      const int truthId = this->SliceTruthId(puritiesMap);
      matchedCompleteness = this->SliceCompleteness(e, sliceHits, allHits, truthId);

      const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata> > pfpMetaVec = pfpMetadataAssoc.at(primary.key());
      const std::vector<art::Ptr<sbn::SimpleFlashMatch> > pfpFMVec = pfpFMAssoc.at(primary.key());
      const std::vector<art::Ptr<anab::T0> > sliceCRTTrackT0s = this->GetCRTTrackT0s(e, primary);
      const std::vector<art::Ptr<anab::T0> > sliceCRTHitT0s = this->GetCRTHitT0s(e, primary);

      if(fVerbose) 
	{
	  std::cout << " matches to MCTruth of origin " << fGenTypeMap[truthId] << " with purity " << puritiesMap[truthId] << std::endl;
	  std::cout << "Number of CRT Track T0 matches: " << sliceCRTTrackT0s.size() << std::endl;
	  std::cout << "Number of CRT Hit T0 matches: " << sliceCRTHitT0s.size() << std::endl;
	}

      if (pfpMetaVec.size() != 1){
	std::cout << "ERROR: ----- Cannot get PFPMetadata" << std::endl;
      }

      if (pfpFMVec.size() != 1){
	std::cout << "ERROR: ----- Cannot find single flash match" << std::endl;
      }

      if (sliceCRTTrackT0s.size() > 1){
	std::cout << "ERROR: ----- Multiple CRT Track t0s" << std::endl;
      }
      else if(sliceCRTTrackT0s.size() == 1) {
	const art::Ptr<anab::T0> crttrackmatcht0 = sliceCRTTrackT0s.front();
	crt_TrackScore = crttrackmatcht0->TriggerConfidence();
      }


      if (sliceCRTHitT0s.size() > 1){
	std::cout << "ERROR: ----- Multiple CRT Hit t0s" << std::endl;
      }
      else if(sliceCRTHitT0s.size() == 1) {
	const art::Ptr<anab::T0> crthitmatcht0 = sliceCRTHitT0s.front();
	crt_HitScore = crthitmatcht0->TriggerConfidence();
      }

      const art::Ptr<larpandoraobj::PFParticleMetadata> pfpMeta = pfpMetaVec.front();
      auto propertiesMap = pfpMeta->GetPropertiesMap();

      const art::Ptr<sbn::SimpleFlashMatch> flashmatch = pfpFMVec.front();
      
      if(fVerbose)
	{
	  for (auto const& propertiesMapIter : propertiesMap)
	    std::cout << propertiesMapIter.first << ": " << propertiesMapIter.second << std::endl;

	  for (auto const& t0 : pfpFMVec)
	    {
	      std::string hasMatch = t0->present ? "yes" : "no";
	      std::cout << "------- A FM Object -------\n"
			<< "Match? " <<  hasMatch << '\n'
			<< "Time: " << t0->time << '\n'
			    << "PE: " << t0->light.pe << '\n'
			    << "Score: " << t0->score.total << '\n' << std::endl;
	    }
	}

      auto propertiesMapIter = propertiesMap.find("NuScore");
      if (propertiesMapIter == propertiesMap.end()){
	std::cout << "Error finding variable" << std::endl;
	continue;
      }
      tpc_NuScore = propertiesMapIter->second;

      propertiesMapIter = propertiesMap.find("CRFracHitsInLongestTrack");
      if (propertiesMapIter == propertiesMap.end()){
	std::cout << "Error finding variable" << std::endl;
	continue;
      }
      tpc_CRFracHitsInLongestTrack = propertiesMapIter->second;

      propertiesMapIter = propertiesMap.find("CRLongestTrackDeflection");
      if (propertiesMapIter == propertiesMap.end()){
	std::cout << "Error finding variable" << std::endl;
	continue;
      }
      tpc_CRLongestTrackDeflection = propertiesMapIter->second;

      propertiesMapIter = propertiesMap.find("CRLongestTrackDirY");
      if (propertiesMapIter == propertiesMap.end()){
	std::cout << "Error finding variable" << std::endl;
	continue;
      }
      tpc_CRLongestTrackDirY = propertiesMapIter->second;

      propertiesMapIter = propertiesMap.find("CRNHitsMax");
      if (propertiesMapIter == propertiesMap.end()){
	std::cout << "Error finding variable" << std::endl;
	continue;
      }
      tpc_CRNHitsMax = propertiesMapIter->second;

      propertiesMapIter = propertiesMap.find("NuEigenRatioInSphere");
      if (propertiesMapIter == propertiesMap.end()){
	std::cout << "Error finding variable" << std::endl;
	continue;
      }
      tpc_NuEigenRatioInSphere = propertiesMapIter->second;

      propertiesMapIter = propertiesMap.find("NuNFinalStatePfos");
      if (propertiesMapIter == propertiesMap.end()){
	std::cout << "Error finding variable" << std::endl;
	continue;
      }
      tpc_NuNFinalStatePfos = propertiesMapIter->second;

      propertiesMapIter = propertiesMap.find("NuNHitsTotal");
      if (propertiesMapIter == propertiesMap.end()){
	std::cout << "Error finding variable" << std::endl;
	continue;
      }
      tpc_NuNHitsTotal = propertiesMapIter->second;

      propertiesMapIter = propertiesMap.find("NuNSpacePointsInSphere");
      if (propertiesMapIter == propertiesMap.end()){
	std::cout << "Error finding variable" << std::endl;
	continue;
      }
      tpc_NuNSpacePointsInSphere = propertiesMapIter->second;

      propertiesMapIter = propertiesMap.find("NuVertexY");
      if (propertiesMapIter == propertiesMap.end()){
	std::cout << "Error finding variable" << std::endl;
	continue;
      }
      tpc_NuVertexY = propertiesMapIter->second;

      propertiesMapIter = propertiesMap.find("NuWeightedDirZ");
      if (propertiesMapIter == propertiesMap.end()){
	std::cout << "Error finding variable" << std::endl;
	continue;
      }
      tpc_NuWeightedDirZ = propertiesMapIter->second;

      pds_FMTotalScore = flashmatch->score.total;
      pds_FMYScore = flashmatch->score.y;
      pds_FMZScore = flashmatch->score.z;
      pds_FMRRScore = flashmatch->score.rr;
      pds_FMRatioScore = flashmatch->score.ratio;
      pds_FMPE = flashmatch->light.pe;

      slicePDG = primary->PdgCode();
      sliceIndex = nSlices;
      matchedIndex = truthId;
      matchedType = fGenTypeMap[truthId];
      matchedPurity = puritiesMap[truthId];
      
      if(matchedType == "Nu" && matchedPurity > 0.8)
	fNuSliceTree->Fill();
      else
	fNotNuSliceTree->Fill();
      
      ++nSlices;
    }

  if(fProcessUnambiguousSlices)
    {
      for (auto const& primary : unambigPrimariesVector)
	{
	  this->ResetVars();
      
	  if(fVerbose)
	    std::cout << "\nCosmic Slice" << std::endl;

	  std::vector<art::Ptr<recob::Hit> > sliceHits = this->GetAllSliceHits(e, primary);
	  std::map<int, float> puritiesMap = this->SlicePurity(e, sliceHits);
	  int truthId = this->SliceTruthId(puritiesMap);
	  matchedCompleteness = this->SliceCompleteness(e, sliceHits, allHits, truthId);

	  const std::vector<art::Ptr<anab::T0> > sliceCRTrackT0s = this->GetCRTTrackT0s(e, primary);
	  const std::vector<art::Ptr<anab::T0> > sliceCRHitT0s = this->GetCRTHitT0s(e, primary);

	  if(fVerbose)
	    {
	      std::cout << " matches to MCTruth of origin " << fGenTypeMap[truthId] << " with purity " << puritiesMap[truthId] << std::endl;
	      std::cout << "Number of CRT Track T0 matches: " << sliceCRTrackT0s.size() << std::endl;
	      std::cout << "Number of CRT Hit T0 matches: " << sliceCRHitT0s.size() << std::endl;
	    }

	  const std::vector<art::Ptr<sbn::SimpleFlashMatch> > pfpFMVec = pfpFMAssoc.at(primary.key());

	  if (pfpFMVec.size() != 1){
	    std::cout << "Cannot find single flash match" << std::endl;
	  }

	  if(fVerbose)
	    {
 	      for (auto const& t0 : pfpFMVec)
		{
		  std::string hasMatch = t0->present ? "yes" : "no";
		  std::cout << "------- A FM Object -------\n"
			    << "Match? " <<  hasMatch << '\n'
			    << "Time: " << t0->time << '\n'
			    << "PE: " << t0->light.pe << '\n'
			    << "Score: " << t0->score.total << '\n' << std::endl;
		}
	    }

	  slicePDG = primary->PdgCode();
	  sliceIndex = nSlices;
	  matchedIndex = truthId;
	  matchedType = fGenTypeMap[truthId];
	  matchedPurity = puritiesMap[truthId];

	  if(matchedType == "Nu" && matchedPurity > 0.8)
	    fNuSliceTree->Fill();
	  else
	    fNotNuSliceTree->Fill();

	  ++nSlices;
	}
    }

}

std::vector<art::Ptr<recob::Hit> > Crumbs::GetAllSliceHits(art::Event const& e, const art::Ptr<recob::PFParticle> pPrimary)
{
  art::Handle<std::vector<recob::PFParticle> > handlePFPs;
  e.getByLabel(fPFParticlesModuleLabel, handlePFPs);

  art::Handle<std::vector<recob::Slice> > handleSlices;
  e.getByLabel(fSliceModuleLabel, handleSlices);

  art::FindManyP<recob::Slice> pfpSliceAssn(handlePFPs,e,fSliceModuleLabel);
  art::FindManyP<recob::Hit> sliceHitAssn(handleSlices,e,fSliceModuleLabel);

  const std::vector<art::Ptr<recob::Slice> > slices = pfpSliceAssn.at(pPrimary.key());
  
  if(slices.size() != 1) std::cout << "Not 1 slice (" << slices.size() << ")" << std::endl;
  
  const std::vector<art::Ptr<recob::Hit> > sliceHits = sliceHitAssn.at(slices.at(0).key());

  return sliceHits;
}

std::map<int, float> Crumbs::SlicePurity(art::Event const& e, std::vector<art::Ptr<recob::Hit> > sliceHits)
{
  std::map<int, int> sliceHitMap;
  std::map<int, float> slicePurityMap;

  auto clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  for (auto const& hit : sliceHits)
    {
      ++sliceHitMap[fTrackToGenMap[TruthMatchUtils::TrueParticleID(clockData,hit,true)]];
    }

  for (auto const& [id, nHits] : sliceHitMap)
    {
      slicePurityMap[id] = (float) nHits / (float) sliceHits.size();
    }

  return slicePurityMap;
}

float Crumbs::SliceCompleteness(art::Event const& e, std::vector<art::Ptr<recob::Hit> > sliceHits, std::vector<art::Ptr<recob::Hit> > allHits, const int matchedGenID)
{
  int nSliceHits(0), nHits(0);
  auto clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  for (auto const& hit : sliceHits)
    {
      if(fTrackToGenMap[TruthMatchUtils::TrueParticleID(clockData,hit,true)] == matchedGenID)
	++nSliceHits;
    }

  for (auto const& hit : allHits)
    {
      if(fTrackToGenMap[TruthMatchUtils::TrueParticleID(clockData,hit,true)] == matchedGenID)
	++nHits;
    }
  
  if(nHits == 0) 
    return 0;

  return (float) nSliceHits / (float) nHits;
}

int Crumbs::SliceTruthId(std::map<int, float> purities)
{
  float maxPur = -1;
  int retId = -999999;

  for (auto const& [id, purity] : purities)
    {
      if(purity > maxPur) 
	{
	  retId = id;
	  maxPur = purity;
	}
    }

  return retId;
}

std::vector<art::Ptr<anab::T0> > Crumbs::GetCRTTrackT0s(art::Event const& e, const art::Ptr<recob::PFParticle> pPrimary)
{
  std::vector<art::Ptr<anab::T0> > t0Vec;

  art::Handle<std::vector<recob::PFParticle> > handlePFPs;
  e.getByLabel(fPFParticlesModuleLabel, handlePFPs);

  art::Handle<std::vector<recob::Slice> > handleSlices;
  e.getByLabel(fSliceModuleLabel, handleSlices);

  art::Handle<std::vector<recob::Track> > handleTracks;
  e.getByLabel(fTrackModuleLabel, handleTracks);

  art::FindOneP<recob::Slice> pfpSliceAssn(handlePFPs,e,fSliceModuleLabel);
  art::FindManyP<recob::PFParticle> slicePFPAssn(handleSlices,e,fSliceModuleLabel);
  art::FindManyP<recob::Track> pfpTrackAssn(handlePFPs,e,fTrackModuleLabel);
  art::FindManyP<anab::T0> trackT0Assn(handleTracks,e,fCRTTrackMatchModuleLabel);

  const art::Ptr<recob::Slice> slice = pfpSliceAssn.at(pPrimary.key());
  //  if(slices.size() != 1) std::cout << "Not 1 slice (" << slices.size() << ")" << std::endl;
  //  const art::Ptr<recob::Slice> slice = slices.front();

  const std::vector<art::Ptr<recob::PFParticle> > pfps = slicePFPAssn.at(slice.key());
  
  for(auto const& pfp : pfps)
    {
      const std::vector<art::Ptr<recob::Track> > tracks = pfpTrackAssn.at(pfp.key());
      if(tracks.size() == 0) continue;
      else if(tracks.size() > 1){
	std::cout << "Multiple tracks associated with PFP" << std::endl;
	continue;
      }
      
      const art::Ptr<recob::Track> track = tracks.front();
      const std::vector<art::Ptr<anab::T0> > t0s = trackT0Assn.at(track.key());
      t0Vec.insert(t0Vec.end(), t0s.begin(), t0s.end());
    }
  
  return t0Vec;
}

std::vector<art::Ptr<anab::T0> > Crumbs::GetCRTHitT0s(art::Event const& e, const art::Ptr<recob::PFParticle> pPrimary)
{
  std::vector<art::Ptr<anab::T0> > t0Vec;

  art::Handle<std::vector<recob::PFParticle> > handlePFPs;
  e.getByLabel(fPFParticlesModuleLabel, handlePFPs);

  art::Handle<std::vector<recob::Slice> > handleSlices;
  e.getByLabel(fSliceModuleLabel, handleSlices);

  art::Handle<std::vector<recob::Track> > handleTracks;
  e.getByLabel(fTrackModuleLabel, handleTracks);

  art::FindOneP<recob::Slice> pfpSliceAssn(handlePFPs,e,fSliceModuleLabel);
  art::FindManyP<recob::PFParticle> slicePFPAssn(handleSlices,e,fSliceModuleLabel);
  art::FindManyP<recob::Track> pfpTrackAssn(handlePFPs,e,fTrackModuleLabel);
  art::FindManyP<anab::T0> trackT0Assn(handleTracks,e,fCRTHitMatchModuleLabel);

  const art::Ptr<recob::Slice> slice = pfpSliceAssn.at(pPrimary.key());
  //  if(slices.size() != 1) std::cout << "Not 1 slice (" << slices.size() << ")" << std::endl;
  //  const art::Ptr<recob::Slice> slice = slices.front();

  const std::vector<art::Ptr<recob::PFParticle> > pfps = slicePFPAssn.at(slice.key());
  
  for(auto const& pfp : pfps)
    {
      const std::vector<art::Ptr<recob::Track> > tracks = pfpTrackAssn.at(pfp.key());
      if(tracks.size() == 0) continue;
      else if(tracks.size() > 1){
	std::cout << "Multiple tracks associated with PFP" << std::endl;
	continue;
      }
      
      const art::Ptr<recob::Track> track = tracks.front();
      const std::vector<art::Ptr<anab::T0> > t0s = trackT0Assn.at(track.key());
      t0Vec.insert(t0Vec.end(), t0s.begin(), t0s.end());
    }
  
  return t0Vec;
}

DEFINE_ART_MODULE(Crumbs)

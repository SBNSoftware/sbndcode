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
#include "nusimdata/SimulationBase/MCTruth.h"
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "art_root_io/TFileService.h"
#include "TTree.h"
#include "lardataobj/AnalysisBase/T0.h"

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
  void beginJob() override;
  void endJob() override;

  void ClearMaps();
  void ResetVars();
  void SetupMaps(art::Event const& e);
  std::vector<art::Ptr<recob::Hit> > GetAllSliceHits(art::Event const& e, const art::Ptr<recob::PFParticle> pPrimary);
  std::map<int,float> SlicePurity(art::Event const& e, std::vector<art::Ptr<recob::Hit> > sliceHits);
  int SliceTruthId(std::map<int, float> purities);


private:

  // Declare member data here.

  bool fVerbose, fProcessUnambiguousSlices;

  std::string fMCParticlesModuleLabel, fPFParticlesModuleLabel, fHitsModuleLabel, fTrackModuleLabel, fShowerModuleLabel, fSliceModuleLabel,
    fGeneratorModuleLabel, fCosmicModuleLabel, fFlashMatchModuleLabel;

  std::map<int, int> fTrackToGenMap;
  std::map<int, std::string> fGenTypeMap;

  sbnd::TPCGeoAlg fTpcGeo;

  TTree *fSliceTree;

  double tpc_NuScore, tpc_CRFracHitsInLongestTrack, tpc_CRLongestTrackDeflection, tpc_CRLongestTrackDirY, tpc_CRNHitsMax,
    tpc_NuEigenRatioInSphere, tpc_NuNFinalStatePfos, tpc_NuNHitsTotal, tpc_NuNSpacePointsInSphere, tpc_NuVertexY, tpc_NuWeightedDirZ;

  unsigned eventID, subRunID, runID, slicePDG, sliceIndex, matchedIndex;
  std::string matchedType;
  double matchedPurity;
};


Crumbs::Crumbs(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fVerbose                    (p.get<bool>("Verbose",false)),
  fProcessUnambiguousSlices   (p.get<bool>("ProcessUnambiguousSlices",false)),
  fMCParticlesModuleLabel     (p.get<std::string>("MCParticlesModuleLabel")),
  fPFParticlesModuleLabel     (p.get<std::string>("PFParticlesModuleLabel")),
  fHitsModuleLabel            (p.get<std::string>("HitsModuleLabel")),
  fTrackModuleLabel           (p.get<std::string>("TrackModuleLabel")),
  fShowerModuleLabel          (p.get<std::string>("ShowerModuleLabel")),
  fSliceModuleLabel           (p.get<std::string>("SliceModuleLabel")),
  fGeneratorModuleLabel       (p.get<std::string>("GeneratorModuleLabel")),
  fCosmicModuleLabel          (p.get<std::string>("CosmicModuleLabel")),
  fFlashMatchModuleLabel      (p.get<std::string>("FlashMatchModuleLabel"))

{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  art::ServiceHandle<art::TFileService> tfs;
  fSliceTree = tfs->make<TTree>("SliceTree","Slice data TTree");

  fSliceTree->Branch("tpc_NuScore",&tpc_NuScore);
  fSliceTree->Branch("tpc_CRFracHitsInLongestTrack",&tpc_CRFracHitsInLongestTrack);
  fSliceTree->Branch("tpc_CRLongestTrackDeflection",&tpc_CRLongestTrackDeflection);
  fSliceTree->Branch("tpc_CRLongestTrackDirY",&tpc_CRLongestTrackDirY);
  fSliceTree->Branch("tpc_CRNHitsMax",&tpc_CRNHitsMax);
  fSliceTree->Branch("tpc_NuEigenRatioInSphere",&tpc_NuEigenRatioInSphere);
  fSliceTree->Branch("tpc_NuNFinalStatePfos",&tpc_NuNFinalStatePfos);
  fSliceTree->Branch("tpc_NuNHitsTotal",&tpc_NuNHitsTotal);
  fSliceTree->Branch("tpc_NuNSpacePointsInSphere",&tpc_NuNSpacePointsInSphere);
  fSliceTree->Branch("tpc_NuVertexY",&tpc_NuVertexY);
  fSliceTree->Branch("tpc_NuWeightedDirZ",tpc_NuWeightedDirZ);

  fSliceTree->Branch("eventID",&eventID);
  fSliceTree->Branch("subRunID",&subRunID);
  fSliceTree->Branch("runID",&runID);
  fSliceTree->Branch("slicePDG",&slicePDG);
  fSliceTree->Branch("sliceIndex",&sliceIndex);
  fSliceTree->Branch("matchedIndex",&matchedIndex);
  fSliceTree->Branch("matchedType",&matchedType);
  fSliceTree->Branch("matchedPurity",&matchedPurity);
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

  slicePDG = 999999; sliceIndex = 999999; matchedIndex = 999999;
  matchedType = "";
  matchedPurity = -999999.;
}

void Crumbs::SetupMaps(art::Event const& e)
{
  art::Handle<std::vector<simb::MCTruth> > handleMCTruthNu;
  e.getByLabel(fGeneratorModuleLabel, handleMCTruthNu);

  art::Handle<std::vector<simb::MCTruth> > handleMCTruthCosmic;
  e.getByLabel(fCosmicModuleLabel, handleMCTruthCosmic);

  art::FindManyP<simb::MCParticle> truthNuMCPAssn(handleMCTruthNu,e,fMCParticlesModuleLabel);
  art::FindManyP<simb::MCParticle> truthCosmicMCPAssn(handleMCTruthCosmic,e,fMCParticlesModuleLabel);

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
  }

  for (unsigned int i = 0; i < handleMCTruthCosmic->size(); ++i){
    const art::Ptr<simb::MCTruth> mcTruth(handleMCTruthCosmic, i);

    fGenTypeMap[i + handleMCTruthNu->size()] = "Cosmic";
    
    const std::vector<art::Ptr<simb::MCParticle> > particles = truthCosmicMCPAssn.at(mcTruth.key());
    
    for (auto const& particle : particles)
      {
	fTrackToGenMap[particle->TrackId()] = i + handleMCTruthNu->size();
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

  art::Handle<std::vector<anab::T0> > handleT0s;
  e.getByLabel(fFlashMatchModuleLabel, handleT0s);

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
  art::FindManyP<recob::Slice> pfpSliceAssoc(handlePFPs, e, fSliceModuleLabel);
  art::FindManyP<anab::T0> sliceT0Assoc(handleSlices, e, fFlashMatchModuleLabel);

  unsigned nSlices(0);
  
  for (auto const& primary : nuPrimariesVector)
    {
      this->ResetVars();

      if(fVerbose)
	std::cout << "\nNeutrino Slice" << std::endl;

      std::vector<art::Ptr<recob::Hit> > sliceHits = this->GetAllSliceHits(e, primary);
      std::map<int, float> puritiesMap = this->SlicePurity(e, sliceHits);
      int truthId = this->SliceTruthId(puritiesMap);

      if(fVerbose)
	std::cout << " matches to MCTruth of origin " << fGenTypeMap[truthId] << " with purity " << puritiesMap[truthId] << std::endl;

      const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata> > pfpMetaVec = pfpMetadataAssoc.at(primary.key());
      const std::vector<art::Ptr<recob::Slice> > pfpSliceVec = pfpSliceAssoc.at(primary.key());

      if (pfpMetaVec.size() != 1){
	std::cout << "Cannot get PFPMetadata" << std::endl;
      }

      if (pfpSliceVec.size() != 1){
	std::cout << "Cannot find single slice" << std::endl;
      }

      const art::Ptr<larpandoraobj::PFParticleMetadata> pfpMeta = pfpMetaVec.front();
      auto propertiesMap = pfpMeta->GetPropertiesMap();

      const art::Ptr<recob::Slice> slice = pfpSliceVec.front();
      const std::vector<art::Ptr<anab::T0> > sliceT0Vec = sliceT0Assoc.at(slice.key());

      if(fVerbose)
	{
	  for (auto const& propertiesMapIter : propertiesMap)
	    std::cout << propertiesMapIter.first << ": " << propertiesMapIter.second << std::endl;

	  std::cout << "Size of T0 Vector: " << sliceT0Vec.size() << std::endl;

	  for (auto const& t0 : sliceT0Vec)
	    {
	      std::cout << "------- A T0 Object -------\n"
			<< "Time: " << t0->Time() << '\n'
			<< "TriggerType: " << t0->TriggerType() << '\n'
			<< "TriggerBits: " << t0->TriggerBits() << '\n'
			<< "ID: " << t0->ID() << '\n'
			<< "TriggerConfidence: " << t0->TriggerConfidence() << '\n' << std::endl;
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

      slicePDG = primary->PdgCode();
      sliceIndex = nSlices;
      matchedIndex = truthId;
      matchedType = fGenTypeMap[truthId];
      matchedPurity = puritiesMap[truthId];

      fSliceTree->Fill();
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

	  if(fVerbose)
	    std::cout << " matches to MCTruth of origin " << fGenTypeMap[truthId] << " with purity " << puritiesMap[truthId] << std::endl;

	  const std::vector<art::Ptr<recob::Slice> > pfpSliceVec = pfpSliceAssoc.at(primary.key());

	  if (pfpSliceVec.size() != 1){
	    std::cout << "Cannot find single slice" << std::endl;
	  }

	  const art::Ptr<recob::Slice> slice = pfpSliceVec.front();
	  const std::vector<art::Ptr<anab::T0> > sliceT0Vec = sliceT0Assoc.at(slice.key());

	  if(fVerbose)
	    {
	      std::cout << "Size of T0 Vector: " << sliceT0Vec.size() << std::endl;
	      
	      for (auto const& t0 : sliceT0Vec)
		{
		  std::cout << "------- A T0 Object -------\n"
			    << "Time: " << t0->Time() << '\n'
			    << "TriggerType: " << t0->TriggerType() << '\n'
			    << "TriggerBits: " << t0->TriggerBits() << '\n'
			    << "ID: " << t0->ID() << '\n'
			    << "TriggerConfidence: " << t0->TriggerConfidence() << '\n' << std::endl;
		}
	    }

	  slicePDG = primary->PdgCode();
	  sliceIndex = nSlices;
	  matchedIndex = truthId;
	  matchedType = fGenTypeMap[truthId];
	  matchedPurity = puritiesMap[truthId];

	  fSliceTree->Fill();
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

void Crumbs::beginJob()
{
  // Implementation of optional member function here.
}

void Crumbs::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(Crumbs)

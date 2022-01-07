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
  void SetupMaps(art::Event const& e, art::Handle<std::vector<simb::MCParticle> > handleMCParticles, art::Handle<std::vector<recob::Hit> > handleHits);
  std::vector<art::Ptr<recob::Hit> > GetAllSliceHits(art::Event const& e, const art::Ptr<recob::PFParticle> pPrimary);
  std::map<int,float> SlicePurity(art::Event const& e, std::vector<art::Ptr<recob::Hit> > sliceHits);
  int SliceTruthId(std::map<int, float> purities);


private:

  // Declare member data here.

  std::string fMCParticlesModuleLabel, fPFParticlesModuleLabel, fHitsModuleLabel, fTrackModuleLabel, fShowerModuleLabel, fSliceModuleLabel,
    fGeneratorModuleLabel, fCosmicModuleLabel;

  std::map<int, art::Ptr<simb::MCParticle> > trueParticleIdMap;
  std::map<int, int> truePrimariesMap, nHitsMap, fTrackToGenMap;
  std::map<int, std::string> fGenTypeMap;

  sbnd::TPCGeoAlg fTpcGeo;
};


Crumbs::Crumbs(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fMCParticlesModuleLabel   (p.get<std::string>("MCParticlesModuleLabel")),
  fPFParticlesModuleLabel   (p.get<std::string>("PFParticlesModuleLabel")),
  fHitsModuleLabel          (p.get<std::string>("HitsModuleLabel")),
  fTrackModuleLabel         (p.get<std::string>("TrackModuleLabel")),
  fShowerModuleLabel        (p.get<std::string>("ShowerModuleLabel")),
  fSliceModuleLabel         (p.get<std::string>("SliceModuleLabel")),
  fGeneratorModuleLabel     (p.get<std::string>("GeneratorModuleLabel")),
  fCosmicModuleLabel        (p.get<std::string>("CosmicModuleLabel"))

{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void Crumbs::ClearMaps() 
{
  trueParticleIdMap.clear();
  truePrimariesMap.clear();
  nHitsMap.clear();
  fTrackToGenMap.clear();
  fGenTypeMap.clear();
}

void Crumbs::SetupMaps(art::Event const& e, art::Handle<std::vector<simb::MCParticle> > handleMCParticles, art::Handle<std::vector<recob::Hit> > handleHits)
{
  for (unsigned int i = 0; i < handleMCParticles->size(); ++i){
    const art::Ptr<simb::MCParticle> pParticle(handleMCParticles, i);

    trueParticleIdMap[pParticle->TrackId()] = pParticle;
  }

  for (unsigned int i = 0; i < handleMCParticles->size(); ++i){
    const art::Ptr<simb::MCParticle> pParticle(handleMCParticles, i);

    int id = pParticle->TrackId();
    int motherId = pParticle->Mother();
    while(trueParticleIdMap.count(motherId) != 0) {
      id = motherId;
      motherId = trueParticleIdMap.at(id)->Mother();
    }

    truePrimariesMap[pParticle->TrackId()] = id;
  }

  auto clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  for (unsigned int i = 0; i < handleHits->size(); ++i){
    const art::Ptr<recob::Hit> pHit(handleHits, i);

    int primaryId = truePrimariesMap[TruthMatchUtils::TrueParticleID(clockData,pHit,true)];
    ++nHitsMap[primaryId];
  }

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
  
  for(auto const& [id, type] : fGenTypeMap)
    std::cout << "MCTruth " << id << " " << type << std::endl;
}

void Crumbs::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  this->ClearMaps();

  art::Handle<std::vector<simb::MCParticle> > handleMCParticles;
  e.getByLabel(fMCParticlesModuleLabel, handleMCParticles);

  art::Handle<std::vector<recob::Hit> > handleHits;
  e.getByLabel(fHitsModuleLabel, handleHits);

  this->SetupMaps(e, handleMCParticles, handleHits);

  art::Handle<std::vector<recob::PFParticle> > handlePFPs;
  e.getByLabel(fPFParticlesModuleLabel, handlePFPs);

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
  
  for (auto const& primary : nuPrimariesVector)
    {
      std::cout << "\nNeutrino Slice" << std::endl;
      std::vector<art::Ptr<recob::Hit> > sliceHits = this->GetAllSliceHits(e, primary);
      std::map<int, float> puritiesMap = this->SlicePurity(e, sliceHits);
      int truthId = this->SliceTruthId(puritiesMap);
      std::cout << " matches to MCTruth of origin " << fGenTypeMap[truthId] << " with purity " << puritiesMap[truthId] << std::endl;

      const std::vector<art::Ptr<larpandoraobj::PFParticleMetadata> > pfpMetaVec = pfpMetadataAssoc.at(primary->Self());

      if (pfpMetaVec.size() !=1){
	std::cout<<"Cannot get PFPMetadata"<<std::endl;
      }

      art::Ptr<larpandoraobj::PFParticleMetadata> pfpMeta = pfpMetaVec.front();
      auto propertiesMap = pfpMeta->GetPropertiesMap();

      for (auto const& propertiesMapIter : propertiesMap)
	{
	  std::cout << propertiesMapIter.first << ": " << propertiesMapIter.second << std::endl;
	}

    }

  for (auto const& primary : unambigPrimariesVector)
    {
      std::cout << "\nCosmic Slice" << std::endl;
      std::vector<art::Ptr<recob::Hit> > sliceHits = this->GetAllSliceHits(e, primary);
      std::map<int, float> puritiesMap = this->SlicePurity(e, sliceHits);
      int truthId = this->SliceTruthId(puritiesMap);
      std::cout << " matches to MCTruth of origin " << fGenTypeMap[truthId] << " with purity " << puritiesMap[truthId] << std::endl;
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

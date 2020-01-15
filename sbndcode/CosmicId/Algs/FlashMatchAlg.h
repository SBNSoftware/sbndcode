#ifndef FLASHMATCHALG_H_SEEN
#define FLASHMATCHALG_H_SEEN


///////////////////////////////////////////////
// FlashMatchAlg.h
//
// Functions for simple flash matching
///////////////////////////////////////////////

// sbndcode
#include "sbndcode/Geometry/GeometryWrappers/TPCGeoAlg.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.h"

// framework
#include "fhiclcpp/ParameterSet.h" 
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

// ROOT
#include "TFile.h"
#include "TH1.h"

// c++
#include <vector>
#include <utility>


namespace sbnd{

  class FlashMatchAlg {
  public:

    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<art::InputTag> SpacePointLabel { Name("SpacePointLabel") };
      fhicl::Atom<art::InputTag> PandoraLabel { Name("PandoraLabel") };
      fhicl::Atom<std::string> InputFile { Name("InputFile") };

      fhicl::Atom<double> BeamFlashMin { Name("BeamFlashMin") };
      fhicl::Atom<double> BeamFlashMax { Name("BeamFlashMax") };
      fhicl::Atom<double> LightWindowMin { Name("LightWindowMin") };
      fhicl::Atom<double> LightWindowMax { Name("LightWindowMax") };
      fhicl::Atom<double> EventMin { Name("EventMin") };
      fhicl::Atom<double> EventMax { Name("EventMax") };
      fhicl::Atom<double> TimeResolution { Name("TimeResolution") };
      fhicl::Atom<double> FlashThreshold { Name("FlashThreshold") };
      fhicl::Atom<double> FlashTimeCorrection { Name("FlashTimeCorrection") };
      fhicl::Atom<double> FlashMatchCut { Name("FlashMatchCut") };

    };

    FlashMatchAlg(const Config& config);

    FlashMatchAlg(const fhicl::ParameterSet& pset) :
      FlashMatchAlg(fhicl::Table<Config>(pset, {})()) {}

    FlashMatchAlg();

    ~FlashMatchAlg();

    void reconfigure(const Config& config);

    // Find the biggest optical flash in the beam window
    std::pair<double, double> BiggestBeamFlash(std::vector<recob::OpHit> ophits);

    // Find all optical flash times in the event
    std::vector<double> OpFlashes(std::vector<recob::OpHit> ophits);

    // Find optical flashes per TPC in event
    std::pair<std::vector<double>, std::vector<double>> OpFlashes(art::ValidHandle<std::vector<recob::OpHit>> pdsHandle);

    // Determine if there is a PDS flash per TPC in time with the neutrino beam
    std::pair<bool, bool> BeamFlash(art::ValidHandle<std::vector<recob::OpHit>> pdsHandle);

    // Calculate the flash matching variables for optical hits in TPC and time window
    std::vector<double> OpVariables(std::vector<recob::OpHit> ophits, int tpc, double start_t, double end_t);

    // Calculate the flash matching score from charge weighted mean and optical hit variables
    double DyScore(double x, double y, std::vector<double> variables);
    double DzScore(double x, double z, std::vector<double> variables);
    double RrScore(double x, std::vector<double> variables);
    double PeScore(double x, std::vector<double> variables);
    double FlashScore(double x, double y, double z, std::vector<double> variables, double w1=1, double w2=1, double w3=1, double w4=0);

    bool FlashMatch(recob::PFParticle pfparticle, std::map< size_t, art::Ptr<recob::PFParticle> > pfParticleMap, const art::Event& event, art::ValidHandle<std::vector<recob::OpHit>> pdsHandle);

  private:

    art::InputTag fSpacePointLabel;
    art::InputTag fPandoraLabel;
    std::string fInputFile;

    double fBeamFlashMin;
    double fBeamFlashMax; 
    double fLightWindowMin; 
    double fLightWindowMax;
    double fEventMin;
    double fEventMax;
    double fTimeResolution;
    double fFlashThreshold;
    double fFlashTimeCorrection;
    double fFlashMatchCut;

    std::vector<float> dysp, dzsp, rrsp, pesp, dymean, dzmean, rrmean, pemean;
    int rr_nbins, dy_nbins, dz_nbins, pe_nbins;
    
    geo::GeometryCore const* fGeometryService;
    TPCGeoAlg fTpcGeo;
    
    opdet::sbndPDMapAlg fChannelMap; //map for photon detector types

  };

}

#endif

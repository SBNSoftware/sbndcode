////////////////////////////////////////////////////////////////////////
// Class:       AnalyzeEventsMin
// Plugin Type: analyzer (Unknown Unknown)
// File:        AnalyzeEventsMin_module.cc
//
// Generated at Thu Jan 26 04:44:06 2023 by Luis Pelegrina gutierrez using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////
// Additional Framework includes
#include "art_root_io/TFileService.h"

// Additional LArSoft includes
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/Geometry/WireReadout.h"

#include "lardataobj/MCBase/MCTrack.h"

// Additional LArSoft includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"

#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"


#include "sbnobj/Common/Reco/CRUMBSResult.h"
#include "sbnobj/Common/Reco/MVAPID.h"
#include "sbnobj/Common/Reco/OpT0FinderResult.h"

#include "sbnobj/Common/Reco/MVAPID.h"
#include "sbnobj/Common/Reco/RangeP.h"
#include "sbnobj/Common/Reco/ScatterClosestApproach.h"
#include "sbnobj/Common/Reco/StoppingChi2Fit.h"

#include "sbnobj/Common/Reco/TPCPMTBarycenterMatch.h"
#include "sbnobj/Common/Reco/Stub.h"


// ROOT includes                                                                                                                                                                        
#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TInterpreter.h"
#include "TTimeStamp.h"
#include <larcorealg/Geometry/Exceptions.h>

#include <vector>
#include <limits>
#include <map>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>

//Truth Matching
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"

constexpr int def_int = std::numeric_limits<int>::min();

namespace test {
  class AnalyzeEventsMin;
}


class test::AnalyzeEventsMin : public art::EDAnalyzer {
public:
  explicit AnalyzeEventsMin(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  AnalyzeEventsMin(AnalyzeEventsMin const &) = delete;

  AnalyzeEventsMin(AnalyzeEventsMin &&) = delete;

  AnalyzeEventsMin &operator=(AnalyzeEventsMin const &) = delete;

  AnalyzeEventsMin &operator=(AnalyzeEventsMin &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

  //Declaration of fhicl parameters
  const std::string slice_label;
  const std::string pfp_label;
  const std::string track_label;
  const std::string slice_to_pfp_label;

  const std::string pfp_to_vertex_label;
  const std::string pfp_to_track_label;
  const std::string track_to_chi2_label;
  const std::string pfp_to_pfpmetadata_label;
};

//Initialize some variables usingn the fhicl file "analysisConfigTruth"
test::AnalyzeEventsMin::AnalyzeEventsMin(fhicl::ParameterSet const &p)
  : EDAnalyzer{p},
  // More initializers here.
    slice_label(p.get<std::string>("slice_label")),
    pfp_label(p.get<std::string>("pfp_label")),
    track_label(p.get<std::string>("track_label")),

    slice_to_pfp_label(p.get<std::string>("slice_to_pfp_label")),
    pfp_to_vertex_label(p.get<std::string>("pfp_to_vertex_label")),
    pfp_to_track_label(p.get<std::string>("pfp_to_track_label")),
    track_to_chi2_label(p.get<std::string>("track_to_chi2_label")),
    pfp_to_pfpmetadata_label(p.get<std::string>("pfp_to_pfpmetadata_label"))
 {
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  }

//Bulk of the analizer
void test::AnalyzeEventsMin::analyze(art::Event const &e) {

  //Fill the event ID related information
  std::cout << "EvtID: " << e.id().event() << " Run: " << e.id().run() << " Subrun: " << e.id().subRun() << std::endl;

  //Access to Slice
  art::Handle <std::vector<recob::Slice>> slice_handle;
  std::vector <art::Ptr<recob::Slice>> slice_vec;
  if (e.getByLabel(slice_label, slice_handle))
    art::fill_ptr_vector(slice_vec, slice_handle);

  //Access to tracks
  art::Handle< std::vector<recob::Track> > track_handle;
  std::vector< art::Ptr<recob::Track> > track_vec;
  if(e.getByLabel(track_label, track_handle))
    art::fill_ptr_vector(track_vec, track_handle);

  //Access to pfpParticles
  art::Handle< std::vector<recob::PFParticle> > pfp_handle;
  std::vector< art::Ptr<recob::PFParticle> > pfp_vec;
  if(e.getByLabel(pfp_label, pfp_handle))
    art::fill_ptr_vector(pfp_vec, pfp_handle);

  //Get all the associations
  art::FindManyP <recob::PFParticle> slice_pfp_assns(slice_handle, e,
                                                     slice_to_pfp_label); //Slice associations to pfp particles
  art::FindManyP <larpandoraobj::PFParticleMetadata> pfp_pfpmeta_assns(pfp_handle, e,
                                                                       pfp_to_pfpmetadata_label); //pfp asociation to MetaData
  art::FindManyP <recob::Vertex> pfp_vertex_assns(pfp_handle, e, pfp_to_vertex_label); //pfparticle asociation to vertex

  for (const art::Ptr <recob::Slice> &slice: slice_vec) {
    std::cout << std::endl << "SliceID: " << slice->ID() << std::endl;

    //Get pfparticles and hits
    std::vector <art::Ptr<recob::PFParticle>> slice_pfp_vec = slice_pfp_assns.at(slice.key());


    //Fill the nu_score of the slice
    double nu_score = -1.;
    for (const art::Ptr <recob::PFParticle> &pfp: slice_pfp_vec) {
      if (pfp->IsPrimary() && (std::abs(pfp->PdgCode()) == 12 || std::abs(pfp->PdgCode()) == 14)) {
        const std::vector <art::Ptr<larpandoraobj::PFParticleMetadata>> pfp_meta_vec = pfp_pfpmeta_assns.at(pfp.key());

        art::Ptr <larpandoraobj::PFParticleMetadata> nu_meta = pfp_meta_vec[0];
        std::map<std::string, float> nu_prop_map = nu_meta->GetPropertiesMap();

        if (nu_prop_map.find("NuScore") != nu_prop_map.end()) {
          nu_score = nu_prop_map["NuScore"];
        }
      }
    }
    std::cout << "NuScore: " << nu_score << std::endl;

    //Get vtx
    for (const art::Ptr <recob::PFParticle> &pfp: slice_pfp_vec) {
      //Only continue if it is a neutrino and a primary particle
      if (!pfp->IsPrimary()) continue;
      std::vector <art::Ptr<recob::Vertex>> pfp_vertex_vec = pfp_vertex_assns.at(pfp.key());
      for (const art::Ptr <recob::Vertex> &vertex: pfp_vertex_vec) {
        std::cout << "vtx: " << vertex->position().X() << " " << vertex->position().Y() << " " << vertex->position().Z()
                  << std::endl;
      }
    }


    //Loop trhough all the particles
    std::cout << "pfps:" << std::endl;
    for (const art::Ptr <recob::PFParticle> &pfp: slice_pfp_vec) {
      //Load associations

      std::cout << "pfpID: " << pfp->Self() << std::endl;

      art::FindManyP <larpandoraobj::PFParticleMetadata> pfp_pfpmeta_assns(pfp_handle, e,
                                                                           pfp_to_pfpmetadata_label); //pfp asociation to MetaData
      art::FindManyP <recob::Track> pfp_track_assns(pfp_handle, e, pfp_to_track_label); //pfp association to tracks
      art::FindManyP <anab::ParticleID> track_chi2_assn(track_handle, e,
                                                        track_to_chi2_label); // Track association to Chi2

      //Get the track score
      std::vector <art::Ptr<larpandoraobj::PFParticleMetadata>> pfp_meta_vec = pfp_pfpmeta_assns.at(pfp.key());;
      std::map<std::string, float> pfp_prop_map = pfp_meta_vec.at(0)->GetPropertiesMap();

      double track_score = -1;
      auto propertiesMapIter = pfp_prop_map.find("TrackScore");
      if (propertiesMapIter != pfp_prop_map.end()) track_score = propertiesMapIter->second;

      std::vector <art::Ptr<recob::Track>> pfp_track_vec = pfp_track_assns.at(pfp.key());


      //START the track and shower analysis
      if (pfp_track_vec.size() != 0) {
        art::Ptr <recob::Track> track = pfp_track_vec[0];

        double chi2mu[3];
        double chi2p[3];

        //Get Chi2
        const std::vector <art::Ptr<anab::ParticleID>> chi2_vec = track_chi2_assn.at(track.key());

          for (unsigned i = 0; i < chi2_vec.size(); i++) {
            const anab::ParticleID &chi2_pid = *chi2_vec[i];
            if (chi2_pid.PlaneID()) {
              int plane_id = chi2_pid.PlaneID().Plane;
              assert(plane_id < 3);

              std::vector <anab::sParticleIDAlgScores> alg_score_vector = chi2_pid.ParticleIDAlgScores();
              for (size_t i_as = 0; i_as < alg_score_vector.size(); i_as++) {
                const anab::sParticleIDAlgScores alg_score = alg_score_vector.at(i_as);
                if (alg_score.fAlgName == "Chi2") {
                  if (TMath::Abs(alg_score.fAssumedPdg) == 13) { // chi2mu
                    chi2mu[plane_id] = alg_score.fValue;
                  } else if (TMath::Abs(alg_score.fAssumedPdg) == 2212) { // chi2pr
                    chi2p[plane_id] = alg_score.fValue;
                  }

                }

              }
            }
          }
        std::cout << " chi2mu0: " << chi2mu[0] << " chi2p0: " << chi2p[0] << std::endl;
        std::cout << " chi2mu1: " << chi2mu[1] << " chi2p1: " << chi2p[1] << std::endl;
        std::cout << " chi2mu2: " << chi2mu[2] << " chi2p2: " << chi2p[2] << std::endl;
        std::cout << "TL: " << track->Length() << " TS: " << track_score << std::endl;
      }

    }

  }
}



//Create a tree when the code is initialized
void test::AnalyzeEventsMin::beginJob() {

}


void test::AnalyzeEventsMin::endJob() {
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::AnalyzeEventsMin)
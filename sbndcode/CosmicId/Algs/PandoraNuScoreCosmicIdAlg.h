#ifndef PANDORANUSCORECOSMICIDALG_H_SEEN
#define PANDORANUSCORECOSMICIDALG_H_SEEN


///////////////////////////////////////////////
// PandoraNuScoreCosmicIdAlg.h
//
// Alg to check pandora pfparticle metadata
// to see if it was tagged as a cosmic
// Ed Tyley, Jan 2020
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
// c++
#include <vector>
#include <iostream>

namespace sbnd{

  class PandoraNuScoreCosmicIdAlg {
    public:

      struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<art::InputTag> PandoraLabel {
          Name("PandoraLabel"),
            Comment("")
        };

        fhicl::Atom<art::InputTag> TpcTrackModuleLabel {
          Name("TpcTrackModuleLabel"),
            Comment("")
        };

        fhicl::Atom<float> NuScoreCut {
          Name("NuScoreCut"),
            Comment("")
        };

      };

      PandoraNuScoreCosmicIdAlg(const Config& config);

      PandoraNuScoreCosmicIdAlg(const fhicl::ParameterSet& pset) :
        PandoraNuScoreCosmicIdAlg(fhicl::Table<Config>(pset, {})()) {}

      PandoraNuScoreCosmicIdAlg();

      ~PandoraNuScoreCosmicIdAlg();

      void reconfigure(const Config& config);

      // Finds any t0s associated with track by pandora, tags if outside beam
      bool PandoraNuScoreCosmicId(recob::Track track, const art::Event& event);

      // Finds any t0s associated with pfparticle by pandora, tags if outside beam
      bool PandoraNuScoreCosmicId(recob::PFParticle pfparticle, std::map< size_t, art::Ptr<recob::PFParticle> > pfParticleMap, const art::Event& event);

      recob::PFParticle GetPFPNeutrino(recob::PFParticle pfparticle, std::map< size_t, art::Ptr<recob::PFParticle> >& pfParticleMap);

      recob::PFParticle GetPFPNeutrino(recob::PFParticle pfp, const std::vector<recob::PFParticle>& pfpVec);

    private:

      art::InputTag fPandoraLabel;
      art::InputTag fTpcTrackModuleLabel;
      float fNuScoreCut;

  };

}

#endif

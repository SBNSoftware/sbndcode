/**
 * @file    MCAssociations.h
 * @brief   This algorithm attempts to decode Track and Hit <--> MCParticle assocations
 * @author  Tracy Usher (usher@slac.stanford.edu)
 * @date    October 25, 2017
 * @see     galleryAnalysis.cpp
 * 
 */

#ifndef MCAssociations_H
#define MCAssociations_H

// LArSoft libraries
#include "gallery/Event.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorProperties.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Track.h"

// canvas libraries
#include "fhiclcpp/ParameterSet.h"

// ROOT libraries
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

// C/C++ standard libraries
#include <vector>
#include <memory> // std::unique_ptr<>


/**
 * @brief Track analysis example.
 * 
 * Configuration
 * --------------
 * 
 * * *MinLength* (real, default: 0): minimum track length, in centimetres
 * 
 */
class MCAssociations
{
public:
  
    MCAssociations(fhicl::ParameterSet const& config);
  
    void setup(const geo::GeometryCore&,
               const detinfo::DetectorProperties&,
               TDirectory*);
  
    void prepare();
  
    void doTrackHitMCAssociations(gallery::Event&);
  
    void finish();
    
private:
    double length(const recob::Track*) const;
    double length(const simb::MCParticle& part, double dx,
                  TVector3& start, TVector3& end, TVector3& startmom, TVector3& endmom,
                  unsigned int tpc = 0, unsigned int cstat = 0) const;
    
    art::InputTag            fHitProducerLabel;
    art::InputTag            fMCTruthProducerLabel;
    art::InputTag            fAssnsProducerLabel;
    art::InputTag            fTrackProducerLabel;
    std::string              fLocalDirName;
    
    geo::GeometryCore const*           fGeometry           = nullptr;
    const detinfo::DetectorProperties* fDetectorProperties = nullptr;   ///< Detector properties service
    TDirectory*                        fDir                = nullptr;
    
    std::unique_ptr<TH1>      fNTracks;
    std::unique_ptr<TH1>      fNHitsPerTrack;
    std::unique_ptr<TH1>      fTrackLength;
    std::unique_ptr<TH2>      fTrackLenVsHits;
    
    std::unique_ptr<TH1>      fNHitsPerPrimary;
    std::unique_ptr<TH1>      fPrimaryLength;
    std::unique_ptr<TH2>      fPrimaryLenVsHits;
    
    std::unique_ptr<TH1>      fNHitsPerReco;
    std::unique_ptr<TH1>      fDeltaNHits;
    std::unique_ptr<TH1>      fPrimaryRecoLength;
    std::unique_ptr<TH1>      fDeltaTrackLen;
    
    std::unique_ptr<TH1>      fPrimaryEfficiency;
    std::unique_ptr<TH1>      fPrimaryCompleteness;
    std::unique_ptr<TH1>      fPrimaryPurity;
    
    std::unique_ptr<TProfile> fPrimaryEffVsMom;
    std::unique_ptr<TProfile> fPrimaryCompVsMom;
    std::unique_ptr<TProfile> fPrimaryPurityVsMom;
    
    std::unique_ptr<TProfile> fPrimaryEffVsLen;
    std::unique_ptr<TProfile> fPrimaryCompVsLen;
    std::unique_ptr<TProfile> fPrimaryPurityVsLen;

    std::unique_ptr<TProfile> fPrimaryEffVsHits;
    std::unique_ptr<TProfile> fPrimaryCompVsHits;
    std::unique_ptr<TProfile> fPrimaryPurityVsHits;
    std::unique_ptr<TProfile> fPrimaryEffVsLogHits;
    std::unique_ptr<TProfile> fPrimaryCompVsLogHits;
    std::unique_ptr<TProfile> fPrimaryPurityVsLogHits;
}; // class MCAssociations

#endif // MCAssociations_H

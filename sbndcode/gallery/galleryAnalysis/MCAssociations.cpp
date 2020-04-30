/**
 * @file    MCAssociations.cpp
 * @brief   Does something with the tracks (implementation file).
 * @author  Tracy Usher (usher@slac.stanford.edu
 * @date    October 24, 2017
 * @see     MCAssociations.h
 * 
 */

#include "MCAssociations.h"

// LArSoft libraries
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/Hit.h"

// canvas libraries
#include "canvas/Persistency/Common/FindMany.h"

// ROOT libraries
#include "TVector3.h"

// C/C++ standard libraries
#include <algorithm> // std::count_if()


MCAssociations::MCAssociations(fhicl::ParameterSet const& config)
    : fHitProducerLabel(config.get<art::InputTag>("HitProducerLabel", "")),
      fMCTruthProducerLabel(config.get<art::InputTag>("MCTruthProducerLabel", "")),
      fAssnsProducerLabel(config.get<art::InputTag>("MCTrackAssnsProducerLabel","")),
      fTrackProducerLabel(config.get<art::InputTag>("TrackProducerLabel","")),
      fLocalDirName(config.get<std::string>("LocalDirName", ""))
  {}

void MCAssociations::setup(const geo::GeometryCore&           geometry,
                           const detinfo::DetectorProperties& detectorProperties,
                           TDirectory*                        outDir)
{
    fGeometry           = &geometry;
    fDetectorProperties = &detectorProperties;
    fDir                = outDir->mkdir(fLocalDirName.c_str());
}

void MCAssociations::prepare()
{
    if (fDir)
    {
        fNTracks = std::make_unique<TH1F>("NTrueTracks", "Number of tracks;number of tracks;events", 100, 0., 200.);
        fNTracks->SetDirectory(fDir);
        fNHitsPerTrack = std::make_unique<TH1F>("NHitsPerTrack", "Number Hits/Track", 250, 0., 250.);
        fNHitsPerTrack->SetDirectory(fDir);
        fTrackLength = std::make_unique<TH1F>("TrackLength", "Length; track length(cm)", 200, 0., 200.);
        fTrackLength->SetDirectory(fDir);
        fTrackLenVsHits = std::make_unique<TH2F>("TrackLenVsHits", "Length;Hits", 200, 0., 200., 250, 0., 250.);
        
        fNHitsPerPrimary = std::make_unique<TH1F>("NHitsPerPrimary", "Number Hits/Primary;log(# hits)", 15, 0.5, 3.5);
        fNHitsPerPrimary->SetDirectory(fDir);
        fPrimaryLength = std::make_unique<TH1F>("PrimaryLength", "Length of Primary;length(cm)", 50, 0., 350.);
        fPrimaryLength->SetDirectory(fDir);
        fPrimaryLenVsHits = std::make_unique<TH2F>("PrimaryLenVsHits", "Length;Hits", 50, 0., 350., 250, 0., 2500.);
        fPrimaryLenVsHits->SetDirectory(fDir);
        
        fNHitsPerReco = std::make_unique<TH1F>("PrimaryRecoNHits", "Number Hits/Track;log(# hits)", 15, 0.5, 3.5);
        fNHitsPerReco->SetDirectory(fDir);
        fDeltaNHits = std::make_unique<TH1F>("DeltaNHits", "Delta Number Hits", 100, -200., 200.);
        fDeltaNHits->SetDirectory(fDir);
        fPrimaryRecoLength = std::make_unique<TH1F>("PrimaryRecLength", "Length of Reco; length(cm)", 50, 0., 350.);
        fPrimaryRecoLength->SetDirectory(fDir);
        fDeltaTrackLen = std::make_unique<TH1F>("DeltaTrackLen", "Delta Track Length", 100, -100., 100.);
        fDeltaTrackLen->SetDirectory(fDir);

        fPrimaryEfficiency = std::make_unique<TH1F>("PrimaryEfficiency", "Efficiency", 101, 0., 1.01);
        fPrimaryEfficiency->SetDirectory(fDir);
        fPrimaryCompleteness = std::make_unique<TH1F>("PrimaryCompleteness", "Completeness", 101, 0., 1.01);
        fPrimaryCompleteness->SetDirectory(fDir);
        fPrimaryPurity = std::make_unique<TH1F>("PrimaryPurity", "Purity", 101, 0., 1.01);
        fPrimaryPurity->SetDirectory(fDir);
        
        fPrimaryEffVsMom = std::make_unique<TProfile>("PrimaryEffVsMom", "Efficiency vs Momentum;Momentum(GeV/c);Efficiency", 25, 0.1, 1.10, 0., 1.1);
        fPrimaryEffVsMom->SetDirectory(fDir);
        fPrimaryCompVsMom = std::make_unique<TProfile>("PrimaryCompVsMom", "Completeness vs Momentum;Momentum(GeV/c);Completeness", 25, 0.1, 1.10, 0., 1.1);
        fPrimaryCompVsMom->SetDirectory(fDir);
        fPrimaryPurityVsMom = std::make_unique<TProfile>("PrimaryPurVsMom", "Purity vs Momentum;Momentum(GeV/c);Purity", 25, 0.1, 1.10, 0., 1.1);
        fPrimaryPurityVsMom->SetDirectory(fDir);
        
        fPrimaryEffVsLen = std::make_unique<TProfile>("PrimaryEffVsLen", "Efficiency vs Length; length; Efficiency", 30, 0., 300., 0., 1.1);
        fPrimaryEffVsLen->SetDirectory(fDir);
        fPrimaryCompVsLen = std::make_unique<TProfile>("PrimaryCompVsLen", "Completeness vs Length; length; Completeness", 30, 0., 300., 0., 1.1);
        fPrimaryCompVsLen->SetDirectory(fDir);
        fPrimaryPurityVsLen = std::make_unique<TProfile>("PrimaryPurVsLen", "Purity vs Length; length; Purity", 30, 0., 300., 0., 1.1);
        fPrimaryPurityVsLen->SetDirectory(fDir);

        fPrimaryEffVsHits = std::make_unique<TProfile>("PrimaryEffVsHits", "Efficiency vs # Hits", 50, 0., 2000., 0., 1.1);
        fPrimaryEffVsHits->SetDirectory(fDir);
        fPrimaryCompVsHits = std::make_unique<TProfile>("PrimaryCompVsHits", "Completeness vs # Hits", 50, 0., 2000., 0., 1.1);
        fPrimaryCompVsHits->SetDirectory(fDir);
        fPrimaryPurityVsHits = std::make_unique<TProfile>("PrimaryPurityVsHits", "Purity vs # Hits", 50, 0., 2000., 0., 1.1);
        fPrimaryPurityVsHits->SetDirectory(fDir);
        
        fPrimaryEffVsLogHits = std::make_unique<TProfile>("PrimaryEffVsLogHits", "Efficiency vs log(# Hits);log(# Hits);Efficiency", 15, 0.5, 3.5, 0., 1.1);
        fPrimaryEffVsLogHits->SetDirectory(fDir);
        fPrimaryCompVsLogHits = std::make_unique<TProfile>("PrimaryCompVsLogHits", "Completeness vs log(# Hits);log(# Hits);Completeness)", 15, 0.5, 3.5, 0., 1.1);
        fPrimaryCompVsLogHits->SetDirectory(fDir);
        fPrimaryPurityVsLogHits = std::make_unique<TProfile>("PrimaryPurityVsLogHits", "Purity vs log(# Hits);log(# Hits);Purity)", 15, 0.5, 3.5, 0., 1.1);
        fPrimaryPurityVsLogHits->SetDirectory(fDir);
    }
    else
    {
        fNTracks.reset();
        fNHitsPerTrack.reset();
        fNHitsPerPrimary.reset();
    }
    
    return;
} // MCAssociations::prepare()

void MCAssociations::doTrackHitMCAssociations(gallery::Event& event)
{
    // First step is to recover the MCTruth object vector...
    const auto& mcParticleHandle = event.getValidHandle<std::vector<simb::MCParticle>>(fMCTruthProducerLabel);
    
    // We also need to recover the hit producer info
    const auto& hitHandle = event.getValidHandle<std::vector<recob::Hit>>(fHitProducerLabel);
    
    // Now see how many reco hits might be associated to this particle
    art::FindMany<recob::Hit, anab::BackTrackerHitMatchingData> hitsPerMCParticle(mcParticleHandle, event, fAssnsProducerLabel);
    
    // Make maps, we really like to make maps
    using HitToPartVecMap = std::map<const recob::Hit*,std::set<const simb::MCParticle*>>;
    using PartToHitVecMap = std::map<const simb::MCParticle*, std::set<const recob::Hit*>>;
    
    HitToPartVecMap hitToPartVecMap;
    PartToHitVecMap partToHitVecMap;

    // Loop through the particles
    for(int mcIdx = 0; mcIdx < mcParticleHandle->size(); mcIdx++)
    {
        try
        {
            const simb::MCParticle& mcParticle = mcParticleHandle->at(mcIdx);
            
            std::vector<const recob::Hit*> hitsVec = hitsPerMCParticle.at(mcIdx);
            
            for(const auto& hit : hitsVec)
            {
                hitToPartVecMap[hit].insert(&mcParticle);
                partToHitVecMap[&mcParticle].insert(hit);
            }
        }
        catch(...) {break;}
    }
    
    // In this section try looking at tracking. Eventually we want to move this out of here...
    // First step is to recover the MCTruth object vector...
    const auto& trackHandle = event.getValidHandle<std::vector<recob::Track>>(fTrackProducerLabel);
    
    // Now see how many reco hits might be associated to this particle
    art::FindMany<recob::Hit> hitsPerTrack(trackHandle, event, fTrackProducerLabel);
    
    // Define a mapping between tracks and associated MCParticles
    using TrackToParticleSetMap  = std::map<const recob::Track*, std::set<const simb::MCParticle*>>;
    using ParticleToTrackSetMap  = std::map<const simb::MCParticle*, std::set<const recob::Track*>>;
    using ParticleToHitSetMap    = std::map<const simb::MCParticle*, std::set<const recob::Hit*>>;
    using TrackToPartHitSetMap   = std::map<const recob::Track*, ParticleToHitSetMap>;
    using TrackToHitsVecMap      = std::map<const recob::Track*, std::vector<const recob::Hit*>>;
    
    TrackToParticleSetMap  trackToParticleSetMap;
    ParticleToTrackSetMap  particleToTrackSetMap;
    
    TrackToPartHitSetMap   trackToPartHitSetMap;
    TrackToHitsVecMap      trackToHitsVecMap;

    // Loop through the tracks and associate via the hits to MCParticles
    for(int trkIdx = 0; trkIdx < trackHandle->size(); trkIdx++)
    {
        const recob::Track& track = trackHandle->at(trkIdx);
        
        std::vector<const recob::Hit*> hitsVec = hitsPerTrack.at(trkIdx);
        
        trackToHitsVecMap[&track] = hitsVec;
        
        ParticleToHitSetMap& particleToHitSetMap = trackToPartHitSetMap[&track];
        
        for(const auto& hit : hitsVec)
        {
            HitToPartVecMap::iterator partItr = hitToPartVecMap.find(hit);
            
            if (partItr != hitToPartVecMap.end())
            {
                for(const auto& particle : partItr->second)
                {
                    particleToHitSetMap[particle].insert(hit);
                    trackToParticleSetMap[&track].insert(particle);
                    particleToTrackSetMap[particle].insert(&track);
                }
            }
        }
    }
    
    // *****************************************************************************************
    // The bits below here should eventually be moved into their own analyzer algorithm
    // but we are in a hurry now so do it all here...
    // Ok, at this point we should be able to relate MCParticles to tracks and hits
    // Let's start by just looking at the primary particle
    const simb::MCParticle& primaryParticle = mcParticleHandle->at(0);
    
    ParticleToTrackSetMap::iterator partTrackItr = particleToTrackSetMap.find(&primaryParticle);
    
    // Define the parameters we want...
    int   numPrimaryHitsTotal = partToHitVecMap[&primaryParticle].size();
    
    // If there are NO reconstructed hits associated to this particle then we don't count
    // But this should really be a check on fiducial volume I think...
    if (numPrimaryHitsTotal > 0)
    {
        const recob::Track* bestTrack(0);
        float               efficiency(0.);
        float               completeness(0.);
        float               purity(0.);
        int                 numTrackHits(0);
        
        // Here we find the best matched track to the MCParticle.
        // Nothing exciting, most hits wins sort of thing...
        if (partTrackItr != particleToTrackSetMap.end())
        {
            // Recover the longest track...
            for(const auto& track : partTrackItr->second)
            {
                if (trackToPartHitSetMap[track][partTrackItr->first].size() > numTrackHits)
                {
                    bestTrack    = track;
                    numTrackHits = trackToPartHitSetMap[track][partTrackItr->first].size();
                }
            }
        
            if (bestTrack)
            {
                int numPrimaryHitsMatch = numTrackHits;
                int numTrackHitsTotal   = trackToHitsVecMap[bestTrack].size();
        
                completeness = float(numPrimaryHitsMatch) / float(numPrimaryHitsTotal);
                purity       = float(numPrimaryHitsMatch) / float(numTrackHitsTotal);
                
                if (completeness > 0.2) efficiency = 1.;
            }
        }
    
        // Calculate the length of this mc particle inside the fiducial volume.
        TVector3 mcstart;
        TVector3 mcend;
        TVector3 mcstartmom;
        TVector3 mcendmom;
    
        double xOffset(0.);
    
        double mcTrackLen = length(primaryParticle, xOffset, mcstart, mcend, mcstartmom, mcendmom);
        double trackLen   = 0.;
        
        if (bestTrack) trackLen = length(bestTrack);
    
        fNTracks->Fill(partToHitVecMap.size(), 1.);
        fNHitsPerPrimary->Fill(std::log10(double(partToHitVecMap[&primaryParticle].size())), 1.);
        fPrimaryLength->Fill(mcTrackLen, 1.);
        fPrimaryLenVsHits->Fill(mcTrackLen, partToHitVecMap[&primaryParticle].size(), 1.);
        
        fPrimaryRecoLength->Fill(trackLen, 1.);
        fDeltaTrackLen->Fill(trackLen-mcTrackLen, 1.);
        
        fNHitsPerReco->Fill(std::log10(numTrackHits), 1.);
        fDeltaNHits->Fill(numTrackHits - int(partToHitVecMap[&primaryParticle].size()), 1.);

        // Loop through the particles again to histogram some secondary info...
        for(int mcIdx = 0; mcIdx < mcParticleHandle->size(); mcIdx++)
        {
            try
            {
                const simb::MCParticle& mcParticle = mcParticleHandle->at(mcIdx);
            
                if (!partToHitVecMap[&mcParticle].empty())
                {
                    // Calculate the length of this mc particle inside the fiducial volume.
                    double secTrackLen = length(mcParticle, xOffset, mcstart, mcend, mcstartmom, mcendmom);
                
                    fNHitsPerTrack->Fill(partToHitVecMap[&mcParticle].size(), 1.);
                    fTrackLength->Fill(secTrackLen, 1.);
                    fTrackLenVsHits->Fill(secTrackLen, partToHitVecMap[&mcParticle].size(), 1.);
                }
            }
            catch(...) {break;}
        }
    
        // Final sets of plots
        fPrimaryEfficiency->Fill(efficiency, 1.);
        fPrimaryCompleteness->Fill(completeness, 1.);
        fPrimaryPurity->Fill(purity, 1.);
        
        fPrimaryEffVsHits->Fill(numPrimaryHitsTotal, efficiency, 1.);
        fPrimaryCompVsHits->Fill(numPrimaryHitsTotal, completeness, 1.);
        fPrimaryPurityVsHits->Fill(numPrimaryHitsTotal, purity, 1.);
        
        double partMom = primaryParticle.P();
        fPrimaryEffVsMom->Fill(partMom, efficiency, 1.);
        fPrimaryCompVsMom->Fill(partMom, completeness, 1.);
        fPrimaryPurityVsMom->Fill(partMom, purity, 1.);
        
        fPrimaryEffVsLen->Fill(mcTrackLen, efficiency, 1.);
        fPrimaryCompVsLen->Fill(mcTrackLen, completeness, 1.);
        fPrimaryPurityVsLen->Fill(mcTrackLen, purity, 1.);

        double logNumHits = std::log10(numPrimaryHitsTotal);
        fPrimaryEffVsLogHits->Fill(logNumHits, efficiency, 1.);
        fPrimaryCompVsLogHits->Fill(logNumHits, completeness, 1.);
        fPrimaryPurityVsLogHits->Fill(logNumHits, purity, 1.);
    }

    return;
} // MCAssociations::processTracks()

void MCAssociations::finish()
{
    if (fNTracks)
    {
        fDir->cd();
        fNTracks->Write();
        fNHitsPerTrack->Write();
        fTrackLength->Write();
        fTrackLenVsHits->Write();
        fNHitsPerPrimary->Write();
        fPrimaryLength->Write();
        fPrimaryLenVsHits->Write();
        fPrimaryEfficiency->Write();
        fPrimaryCompleteness->Write();
        fPrimaryPurity->Write();
        fPrimaryEffVsHits->Write();
        fPrimaryCompVsHits->Write();
        fPrimaryPurityVsHits->Write();
        fPrimaryEffVsMom->Write();
        fPrimaryCompVsMom->Write();
        fPrimaryPurityVsMom->Write();
        fPrimaryEffVsLen->Write();
        fPrimaryCompVsLen->Write();
        fPrimaryPurityVsLen->Write();
        fPrimaryEffVsLogHits->Write();
        fPrimaryCompVsLogHits->Write();
        fPrimaryPurityVsLogHits->Write();
        fPrimaryRecoLength->Write();
        fDeltaTrackLen->Write();
        fNHitsPerReco->Write();
        fDeltaNHits->Write();
    }
} // MCAssociations::finish()

// Length of reconstructed track.
//----------------------------------------------------------------------------
double MCAssociations::length(const recob::Track* track) const
{
    return track->Length();
}

// Length of MC particle.
//----------------------------------------------------------------------------
double MCAssociations::length(const simb::MCParticle& part, double dx,
                              TVector3& start, TVector3& end, TVector3& startmom, TVector3& endmom,
                              unsigned int tpc, unsigned int cstat) const
{
    // Need the readout window size...
    double readOutWindowSize = fDetectorProperties->ReadOutWindowSize();
    
    double result = 0.;
    TVector3 disp;
    int n = part.NumberTrajectoryPoints();
    bool first = true;
    
    // The following for debugging purposes
    int findTrackID(-1);
    
    if (part.TrackId() == findTrackID) std::cout << ">>> length, mcpart: " << part << std::endl;
    
    // Loop over the complete collection of trajectory points
    for(int i = 0; i < n; ++i)
    {
        TVector3 posInTPC(0.,0.,0.);
        TVector3 posVec = part.Position(i).Vect();
        double   pos[]  = {posVec.X(),posVec.Y(),posVec.Z()};
        
        // Try identifying the TPC we are in
        unsigned int tpc(0);
        unsigned int cstat(0);

        // Need to make sure this position is in an active region of the TPC
        // If the particle is not in the cryostat then we skip
        try
        {
            const geo::TPCGeo& tpcGeo = fGeometry->PositionToTPC(pos, tpc, cstat);
            
            TVector3 activePos = posVec - tpcGeo.GetActiveVolumeCenter();
            
            if (part.TrackId() == findTrackID)
                std::cout << "   --> traj point: " << i << ", pos: " << posVec.X() << "/" << posVec.Y() << "/" << posVec.Z() << ", active pos: " << activePos.X() << "/" << activePos.Y() << "/" << activePos.Z() << std::endl;
            
            if (std::fabs(activePos.X()) > tpcGeo.ActiveHalfWidth() || std::fabs(activePos.Y()) > tpcGeo.ActiveHalfHeight() || std::fabs(activePos.Z()) > 0.5 * tpcGeo.ActiveLength()) continue;
            
            posInTPC = TVector3(activePos.X() + tpcGeo.ActiveHalfWidth(), activePos.Y(), activePos.Z());
        } catch(...) {continue;}
        
        // Make fiducial cuts.
        // There are two sets here:
        // 1) We check the original x,y,z position of the trajectory points and require they be
        //    within the confines of the physical TPC
        // 2) We then check the timing of the presumed hit and make sure it lies within the
        //    readout window for this simulation
        pos[0] += dx;
        double ticks = fDetectorProperties->ConvertXToTicks(posInTPC.X(), 0, tpc, cstat);
        
        if (part.TrackId() == findTrackID)
            std::cout << "   ==> tpc: " << tpc << ", cstat: " << cstat << ", ticks: " << ticks << std::endl;

        // Currently it appears that the detector properties are not getting initialized properly and the returned
        // number of ticks is garbage so we need to skip this for now. Will be important when we get CR's
//        if(ticks >= 0. && ticks < readOutWindowSize)
        {
            if(first)
            {
                start = pos;
                startmom = part.Momentum(i).Vect();
            }
            else
            {
                disp -= pos;
                result += disp.Mag();
            }
            first = false;
            disp = pos;
            end = pos;
            endmom = part.Momentum(i).Vect();
        }
        
        if (part.TrackId() == findTrackID)
        {
            try
            {
                std::cout << ">>> Track #" << findTrackID << ", pos: " << posVec.X() << ", " << posVec.Y() << ", " << posVec.Z() << ", ticks: " << ticks << ", nearest Y wire: ";
                geo::WireID wireID = fGeometry->NearestWireID(pos, 2);
                std::cout << wireID << std::endl;
            }
            catch(...) {}
        }
    }
    
    return result;
}

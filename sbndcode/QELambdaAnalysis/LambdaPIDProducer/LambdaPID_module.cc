////////////////////////////////////////////////////////////////////////
// Class:       LambdaPID
// Plugin Type: Producer (art v3_05_01)
// File:        LambdaPID_module.cc
////////////////////////////////////////////////////////////////////////

#include "sbndcode/QELambdaAnalysis/LambdaPIDProducer/LambdaPID_module.hh"

#include <cmath>
#include <limits>
#include <map>
#include <optional>


namespace {
  template <typename T>
  bool fetchHandle(art::Event const& e, std::string const& label, art::Handle<T>& handle)
  {
    handle = e.getHandle<T>(label);
    if (!handle.isValid()) {
      mf::LogWarning("LambdaPID")
        << "Product not found (label '" << label << "'), event "
        << e.id().run() << "_" << e.id().subRun() << "_" << e.id().event();
      return false;
    }
    return true;
  }
}


// -------- Constructor
qelambda::LambdaPID::LambdaPID(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , fReco2Label( p.get< std::string >("Reco2Label") )
  , fTrackLabel(p.get<std::string>("TrackLabel"))
  , fCalorimetryLabel(p.get<std::string>("CalorimetryLabel"))
  , fParticleIDLabel(p.get<std::string>("ParticleIDLabel"))

  
{
  // PID Lambda products
  produces<qelambda::QELambdaPID>();
}

std::vector<art::Ptr<recob::PFParticle>> qelambda::LambdaPID::getSlicePrimPfp(art::Event const &e,std::vector<art::Ptr<recob::PFParticle>> slcPfpVec)
{

  std::vector<art::Ptr<recob::PFParticle>> primPfpVec;

  // Loop over pfp in slice to get prim information
  for(const art::Ptr<recob::PFParticle> &pfp : slcPfpVec) {
    // Only look for primary pfps
    if (pfp->IsPrimary()) {
      primPfpVec.push_back(pfp);
    }
  }

  return primPfpVec;

}

// -------- Main function
void qelambda::LambdaPID::produce(art::Event &e)
{

  fLambdaPID.SetDefault();
  fProtonMom.SetXYZ(0,0,0);
  fPionMom.SetXYZ(0,0,0);
  fRecoStatus = false;

  auto putProduct = [&]() {
    e.put(std::make_unique<qelambda::QELambdaPID>(fLambdaPID));
  };

  constexpr double kNoScore = std::numeric_limits<double>::max();

  const bool hasReco2 = !fReco2Label.empty();
  const std::string pfpLabel = hasReco2 ? fReco2Label : fTrackLabel;

  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  art::Handle<std::vector<recob::Track>>      trackHandle;
  art::Handle<std::vector<recob::Slice>>      sliceHandle;     

  // Evaluate both (no short-circuit) so every missing product gets logged
  bool handlesOk = true;
  handlesOk &= fetchHandle(e, pfpLabel,    pfpHandle);
  handlesOk &= fetchHandle(e, fTrackLabel, trackHandle);

  if (!handlesOk) {
    mf::LogWarning("LambdaPID") << "  Missing required products, writing default LambdaPID";
    putProduct();
    return;
  }

  if (hasReco2) {
    sliceHandle = e.getHandle<std::vector<recob::Slice>>(fReco2Label);
  }

  if (!sliceHandle.isValid()) {
    mf::LogWarning("LambdaPID") << "Slice or Cluster handle not valid, writing default LambdaPID";
    putProduct();
    return;
  }

  // --- Associations
  // PFParticle to track
  art::FindManyP<recob::Track> pfp_track_assns (pfpHandle, e, fTrackLabel);
  // Track to calorimetry
  art::FindManyP<anab::Calorimetry> track_to_calo_assns(trackHandle, e, fCalorimetryLabel);
  // Track to PID
  art::FindManyP<anab::ParticleID> track_to_pid_assns(trackHandle, e, fParticleIDLabel);
  // PFParticle to vertex
  art::FindManyP<recob::Vertex> pfp_vertex_assns(pfpHandle, e, pfpLabel);

  if (!pfp_track_assns.isValid() || !track_to_calo_assns.isValid() || !track_to_pid_assns.isValid()) {
    mf::LogWarning("LambdaPID") << "PFParticle->Track, Track->Calorimetry or Track->ParticleID assns not available, "
                                << "writing default LambdaPID";
    putProduct();
    return;
  }

  // --- Variables initialization
  // Direction of the low/high ionizing tracks
  TVector3 lowIonizingDir(0,0,0); // (0,0,0) is a dummy value, it will be overwritten
  TVector3 highIonizingDir(0,0,0); // (0,0,0) is a dummy value, it will be overwritten
  // Minimum distance to the vertex for low/high ionizing tracks
  double minDistLowIonizing = 1e9;
  double minDistHighIonizing = 1e9;

  // --- Get candidate PFParticles
  std::vector<art::Ptr<recob::PFParticle>> pfpVect;
  std::vector<art::Ptr<recob::Slice>> sliceVect;

  if (sliceHandle.isValid()) art::fill_ptr_vector(sliceVect, sliceHandle);

  bool pfpFromSlice = false;
  if (!sliceVect.empty()) {

    art::FindManyP<recob::PFParticle> slice_pfp_assns(sliceHandle, e, fReco2Label);
    if (slice_pfp_assns.isValid()) {
      pfpVect = slice_pfp_assns.at(sliceVect[0].key());
      pfpFromSlice = true;
    } else {
      mf::LogWarning("LambdaPID") << "Slice->PFParticle assns not available, using all PFParticles";
    }
  }
  if (!pfpFromSlice) {
    // No slices: use all PFParticles
    art::fill_ptr_vector(pfpVect, pfpHandle);
  }

  // Primary PFParticle and its vertex
  art::Ptr<recob::PFParticle> prim;
  for (auto const& p : pfpVect) {
    if (p->IsPrimary()) { prim = p; break; }
  }

  bool hasVertex = false;
  if (prim.isNonnull()) {
    auto const& vtxVec = pfp_vertex_assns.at(prim.key());
    if (!vtxVec.empty()) {
      const auto pos = vtxVec.front()->position();
      fLambdaPID.VtxX = pos.X();
      fLambdaPID.VtxY = pos.Y();
      fLambdaPID.VtxZ = pos.Z();
      hasVertex = true;
    }
  }
 
  if (!hasVertex) {
    mf::LogWarning("LambdaPID") << "No primary PFParticle vertex found, writing default LambdaPID";
    putProduct();
    return;
  }

  // --- Loop over the PFParticles
  for(const art::Ptr<recob::PFParticle> &pfp : pfpVect){

    // std::cout<<"   ** PFParticle: "<<pfp->Self()<<"      PDG:"<<pfp->PdgCode()<<"  Primary="<<pfp->IsPrimary()<<" Mother="<<pfp->Parent()<<std::endl;

    // Get the associated tracks
    std::vector<art::Ptr<recob::Track>> trackVector = pfp_track_assns.at(pfp.key());
    if (trackVector.empty()) continue;

    // Loop over the tracks
    for (size_t i = 0; i < trackVector.size(); ++i) {
      art::Ptr<recob::Track> track = trackVector[i];
      if (track.isNull()) continue;
      
      // --- Get the associated calorimetry and PID objects
      std::vector<art::Ptr<anab::Calorimetry>> caloV = track_to_calo_assns.at(track.key());
      std::vector<art::Ptr<anab::ParticleID>> pidV = track_to_pid_assns.at(track.key());

       // --- Calorimetry: use the plane with most hits (collection plane has priority)
      double KE = 0.0;
      size_t maxNHit = 0;
      int bestPlane = -1;
      bool foundCalo = false;
 
      for (const auto& caloPtr : caloV) {
        if (caloPtr.isNull()) continue;
 
        const anab::Calorimetry& calo = *caloPtr;
        const int planeID = calo.PlaneID().Plane;
        const size_t nHits = calo.dEdx().size();
 
        if (!foundCalo ||
            nHits > maxNHit ||
            (nHits == maxNHit && planeID == 2 && bestPlane != 2)) {
          foundCalo = true;
          maxNHit = nHits;
          bestPlane = planeID;
          KE = calo.KineticEnergy();
        }
      }
 
      if (!foundCalo || !(KE > 0.0) || bestPlane < 0) continue;

      // --- Loop over the PID objects
      double scorePion = 0;
      double scoreproton = 0;
      int nPionScores = 0;
      int nprotonScores = 0;

      for (size_t j = 0; j < pidV.size(); ++j) {
        if (pidV[j].isNull()) continue;
        const anab::ParticleID& pid = *pidV[j]; 

        if (pid.PlaneID().Plane != static_cast<unsigned int>(bestPlane)) continue;

        for(size_t k=0; k<pid.ParticleIDAlgScores().size(); k++){
          const auto& algScore = pid.ParticleIDAlgScores()[k];

          // skip 0 and 321 PDG, only possible templates are 211, 13 and 2212s
          if(algScore.fAssumedPdg == 0 || algScore.fAssumedPdg == 321) continue;

          if(algScore.fAssumedPdg == 211 || algScore.fAssumedPdg == 13){
            scorePion = scorePion + algScore.fValue;
            nPionScores++;
          }
          else if(algScore.fAssumedPdg == 2212){
            scoreproton = scoreproton + algScore.fValue;
            nprotonScores++;
          }
          
        }

      } // end of PID loop

      // Average the scores
      scorePion = nPionScores > 0 ? scorePion / nPionScores : kNoScore;
      scoreproton = nprotonScores > 0 ? scoreproton / nprotonScores : kNoScore;

      // --- Fill Lambda PID struct
      if (nPionScores == 0 && nprotonScores == 0) {
        mf::LogWarning("LambdaPID") << "     --  Track without PID scores in the selected views, skipping";
        continue;
      }

      // Select the best PID
      // double minScore = std::min(scorePion, scoreproton);
      int bestPDG = scorePion<scoreproton ? 211 : 2212;

      // Print results
      // std::cout << "     --  KE: "<<KE<<"  Best PID: "<<bestPDG<<"  Min Score: "<<minScore<<std::endl;

      // --- Check if the track is closest to the neutrino vertex
      double dist = sqrt( pow(track->Vertex().X()-fLambdaPID.VtxX, 2) + pow(track->Vertex().Y()-fLambdaPID.VtxY, 2) + pow(track->Vertex().Z()-fLambdaPID.VtxZ, 2) );
        
      if(bestPDG == 2212 && dist<minDistHighIonizing){
        fLambdaPID.ProtonKE = KE;
        highIonizingDir.SetXYZ(track->StartDirection().X(), track->StartDirection().Y(), track->StartDirection().Z());
        minDistHighIonizing = dist;
      }
      else if( (bestPDG == 211 || bestPDG == 13) && dist<minDistLowIonizing){
        fLambdaPID.PionKE = KE;
        lowIonizingDir.SetXYZ(track->StartDirection().X(), track->StartDirection().Y(), track->StartDirection().Z());
        minDistLowIonizing = dist;
      }        
    } // end of track loop
      
  }// end of PFParticle loop

  if (fLambdaPID.PionKE>0) {
    double pPionMag = sqrt( fLambdaPID.PionKE * (fLambdaPID.PionKE + 2 * PionMass) );
    fPionMom.SetXYZ(pPionMag*lowIonizingDir.X(), pPionMag*lowIonizingDir.Y(), pPionMag*lowIonizingDir.Z());
    fLambdaPID.PionPX = fPionMom.X();
    fLambdaPID.PionPY = fPionMom.Y();
    fLambdaPID.PionPZ = fPionMom.Z();
  }

  if (fLambdaPID.ProtonKE>0) {
    double pProtonMag = sqrt( fLambdaPID.ProtonKE * (fLambdaPID.ProtonKE + 2 * ProtonMass) );
    fProtonMom.SetXYZ(pProtonMag*highIonizingDir.X(), pProtonMag*highIonizingDir.Y(), pProtonMag*highIonizingDir.Z());
    fLambdaPID.ProtonPX = fProtonMom.X();
    fLambdaPID.ProtonPY = fProtonMom.Y();
    fLambdaPID.ProtonPZ = fProtonMom.Z();
  }

  putProduct();

}
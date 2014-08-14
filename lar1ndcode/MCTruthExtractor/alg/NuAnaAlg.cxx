
#include "NuAnaAlg.h"

namespace lar1nd{

  NuAnaAlg::NuAnaAlg(){
    std::cout << "This is the fucking default ctor.  Deal with it.\n";
  }


  void NuAnaAlg::configureGeometry(art::ServiceHandle<geo::Geometry> geom){
  
    xlow  = - geom -> DetHalfWidth();
    xhigh = geom -> DetHalfWidth();
    ylow  = - geom -> DetHalfHeight();
    yhigh = geom -> DetHalfHeight();
    zlow  = 0.0;
    zhigh = geom -> DetLength();

    return;
  }

  // get the basic neutrino info:
  void NuAnaAlg::packNeutrinoInfo(simb::MCNeutrino * neutrino, 
                                  int& nuchan,
                                  int& inno,
                                  double& enugen,
                                  int& isCC,
                                  int& mode,
                                  double& thetaLep,
                                  double& phiLep,
                                  double& Elep,
                                  TLorentzVector& neutMom,
                                  TVector3& vertex){
    // Fill in all the variables that need to be filled!

    // 1000 is the offset value used for NUANCE
    nuchan = (neutrino->InteractionType() - 1000); 
    inno = neutrino->Nu().PdgCode();
    enugen = neutrino->Nu().E();
    // Is it a CC or NC interaction?
    if (neutrino->CCNC()==simb::kCC) isCC = 1;
    else if (neutrino->CCNC()==simb::kNC)  isCC = 0;
    mode = neutrino->Mode();

    simb::MCParticle lepton = neutrino->Lepton();
    // thetaLep is the angle the lepton makes with the beam
    double tempNumerator = lepton.Px()*lepton.Px() + lepton.Py()*lepton.Py();
    thetaLep = atan(sqrt(tempNumerator)/lepton.Pz());
    phiLep = atan(lepton.Px()/lepton.Py());
    Elep = lepton.E();
    vertex = neutrino->Nu().Position().Vect();

    neutMom = neutrino->Nu().Momentum();

  }

  void NuAnaAlg::packFluxInfo(art::Ptr<simb::MCFlux > flux, 
                              int& ptype, int& tptype, int& ndecay,
                              TVector3& ParentVertex,
                              TVector3& nuParentMomAtDecay,
                              TVector3& nuParentMomAtProd,
                              TVector3& nuParentMomTargetExit){

    ptype = flux -> fptype;
    tptype = flux->ftptype;
    ndecay = flux->fndecay;

    ParentVertex.SetXYZ         (flux -> fvx,   flux -> fvy,   flux -> fvz);
    nuParentMomAtDecay.SetXYZ   (flux -> fpdpx, flux -> fpdpy, flux -> fpdpz);
    nuParentMomAtProd.SetXYZ    (flux -> fppdxdz, flux -> fppdydz, flux -> fpppz);
    nuParentMomTargetExit.SetXYZ(flux -> ftpx,  flux -> ftpy,  flux -> ftpz);

  }

  // Pack up the genie info:
  void NuAnaAlg::packGenieInfo(art::Ptr<simb::MCTruth>  truth,
                               std::vector<int> & GeniePDG,
                               std::vector<TLorentzVector>& GenieMomentum,
                               std::vector<std::string>& GenieProc,
                               int& NPi0FinalState,
                               int& NGamma){
    int i = 0;
    while( i < truth -> NParticles()){
      auto part = truth -> GetParticle(i);
      if (part.StatusCode() == 1){
        GeniePDG.push_back(part.PdgCode());
        GenieMomentum.push_back(part.Momentum());
        GenieProc.push_back(part.Process());
        if (part.PdgCode() == 22) NGamma ++;
        if (part.PdgCode() == 111) NPi0FinalState ++;
        // if (abs(part.PdgCode()) == 211) NChargedPions ++;
      }
      i++;
    }
  }

  simb::MCParticle NuAnaAlg::getParticleByID(
          art::Handle< std::vector<simb::MCParticle> > & mclistLARG4, int TrackId) const{
    for(unsigned int i = 0; i < mclistLARG4 -> size(); i ++){
      if ( mclistLARG4 -> at(i).TrackId() == TrackId){
        return mclistLARG4 -> at(i);
      }
    }
    simb::MCParticle part;
    return part;
  }
  simb::MCParticle NuAnaAlg::getParticleByPDG(
          art::Handle< std::vector<simb::MCParticle> > & mclistLARG4, int PDG) const{
    
    for(unsigned int i = 0; i < mclistLARG4 -> size(); i ++){
      if ( mclistLARG4 -> at(i).PdgCode() == PDG){
        return mclistLARG4 -> at(i);
      }
    }
    simb::MCParticle part;
    return part;
  }
  bool NuAnaAlg::isInTPC(TVector3 & v) const{
    if (v.X() > xhigh || v.X() < xlow) return false;
    if (v.Y() > yhigh || v.Y() < ylow) return false;
    if (v.Z() > zhigh || v.Z() < zlow) return false;
    return true;
  }

  // Method to take in a photon and determine where it started converting.
  // Looks at the photon's energy at each step.  
  void NuAnaAlg::GetPhotonConversionInfo(simb::MCParticle& photon,
                                      TLorentzVector& ConversionPos,
                                      TLorentzVector& ConversionMom){
  // Just loop over the photons Trajectory points in momentum space 
  // and wait for a change. If it gets through the whole trajectory,
  // the conversion point must be the end of the trajectory.
  // When it changed, return the trajectory point immediately 
  // before to get the last point of the photon that's not converted.
  
  // Going to be watching for changes in energy, 
  // so we'd better know the start energy.
  double E = photon.E(0);
  // std::cout << "Starting Energy is " << E << std::endl;
  // std::cout << "Looking at photon 4 vector, momentum currently starts \n\t("
  //     << photon.Momentum().X() << ", "
  //     << photon.Momentum().Y() << ", "
  //     << photon.Momentum().Z() << ", "
  //     << photon.Momentum().T() << ") ";
  // std::cout << " At position \n\t( "
  //     << photon.Position().X() << ", "
  //     << photon.Position().Y() << ", "
  //     << photon.Position().Z() << ", "
  //     << photon.Position().T() << ") " << std::endl;
  //Catch some special cases first. 
  if (photon.NumberTrajectoryPoints() == 0) return;   
  if (photon.NumberTrajectoryPoints() == 1){
    ConversionPos = photon.EndPosition();
    ConversionMom = photon.EndMomentum();
    return;
  }

  for(unsigned int i = 1; i < photon.NumberTrajectoryPoints(); i ++){
    //ok, check if the energy has changed.
    // std::cout << "Looking at photon 4 vector, momentum currently is \n\t("
    //     << photon.Momentum(i).X() << ", "
    //     << photon.Momentum(i).Y() << ", "
    //     << photon.Momentum(i).Z() << ", "
    //     << photon.Momentum(i).T() << ") ";
    // std::cout << " At position \n\t("
    //     << photon.Position(i).X() << ", "
    //     << photon.Position(i).Y() << ", "
    //     << photon.Position(i).Z() << ", "
    //     << photon.Position(i).T() << ") " << std::endl;
    if (photon.E(i) != E) {
      //then the energy is different.  Must be scattering or something.
      //set the conversion points and bail!
      ConversionPos = photon.Position(i); 
      ConversionMom = photon.Momentum(i-1); //<- This should be OK since i starts at 1.
//      std::cout << "Conversion Pos is " << ConversionPos << std::endl;
      return;
    }
  } 

  // If we made it out here, the photon must have ended without changing energy.
  // So send back end position, momentum.
  ConversionPos = photon.EndPosition();
  ConversionMom = photon.EndMomentum();

  return;
  }

}








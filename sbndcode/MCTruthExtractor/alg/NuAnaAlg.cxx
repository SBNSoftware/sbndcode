#include "sbndcode/MCTruthExtractor/alg/NuAnaAlg.h"
////#define CUSTOM_NUTOOLS

namespace sbnd{

  NuAnaAlg::NuAnaAlg(){
  }

  // NuAnaAlg::~NuAnaAlg(){
  //   // if (reweight) delete reweight;
  //   // reweight = NULL;
  // }

  void NuAnaAlg::configureGeometry(art::ServiceHandle<geo::Geometry> geom){

    xlow  = - geom -> DetHalfWidth();
    xhigh = geom -> DetHalfWidth();
    ylow  = - geom -> DetHalfHeight();
    yhigh = geom -> DetHalfHeight();
    zlow  = 0.0;
    zhigh = geom -> DetLength();

    if (geom -> DetectorName() == "uboone_basic"){
      xlow = 0.0;
      xhigh = 2*xhigh;
    }

    std::cout << "The dimensions of this detector are read to be:\n"
              << "  x: " << xlow << " to " << xhigh << "\n"
              << "  y: " << ylow << " to " << yhigh << "\n"
              << "  z: " << zlow << " to " << zhigh << "\n";

    return;
  }

/**
*
* This function configures the reweighting machinery.  I am guess that,
* due to the overhead in configuring the genie reweighting, it's faster
* to make one reweight object for each weight needed.
* This will also make a "total" reweight object that will set all the
* switches on for all weighting parameters.
*/
  void NuAnaAlg::configureReWeight(const std::vector<reweight> & weights,
                         const std::vector<std::vector<float>>& reweightingSigmas){

    // Expand the vector to the correct number of physical knobs plus 1
    reweightVector.resize(weights.size()+1);

    // if (weights.size()+1 != reweightingSigmas.size()){
    //   std::cerr << "Error configuring the reweights, the number of weights must be"
    //             << " equal to the number of boundaries (range).\n";
    //   exit(-1);
    // }

    // don't forget the last vector with all of the weights
    reweightVector.back().resize(reweightingSigmas.front().size());
    for (auto & ptr : reweightVector.back()) ptr = new rwgt::NuReweight;

    for (unsigned int i_weight = 0; i_weight < reweightingSigmas.front().size(); ++i_weight)
    {
      // reweightVector.back().at(i_weight)
      //     -> ReweightNCEL(reweightingSigmas[kNCEL][i_weight]);
      reweightVector.back().at(i_weight)
          -> ReweightQEMA(reweightingSigmas[kQEMA][i_weight]);
      // reweightVector.back().at(i_weight)
      //     -> ReweightQEVec(reweightingSigmas[kQEVec][i_weight]);
      reweightVector.back().at(i_weight)
          -> ReweightResGanged(reweightingSigmas[kResGanged][i_weight]);
      reweightVector.back().at(i_weight)
          -> ReweightCCRes(reweightingSigmas[kCCRes][i_weight]);
      reweightVector.back().at(i_weight)
          -> ReweightNCRes(reweightingSigmas[kNCRes][i_weight]);
      // reweightVector.back().at(i_weight)
          // -> ReweightCoh(reweightingSigmas[kCoh][i_weight]);
      reweightVector.back().at(i_weight)
          -> ReweightNonResRvp1pi(reweightingSigmas[kNonResRvp1pi][i_weight]);
      reweightVector.back().at(i_weight)
          -> ReweightNonResRvbarp1pi(reweightingSigmas[kNonResRvbarp1pi][i_weight]);
      reweightVector.back().at(i_weight)
          -> ReweightNonResRvp2pi(reweightingSigmas[kNonResRvp2pi][i_weight]);
      reweightVector.back().at(i_weight)
          -> ReweightNonResRvbarp2pi(reweightingSigmas[kNonResRvbarp2pi][i_weight]);
      // reweightVector.back().at(i_weight)
          // -> ReweightResDecay(reweightingSigmas[kResDecay][i_weight]);
      reweightVector.back().at(i_weight)
          -> ReweightNC(reweightingSigmas[kNC][i_weight]);
      // reweightVector.back().at(i_weight)
          // -> ReweightDIS(reweightingSigmas[kDIS][i_weight]);
      reweightVector.back().at(i_weight)
          -> ReweightDISnucl(reweightingSigmas[kDISnucl][i_weight]);
      // reweightVector.back().at(i_weight)
      //     -> ReweightAGKY(reweightingSigmas[kAGKY][i_weight]);
    }


    // loop over the physical knobs and expand to the correct number of weights
    for (unsigned int i_reweightingKnob = 0;
         i_reweightingKnob < reweightVector.size()-1;
         i_reweightingKnob++)
    {

      // resize this row to accomodate all of the weight points
      reweightVector[i_reweightingKnob].resize(
            reweightingSigmas[i_reweightingKnob].size());

      for (unsigned int weight_point = 0;
           weight_point < reweightingSigmas[i_reweightingKnob].size();
           weight_point++){

        // Figure out what is the value going in to this reweight
        // double stepSize = (reweightingSigmas[i_reweightingKnob]
        //                - rangeLow[i_reweightingKnob])/(nWeights-1);
        // double reweightingValue = rangeLow[i_reweightingKnob]
        //                         + weight_point*stepSize;

        reweightVector[i_reweightingKnob][weight_point] = new rwgt::NuReweight;

        switch (weights[i_reweightingKnob]){
          case kNCEL:
            // reweightVector[i_reweightingKnob][weight_point]
            //   -> ReweightNCEL(reweightingSigmas[kNCEL][weight_point]);
            break;
          case kQEMA:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightQEMA(reweightingSigmas[kQEMA][weight_point]);
            break;
          case kQEVec:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightQEVec(reweightingSigmas[kQEVec][weight_point]);
            break;
          case kResGanged:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightResGanged(reweightingSigmas[kResGanged][weight_point]);
            break;
          case kCCRes:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightCCRes(reweightingSigmas[kCCRes][weight_point]);
            break;
          case kNCRes:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightNCRes(reweightingSigmas[kNCRes][weight_point]);
            break;
          case kCoh:
          //   reweightVector[i_reweightingKnob][weight_point]
          //     -> ReweightCoh(reweightingSigmas[kCoh][weight_point]);
            break;
          case kNonResRvp1pi:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightNonResRvp1pi(reweightingSigmas[kNonResRvp1pi][weight_point]);
            break;
          case kNonResRvbarp1pi:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightNonResRvbarp1pi(reweightingSigmas[kNonResRvbarp1pi][weight_point]);
            break;
          case kNonResRvp2pi:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightNonResRvp2pi(reweightingSigmas[kNonResRvp2pi][weight_point]);
            break;
          case kNonResRvbarp2pi:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightNonResRvbarp2pi(reweightingSigmas[kNonResRvbarp2pi][weight_point]);
            break;
          case kResDecay:
            // reweightVector[i_reweightingKnob][weight_point]
            //   -> ReweightResDecay(reweightingSigmas[kResDecay][weight_point]);
            break;
          case kNC:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightNC(reweightingSigmas[kNC][weight_point]);
            break;
          case kDIS:
            // reweightVector[i_reweightingKnob][weight_point]
            //   -> ReweightDIS(reweightingSigmas[kDIS][weight_point]);
            break;
          case kDISnucl:
            reweightVector[i_reweightingKnob][weight_point]
              -> ReweightDISnucl(reweightingSigmas[kDISnucl][weight_point]);
            break;
          case kAGKY:
            // reweightVector[i_reweightingKnob][weight_point]
            //   -> ReweightAGKY(reweightingSigmas[kAGKY][weight_point]);
            break;
          case kNReWeights:
            break;
        }

      } //loop over nWeights
    } // loop over physical knobs







    // std::cout << "\n\n\nsetup finished, running"
    //           << " configure on each weight......\n\n\n";

    // // Tell all of the reweight drivers to configure themselves:
    // for(auto & vec : reweightVector){
    //   for (auto & driver : vec){
    //     driver -> Configure();
    //   }
    // }

    return;

  }

  unsigned int NuAnaAlg::prepareSigmas(int NWeights,
                           unsigned int RandSeed,
                           std::vector<std::vector<float> > & reweightingSigmas)
  {

    TRandom rand;
    rand.SetSeed(RandSeed);

    reweightingSigmas.resize(kNReWeights);
    for (unsigned int i = 0; i < reweightingSigmas.size(); ++i)
    {
      reweightingSigmas[i].resize(NWeights);
      for (int j = 0; j < NWeights; j ++)
        reweightingSigmas[i][j] = rand.Gaus(0,1);
    }
    return rand.GetSeed();
  }

  void NuAnaAlg::parseWeights(const std::vector<std::string> & string_weights,
                              std::vector<reweight> & enum_weights){

    enum_weights.reserve(string_weights.size());
    for( auto & s : string_weights){
      if (s == "NCEL") enum_weights.push_back(kNCEL);
      else if (s == "QEMA") enum_weights.push_back(kQEMA);
      else if (s == "QEVec") enum_weights.push_back(kQEVec);
      else if (s == "ResGanged") enum_weights.push_back(kResGanged);
      else if (s == "CCRes") enum_weights.push_back(kCCRes);
      else if (s == "NCRes") enum_weights.push_back(kNCRes);
      else if (s == "Coh") enum_weights.push_back(kCoh);
      else if (s == "NonResRvp1pi") enum_weights.push_back(kNonResRvp1pi);
      else if (s == "NonResRvbarp1pi") enum_weights.push_back(kNonResRvbarp1pi);
      else if (s == "NonResRvp2pi") enum_weights.push_back(kNonResRvp2pi);
      else if (s == "NonResRvbarp2pi") enum_weights.push_back(kNonResRvbarp2pi);
      else if (s == "ResDecay") enum_weights.push_back(kResDecay);
      else if (s == "NC") enum_weights.push_back(kNC);
      else if (s == "DIS") enum_weights.push_back(kDIS);
      else if (s == "DISnucl") enum_weights.push_back(kDISnucl);
      else if (s == "AGKY") enum_weights.push_back(kAGKY);
    }



  }


  void NuAnaAlg::calcWeight(simb::MCTruth const& mctruth,
                            simb::GTruth const& gtruth,
                            std::vector<std::vector<float>>& weights){
    // return reweight.CalcWeight(*mctruth,*gtruth);
    // if (weights.size() == 0) weights.resize(1);
    // weights.front().push_back( reweight -> CalcWeight(*mctruth,*gtruth) );

    // weights needs to be the size of the reweighting vector
    if (weights.size() != reweightVector.size())
      weights.resize(reweightVector.size());

    for (unsigned int i_weight = 0; i_weight < reweightVector.size(); i_weight ++){
      if (weights[i_weight].size() != reweightVector[i_weight].size()){
        weights[i_weight].resize(reweightVector[i_weight].size());
      }
      for (unsigned int i_reweightingKnob = 0;
           i_reweightingKnob < reweightVector[i_weight].size();
           i_reweightingKnob ++)
      {
        weights[i_weight][i_reweightingKnob]
          = reweightVector[i_weight][i_reweightingKnob]
            -> CalcWeight(mctruth,gtruth);
      }
    }

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
                                  std::vector<float>& neutMom,
                                  std::vector<float>& vertex){
    // Fill in all the variables that need to be filled!

    // prepare the vectors:
    neutMom.clear();
    neutMom.resize(4);
    vertex.clear();
    vertex.resize(3);

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

    vertex[0]  = neutrino->Nu().Position().X();
    vertex[1]  = neutrino->Nu().Position().Y();
    vertex[2]  = neutrino->Nu().Position().Z();
    neutMom[0] = neutrino->Nu().Momentum().E();
    neutMom[1] = neutrino->Nu().Momentum().X();
    neutMom[2] = neutrino->Nu().Momentum().Y();
    neutMom[3] = neutrino->Nu().Momentum().Z();


  }

  void NuAnaAlg::packFluxInfo(simb::MCFlux const& flux,
                              int& ptype, int& tptype, int& ndecay,
                              std::vector<float>& neutVertexInWindow,
                              std::vector<float>& ParentVertex,
                              std::vector<float>& nuParentMomAtDecay,
                              std::vector<float>& nuParentMomAtProd,
                              std::vector<float>& nuParentMomTargetExit){

    ptype  = flux.fptype;
    tptype = flux.ftptype;
    ndecay = flux.fndecay;

    neutVertexInWindow.clear();
    ParentVertex.clear();
    nuParentMomAtDecay.clear();
    nuParentMomAtProd.clear();
    nuParentMomTargetExit.clear();

    neutVertexInWindow.push_back(flux.fgenx);
    neutVertexInWindow.push_back(flux.fgeny);
    neutVertexInWindow.push_back(flux.fgenz);
    ParentVertex.push_back(flux.fvx);
    ParentVertex.push_back(flux.fvy);
    ParentVertex.push_back(flux.fvz);
    nuParentMomAtDecay.push_back(flux.fpdpx);
    nuParentMomAtDecay.push_back(flux.fpdpy);
    nuParentMomAtDecay.push_back(flux.fpdpz);
    nuParentMomAtProd.push_back(flux.fppdxdz);
    nuParentMomAtProd.push_back(flux.fppdydz);
    nuParentMomAtProd.push_back(flux.fpppz);
    nuParentMomTargetExit.push_back(flux.ftpx);
    nuParentMomTargetExit.push_back(flux.ftpy);
    nuParentMomTargetExit.push_back(flux.ftpz);

    return;
  }

  // Pack up the genie info:
  void NuAnaAlg::packGenieInfo(simb::MCTruth const&  truth,
                               std::vector<int> & GeniePDG,
                               std::vector<std::vector<float>>& GenieMomentum,
                               std::vector<std::string>& GenieProc,
                               int& NPi0FinalState,
                               int& NGamma,
                               int& NChargedPions){
    int i = 0;
    std::vector<float> tempMomentum;
    tempMomentum.resize(4);
    while( i < truth.NParticles()){
      auto part = truth.GetParticle(i);
      if (part.StatusCode() == 1){
        GeniePDG.push_back(part.PdgCode());
        GenieMomentum.push_back(tempMomentum);
        pack4Vector(part.Momentum(),GenieMomentum.back());
        GenieProc.push_back(part.Process());
        if (part.PdgCode() == 22)  NGamma ++;
        if (part.PdgCode() == 111) NPi0FinalState ++;
        if (abs(part.PdgCode()) == 211) NChargedPions ++;
      }
      i++;
    }
  }

  simb::MCParticle const* NuAnaAlg::getParticleByID(
          std::vector<simb::MCParticle> const& mclistLARG4,
          int TrackId) const
  {
    for(simb::MCParticle const& particle : mclistLARG4) {
      if ( particle.TrackId() == TrackId){
        return &particle;
      }
    }
    return nullptr;
  }

  simb::MCParticle const& NuAnaAlg::getParticleByIDStrict(
          std::vector<simb::MCParticle> const& mclistLARG4,
          int TrackId) const
  {
    if (auto particle = getParticleByID(mclistLARG4, TrackId)) {
      return *particle;
    }
    throw cet::exception("NuAnaAlg") << "Particle with track ID: " << TrackId << " not found.";
  }

  bool NuAnaAlg::isInTPC(const TVector3 & v) const{
    if (v.X() > xhigh || v.X() < xlow) return false;
    if (v.Y() > yhigh || v.Y() < ylow) return false;
    if (v.Z() > zhigh || v.Z() < zlow) return false;
    return true;
  }

  // Method to take in a photon and determine where it started converting.
  // Looks at the photon's energy at each step.
  void NuAnaAlg::GetPhotonConversionInfo(simb::MCParticle const& photon,
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
        ConversionMom = photon.Momentum(i-1); //<- This i OK since i starts at 1.
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

  void NuAnaAlg::packLarg4Info( std::vector<simb::MCParticle> const& mclarg4, int isCC,
                                int NPi0FinalState, int NGamma, int NChargedPions,
                                std::vector<std::vector<float>> & leptonPos,
                                std::vector<std::vector<float>> & leptonMom,
                                std::vector<std::vector<float>> & p1PhotonConversionPos,
                                std::vector<std::vector<float>> & p1PhotonConversionMom,
                                std::vector<std::vector<float>> & p2PhotonConversionPos,
                                std::vector<std::vector<float>> & p2PhotonConversionMom,
                                std::vector<std::vector<float>> & miscPhotonConversionPos,
                                std::vector<std::vector<float>> & miscPhotonConversionMom,
                                std::vector<std::vector<float>> & pionPos,
                                std::vector<std::vector<float>> & pionMom,
                                std::vector<std::vector<std::vector<float>>> &chargedPionPos,
                                std::vector<std::vector<std::vector<float>>> &chargedPionMom,
                                std::vector<int> & chargedPionSign)
  {

      // This function loops over the larg4 object and packs all the info
      // up into the vectors provided (which are written to ntuple)
      // We know from the genie function how many pi0 and gammas there are already
      std::vector<float> empty4Vector;
      empty4Vector.resize(4);

      // prepare the pi0 and photon vectors:
      p1PhotonConversionPos.reserve(NPi0FinalState);
      p1PhotonConversionMom.reserve(NPi0FinalState);
      p2PhotonConversionPos.reserve(NPi0FinalState);
      p2PhotonConversionMom.reserve(NPi0FinalState);
      pionPos.reserve(NPi0FinalState);
      pionMom.reserve(NPi0FinalState);
      miscPhotonConversionPos.reserve(NGamma);
      miscPhotonConversionMom.reserve(NGamma);

      // some counting variables to ensure everything is done correctly
      int nPrimaryLepton(0);  // should be == 1 at the end
      int nPrimaryGamma(0);   // should be == NGamma at the end
      int nPrimaryPi0(0);     // should be == NPi0FinalState

      chargedPionPos.resize(NChargedPions);
      chargedPionMom.resize(NChargedPions);


      // loop over the particles, extending the pion vectors as it goes.
      for (simb::MCParticle const& particle : mclarg4) {

        if (particle.Mother() == 0 ){ // then this is a primary
          // std::cout << "On particle " << particle.TrackId()
          //           << " with PDG " << particle.PdgCode() << std::endl;

          // For older files, set the nPrimaryLepton up since it didn't track neutrinos
          if (isCC == 0) nPrimaryLepton ++;

          if (abs(particle.PdgCode()) == 11 ||
              abs(particle.PdgCode()) == 12 ||
              abs(particle.PdgCode()) == 13 ||
              abs(particle.PdgCode()) == 14 )
          {

            // then this is definitely the lepton.
            nPrimaryLepton ++;
            int nTrajectoryPoints = particle.NumberTrajectoryPoints();
            leptonPos.reserve(nTrajectoryPoints);
            leptonMom.reserve(nTrajectoryPoints);
            for (int traj_point = 0; traj_point < nTrajectoryPoints; ++traj_point){
              leptonPos.push_back(empty4Vector);
              pack4Vector(particle.Position(traj_point),leptonPos.back());
              leptonMom.push_back(empty4Vector);
              pack4Vector(particle.Momentum(traj_point),leptonMom.back());
              if (!isInTPC(particle.Position(traj_point).Vect())) break;
            }
          } // end lepton if block

          if (particle.PdgCode() == 22)  //  misc gamma
          {
            nPrimaryGamma ++;
            TLorentzVector conversionPoint, conversionMom;
            GetPhotonConversionInfo(particle, conversionPoint,conversionMom);
            miscPhotonConversionPos.push_back(empty4Vector);
            miscPhotonConversionMom.push_back(empty4Vector);
            pack4Vector(conversionPoint,miscPhotonConversionPos.back());
            pack4Vector(conversionMom,  miscPhotonConversionMom.back());
          }

          if (particle.PdgCode() == 111)  // neutral pion
          {
            pionPos.push_back(empty4Vector);
            pionMom.push_back(empty4Vector);
            pack4Vector(particle.EndPosition(),pionPos.back());
            pack4Vector(particle.Momentum(),pionMom.back());
            nPrimaryPi0 ++;
            // get the conversion point of each photon
            if ( particle.NumberDaughters() == 2){
              // Get the info:
              TLorentzVector conversionPoint, conversionMom;
              simb::MCParticle const& daughter0
                  = getParticleByIDStrict(mclarg4,particle.Daughter(0));
              GetPhotonConversionInfo(daughter0, conversionPoint,conversionMom);
              // Pack it into the vectors
              p1PhotonConversionPos.push_back(empty4Vector);
              p1PhotonConversionMom.push_back(empty4Vector);
              pack4Vector(conversionPoint,p1PhotonConversionPos.back());
              pack4Vector(conversionMom,  p1PhotonConversionMom.back());

              // get photon2:
              simb::MCParticle const& daughter1
                  = getParticleByIDStrict(mclarg4,particle.Daughter(1));
              GetPhotonConversionInfo(daughter1, conversionPoint,conversionMom);
              // Pack it into the vectors:
              p2PhotonConversionPos.push_back(empty4Vector);
              p2PhotonConversionMom.push_back(empty4Vector);
              pack4Vector(conversionPoint,p2PhotonConversionPos.back());
              pack4Vector(conversionMom,  p2PhotonConversionMom.back());
            }
            else{ // this is the "dalitz" decay
              // Only take the photon, but be sure to find it
              TLorentzVector conversionPoint, conversionMom;
              p2PhotonConversionPos.push_back(empty4Vector);
              p2PhotonConversionMom.push_back(empty4Vector);
              for (int daughter = 0; daughter < 3; daughter ++ ){
                simb::MCParticle const* daughterParticle
                  = getParticleByID(mclarg4,particle.Daughter(daughter));
                if (daughterParticle and daughterParticle->PdgCode() == 22) {
                  GetPhotonConversionInfo(*daughterParticle,
                                          conversionPoint,
                                          conversionMom);
                  p1PhotonConversionPos.push_back(empty4Vector);
                  p1PhotonConversionMom.push_back(empty4Vector);
                  pack4Vector(conversionPoint,p1PhotonConversionPos.back());
                  pack4Vector(conversionMom,  p1PhotonConversionMom.back());
                  break;
                }
              }
            }

          }

          if (abs(particle.PdgCode()) == 211)  // charged pion
          {
            chargedPionSign.push_back(particle.PdgCode() / 211);
            int nTrajectoryPoints = particle.NumberTrajectoryPoints();
            unsigned int index = chargedPionSign.size() - 1;
            chargedPionPos[index].reserve(nTrajectoryPoints);
            chargedPionMom[index].reserve(nTrajectoryPoints);
            for (int traj_point = 0; traj_point < nTrajectoryPoints; ++ traj_point){
              chargedPionPos[index].push_back(empty4Vector);
              pack4Vector(particle.Position(traj_point),
                          chargedPionPos[index].back());
              chargedPionMom[index].push_back(empty4Vector);
              pack4Vector(particle.Momentum(traj_point),
                          chargedPionMom[index].back());
            }
          }
        } // end of if primary
      } // end of loop over particles

    // if (nPrimaryLepton < 1){
    //   std::cerr << "Major problems here ... nPrimaryLepton\n"
    //             << "  Should be " << 1
    //             << " but is " << nPrimaryLepton << "\n";
    //   exit(-1);
    // }
    if (nPrimaryGamma != NGamma){
      std::cerr << "Major problems here ... nPrimaryGamma\n"
                << "  Should be " << NGamma
                << " but is " << nPrimaryGamma << "\n";
      exit(-1);
    }
    if (nPrimaryPi0 != NPi0FinalState){
      std::cerr << "Major problems here ... nPrimaryPi0\n"
                << "  Should be " << NPi0FinalState
                << " but is " << nPrimaryPi0 << "\n";
      exit(-1);
    }
  } // end of packLarg4Info

  void NuAnaAlg::pack4Vector(const TLorentzVector& input,
                             std::vector<float>& output) const{
    if (output.size() != 4) output.resize(4);
    output[0] = input.E();
    output[1] = input.X();
    output[2] = input.Y();
    output[3] = input.Z();
    return;
  }
  void NuAnaAlg::pack3Vector(const TVector3& input,
                             std::vector<float>& output) const{
    if (output.size() != 3) output.resize(3);
    output[0] = input.X();
    output[1] = input.Y();
    output[2] = input.Z();
    return;
  }


#ifdef CUSTOM_NUTOOLS
    void NuAnaAlg::packFluxWeight(simb::MCFlux const& flux,
                            std::vector<std::vector<float>>& eventReweight)
    {
      // This is getting the flux weights from the flux object.
      // It's a total hack, it's hardcoded, and requires a custom version of
      // the nutools software.
      eventReweight.clear();
      eventReweight.resize(7);
      for (int i = 0; i < 7; i ++)
      {
        eventReweight[i].resize(1000);
        for (int j = 0; j < 1000; ++j)
        {
          eventReweight[i][j]=flux.eventReweight[i*1000+j];
        }
      }

      return;
    }

#else
     void NuAnaAlg::packFluxWeight(simb::MCFlux const& flux,
                             std::vector<std::vector<float>>&)
    {

      std::cerr << "You are asking to pack flux weights but the version "
                << "of nutools you use does not appear to support it.\n";
      return;
    }

#endif

} // end of namespace sbnd

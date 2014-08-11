#ifndef LAR1ND_NUANA
#define LAR1ND_NUANA value



/// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

// Larsoft includes
#include "Geometry/Geometry.h"

// nusoft includes
#include "SimulationBase/MCFlux.h"
#include "SimulationBase/MCNeutrino.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"
#include "SummaryData/POTSummary.h"

// root includes
#include <TTree.h>
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// c++ includes
#include <vector>
#include <math.h>
#include <iomanip>
#include <iostream>

// LArSoft includes                                                                                         
#include "NuAnaAlg.h"

namespace lar1nd{
  /// A module to extract the info from the Monte Carlo generator GENIE
  class NuAna : public art::EDAnalyzer {
  public:
    explicit NuAna(fhicl::ParameterSet const& pset);
    virtual ~NuAna(){};
    void beginJob();    
    void reconfigure(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& evt);
    void beginSubRun (const art::SubRun& subrun);
    void reset();

  private:

    // This alg does all of the computing, while this module interfaces with
    // art and packs all of the variables into the ttree
    NuAnaAlg fNuAnaAlg;

    //Variables needed from the .fcl file to get things out of the event
    std::string fGenieModuleLabel;
    std::string fLarg4ModuleLabel;

    // Other variables from the fcl file
    std::string fMode;      // mode is beam mode, nu or nubar
    double fBaseline;       // baseline is 100, 470, or 600 meters
    bool fFullOscTrue;      // fullosctrue is obvious.

    // #---------------------------------------------------------------
    // #This is the list of analysis variables needed for the ttree:
    // #---------------------------------------------------------------
    // These are the ttrees themselves:
    TTree* fTreeTot;    //This tree stores all the important information from an event
    TTree* PoTTree;     //This tree only stores the POT of the file, nothing else.

    std::vector<float> eventWeights;  

    TLorentzVector     neutMom;      // Neutrino Momentum
    
    //Track the muon, if there is one.
    //If there is no muon, these vectors will be size 0
    std::vector< TLorentzVector > muonPos;
    std::vector< TLorentzVector > muonMom;
    
    //Track the electron, if there is one.
    //If there is no muon, these vectors will be size 0
    std::vector< TLorentzVector > electronPos;
    std::vector< TLorentzVector > electronMom;
    
    int NPi0FinalState;     // Number of neutral pions in Final State Particles
    int NPi0;               // Number of neutral pions in larg4, which could be different
    int NGamma;             // Number of gammas in the Final State Particles
    int NChargedPions;      // Number of charged pions in the Final State Particles
    
    bool foundAllPhotons;   // True by default, goes to false ONLY if there
                            // are final state photons that the alg didn't
                            // find in the larg4 particle list
                             
    // Keep track of where photons convert and what their energy is.
    // Except in categories 1, 3, 5, these are probably not filled.
    // The vectors here labeled p1 and p2 keep track only of photons from neutral pion decay
    std::vector< TLorentzVector > p1PhotonConversionPos;  // Store the position of conversion
    std::vector< TLorentzVector > p1PhotonConversionMom;  // Store the momentum 4vector at conversion
    std::vector< TLorentzVector > p2PhotonConversionPos;  // Store the position of conversion
    std::vector< TLorentzVector > p2PhotonConversionMom;  // Store the momentum 4vector at conversion
    
    // These vectors keep track of any other photon conversions we need to look at.
    std::vector< TLorentzVector > miscPhotonConversionPos;  // Store the position of conversion
    std::vector< TLorentzVector > miscPhotonConversionMom;    // Store the momentum 4vector at conversion
    
    std::vector< std::vector < TLorentzVector > > chargedPionPos; // store position of charged pions
    std::vector< std::vector < TLorentzVector > > chargedPionMom; // store momentum of charged pions
    std::vector< int > chargePionSign;  // Sign (+/-) of the charged pions, one to one with above vectors
    
    // Track any pions
    // Only keeping the position of decay and the momentum at decay
    std::vector< TLorentzVector > pionPos;  // Position at decay
    std::vector< TLorentzVector > pionMom;  // Momentum at decay                            
    
    //Info from the genie truth:
    std::vector<int>    GeniePDG;       // Contains the pdg of the FSP in genie generated
    std::vector<TLorentzVector> GenieMomentum;
    // std::vector<double> GenieE;         // Energy of particles in genie FSP
    // std::vector<double> GeniePx;        // Contains the x momentum of the particles in genie
    // std::vector<double> GeniePy;        // Contains the y momentum of the particles in genie
    // std::vector<double> GeniePz;        // Contains the z momentum of the particles in genie   
    std::vector<std::string> GenieProc; // Contains the process information 
    

    int iflux;                      // represents the sample, 0 = nu, 1 = nu_fosc, 2 = nubar, 3 = nubar_fosc
    int ibkg;                       // A bit outdated, represents what type of background an event might be
    int nuchan;                     // Type of neutrino interaction
    int inno;                       // Type of neutrino, PDG code.  Not sure why it's call inno...
    int parid;                      // ID of the parent, for flux reweighing
    int ndecay;                     // Type of decay, for flux reweighing
    int isCC;                       // isCC event? isCC == 1 means CC, isCC == 0 means NC
    int mode;                       // beam mode
    double enugen,energy;           // Energy of the neutrino (both)
    double nuleng;                  // Length the neutrino traveled. 
    double wgt;                     // Some weighting function, not filled here (but used in ntuples)
    TVector3 vertex;                // Vertex location
    double Vdist;                   // Not sure why this is here...
    double ParVx, ParVy, ParVz;     // Parent Vertex (not in detector)
    double ParPx, ParPy, ParPz;     // Parent Momentum
    double pdpx, pdpy, pdpz;        // nu parent momentum at the time of decay
    double pppx, pppy, pppz;        // nu parent momentum at production point
    double tpx, tpy, tpz;           // momentum of parent off the target
    int    ptype, tptype;
    double LepPx, LepPy, LepPz;     // Lepton Momentum, initial
    double Elep, ElepSmeared;       // Energy of the produced lepton
    double thetaLep, phiLep, thetaLepSmeared, phiLepSmeared; //angles in the detector of the produced lepton
    double POT;                     //POT of the whole file.

    /*
    For more information on some of the neutrino parentage information, start here
    http://genie.hepforge.org/doxygen/html/classgenie_1_1flux_1_1GSimpleNtpNuMI.html
     */
    
    // #---------------------------------------------------------------
    // # End of the list of variables for the tree
    // #---------------------------------------------------------------
    

  };
}

namespace lar1nd{
  
  NuAna::NuAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {
    this -> reconfigure(pset);

    std::cout << "\nThe mode is " << fMode << std::endl;
    std::cout << "This set of events is ";
    if (fFullOscTrue) std::cout << "fullosc" << std::endl;
    else std::cout << "not fullosc" << std::endl;
    std::cout << "The baseline for this detector is " << fBaseline << "m." << std::endl << std::endl;

    return;

  }

  void NuAna::reconfigure(fhicl::ParameterSet const& pset){
    fMode             = pset.get< std::string > ("Mode");
    fFullOscTrue      = pset.get< bool        > ("FullOsc");
    fBaseline         = pset.get< double      > ("Baseline");
    fGenieModuleLabel = pset.get< std::string > ("GenieModuleLabel");
    fLarg4ModuleLabel = pset.get< std::string > ("LArG4ModuleLabel");
    return;
  }

  void NuAna::beginJob(){

    // This function sets up the ttrees
    // get access to the TFile service  
    art::ServiceHandle<art::TFileService> tfs;

    PoTTree  = tfs->make<TTree>("POT", "POT");
    PoTTree ->Branch("POT",&POT,"POT/D");

    fTreeTot = tfs->make<TTree>("EventsTot", "Event info for ALL types");

    fTreeTot->Branch("iflux",          &iflux,          "iflux/I");
    fTreeTot->Branch("ibkg",           &ibkg,           "ibkg/I");
    fTreeTot->Branch("nuchan",         &nuchan,         "nuchan/I");
    fTreeTot->Branch("inno",           &inno,           "inno/I");
    fTreeTot->Branch("enugen",         &enugen,         "enugen/D");
    fTreeTot->Branch("energy",         &energy,         "energy/D");
    fTreeTot->Branch("nuleng",         &nuleng,         "nuleng/D");
    fTreeTot->Branch("parid",          &parid,          "parid/I");
    fTreeTot->Branch("ndecay",         &ndecay,         "ndecay/I");
    fTreeTot->Branch("wgt",            &wgt,            "wgt/D");
    fTreeTot->Branch("NPi0",           &NPi0,           "NPi0/I");
    fTreeTot->Branch("NPi0FinalState", &NPi0FinalState, "NPi0FinalState/I");
    fTreeTot->Branch("NGamma",         &NGamma,         "NGamma/I");
    fTreeTot->Branch("FoundPhotons",   &foundAllPhotons,"FoundAllPhotons/B");
    fTreeTot->Branch("GeniePDG",       &GeniePDG);
    fTreeTot->Branch("GenieMomentum",  &GenieMomentum);  
    fTreeTot->Branch("GenieProc",      &GenieProc);

    fTreeTot->Branch("isCC",  &isCC,"isCC/I");
    fTreeTot->Branch("mode",  &mode,"mode/I");
    fTreeTot->Branch("vertex",&vertex);
    fTreeTot->Branch("ParVx", &ParVx,"ParVx/D");
    fTreeTot->Branch("ParVy", &ParVy,"ParVy/D");
    fTreeTot->Branch("ParVz", &ParVz,"ParVz/D");
    fTreeTot->Branch("ParPx", &ParPx,"ParPx/D");
    fTreeTot->Branch("ParPy", &ParPy,"ParPy/D");
    fTreeTot->Branch("ParPz", &ParPz,"ParPz/D");
    fTreeTot->Branch("LepPx", &LepPx,"LepPx/D");
    fTreeTot->Branch("LepPy", &LepPy,"LepPy/D");
    fTreeTot->Branch("LepPz", &LepPz,"LepPz/D");

    fTreeTot->Branch("pdpx",  &pdpx,"pdpx/D");      // nu parent momentum at the time of decay
    fTreeTot->Branch("pdpy",  &pdpy,"pdpy/D");
    fTreeTot->Branch("pdpz",  &pdpz,"pdpz/D");        
    fTreeTot->Branch("pppx",  &pppx,"pppx/D");      // nu parent momentum at production point
    fTreeTot->Branch("pppy",  &pppy,"pppy/D");
    fTreeTot->Branch("pppz",  &pppz,"pppz/D");        
    fTreeTot->Branch("tpx",   &tpx,"tpx/D");         // momentum of parent off the target
    fTreeTot->Branch("tpy",   &tpy,"tpy/D");
    fTreeTot->Branch("tpx",   &tpx,"tpx/D");           
    
    
    fTreeTot->Branch("ptype", &ptype,"ptype/I");
    fTreeTot->Branch("tptype",&tptype,"tptype/I");
    fTreeTot->Branch("neutMom", "neutMom", &neutMom, 32000, 0);


    fTreeTot->Branch("ThetaLep", &thetaLep, "ThetaLep/D");
    fTreeTot->Branch("PhiLep", &phiLep, "PhiLep/D");
    fTreeTot->Branch("ThetaLepSmeared", &thetaLepSmeared, "ThetaLepSmeared/D");
    fTreeTot->Branch("PhiLepSmeared", &phiLepSmeared, "PhiLepSmeared/D");
    fTreeTot->Branch("Elep", &Elep, "Elep/D");
    fTreeTot->Branch("ElepSmeared", &ElepSmeared, "ElepSmeared/D");
    
    fTreeTot->Branch("MuonPos","MuonPos", &muonPos, 32000, 0);
    fTreeTot->Branch("MuonMom","MuonMom", &muonMom, 32000, 0);
    fTreeTot->Branch("ElectronPos","ElectronPos", &electronPos, 32000, 0);
    fTreeTot->Branch("ElectronMom","ElectronMom", &electronMom, 32000, 0);
    fTreeTot->Branch("p1PhotonConversionPos", "p1PhotonConversionPos", &p1PhotonConversionPos, 32000, 0);
    fTreeTot->Branch("p1PhotonConversionMom","p1PhotonConversionMom", &p1PhotonConversionMom, 32000, 0);
    fTreeTot->Branch("p2PhotonConversionPos","p2PhotonConversionPos", &p2PhotonConversionPos, 32000, 0);
    fTreeTot->Branch("p2PhotonConversionMom","p2PhotonConversionMom", &p2PhotonConversionMom, 32000, 0);
    fTreeTot->Branch("miscPhotonConversionPos","miscPhotonConversionPos", &miscPhotonConversionPos, 32000, 0);
    fTreeTot->Branch("miscPhotonConversionMom","miscPhotonConversionMom", &miscPhotonConversionMom, 32000, 0);
    fTreeTot->Branch("PionPos","PionPos", &pionPos, 32000, 0);
    fTreeTot->Branch("PionMom","PionMom", &pionMom, 32000, 0);
    fTreeTot->Branch("ChargedPionPos","ChargedPionPos",&chargedPionPos, 32000,0);
    fTreeTot->Branch("ChargedPionMom","ChargedPionMom",&chargedPionMom, 32000,0);
    fTreeTot->Branch("ChargePionSign","ChargePionSign",&chargePionSign, 32000,0);
    

    // fTreeTot->Branch("MultiWeight","MultiWeight",&eventReweight,32000,0);
    art::ServiceHandle<geo::Geometry> geom;
    // configure the geometry in the worker function:
    fNuAnaAlg.configureGeometry(geom);

    return;
  }

  void NuAna::analyze(const art::Event& evt){
    return;
  }
  void NuAna::beginSubRun(const art::SubRun& subrun){
    return;
  }
  void NuAna::reset(){
    return;
  }
}


namespace lar1nd {

  // A macro required for a JobControl module.                                                              
  DEFINE_ART_MODULE(NuAna)

} // namespace lar1nd                                                                                       

#endif
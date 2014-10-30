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
#include "SimulationBase/GTruth.h"
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
#include "TROOT.h"

// c++ includes
#include <vector>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>

// LArSoft includes                                                                                         
#include "alg/NuAnaAlg.h"




namespace lar1nd{
  /// A module to extract the info from the Monte Carlo generator GENIE
  class NuAna : public art::EDAnalyzer {
  public:
    explicit NuAna(fhicl::ParameterSet const& pset);
    virtual ~NuAna(){};

    void beginJob();    
    // void reconfigure(fhicl::ParameterSet const& pset);
    void analyze(const art::Event& evt);
    void beginSubRun (const art::SubRun& subrun);
    void reset();

  private:

    // This alg does all of the computing, while this module interfaces with
    // art and packs all of the variables into the ttree
    NuAnaAlg fNuAnaAlg;



    // Other variables from the fcl file
    std::string fMode;      // mode is beam mode, nu or nubar
    double fBaseline;       // baseline is 100, 470, or 600 meters
    bool fFullOscTrue;      // fullosctrue is obvious.
    
    //Variables needed from the .fcl file to get things out of the event
    std::string fGenieModuleLabel;
    std::string fLarg4ModuleLabel;


    bool fFluxReweight;

    bool fXSecReweight;
    std::vector<std::string> fWeights;
    std::vector<float> fWeightRangeSigma;
    unsigned int fRandSeed;
    int fNWeights;

    unsigned int FinalRandSeed;

    std::vector< std::vector<float> > reweightingSigmas;

    // #---------------------------------------------------------------
    // #This is the list of analysis variables needed for the ttree:
    // #---------------------------------------------------------------
    // These are the ttrees themselves:
    TTree* fTreeTot;    //This tree stores all the important information from an event
    // TTree* PoTTree;     //This tree only stores the POT of the file, nothing else.

    std::vector<std::vector<float> > eventWeights;  

    std::vector<float>     neutMom;      // Neutrino Momentum
    
    // Information about the lepton:
    // muons are tracked till they exit, electrons are just the initial
    // particle's position and momemtum (since they are showering)
    std::vector< std::vector<float> > leptonPos; // pos of lepton till it exits
    std::vector< std::vector<float> > leptonMom; // mom of lepton till it exits
    double Elep;                // Energy of the produced lepton
    double thetaLep, phiLep;    // angles in the detector of the produced lepton
    
    
    int NPi0FinalState;     // Number of neutral pions in Final State Particles
    int NGamma;             // Number of gammas in the Final State Particles
    int NChargedPions;      // Number of charged pions in the Final State Particles
    
    bool foundAllPhotons;   // True by default, goes to false ONLY if there
                            // are final state photons that the alg didn't
                            // find in the larg4 particle list
                             
    // Keep track of where photons convert and what their energy is.
    // Except in categories 1, 3, 5, these are probably not filled.
    // The vectors here labeled p1 and p2 keep track only of photons from neutral pion decay
    std::vector< std::vector<float> > p1PhotonConversionPos;  // Store the position of conversion
    std::vector< std::vector<float> > p1PhotonConversionMom;  // Store the momentum 4vector at conversion
    std::vector< std::vector<float> > p2PhotonConversionPos;  // Store the position of conversion
    std::vector< std::vector<float> > p2PhotonConversionMom;  // Store the momentum 4vector at conversion
    
    // These vectors keep track of any other photon conversions we need to look at.
    std::vector< std::vector<float> > miscPhotonConversionPos;  // Store the position of conversion
    std::vector< std::vector<float> > miscPhotonConversionMom;    // Store the momentum 4vector at conversion
    
    std::vector< std::vector < std::vector<float> > > chargedPionPos; // store position of charged pions
    std::vector< std::vector < std::vector<float> > > chargedPionMom; // store momentum of charged pions
    std::vector< int > chargedPionSign;  // Sign (+/-) of the charged pions, one to one with above vectors
    
    // Track any pions
    // Only keeping the position of decay and the momentum at decay
    std::vector< std::vector<float> > pionPos;  // Position at decay
    std::vector< std::vector<float> > pionMom;  // Momentum at decay                            
    
    //Info from the genie truth:
    std::vector<int>    GeniePDG;       // Contains the pdg of the FSP in genie generated
    std::vector<std::vector<float>> GenieMomentum;
    std::vector<std::string> GenieProc; // Contains the process information 
    // Genie Reweight vectors:
    std::vector< std::vector< float > > genieReweights;

    int iflux;                      // represents the sample, 0 = nu, 1 = nu_fosc, 2 = nubar, 3 = nubar_fosc
    int nuchan;                     // Type of neutrino interaction
    int inno;                       // Type of neutrino, PDG code.  Not sure why it's call inno...
    int parid;                      // ID of the parent, for flux reweighing
    int ndecay;                     // Type of decay, for flux reweighing
    int isCC;                       // isCC event? isCC == 1 means CC, isCC == 0 means NC
    int mode;                       // beam mode
    double enugen;                  // Energy of the neutrino
    double nuleng;                  // Length the neutrino traveled. 
    std::vector<float> vertex;                // Vertex location
    std::vector<float> neutVertexInWindow;    // origin of neutrino in flux window
    std::vector<float> ParentVertex;          // Parent Vertex (not in detector)
    std::vector<float> nuParentMomAtDecay;    // nu parent momentum at the time of decay
    std::vector<float> nuParentMomAtProd;     // nu parent momentum at production point
    std::vector<float> nuParentMomTargetExit; // parent particle moment at target exit
    int    ptype, tptype;
    double POT;                     // POT of the whole file.

    /*
    For more information on some of the neutrino parentage information, start here
    http://genie.hepforge.org/doxygen/html/classgenie_1_1flux_1_1GSimpleNtpNuMI.html
     */
    
    // #---------------------------------------------------------------
    // # End of the list of variables for the tree
    // #---------------------------------------------------------------
    

  }; // end of class NuAna
// }

// namespace lar1nd{
  
  NuAna::NuAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    , fMode             (pset.get< std::string >              ("Mode"))
    , fBaseline         (pset.get< double >                   ("Baseline"))
    , fFullOscTrue      (pset.get< bool >                     ("FullOsc"))
    , fGenieModuleLabel (pset.get< std::string >              ("GenieModuleLabel"))
    , fLarg4ModuleLabel (pset.get< std::string >              ("LArG4ModuleLabel"))
    , fFluxReweight     (pset.get< bool >                     ("FluxReweight"))
    , fXSecReweight     (pset.get< bool >                     ("XSecReweight"))
    , fWeights          (pset.get< std::vector<std::string> > ("Weights"))
    , fRandSeed         (pset.get< unsigned int >             ("RandSeed"))
    , fNWeights         (pset.get< int >                      ("NWeights"))
  {

    // This function sets up the ttrees
    // get access to the TFile service  
    art::ServiceHandle<art::TFileService> tfs;

    // PoTTree  = tfs->make<TTree>("POT", "POT");

    std::ofstream libraryFile;
    libraryFile.open("loadLibs.C");
    libraryFile << "#include <vector>\n";
    libraryFile << "#ifdef __MAKECINT__\n";
    libraryFile << "#pragma link C++ class std::vector<std::vector<float> >+;\n";
    libraryFile << "#pragma link C++ class std::vector<std::vector<std::vector<float> > >+;\n";
    libraryFile << "#endif\n";
    libraryFile.close();
    gROOT -> ProcessLine(".L loadLibs.C+");

    fTreeTot = tfs->make<TTree>("EventsTot", "Event info for ALL types");
    
    fTreeTot->Branch("POT",&POT,"POT/D");
    
    // Neutrino/event variables:
    fTreeTot->Branch("iflux",    &iflux,    "iflux/I");
    fTreeTot->Branch("nuchan",   &nuchan,   "nuchan/I");
    fTreeTot->Branch("inno",     &inno,     "inno/I");
    fTreeTot->Branch("enugen",   &enugen,   "enugen/D");
    fTreeTot->Branch("nuleng",   &nuleng,   "nuleng/D");
    fTreeTot->Branch("isCC",     &isCC,     "isCC/I");
    fTreeTot->Branch("mode",     &mode,     "mode/I");
    fTreeTot->Branch("ThetaLep", &thetaLep, "ThetaLep/D");
    fTreeTot->Branch("PhiLep",   &phiLep,   "PhiLep/D");
    fTreeTot->Branch("Elep",     &Elep,     "Elep/D");
    fTreeTot->Branch("neutMom", "neutMom", &neutMom, 32000, 0);
    fTreeTot->Branch("vertex",  "vertex",  &vertex, 32000, 0);

    // Genie Variables
    fTreeTot->Branch("GeniePDG",       &GeniePDG);
    fTreeTot->Branch("GenieMomentum", "std::vector<std::vector<float> >",  &GenieMomentum);  
    fTreeTot->Branch("GenieProc",      &GenieProc);
          // "var","std::vector<std::vector<float> >",&vec
    
    // Flux variables:
    fTreeTot->Branch("ptype", &ptype, "ptype/I");
    fTreeTot->Branch("tptype",&tptype,"tptype/I");
    fTreeTot->Branch("ndecay",&ndecay,"ndecay/I");
    fTreeTot->Branch("neutVertexInWindow", "neutVertexInWindow", &neutVertexInWindow, 32000, 0);
    fTreeTot->Branch("ParentVertex", "ParentVertex", &ParentVertex, 32000, 0);
    fTreeTot->Branch("nuParentMomAtDecay", "nuParentMomAtDecay", &nuParentMomAtDecay, 32000, 0);
    fTreeTot->Branch("nuParentMomAtProd",  "nuParentMomAtProd",  &nuParentMomAtProd, 32000, 0);
    fTreeTot->Branch("nuParentMomTargetExit", "nuParentMomTargetExit", &nuParentMomTargetExit, 32000, 0);

    // larg4 info
    fTreeTot->Branch("leptonPos","std::vector< std::vector<float> > ", &leptonPos, 32000, 0);
    fTreeTot->Branch("leptonMom","std::vector< std::vector<float> > ", &leptonMom, 32000, 0);
    fTreeTot->Branch("NPi0FinalState", &NPi0FinalState, "NPi0FinalState/I");
    fTreeTot->Branch("NGamma",         &NGamma,         "NGamma/I");
    fTreeTot->Branch("FoundPhotons",   &foundAllPhotons,"FoundAllPhotons/B");
    fTreeTot->Branch("p1PhotonConversionPos","std::vector< std::vector<float> > ", &p1PhotonConversionPos, 32000, 0);
    fTreeTot->Branch("p1PhotonConversionMom","std::vector< std::vector<float> > ", &p1PhotonConversionMom, 32000, 0);
    fTreeTot->Branch("p2PhotonConversionPos","std::vector< std::vector<float> > ", &p2PhotonConversionPos, 32000, 0);
    fTreeTot->Branch("p2PhotonConversionMom","std::vector< std::vector<float> > ", &p2PhotonConversionMom, 32000, 0);
    fTreeTot->Branch("miscPhotonConversionPos","std::vector< std::vector<float> > ", &miscPhotonConversionPos, 32000, 0);
    fTreeTot->Branch("miscPhotonConversionMom","std::vector< std::vector<float> > ", &miscPhotonConversionMom, 32000, 0);
    fTreeTot->Branch("PionPos","std::vector< std::vector<float> >", &pionPos, 32000, 0);
    fTreeTot->Branch("PionMom","std::vector< std::vector<float> >", &pionMom, 32000, 0);
    fTreeTot->Branch("ChargedPionPos","std::vector<std::vector< std::vector<float> > >",&chargedPionPos, 32000,0);
    fTreeTot->Branch("ChargedPionMom","std::vector<std::vector< std::vector<float> > >",&chargedPionMom, 32000,0);
    fTreeTot->Branch("ChargedPionSign","ChargedPionSign",&chargedPionSign, 32000,0);
    

    if (fXSecReweight){
      fTreeTot->Branch("MultiWeight","MultiWeight",&eventWeights,32000,0);
      fTreeTot->Branch("FinalRandSeed", &FinalRandSeed, "FinalRandSeed/I");
    }
    if (fFluxReweight){
      fTreeTot->Branch("MultiWeight","MultiWeight",&eventWeights,32000,0);
    }

  }


  void NuAna::beginJob(){

    std::cout << "\nThe mode is " << fMode << std::endl;
    std::cout << "This set of events is ";
    if (fFullOscTrue) std::cout << "fullosc" << std::endl;
    else std::cout << "not fullosc" << std::endl;
    std::cout << "The baseline for this detector is " 
              << fBaseline << "m." << std::endl << std::endl;

    reset();
    // gROOT->ProcessLine(".L loadDictionaries.C+");

    art::ServiceHandle<geo::Geometry> geom;
    // configure the geometry in the worker function:
    fNuAnaAlg.configureGeometry(geom);

    if (fXSecReweight){
      std::vector<reweight> reweights;
      fNuAnaAlg.parseWeights(fWeights, reweights);
      FinalRandSeed = fNuAnaAlg.prepareSigmas(fNWeights,
                                              fRandSeed, reweightingSigmas);
      fNuAnaAlg.configureReWeight(reweights, reweightingSigmas);
    }

    return;
  }


  void NuAna::beginSubRun(const art::SubRun& subrun){
    
    // Go through POTSummary objects 
    art::Handle< sumdata::POTSummary > potHandle;
    subrun.getByLabel(fGenieModuleLabel, potHandle);
    const sumdata::POTSummary& potSum = (*potHandle);
    POT = potSum.totpot;
    // PoTTree->Fill();
    return;
  }
  void NuAna::reset(){
    
    // This method takes ANYTHING that goes into the ntuple and sets it to default.
    
    // Start by making sure the appropriate vectors are cleared and emptied:
    eventWeights.clear();
    leptonPos.clear();
    leptonMom.clear();
    p1PhotonConversionPos.clear();
    p1PhotonConversionMom.clear();
    p2PhotonConversionPos.clear();
    p2PhotonConversionMom.clear();
    miscPhotonConversionPos.clear();
    miscPhotonConversionMom.clear();
    pionPos.clear();
    pionMom.clear();
    chargedPionPos.clear();
    chargedPionMom.clear();
    chargedPionSign.clear();
    
    GeniePDG.clear();
    GenieMomentum.clear();
    GenieProc.clear();
    
    NPi0FinalState  = 0;
    NChargedPions   = 0;
    NGamma          = 0;
    foundAllPhotons = true;

    ptype  = 0;
    tptype = 0;

    // POT = 0;

    // Clear out the root vectors
    vertex.clear();
    ParentVertex.clear();
    nuParentMomAtDecay.clear();
    nuParentMomAtProd.clear();
    nuParentMomTargetExit.clear();
    neutMom.clear();

    // iflux    = -999;             // represents the sample, 0 = nu, 1 = nu_fosc, 2 = nubar, 3 = nubar_fosc
    nuchan   = -999;             // Type of neutrino interaction
    inno     = -999;             // Type of neutrino, PDG code.  Not sure why it's inno...
    parid    = -999;             // ID of the parent, for flux reweighing
    ndecay   = -999;             // Type of decay, for flux reweighing
    isCC     = -999;             // isCC event? isCC == 1 means CC, isCC == 0 means NC
    mode     = -999;             // beam mode
    enugen   = -999;             // Energy of the neutrino (both)
    nuleng   = -999;             // Length the neutrino traveled.    
    Elep     = -999;
    thetaLep = -999;
    phiLep   = -999;
    // POT      = -999;             //POT of the whole file.

    // eventReweight.clear();
    // eventReweight.resize(7);
    // for (auto & weight : eventReweight)
    // {
    //   weight.clear();
    //   weight.resize(1000); 
    // }

    return;
  }

  void NuAna::analyze(const art::Event& evt){

    this -> reset();
    
    //get the MC generator information out of the event       
    //these are all handles to mc information.
    art::Handle< std::vector<simb::MCTruth> > mclistGENIE;  
    art::Handle< std::vector<simb::MCParticle> > mclistLARG4;
    art::Handle< std::vector<simb::MCFlux> > mcflux;
    art::Handle< std::vector<simb::GTruth> > mcgtruth;

    // Start by getting the POT for this neutrino's file:
    // art::Handle< sumdata::POTSummary > potHandle;
    // evt.getByLabel(fGenieModuleLabel, potHandle);
    // POT = potHandle->totpot;


    //actually go and get the stuff
    evt.getByLabel(fGenieModuleLabel,mclistGENIE);
    evt.getByLabel(fGenieModuleLabel,mcflux);
    evt.getByLabel(fGenieModuleLabel,mcgtruth);
    if (!fFullOscTrue) 
        evt.getByLabel(fLarg4ModuleLabel,mclistLARG4);

    // contains the mctruth object from genie
    art::Ptr<simb::MCTruth> mc(mclistGENIE,0);
    
    // contains the mcflux object
    art::Ptr<simb::MCFlux > flux(mcflux,0);

    // contains the gtruth object
    art::Ptr<simb::GTruth > gtruth(mcgtruth,0);



    // Contains the neutrino info
    simb::MCNeutrino neutrino = mc -> GetNeutrino();
    
    // std::cout << "The interaction info is: \n" 
    //           << "  gtruth->ftgtPDG................." << gtruth->ftgtPDG << "\n"
    //           << "  gtruth->ftgtZ..................." << gtruth->ftgtZ << "\n"
    //           << "  gtruth->ftgtA..................." << gtruth->ftgtA << "\n"
    //           << "  gtruth->fGint..................." << gtruth->fGint << "\n"
    //           << "  gtruth->fGscatter..............." << gtruth->fGscatter << "\n"
    //           << "  gtruth->fweight................." << gtruth->fweight << "\n"
    //           << "  neutrino.Mode()................." << neutrino.Mode() << "\n"
    //           << "  neutrino.InteractionType()......" << neutrino.InteractionType() << "\n"
    //           << "  neutrino.CCNC()................." << neutrino.CCNC() << "\n"
    //           << "  neutrino.Target()..............." << neutrino.Target() << "\n"
    //           << std::endl;

    // Now start packing up the variables to fill the tree
    // In general, have the algorithm do this:
    

    // get the basic neutrino info:
    fNuAnaAlg.packNeutrinoInfo( &neutrino,
                                nuchan,
                                inno,
                                enugen,
                                isCC,
                                mode,
                                thetaLep,
                                phiLep,
                                Elep,
                                neutMom,
                                vertex);


    // Pack up the flux info:
    fNuAnaAlg.packFluxInfo(     flux, 
                                ptype, tptype, ndecay,
                                neutVertexInWindow,  
                                ParentVertex,
                                nuParentMomAtDecay,
                                nuParentMomAtProd,
                                nuParentMomTargetExit);

    // nuleng needs special attention:
    TVector3 length(vertex[0] - ParentVertex[0],
                    vertex[1] - ParentVertex[1],
                    100*fBaseline + vertex[2] - ParentVertex[2]);
    nuleng = length.Mag();

    // Pack up the genie info:
    fNuAnaAlg.packGenieInfo(    mc,
                                GeniePDG,
                                GenieMomentum,
                                GenieProc,
                                NPi0FinalState,
                                NGamma,
                                NChargedPions);
    

    // pack up the larg4 photon info:
    if(!fFullOscTrue)
        fNuAnaAlg.packLarg4Info(mclistLARG4, isCC, NPi0FinalState, 
                                NGamma, NChargedPions,
                                leptonPos,
                                leptonMom,
                                p1PhotonConversionPos,
                                p1PhotonConversionMom,
                                p2PhotonConversionPos,
                                p2PhotonConversionMom,
                                miscPhotonConversionPos,
                                miscPhotonConversionMom,
                                pionPos,
                                pionMom,
                                chargedPionPos,
                                chargedPionMom,
                                chargedPionSign);




    

    // If needed, set up the weights.
    
    if (fFluxReweight){
        fNuAnaAlg.packFluxWeight(flux, eventWeights);
    }
    
    // Find a new weight for this event:
    // std::vector<std::vector<float> > weights;
    if (fXSecReweight){
      fNuAnaAlg.calcWeight(mc, gtruth,eventWeights);
      // std::cout << "Printing weights for this event:\n";
      // for (auto s : fWeights) std::cout << s << "\t";
      //    std::cout << "Total\n";
      //    for (unsigned int row = 0; row < eventWeights.front().size(); row ++){
      //      for (unsigned int column = 0; column < eventWeights.size(); column ++){
      //        std::cout << eventWeights[column][row] << "\t";
      //      }
      //    std::cout << std::endl;
      //    }

    }

    // Fill the ttree:
    fTreeTot->Fill();


    return;
  }

  DEFINE_ART_MODULE(NuAna)


}

#endif
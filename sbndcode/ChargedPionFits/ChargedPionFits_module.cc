////////////////////////////////////////////////////////////////////////
// Class:       ChargedPionFits
// Plugin Type: analyzer (art v2_11_03)
// File:        ChargedPionFits_module.cc
//
// Generated at Tue Sep 25 15:08:32 2018 by Rhiannon Jones using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/Atom.h"

#include <sstream>
#include <cmath>
#include <ctime>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <numeric>
#include <algorithm>

#include "TROOT.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"

#include <typeinfo>

namespace pion {
  class ChargedPionFits;
}


class pion::ChargedPionFits : public art::EDAnalyzer {
public:
  explicit ChargedPionFits(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ChargedPionFits(ChargedPionFits const &) = delete;
  ChargedPionFits(ChargedPionFits &&) = delete;
  ChargedPionFits & operator = (ChargedPionFits const &) = delete;
  ChargedPionFits & operator = (ChargedPionFits &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p);

private:

  // Declare member data here.
  float m_detectorLengthX;
  float m_detectorLengthY;
  float m_detectorLengthZ;
  float m_coordinateOffsetX;
  float m_coordinateOffsetY;
  float m_coordinateOffsetZ;
  float m_selectedBorderX;
  float m_selectedBorderY;
  float m_selectedBorderZ;

  // Handle labels 
  std::string m_generator_label;
  std::string m_geant_label;
 
  // TTrees to write information to
  TTree *mcparticle_tree;

  // MCParticle variables
  int p_pdgcode;
  double p_momentum[3], norm_xz[3], norm_yz[3];
  double p_thetaxz, p_thetayz, p_momentum_mag;
};


pion::ChargedPionFits::ChargedPionFits(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{
  this->reconfigure(p);
}

void pion::ChargedPionFits::analyze(art::Event const & e)
{
  // Get the MCTruth handle
  art::Handle< std::vector< simb::MCTruth > > mct_handle;
  e.getByLabel(m_generator_label, mct_handle);
  int mct_size = mct_handle->size();

  // Loop over the truth, make sure the neutrino is from the beam and 
  // get the associated MCParticles
  if(mct_size != 1) return;
  
  if(mct_handle.isValid()){
    
    art::Ptr< simb::MCTruth > mct(mct_handle, 0);
    simb::MCNeutrino nu = mct->GetNeutrino();

    art::FindManyP< simb::MCParticle > fmcp( mct_handle, e, m_geant_label );
    std::vector< art::Ptr<simb::MCParticle> > mcp_assn = fmcp.at(0);
    
    for(art::Ptr<simb::MCParticle> part : mcp_assn) {
      // Primary charged pions only 
      if(part->Process() == "primary" && (part->PdgCode() == 211 || part->PdgCode() == -211)){

        // Get the pion momentum and the angles in the xz and yz planes
        p_pdgcode             = part->PdgCode();
        p_momentum[0]         = part->Px();
        p_momentum[1]         = part->Py();
        p_momentum[2]         = part->Pz();

        TVector3 dir;
        double mom   = sqrt(pow(p_momentum[0], 2) + pow(p_momentum[1], 2) + pow(p_momentum[2], 2));

        dir[0] = p_momentum[0] * (1. / mom);
        dir[1] = p_momentum[1] * (1. / mom);
        dir[2] = p_momentum[2] * (1. / mom);

        p_momentum_mag = mom;
        p_thetaxz = (180./ M_PI) * std::atan2(dir[0], dir[2]);
        p_thetayz = (180./ M_PI) * std::atan2(dir[1], dir[2]);
        mcparticle_tree->Fill();
      }
    }
  }
}

void pion::ChargedPionFits::beginJob()
{
  // The MCParticle tree
  mcparticle_tree = new TTree("particle_tree", "MCParticle tree: True SBND initial and final state topology information");
  
  mcparticle_tree->Branch("p_pdgcode",      &p_pdgcode,      "p_pdgcode/I");
  mcparticle_tree->Branch("p_thetaxz",      &p_thetaxz,      "p_thetaxz/D");
  mcparticle_tree->Branch("p_thetayz",      &p_thetayz,      "p_thetayz/D");
  mcparticle_tree->Branch("p_momentum",     &p_momentum,     "p_momentum[3]/D");
  mcparticle_tree->Branch("p_momentum_mag", &p_momentum_mag, "p_momentum_mag/D");
  
  mcparticle_tree->SetDirectory(0);

  norm_xz[0] = 0;
  norm_xz[1] = 1;
  norm_xz[2] = 0;
  norm_yz[0] = 1;
  norm_yz[1] = 0;
  norm_yz[2] = 0;
}

void pion::ChargedPionFits::endJob()
{
  std::cout << " Entries : " << mcparticle_tree->GetEntries() << std::endl;

  TFile file("result_output_file.root", "RECREATE");
  mcparticle_tree->Write();
  file.Write();
  file.Close();
  delete mcparticle_tree;
}

void pion::ChargedPionFits::reconfigure(fhicl::ParameterSet const & p)
{
  // Geometry
  m_detectorLengthX     = p.get<float>("DetectorLengthX");
  m_detectorLengthY     = p.get<float>("DetectorLengthY");
  m_detectorLengthZ     = p.get<float>("DetectorLengthZ");
  m_coordinateOffsetX   = p.get<float>("CoordinateOffsetX");
  m_coordinateOffsetY   = p.get<float>("CoordinateOffsetY");
  m_coordinateOffsetZ   = p.get<float>("CoordinateOffsetZ");
  m_selectedBorderX     = p.get<float>("SelectedBorderX");
  m_selectedBorderY     = p.get<float>("SelectedBorderY");
  m_selectedBorderZ     = p.get<float>("SelectedBorderZ");

  // Handle labels 
  m_generator_label              = p.get<std::string>("TruthLabel");
  m_geant_label                  = p.get<std::string>("G4Label");

}

DEFINE_ART_MODULE(pion::ChargedPionFits)

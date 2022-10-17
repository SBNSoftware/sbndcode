#include <algorithm>    // std::min_element, std::max_element
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <math.h>
#include <climits>

#include "IO/InputManager.hh"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TPad.h"
#include "TTree.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TLine.h"
#include "TMath.h"

int main(int argc, char** argv) {


  // --- Setting up the command line flags --- ///
  int opt;
  extern char *optarg;
  extern int   optopt;

  int         nEvent      = -1;
  std::string InFileName  = "";
  std::string InTreeName  = "";
  std::string OutFileName = "";

  while ( (opt = getopt(argc, argv, "i:t:o:n")) != -1) {
    switch (opt) {
    case 'i':
      InFileName = optarg;
      break;
    case 't':
      InTreeName = optarg;
      break;
    case 'o':
      OutFileName = optarg;
      break;
    case 'n':
      nEvent = atoi(optarg);
      break;
    case '?':
      std::cerr << "Unknown option: " << char(optopt) << "!" << std::endl;
      break;
    }
  }

  // --- Reading in the TTree --- //
  FullGeoAnaInputManager *im = new FullGeoAnaInputManager();
  im->SetInputFile(InFileName.c_str());
  im->SetInputTree(InTreeName.c_str());
  im->LoadTree    (          );

  // --- Grabbing the number of events to run over --- //
  int fNEvent = 0;
  if (nEvent == -1) fNEvent = im->GetEntries();
  else              fNEvent = std::min(nEvent, im->GetEntries());

  // --- Defining the plots we're going to make --- //
  TH1D *Angle = new TH1D("Angle", "Angle", 20, 0, 180);

  // --- Stepping through each event --- //
  for (int CurrentEvent=0; CurrentEvent<fNEvent; ++CurrentEvent) {
    im->GetEntry(CurrentEvent);
    std::cout << CurrentEvent << std::endl;
     TVector3 start_vertex  ;
     TVector3 muon_end_vtx  ;
     TVector3 proton_end_vtx;
 
     // Get the proton and the muon
     for (size_t i=0; i<im->True_Bck_Mother->size(); ++i) {
       if ((*im->True_Bck_Mother)[i] == 0) {
 	if ((*im->True_Bck_PDG)[i] == 2212) {
 
 	  TVector3 source((*im->True_Bck_VertX)[i],
 			 (*im->True_Bck_VertY)[i],
 			 (*im->True_Bck_VertZ)[i]);
 
 	  TVector3 p_end((*im->True_Bck_EndX)[i],
 			(*im->True_Bck_EndY)[i],
 			(*im->True_Bck_EndZ)[i]);
 
 	  start_vertex   = source;
 	  proton_end_vtx = p_end;
 	}
       }
 
       if ((*im->True_Bck_Mother)[i] == 0) {
 	if ((*im->True_Bck_PDG)[i] == 13) {
 	  TVector3 m_end((*im->True_Bck_EndX)[i],
 			(*im->True_Bck_EndY)[i],
 			(*im->True_Bck_EndZ)[i]);
 	  muon_end_vtx = m_end;
 	}
       }      
     }
 
     double angle = (muon_end_vtx - start_vertex).Angle(proton_end_vtx - start_vertex) * (180.0 / TMath::Pi());
     Angle->Fill(angle);

  }

  // --- Setting up the canvas for the output --- //
  TCanvas *c = new TCanvas("c", "c");
  c->Print((OutFileName + ".pdf[").c_str());

  Angle->Draw();
  c->Print((OutFileName + ".pdf").c_str());
  c->Print((OutFileName + ".pdf]").c_str());


  return 0;

}

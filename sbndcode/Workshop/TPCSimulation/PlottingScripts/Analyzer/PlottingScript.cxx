// Include anything from C++ you might need below

// Include anything from ROOT you might need below
#include "TCanvas.h"

// Including the Input Manager
#include "IO/InputManager.hh"

int main(int argc, char** argv) {


  // --- Setting up the command line flags --- ///
  int opt;
  extern char *optarg;
  extern int   optopt;

  int         nEvent      = -1; // -n, number of events to run
  std::string InFileName  = ""; // -i, input ROOT file
  std::string InTreeName  = ""; // -t, input TTree name
  std::string OutFileName = ""; // -o, output file name, no extension

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
  im->LoadTree    (                  );

  // --- Grabbing the number of events to run over --- //
  int fNEvent = 0;
  if (nEvent == -1) fNEvent = im->GetEntries();
  else              fNEvent = std::min(nEvent, im->GetEntries());

  // --- Defining the plots we're going to make --- //
  // Put the plotting objects here, for example
  // TH1D *my_hist = new TH1D();
  
  // --- Stepping through each event --- //
  for (int CurrentEvent=0; CurrentEvent<fNEvent; ++CurrentEvent) {
    im->GetEntry(CurrentEvent);

    /*
    | You can get variables from the input manager using this.
    | Suppose we want to check the i-th mother id for the 
    | current event. We would do:
    | std::cout << (*im->True_Bck_Mother)[i] << std::endl;
    */
    
    // Determine which particle is your muon and proton
    
    // Work out the angle between the two particle trajectories
 
    // Fill the information into your plot


  }

  // --- Setting up the canvas for the output --- //

  /*
  | The c->Print() statement with ".pdf[" opens a new pdf file.
  | Each line with ".pdf" put after drawing a plot object
  | prints that plot to a new page of the pdf.
  | The line with ".pdf]" closes the output pdf. without this
  | you'll just have a corrupted pdf file
  */

  TCanvas *c = new TCanvas("c", "c");
  c->Print((OutFileName + ".pdf[").c_str());
  // Draw the plot you've made


  c->Print((OutFileName + ".pdf").c_str());
  c->Print((OutFileName + ".pdf]").c_str());


  return 0;

}

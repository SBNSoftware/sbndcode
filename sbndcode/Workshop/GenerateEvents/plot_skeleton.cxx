#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

// Draw histogram
void DrawHist(TH1D* hist, std::string name, std::string xaxis){

  TCanvas *canvas = new TCanvas(name.c_str(), "", 900, 600);

  hist->GetXaxis()->SetTitle(xaxis.c_str());
  hist->Draw();

  canvas->SaveAs((name+".png").c_str());
}

void plot_skeleton() {

  // open the input file
  TFile f("anatree_output.root");

  // setup a "TTreeReader" to read from the "anatree" Tree
  // Each entry is an art::Event
  TTreeReader anatree_reader("analysistree/anatree", &f);

  // Set up a reader for values in the tree -- use this as an example for retrieving values from the TTree
  // You can find the names by opening the tree in a TBrowser
  TTreeReaderValue<int> n_particles (anatree_reader, "geant_list_size"); // Number of particles simulated by g4
  TTreeReaderArray<float> P (anatree_reader, "P"); // Array of momenta

  // Define histograms
  TH1D* hMomentum = new TH1D("hMomentum", "", 20, 0, 2);

  // Loop over the analysis tree entries
  while (anatree_reader.Next()) {

    // Loop over the number of particles in the event
    for(size_t i = 0; i < *n_particles; i++){

      // Check if its a muon or proton originating from the vertex
      // and store their properties in a variable for use outside
      // this loop
     
      // Fill a histogram with all particle momentums -- just an example
      hMomentum->Fill(P[i]);
    }
      
    // Calculate the angle between the muon and proton 
      
    // Fill a histogram with the angle

  }

  // Save your histograms with the DrawHist funtion here
  DrawHist(hMomentum, "momentum", "P [GeV]");

}

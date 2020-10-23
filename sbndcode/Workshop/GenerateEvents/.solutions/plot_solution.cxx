#include "TFile.h"
#include "TH1F.h"
#include "TVector3.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

// Draw histogram
void DrawHist(TH1D* hist, std::string name, std::string xaxis){

  TCanvas *canvas = new TCanvas(name.c_str(), "", 900, 600);

  hist->GetXaxis()->SetTitle(xaxis.c_str());
  hist->Draw();

  canvas->SaveAs((name+".png").c_str());
}

void plot_solution() {

  // open the input file
  TFile f("anatree_output.root");

  // setup a "TTreeReader" to read from the "anatree" Tree
  // Each entry is an art::Event
  TTreeReader anatree_reader("analysistree/anatree", &f);

  // Set up a reader for values in the tree
  // You can find the names by opening the tree in a TBrowser
  TTreeReaderValue<int> n_particles (anatree_reader, "geant_list_size"); // Number of particles simulated by g4
  TTreeReaderArray<int> pdg (anatree_reader, "pdg"); // Array of momenta
  TTreeReaderArray<float> start_x (anatree_reader, "StartPointx"); // Array of starting x positions
  TTreeReaderArray<float> start_y (anatree_reader, "StartPointy"); // Array of starting y positions
  TTreeReaderArray<float> start_z (anatree_reader, "StartPointz"); // Array of starting z positions
  TTreeReaderArray<float> end_x (anatree_reader, "EndPointx"); // Array of end x positions
  TTreeReaderArray<float> end_y (anatree_reader, "EndPointy"); // Array of end y positions
  TTreeReaderArray<float> end_z (anatree_reader, "EndPointz"); // Array of end z positions

  // Define histograms
  TH1D* hAngle = new TH1D("hAngle", "", 20, 0, 180);

  // The simulated vertex
  TVector3 vertex(-100, 0, 50);

  // Loop over the analysis tree
  while (anatree_reader.Next()) {

    // Want to save the muon and proton end positions
    TVector3 mu_end;
    TVector3 p_end;

    // Loop over the number of particles in the event
    for(size_t i = 0; i < *n_particles; i++){
      TVector3 start(start_x[i], start_y[i], start_z[i]);
      TVector3 end(end_x[i], end_y[i], end_z[i]);

      // Check that particle comes from the vertex
      if(start != vertex) continue;
      
      // Get the muon end position
      if(pdg[i] == 13) mu_end = end;
      // Get the proton end position
      if(pdg[i] == 2212) p_end = end;
    }

    // Calculate the angle between them, convert to degrees
    double angle = (mu_end - vertex).Angle(p_end - vertex) * 180./TMath::Pi();
    // Fill the histogram
    hAngle->Fill(angle);

  }

  // Draw the histogram
  DrawHist(hAngle, "angle", "#theta_{#mu p} [deg]");

}

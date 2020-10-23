#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

void bonus_solution() {

  // open the input file
  TFile f("anatree_output.root");

  // setup a "TTreeReader" to read from the "anatree" Tree
  // Each entry is an art::Event
  TTreeReader anatree_reader("analysistree/anatree", &f);

  // Set up a reader for values in the tree
  // You can find the names by opening the tree in a TBrowser
  TTreeReaderValue<int> n_nu (anatree_reader, "mcevts_truth"); // Number of neutrino interactions
  TTreeReaderArray<int> nu_pdg (anatree_reader, "nuPDG_truth"); // Array of neutrino pdg
  TTreeReaderArray<int> ccnc (anatree_reader, "ccnc_truth"); // Array of neutrino CC or NC
  TTreeReaderArray<float> vertex_x (anatree_reader, "nuvtxx_truth"); // Array of neutrino vertex x positions
  TTreeReaderArray<float> vertex_y (anatree_reader, "nuvtxy_truth"); // Array of neutrino vertex y positions
  TTreeReaderArray<float> vertex_z (anatree_reader, "nuvtxz_truth"); // Array of neutrino vertex z positions
  TTreeReaderArray<int> n_primaries (anatree_reader, "genie_no_primaries"); // Array of number of primary particles produced by GENIE
  TTreeReaderArray<int> primary_pdg (anatree_reader, "genie_primaries_pdg"); // Array of pdg codes for primary particles
  TTreeReaderArray<int> primary_P (anatree_reader, "genie_P"); // Array of momenta for primary particles

  // Loop over the analysis tree
  while (anatree_reader.Next()) {

    // Loop over the number of particles in the event
    for(size_t i = 0; i < *n_nu; i++){
      // Get the neutrino flavour
      if(nu_pdg[i] < 0) std::cout<<"anti-";
      if(std::abs(nu_pdg[i])==12) std::cout<<"nu_e";
      if(std::abs(nu_pdg[i])==14) std::cout<<"nu_mu";
      // Pretty silly numbering convention, makes (slightly) more sense inside LArSoft
      if(ccnc[i] == 0) std::cout<<" CC ";
      if(ccnc[i] == 1) std::cout<<" NC ";

      // Count the number of protons
      int n_pr = 0;
      // The analysis tree is a bit broken here if we have neutrino pile-up!
      for(size_t j = 0; j < n_primaries[i]; j++){
        if(primary_pdg[j] == 2212 && primary_P[j] > 0.3) n_pr++;
      }
      std::cout<<n_pr<<"p:";

      // Get the neutrino vertex
      std::cout<<" vertex = ("<<vertex_x[i]<<", "<<vertex_y[i]<<", "<<vertex_z[i]<<")\n";
    }

  }

}

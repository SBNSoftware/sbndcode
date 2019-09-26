#include <string>
#include <vector>
#include <TChain.h>
#include <TTree.h>

#include <map>
#include <cmath>
#include <iostream>
#include <ctime>
#include <algorithm>

#include <TFile.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLine.h>
#include <TVector3.h>
#include <TLegend.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THStack.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCanvas.h>
#include <sstream>
#include <fstream>
#include <math.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

//TODO Unfolding/systematics

// Global configuration parameters
// File configurations
TString fInputFile;
TString fTreePath;
TString fMetaDataPath;
TString fOutputFile;
// Neutrino interaction configurations
std::vector<int> fNuPdg;
std::vector<int> fIsCC;
std::vector<bool> fContainedLepton;
std::vector<bool> fContainedParticles;
std::vector<double> fFiducial;
bool fPlotByFsi;
std::vector<int> fNumProtons;
std::vector<int> fNumPiPM;
std::vector<int> fNumPi0;
std::vector<int> fInteractionType;
// Plotting variable configurations
TString fStage;
std::vector<TString> fPlotVariables;
std::vector<bool> fShowPlots;
double fPotScale;
// Plotting option configurations
bool fPlotStacked;
TString fStackBy;
std::vector<double> fMinValue;
std::vector<double> fMaxValue;
std::vector<int> fNumBins;
std::vector<std::vector<double>> fBinEdges;
double fMaxError;
bool fPlotXSec;
bool fPlotFilled;
// Optional extras
bool fShowInfo;
bool fShowStatError;
bool fShowErrorBars;
bool fPlotEffPur;
bool fPlotResponse; //TODO
bool fUnfold; //TODO

// Set by functions
double fPot;
double fPotScaleFac = 1;
double fFiducialMass;
double fFlux;
double fTargets;
// Constants
// Integrated flux for each neutrino species
std::map<int, double> fNuFlux = {{14, 7.91}, {-14, 0.5475}, {12, 0.04696}, {-12, 0.004732}};
const std::vector<int> fCols = {46, 33, 38, 42, 40, 30, 49};

// Structure for holding interaction information
class Interaction
{
  public:

  bool selected;
  bool true_selected;
  std::string fsi;
  std::string int_type;
  std::string nu_type;
  std::vector<double> variables;
  std::vector<double> true_variables;

  Interaction(bool s, bool ts, std::string f, std::string i, std::string n, std::vector<double> v, std::vector<double> tv)
  {
    selected = s;
    true_selected = ts;
    fsi = f;
    int_type = i;
    nu_type = n;
    variables = v;
    true_variables = tv;
  }
};

// Structure for holding plot titles
class Titles
{
  public:

  std::vector<TString> hist_titles;
  std::vector<TString> names;
  std::vector<TString> units;
  TString data_type;
  TString part_cont;
  TString lep_cont;
  TString is_cc;
  TString n_pr;
  TString n_pipm;
  TString n_pi0;
  TString int_type;
  TString pot;
  TString mass;

  Titles(std::vector<TString> ht, std::vector<TString> n, std::vector<TString> u, TString dt, TString pc, TString lc, 
        TString ic, TString npr, TString npi, TString npi0, TString it, TString p, TString m)
  {
    hist_titles = ht;
    names = n;
    units = u;
    data_type = dt;
    part_cont = pc;
    lep_cont = lc;
    is_cc = ic;
    n_pr = npr;
    n_pipm = npi;
    n_pi0 = npi0;
    int_type = it;
    pot = p;
    mass = m;
  }

};

// Set some global style configurations here
void SetStyle(){

  Int_t font = 62;
  Double_t font_size = 0.06;
  Double_t line_width = 3.;
  Double_t label_offset = 0.01;

  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(55);
  gStyle->SetMarkerStyle(8);
  // Widths
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameLineWidth(line_width);
  gStyle->SetGridWidth(line_width);
  gStyle->SetLineWidth(line_width);
  // Fonts
  gStyle->SetTitleFont(font, "title");
  gStyle->SetTitleFont(font, "x");
  gStyle->SetTitleFont(font, "y");
  gStyle->SetTitleFont(font, "z");
  gStyle->SetLabelFont(font, "x");
  gStyle->SetLabelFont(font, "y");
  gStyle->SetLabelFont(font, "z");
  gStyle->SetTextFont(font);
  gStyle->SetLegendFont(font);
  // Sizes
  gStyle->SetTitleSize(1.3*font_size, "title");
  gStyle->SetTitleSize(font_size, "x");
  gStyle->SetTitleSize(font_size, "y");
  gStyle->SetTitleSize(font_size, "z");
  gStyle->SetLabelSize(font_size, "x");
  gStyle->SetLabelSize(font_size, "y");
  gStyle->SetLabelSize(font_size, "z");
  gStyle->SetMarkerSize(0.6);
  // Offsets
  gStyle->SetLabelOffset(label_offset, "x");
  gStyle->SetLabelOffset(label_offset, "y");
  gStyle->SetLabelOffset(label_offset, "z");
  // Legend
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);

}

// Convert comma separated string to vector of strings
std::vector<std::string> ToVector(std::string values, std::string delim = ","){

  // Expects comma separated string
  size_t pos = 0;
  std::vector<string> return_values;

  while ((pos = values.find(delim)) != std::string::npos) {
    std::string value = values.substr(0, pos);

    // Get rid of any []
    value.erase(std::remove(value.begin(), value.end(), '['), value.end());
    value.erase(std::remove(value.begin(), value.end(), ']'), value.end());

    return_values.push_back(value);
    values.erase(0, pos + delim.length());
  }
  values.erase(std::remove(values.begin(), values.end(), '['), values.end());
  values.erase(std::remove(values.begin(), values.end(), ']'), values.end());
  return_values.push_back(values);


  return return_values;
}

// Convert string to bool
std::vector<bool> ToBools(std::string values){
  std::vector<bool> return_bools;
  for(auto const& value : ToVector(values)){
    return_bools.push_back(value=="true");
  }
  return return_bools;
}

// Convert string to int
std::vector<int> ToInts(std::string values){
  std::vector<int> return_ints;
  for(auto const& value : ToVector(values)){
    return_ints.push_back(stoi(value));
  }
  return return_ints;
}

// Convert string to double
std::vector<double> ToDoubles(std::string values){
  std::vector<double> return_doubles;
  for(auto const& value : ToVector(values)){
    return_doubles.push_back(stod(value));
  }
  return return_doubles;
}

// Convert string to TString
std::vector<TString> ToTStrings(std::string values){
  std::vector<TString> return_doubles;
  for(auto const& value : ToVector(values)){
    return_doubles.push_back(TString(value));
  }
  return return_doubles;
}

// Read in the plotting configuration
void Configure(const std::string config_filename) {

  std::ifstream config_file;
  config_file.open(config_filename);
  if(!config_file){
    std::cerr << "Unable to open configuration file!\n";
    exit(1);
  }

  std::string delim = ":";
  std::string line;
  while(std::getline(config_file, line)){
    std::string key = line.substr(0, line.find(delim));
    std::string value = line.erase(0, line.find(delim)+delim.length());
    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
    // File configurations
    if(key.find("InputFile") != std::string::npos)    fInputFile = TString(value);
    if(key.find("TreePath") != std::string::npos)     fTreePath = TString(value);
    if(key.find("MetaDataPath") != std::string::npos) fMetaDataPath = TString(value);
    if(key.find("OutputFile") != std::string::npos)   fOutputFile = TString(value);
    // Neutrino configurations
    if(key.find("NuPdg") != std::string::npos)              fNuPdg = ToInts(value);
    if(key.find("IsCC") != std::string::npos)               fIsCC = ToInts(value);
    if(key.find("ContainedLepton") != std::string::npos)    fContainedLepton = ToBools(value);
    if(key.find("ContainedParticles") != std::string::npos) fContainedParticles = ToBools(value);
    if(key.find("Fiducial") != std::string::npos)           fFiducial = ToDoubles(value);
    if(key.find("PlotByFsi") != std::string::npos)          fPlotByFsi = (value=="true");
    if(key.find("NumProtons") != std::string::npos)         fNumProtons = ToInts(value);
    if(key.find("NumPiPM") != std::string::npos)            fNumPiPM = ToInts(value);
    if(key.find("NumPi0") != std::string::npos)             fNumPi0 = ToInts(value);
    if(key.find("InteractionType") != std::string::npos)    fInteractionType = ToInts(value);
    // Plotting variable configurations
    if(key.find("Stage") != std::string::npos)        fStage = TString(value);
    if(key.find("PlotVariable") != std::string::npos) fPlotVariables = ToTStrings(value);
    if(key.find("ShowPlots") != std::string::npos)    fShowPlots = ToBools(value);
    if(key.find("PotScale") != std::string::npos)     fPotScale = stod(value);
    // Plotting option configurations
    if(key.find("PlotStacked") != std::string::npos) fPlotStacked = (value=="true");
    if(key.find("StackBy") != std::string::npos)     fStackBy = TString(value);
    if(key.find("MinValue") != std::string::npos)    fMinValue = ToDoubles(value);
    if(key.find("MaxValue") != std::string::npos)    fMaxValue = ToDoubles(value);
    if(key.find("NumBins") != std::string::npos)     fNumBins = ToInts(value);
    if(key.find("BinEdges") != std::string::npos){
      for(auto const& val : ToVector(value, "],[")) fBinEdges.push_back(ToDoubles(val));
    }
    if(key.find("MaxError") != std::string::npos)    fMaxError = stod(value);
    if(key.find("PlotXSec") != std::string::npos)    fPlotXSec = (value=="true");
    if(key.find("PlotFilled") != std::string::npos)  fPlotFilled = (value=="true");
    // Optional extras
    if(key.find("ShowInfo") != std::string::npos)      fShowInfo = (value=="true");
    if(key.find("ShowStatError") != std::string::npos) fShowStatError = (value=="true");
    if(key.find("ShowErrorBars") != std::string::npos) fShowErrorBars = (value=="true");
    if(key.find("PlotEffPur") != std::string::npos)    fPlotEffPur = (value=="true");
    if(key.find("PlotResponse") != std::string::npos)  fPlotResponse = (value=="true");
    if(key.find("Unfold") != std::string::npos)        fUnfold = (value=="true");
  }

  if(fPlotVariables.size() != fMinValue.size()
     || fPlotVariables.size() != fMaxValue.size()
     || fPlotVariables.size() != fNumBins.size()
     || fPlotVariables.size() != fBinEdges.size()){
    std::cout<<"Must have same number of binning parameters as plotting variables!\n";
    exit(1);
  }
  if(fPlotVariables.size() > 3){
    std::cout<<"Sorry, ROOT doesn't do >3D histograms...\n";
    exit(1);
  }
  if(fPlotVariables.size() < 1){
    std::cout<<"Need something to plot in.\n";
    exit(1);
  }
  if(!(fStackBy=="fsi" || fStackBy=="int" || fStackBy=="nu")){
    std::cout<<"Unknown stack by parameter, not stacking.\n";
    fPlotStacked = false;
  }

}

// Check if an interaction matches the selection criteria
bool IsSelected(int nu_pdg, int cc, bool lep_contained, bool particles_contained, int n_pr, int n_pipm, int n_pi0, int int_type){
  bool selected = true;

  if (std::find(fNuPdg.begin(), fNuPdg.end(), nu_pdg) == fNuPdg.end())
    selected = false;
  if (std::find(fIsCC.begin(), fIsCC.end(), cc) == fIsCC.end())
    selected = false;
  if (std::find(fContainedLepton.begin(), fContainedLepton.end(), lep_contained) == fContainedLepton.end())
    selected = false;
  if (std::find(fContainedParticles.begin(), fContainedParticles.end(), particles_contained) == fContainedParticles.end())
    selected = false;
  if (fPlotByFsi){
    if (std::find(fNumProtons.begin(), fNumProtons.end(), -1) == fNumProtons.end()){
      if (std::find(fNumProtons.begin(), fNumProtons.end(), n_pr) == fNumProtons.end())
        selected = false;
    }
    if (std::find(fNumPiPM.begin(), fNumPiPM.end(), -1) == fNumPiPM.end()){
      if (std::find(fNumPiPM.begin(), fNumPiPM.end(), n_pipm) == fNumPiPM.end())
        selected = false;
    }
    if (std::find(fNumPi0.begin(), fNumPi0.end(), -1) == fNumPi0.end()){
      if (std::find(fNumPi0.begin(), fNumPi0.end(), n_pi0) == fNumPi0.end())
        selected = false;
    }
  }
  else{
    if (std::find(fInteractionType.begin(), fInteractionType.end(), -1) == fInteractionType.end()){
      if (std::find(fInteractionType.begin(), fInteractionType.end(), int_type) == fInteractionType.end())
        selected = false;
    }
  }

  return selected;
}

// Read in true variables
std::vector<Interaction> ReadData(){

  // Open the root tree
  TFile data_file(fInputFile);
  if(!data_file.IsOpen()){
    std::cout<<"Could not read input file!\n";
    exit(1);
  }

  TString prefix = fStage+"_";

  //Read in TTree
  TTreeReader tree_reader(fTreePath, &data_file);

  // True vertex
  TTreeReaderValue<double>       vtx_x(tree_reader, "vtx_x");
  TTreeReaderValue<double>       vtx_y(tree_reader, "vtx_y");
  TTreeReaderValue<double>       vtx_z(tree_reader, "vtx_z");

  // True quantities for stacked hists
  TTreeReaderValue<bool>         true_lep_contained(tree_reader, "true_lep_contained");
  TTreeReaderValue<bool>         true_particles_contained(tree_reader, "true_particles_contained");
  TTreeReaderValue<int>          true_cc(tree_reader, "true_cc");
  TTreeReaderValue<int>          true_nu_pdg(tree_reader, "true_nu_pdg");
  TTreeReaderValue<int>          true_int_type(tree_reader, "true_int_type");
  TTreeReaderValue<unsigned int> true_n_pipm(tree_reader, "true_n_pipm");
  TTreeReaderValue<unsigned int> true_n_pi0(tree_reader, "true_n_pi0");
  TTreeReaderValue<unsigned int> true_n_pr(tree_reader, "true_n_pr");

  // Need true values for efficiency/purity/response
  TTreeReaderValue<double>       true_nu_energy(tree_reader, "true_nu_energy");
  TTreeReaderValue<double>       true_lep_mom(tree_reader, "true_lep_mom");
  TTreeReaderValue<double>       true_lep_theta(tree_reader, "true_lep_theta");
  TTreeReaderValue<double>       true_pr1_mom(tree_reader, "true_pr1_mom");
  TTreeReaderValue<double>       true_pr1_theta(tree_reader, "true_pr1_theta");
  TTreeReaderValue<double>       true_lep_pr1_angle(tree_reader, "true_lep_pr1_angle");
  TTreeReaderValue<double>       true_pipm1_mom(tree_reader, "true_pipm1_mom");
  TTreeReaderValue<double>       true_pipm1_theta(tree_reader, "true_pipm1_theta");
  TTreeReaderValue<double>       true_lep_pipm1_angle(tree_reader, "true_lep_pipm1_angle");
  TTreeReaderValue<double>       true_delta_pt(tree_reader, "true_delta_pt");
  TTreeReaderValue<double>       true_delta_alphat(tree_reader, "true_delta_alphat");
  TTreeReaderValue<double>       true_delta_phit(tree_reader, "true_delta_phit");

  //Associate TTree values with variables
  TTreeReaderValue<bool>         lep_contained(tree_reader, prefix+"lep_contained");
  TTreeReaderValue<bool>         particles_contained(tree_reader, prefix+"particles_contained");
  TTreeReaderValue<int>          cc(tree_reader, prefix+"cc");
  TTreeReaderValue<int>          nu_pdg(tree_reader, prefix+"nu_pdg");
  TTreeReaderValue<int>          int_type(tree_reader, prefix+"int_type");
  TTreeReaderValue<unsigned int> n_pipm(tree_reader, prefix+"n_pipm");
  TTreeReaderValue<unsigned int> n_pi0(tree_reader, prefix+"n_pi0");
  TTreeReaderValue<unsigned int> n_pr(tree_reader, prefix+"n_pr");

  TTreeReaderValue<double>       nu_energy(tree_reader, prefix+"nu_energy");
  TTreeReaderValue<double>       lep_mom(tree_reader, prefix+"lep_mom");
  TTreeReaderValue<double>       lep_theta(tree_reader, prefix+"lep_theta");
  TTreeReaderValue<double>       pr1_mom(tree_reader, prefix+"pr1_mom");
  TTreeReaderValue<double>       pr1_theta(tree_reader, prefix+"pr1_theta");
  TTreeReaderValue<double>       lep_pr1_angle(tree_reader, prefix+"lep_pr1_angle");
  TTreeReaderValue<double>       pipm1_mom(tree_reader, prefix+"pipm1_mom");
  TTreeReaderValue<double>       pipm1_theta(tree_reader, prefix+"pipm1_theta");
  TTreeReaderValue<double>       lep_pipm1_angle(tree_reader, prefix+"lep_pipm1_angle");
  TTreeReaderValue<double>       delta_pt(tree_reader, prefix+"delta_pt");
  TTreeReaderValue<double>       delta_alphat(tree_reader, prefix+"delta_alphat");
  TTreeReaderValue<double>       delta_phit(tree_reader, prefix+"delta_phit");

  // Vector to fill
  std::vector<Interaction> interactions;
 
  // Loop over all the interactions
  while (tree_reader.Next()) {
   
    bool selected = IsSelected(*nu_pdg, *cc, *lep_contained, *particles_contained, 
                               *n_pr, *n_pipm, *n_pi0, *int_type);
    bool true_selected = IsSelected(*true_nu_pdg, *true_cc, *true_lep_contained, 
                                    *true_particles_contained, *true_n_pr, *true_n_pipm, *true_n_pi0, *true_int_type);

    // Check true vertex inside fiducial volume
    if(std::find(fFiducial.begin(), fFiducial.end(), -1) == fFiducial.end() && fFiducial.size() == 6){
      if(*vtx_x < -200+fFiducial[0] || *vtx_x > 200-fFiducial[3] ||
         *vtx_y < -200+fFiducial[1] || *vtx_x > 200-fFiducial[4] ||
         *vtx_z < 0+fFiducial[2] || *vtx_x > 500-fFiducial[5]){ 
        selected = false;
        true_selected = false;
      }
    }

    std::vector<double> variables;
    std::vector<double> true_variables;
    int index = 0;
    for(auto const& var : fPlotVariables){
      bool apply_cos = false;
      TString plot_var = var;
      if(var(0, 4) == "cos_"){ 
        apply_cos = true;
        plot_var = var(4, plot_var.Length());
      }

      if (plot_var == "lep_contained"){ 
        variables.push_back((double)(*lep_contained));
        true_variables.push_back((double)(*true_lep_contained));
      }
      if (plot_var == "particles_contained"){ 
        variables.push_back((double)(*particles_contained));
        true_variables.push_back((double)(*true_particles_contained));
      }
      if (plot_var == "cc"){ 
        variables.push_back((double)(*cc));
        true_variables.push_back((double)(*true_cc));
      }
      if (plot_var == "nu_pdg"){ 
        variables.push_back((double)(*nu_pdg));
        true_variables.push_back((double)(*true_nu_pdg));
      }
      if (plot_var == "int_type"){ 
        variables.push_back((double)(*int_type));
        true_variables.push_back((double)(*true_int_type));
      }
      if (plot_var == "n_pr"){ 
        variables.push_back((double)(*n_pr));
        true_variables.push_back((double)(*true_n_pr));
      }
      if (plot_var == "n_pipm"){ 
        variables.push_back((double)(*n_pipm));
        true_variables.push_back((double)(*true_n_pipm));
      }
      if (plot_var == "n_pi0"){ 
        variables.push_back((double)(*n_pi0));
        true_variables.push_back((double)(*true_n_pi0));
      }
      if (plot_var == "nu_energy"){ 
        variables.push_back(*nu_energy);
        true_variables.push_back(*true_nu_energy);
      }
      if (plot_var == "lep_mom"){ 
        variables.push_back(*lep_mom);
        true_variables.push_back(*true_lep_mom);
      }
      if (plot_var == "lep_theta"){ 
        variables.push_back(*lep_theta);
        true_variables.push_back(*true_lep_theta);
      }
      if (plot_var == "pr1_mom"){ 
        variables.push_back(*pr1_mom);
        true_variables.push_back(*true_pr1_mom);
      }
      if (plot_var == "pr1_theta"){ 
        variables.push_back(*pr1_theta);
        true_variables.push_back(*true_pr1_theta);
      }
      if (plot_var == "lep_pr1_angle"){ 
        variables.push_back(*lep_pr1_angle);
        true_variables.push_back(*true_lep_pr1_angle);
      }
      if (plot_var == "pipm1_mom"){ 
        variables.push_back(*pipm1_mom);
        true_variables.push_back(*true_pipm1_mom);
      }
      if (plot_var == "pipm1_theta"){ 
        variables.push_back(*pipm1_theta);
        true_variables.push_back(*true_pipm1_theta);
      }
      if (plot_var == "lep_pipm1_angle"){ 
        variables.push_back(*lep_pipm1_angle);
        true_variables.push_back(*true_lep_pipm1_angle);
      }
      if (plot_var == "delta_pt"){ 
        variables.push_back(*delta_pt);
        true_variables.push_back(*true_delta_pt);
      }
      if (plot_var == "delta_alphat"){ 
        variables.push_back(*delta_alphat);
        true_variables.push_back(*true_delta_alphat);
      }
      if (plot_var == "delta_phit"){ 
        variables.push_back(*delta_phit);
        true_variables.push_back(*true_delta_phit);
      }

      if(apply_cos) variables[index] = cos(variables[index]);
      index++;
    }

    // FSI: 0pi0p, 0pi1p, 0pi2+p, 1pi, 2+pi, 1+pi0
    std::string fsi_string = "other";
    if(*true_n_pi0 >= 1) fsi_string = "#geq1#pi^{0}";
    else if(*true_n_pipm == 1) fsi_string = "1#pi^{#pm}";
    else if(*true_n_pipm >= 2) fsi_string = "#geq2#pi^{#pm}";
    else if(*true_n_pr == 0) fsi_string = "0#pi0p";
    else if(*true_n_pr == 1) fsi_string = "0#pi1p";
    else if(*true_n_pr >= 2) fsi_string = "0#pi#geq2p";

    // Int: QE, RES, COH, DIS, MEC
    std::string int_string = "other";
    if(*true_int_type == 0) int_string = "QE";
    else if(*true_int_type == 1) int_string = "RES";
    else if(*true_int_type == 2) int_string = "DIS";
    else if(*true_int_type == 3) int_string = "COH";
    else if(*true_int_type == 10) int_string = "MEC";

    std::string nu_string = "other";
    if(*true_nu_pdg == -12 && *true_cc) nu_string = "#bar{#nu}_{e} CC";
    if(*true_nu_pdg == -12 && !*true_cc) nu_string = "#bar{#nu}_{e} NC";
    if(*true_nu_pdg == -14 && *true_cc) nu_string = "#bar{#nu}_{#mu} CC";
    if(*true_nu_pdg == -14 && !*true_cc) nu_string = "#bar{#nu}_{#mu} NC";
    if(*true_nu_pdg == 12 && *true_cc) nu_string = "#nu_{e} CC";
    if(*true_nu_pdg == 12 && !*true_cc) nu_string = "#nu_{e} NC";
    if(*true_nu_pdg == 14 && *true_cc) nu_string = "#nu_{#mu} CC";
    if(*true_nu_pdg == 14 && !*true_cc) nu_string = "#nu_{#mu} NC";

    // Check that the event matches the selection
    if(!selected && !true_selected) continue;
    // Check that all variables are filled FIXME is this right?
    if(std::find(true_variables.begin(), true_variables.end(), -99999) != true_variables.end()) continue;
    
    Interaction interaction(selected, true_selected, fsi_string, int_string, nu_string, variables, true_variables);
    interactions.push_back(interaction);
    
  }

  // Throw exception if nothing to plot
  if (interactions.size() == 0){
    std::cout << "No events match your input parameters." << std::endl;
    throw std::exception();
  }

  return interactions;
 
}

// Get the POT, flux and target number from the metadata
void GetMetaData(){

  // Open the root tree file
  TFile data_file(fInputFile);

  //Read in TTree
  TTreeReader tree_reader(fMetaDataPath, &data_file);
  TTreeReaderValue<double> pot(tree_reader, "pot");
  double pot_count = 0;
  
  while (tree_reader.Next()) {
    pot_count += *pot;
  }
  
  fPot = pot_count;

  if(fPotScale > 0){
    fPotScaleFac = fPotScale/fPot;
  }

  double flux_factor = 0;
  for(auto const& pdg : fNuPdg){
    if(fNuFlux.find(pdg) == fNuFlux.end()){
      std::cout<<"Unknown neutrino PDG code!\n";
      exit(1);
    }
    flux_factor += fNuFlux[pdg];
  }
  fFlux = flux_factor * 1e-6 * fPot * fPotScaleFac / (10000*0.05); // [cm^-2]

  double volume = 400*400*500; // [cm^3]
  if(std::find(fFiducial.begin(), fFiducial.end(), -1) == fFiducial.end() && fFiducial.size() == 6){
    volume = (400-fFiducial[0]-fFiducial[3])*(400-fFiducial[1]-fFiducial[4])*(500-fFiducial[2]-fFiducial[5]); // [cm^3]
  }
  fFiducialMass = 1.3973*volume/1e6; //[tons]
  fTargets = 6.022e23 * fFiducialMass * 1e3 * 40/ (0.03995); // [/nucleon]

}

// Set the minimum bin value
double GetMinBin(std::vector<double> data, int i){
  
  double min_bins = data[0];

  if(fPlotVariables[i]=="particles_contained"||
     fPlotVariables[i]=="lep_contained"||
     fPlotVariables[i]=="cc"){
    min_bins = 0;
  }
  if(fPlotVariables[i]=="nu_pdg") min_bins = -14;
  if(fPlotVariables[i]=="int_type") min_bins = 0;
  if(fPlotVariables[i]=="n_pr") min_bins = 0;
  if(fPlotVariables[i]=="n_pipm") min_bins = 0;
  if(fPlotVariables[i]=="n_pi0") min_bins = 0;

  return (min_bins);
 
}

// Set the maximum bin value
double GetMaxBin(std::vector<double> data, int i){
  
  double max_bins = data[data.size()-1];

  if(fPlotVariables[i]=="particles_contained"||
     fPlotVariables[i]=="lep_contained"||
     fPlotVariables[i]=="cc"){
    max_bins = 2;
  }
  if(fPlotVariables[i]=="nu_pdg") max_bins = 15;
  if(fPlotVariables[i]=="int_type") max_bins = 11;
  if(fPlotVariables[i]=="n_pr") max_bins = 11;
  if(fPlotVariables[i]=="n_pipm") max_bins = 11;
  if(fPlotVariables[i]=="n_pi0") max_bins = 11;

  return (max_bins);
 
}

// Set the default number of bins
int DefaultBins(int data_size, int dims, int i){
  
  int dflt_bins = (int)pow(2*pow(data_size,.33), 1./dims);

  if(fPlotVariables[i]=="particles_contained"||
     fPlotVariables[i]=="lep_contained"||
     fPlotVariables[i]=="cc"){
    dflt_bins = 2;
  }
  if(fPlotVariables[i]=="nu_pdg") dflt_bins = 29;
  if(fPlotVariables[i]=="int_type") dflt_bins = 11;
  if(fPlotVariables[i]=="n_pr") dflt_bins = 11;
  if(fPlotVariables[i]=="n_pipm") dflt_bins = 11;
  if(fPlotVariables[i]=="n_pi0") dflt_bins = 11;

  return (dflt_bins);
 
}


// Recursive function to merge bins to get up to required statistical error
std::vector<double> ChangeBinning(TH1D* hist, double max){

  std::vector<double> bin_edges;
  for(size_t i = 0; i < hist->GetNbinsX(); i++){
    bin_edges.push_back(hist->GetBinLowEdge(i+1));
  }
  bin_edges.push_back(max);

  for(size_t i = 0; i < hist->GetNbinsX(); i++){
    if(hist->GetBinError(i+1)/hist->GetBinContent(i+1) > fMaxError){

      if(i != hist->GetNbinsX()-1){
        bin_edges.erase(bin_edges.begin()+i+1);
        double edges_array[bin_edges.size()];
        std::copy(bin_edges.begin(), bin_edges.end(), edges_array);
        TH1D* new_hist = (TH1D*)hist->Rebin(bin_edges.size()-1, "new", edges_array);
        return ChangeBinning(new_hist, max);
      }

      else{
        bin_edges.erase(bin_edges.begin()+(bin_edges.size()-2));
        return bin_edges;
      }

    }
  }
  return bin_edges;
}


// Function to get or calculate the number of bins, and min and max bins
std::vector<std::vector<double>> GetBinning(std::vector<std::vector<double>> data_v){

  std::map<int, std::vector<double>> data_map;
  for(auto const& data : data_v){
    for(size_t i = 0; i < data.size(); i++){
      data_map[i].push_back(data[i]);
    }
  }
  for(size_t i = 0; i < data_map.size(); i++){
    std::sort(data_map[i].begin(), data_map[i].end());
  }
  std::vector<double> hist_min = fMinValue;
  std::vector<double> hist_max = fMaxValue;
  std::vector<int> hist_bins = fNumBins;
  for(size_t i = 0; i < fPlotVariables.size(); i++){
    if(fMinValue[i] < 0){
      hist_min[i] = GetMinBin(data_map[i], i);
    }
    if(fMaxValue[i] <= 0){
      hist_max[i] = GetMaxBin(data_map[i], i);
    }
    if(fNumBins[i] <= 0){
      int data_size = 0;
      for(auto const& data : data_v){
        if(data[i] >= hist_min[i] && data[i] <= hist_max[i]) data_size++;
      }
      hist_bins[i] = DefaultBins(data_size*fPotScaleFac, data_map.size(), i);
    }
  }
  
  std::vector<std::vector<double>> all_bin_edges;
  for(size_t i = 0; i < fPlotVariables.size(); i++){
    // If bin edges are set by the user
    if(fBinEdges[i].size() > 1){
      all_bin_edges.push_back(fBinEdges[i]);
    }
    else{
      TH1D *temp_hist = new TH1D("temp_hist", "", hist_bins[i], hist_min[i], hist_max[i]);
      std::vector<double> bin_edges;
      for(size_t j = 0; j < temp_hist->GetNbinsX(); j++){
        bin_edges.push_back(temp_hist->GetBinLowEdge(j+1));
      }
      delete temp_hist;
      bin_edges.push_back(hist_max[i]);
      all_bin_edges.push_back(bin_edges);
    }
  }

  // See if bin edges set by user
  //

  // If a maximum bin error is set
  if(fMaxError > 0){
    for(size_t i = 0; i < fPlotVariables.size(); i++){
      // Get the optimal binning for the first variable
      TH1D *temp_hist = new TH1D("temp_hist", "", hist_bins[i], hist_min[i], hist_max[i]);
      // Loop over data
      for(auto const& data : data_v){
        // Fill temporary histogram
        temp_hist->Fill(data[i]);
      }
      // Include scale factor bin by bin as Scale() won't change errors
      for(size_t n = 1; n <= temp_hist->GetNbinsX(); n++){
        temp_hist->SetBinContent(n, temp_hist->GetBinContent(n)*fPotScaleFac);
      }
      // Change the binning so that all bin errors below maximum
      std::vector<double> bin_edges = ChangeBinning(temp_hist, hist_max[i]);
      delete temp_hist;
      all_bin_edges[i] = bin_edges;
    }
  }

  return all_bin_edges;
}

// Set info based on configuration
Titles GetTitles(){
  
  // Set the units and histogram titles based on plotting variable
  std::vector<TString> names;
  std::vector<TString> units;
  std::vector<TString> hist_titles;
  int index = 0;
  for(auto const& var : fPlotVariables){
    bool apply_cos = false;
    TString plot_var = var;
    if(var(0, 4) == "cos_"){ 
      apply_cos = true;
      plot_var = var(4, plot_var.Length());
    }

    if (plot_var == "nu_energy"){
      names.push_back("E_{#nu}");
      units.push_back("GeV");
      hist_titles.push_back("neutrino energy");
    } 
    else if (plot_var == "lep_mom"){
      names.push_back("P_{lep}");
      units.push_back("GeV");
      hist_titles.push_back("lepton momentum");
    }
    else if (plot_var == "pr1_mom"){
      names.push_back("P_{p}");
      units.push_back("GeV");
      hist_titles.push_back("leading proton momentum");
    }
    else if (plot_var == "pipm1_mom"){
      names.push_back("P_{#pi}");
      units.push_back("GeV");
      hist_titles.push_back("leading #pi^{#pm} momentum");  
    }
    else if (plot_var == "lep_theta"){
      names.push_back("#theta_{lep}");
      units.push_back("rad");
      hist_titles.push_back("lepton #theta");
    }
    else if (plot_var == "pr1_theta"){
      names.push_back("#theta_{p}");
      units.push_back("rad");
      hist_titles.push_back("leading proton #theta");
    }
    else if (plot_var == "pipm1_theta"){
      names.push_back("#theta_{#pi}");
      units.push_back("rad");
      hist_titles.push_back("leading #pi^{#pm} #theta");
    }
    else if (plot_var == "lep_pr1_angle"){
      names.push_back("#theta_{lep,p}");
      units.push_back("rad");
      hist_titles.push_back("angle between lepton and proton");
    }
    else if (plot_var == "lep_pipm1_angle"){
      names.push_back("#theta_{lep,#pi}");
      units.push_back("rad");
      hist_titles.push_back("angle between lepton and #pi^{#pm}");
    }
    else if (plot_var == "nu_pdg"){
      names.push_back("PDG code");
      units.push_back("");
      hist_titles.push_back("neutrino PDG");
    }
    else if (plot_var == "lep_contained"){
      names.push_back("lepton contained?");
      units.push_back("");
      hist_titles.push_back("lepton containment");
    }
    else if (plot_var == "particles_contained"){
      names.push_back("particles contained?");
      units.push_back("");
      hist_titles.push_back("secondary particle containment");
    }
    else if (plot_var == "cc"){
      names.push_back("Is CC?");
      units.push_back("");
      hist_titles.push_back("charged or neutral current");
    }
    else if (plot_var == "int_type"){
      names.push_back("interaction code");
      units.push_back("");
      hist_titles.push_back("QE(0) RES(1) DIS(2) COH(3) MEC(10)");
    }
    else if (plot_var == "n_pr"){
      names.push_back("N_{p}");
      units.push_back("");
      hist_titles.push_back("number of protons");
    }
    else if (plot_var == "n_pipm"){
      names.push_back("N_{#pi}");
      units.push_back("");
      hist_titles.push_back("number of #pi^{#pm}");
    }
    else if (plot_var == "n_pi0"){
      names.push_back("N_{#pi^{0}}");
      units.push_back("");
      hist_titles.push_back("number of #pi^{0}");
    }
    else if (plot_var == "delta_pt"){
      names.push_back("#delta p_{T}");
      units.push_back("GeV");
      hist_titles.push_back("#delta p_{T}");
    }
    else if (plot_var == "delta_alphat"){
      names.push_back("#delta #alpha_{T}");
      units.push_back("deg");
      hist_titles.push_back("#delta #alpha_{T}");
    }
    else if (plot_var == "delta_phit"){
      names.push_back("#delta #phi_{T}");
      units.push_back("deg");
      hist_titles.push_back("#delta #phi_{T}");
    }
    else{
      std::cout<<"Invalid plotting variable!\n";
      exit(1);
    }

    // TODO remove units
    if(apply_cos){
      names[index] = "cos"+ names[index];
      units[index] = "";
      hist_titles[index] = "cos "+hist_titles[index];
    }
    index++;
  }

  // True or reco
  TString data_type;
  if (fStage == "reco")
    data_type = "Reconstructed";
  else if (fStage == "true")
    data_type = "Truth";
  else if (fStage == "smeareff")
    data_type = "Smearing + Efficiency";
  else{
    std::cout<<"Unrecognised stage!\n";
    exit(1);
  }

  // Particle containment
  TString part_cont;
  if (fContainedParticles.size() == 1){
    if (fContainedParticles[0])
      part_cont = "Contained Particles";
    else
      part_cont = "Exiting Particles";
  }
  else
    part_cont = "Cont+Exit Particles";

  // Lepton containment
  TString lep_cont;
  if (fContainedLepton.size() == 1){
    if (fContainedLepton[0])
      lep_cont = "Contained Lepton";
    else
      lep_cont = "Exiting Lepton";
  }
  else
    lep_cont = "Cont+Exit Lepton";
  
  // Neutrino pdg
  TString nu_type;
  for(auto const& pdg : fNuPdg){
    if(pdg == 12) nu_type += "#nu_{e} ";
    if(pdg == 14) nu_type += "#nu_{#mu} ";
    if(pdg == -12) nu_type += "#bar{#nu}_{e} ";
    if(pdg == -14) nu_type += "#bar{#nu}_{#mu} ";
  }

  // Charged or neutral current
  TString is_cc;
  if (fIsCC.size() == 1){
    if (fIsCC[0])
      is_cc = "CC";
    else
      is_cc = "NC";
  }
  else
    is_cc = "CC+NC";

  is_cc = nu_type + is_cc;

  // Number of protons
  TString n_pr;
  if(std::find(fNumProtons.begin(), fNumProtons.end(), -1) != fNumProtons.end())
    n_pr = "All ";
  else{
    for(auto const& num : fNumProtons)
      n_pr += std::to_string(num) + " ";
  }
  n_pr += "Proton(s)";

  // Number of charged pions
  TString n_pipm;
  if(std::find(fNumPiPM.begin(), fNumPiPM.end(), -1) != fNumPiPM.end())
    n_pipm = "All ";
  else{
    for(auto const& num : fNumPiPM)
      n_pipm += std::to_string(num) + " ";
  }
  n_pipm += "Charged Pion(s)";

  // Number of neutral pions
  TString n_pi0;
  if(std::find(fNumPi0.begin(), fNumPi0.end(), -1) != fNumPi0.end())
    n_pi0 = "All ";
  else{
    for(auto const& num : fNumPi0)
      n_pi0 += std::to_string(num) + " ";
  }
  n_pi0 += "Neutral Pion(s)";

  // Interaction type
  TString int_type;
  if(std::find(fInteractionType.begin(), fInteractionType.end(), -1) != fInteractionType.end())
    int_type = "All Interactions";
  else{
    for(auto const& i_type : fInteractionType){
      if(i_type == 0) int_type += "QE ";
      if(i_type == 1) int_type += "RES ";
      if(i_type == 2) int_type += "DIS ";
      if(i_type == 3) int_type += "COH ";
      if(i_type == 10) int_type += "MEC ";
    }
  }

  // POT
  std::stringstream pot_stream;
  double pot = fPotScale;
  if(fPotScale <= 0) pot = fPot;
  pot_stream << std::setprecision(3) << "POT = " << pot << "}";
  std::string pot_string = pot_stream.str();
  pot_string.replace(pot_string.find("e"), 1, "#times10^{");
  pot_string.replace(pot_string.find("+"), 1, "");

  // Fiducial mass
  std::stringstream mass_stream;
  mass_stream << std::setprecision(3) << "Fid Mass = " << fFiducialMass << " t";
  std::string mass_string = mass_stream.str();

  Titles titles(hist_titles, names, units, data_type, part_cont, lep_cont, is_cc, n_pr, 
                n_pipm, n_pi0, int_type, TString(pot_string), TString(mass_string));

  return(titles);
  
}

// Draw additional information on to hist
void DrawInfo(Titles titles, double width, double height, double size){

  TLatex *POT        = new TLatex(width, .97*height, titles.pot);
  TLatex *mass       = new TLatex(width, .91*height, titles.mass);
  TLatex *data_type  = new TLatex(width, .85*height, titles.data_type);
  TLatex *is_cc      = new TLatex(width, .79*height, titles.is_cc);
  TLatex *part_cont  = new TLatex(width, .72*height, titles.part_cont);
  TLatex *lep_cont   = new TLatex(width, .66*height, titles.lep_cont);
  TLatex *n_pr       = new TLatex(width, .60*height, titles.n_pr);
  TLatex *n_pipm     = new TLatex(width, .54*height, titles.n_pipm);
  TLatex *n_pi0      = new TLatex(width, .48*height, titles.n_pi0);
  TLatex *int_type   = new TLatex(width, .60*height, titles.int_type);

  // Set the text size
  POT->SetTextSize(size);
  mass->SetTextSize(size);
  data_type->SetTextSize(size);
  part_cont->SetTextSize(size);
  lep_cont->SetTextSize(size);
  is_cc->SetTextSize(size);
  n_pr->SetTextSize(size);
  n_pipm->SetTextSize(size);
  n_pi0->SetTextSize(size);
  int_type->SetTextSize(size);

  // Draw the info text
  POT->Draw("same");
  mass->Draw("same");
  data_type->Draw("same");
  part_cont->Draw("same");
  lep_cont->Draw("same");
  is_cc->Draw("same");
  if(fPlotByFsi){
    n_pr->Draw("same");
    n_pipm->Draw("same");
    n_pi0->Draw("same");
  }
  else{
    int_type->Draw("same");
  }

}

// Get the Y axis title when plotting cross sections
TString GetXSecTitle(Titles titles, int i, int j = -1, int k = -1){

  TString xsec_title = "d#sigma/d"+titles.names[i]+" [10^{-38}#frac{cm^{2}}{"+titles.units[i]+" n}]";

  if(j != -1){
    xsec_title = "d^{2}#sigma/d"+titles.names[i]+"d"+titles.names[j]+" [10^{-38}#frac{cm^{2}}{"+titles.units[i]+" "+titles.units[j]+" n}]";

    if(k != -1){
      xsec_title = "d^{3}#sigma/d"+titles.names[i]+"d"+titles.names[j]+"d"+titles.names[k]+" [10^{-38}#frac{cm^{2}}{"+titles.units[i]+" "+titles.units[j]+" "+titles.units[k]+" n}]";
    }

  }
  return TString(xsec_title);
}

// Plot a 1D stacked hist with statistical errors on the bottom
void Plot1DWithErrors(THStack* hstack, TLegend* legend, TH1D* error_bands, Titles titles, TH1D* total_hist, size_t i, size_t j = -1, size_t k = -1){

  // Create the canvas
  TString name = hstack->GetName();
  TCanvas *canvas = new TCanvas("canvas"+name,"canvas");

  // Split the pad for histogram and error plot
  double pad_split = .3;
  TPad *upper_pad = new TPad("upper_pad", "" , 0., pad_split, 1.0, 1.0);
  upper_pad->SetTopMargin(0.12);
  upper_pad->SetBottomMargin(0.075);
  upper_pad->SetLeftMargin(0.12);
  upper_pad->SetRightMargin(0.05);

  TPad *lower_pad = new TPad("lower_pad", "", 0., 0., 1., pad_split);
  lower_pad->SetTopMargin(0.01);
  lower_pad->SetBottomMargin(0.34);
  lower_pad->SetLeftMargin(0.12);
  lower_pad->SetRightMargin(0.05);

  upper_pad->Draw();
  lower_pad->Draw();

  // Fill the upper pad with histogram, info and legend
  upper_pad->cd();

  // Draw the stacked histogram and legend
  hstack->Draw("HIST");
  if(fShowErrorBars){
    total_hist->SetLineWidth(2);
    total_hist->SetMarkerStyle(1);
    total_hist->Draw("E1 X0 SAME");
  }

  if(fPlotStacked){
    legend->SetNColumns(legend->GetNRows());
    legend->SetFillStyle(0);
    legend->Draw();
  }
  // Set the titles
  if(fPlotXSec){
    hstack->GetYaxis()->SetTitle(GetXSecTitle(titles, i, j, k));
  }
  else if(fMaxError > 0){
    hstack->GetYaxis()->SetTitle("Events (/bin width)");
  }
  else{
    hstack->GetYaxis()->SetTitle("Events");
  }
  // X axis config
  hstack->GetXaxis()->SetLabelOffset(0.1);
  hstack->GetXaxis()->SetTitleOffset(1.8);
  hstack->GetXaxis()->SetTickLength(0.04);
  // Y axis config
  hstack->GetYaxis()->SetTitleOffset(0.8);
  double title_size = 1.1*hstack->GetYaxis()->GetTitleSize();
  if(fPlotXSec && fPlotVariables.size()==1){ 
    title_size = 1.0*hstack->GetYaxis()->GetTitleSize();
    hstack->GetYaxis()->SetTitleOffset(0.9);
  }
  if(fPlotXSec && fPlotVariables.size()==2){ 
    title_size = 0.8*hstack->GetYaxis()->GetTitleSize();
    hstack->GetYaxis()->SetTitleOffset(1.0);
  }
  if(fPlotXSec && fPlotVariables.size()==3){ 
    title_size = 0.6*hstack->GetYaxis()->GetTitleSize();
    hstack->GetYaxis()->SetTitleOffset(1.1);
  }
  hstack->GetYaxis()->SetTitleSize(title_size);
  hstack->GetYaxis()->SetNdivisions(110);
  hstack->GetYaxis()->SetTickLength(0.015);
  canvas->Modified();

  // Info text
  // Text position and content
  double width = 0.7*(hstack->GetXaxis()->GetXmax()-hstack->GetXaxis()->GetXmin())+hstack->GetXaxis()->GetXmin();
  double height = hstack->GetMaximum();
  double upper_text_size = 0.7*hstack->GetYaxis()->GetTitleSize();
  if(fShowInfo) DrawInfo(titles, width, height, upper_text_size);
 
  // Fill the lower pad with percentage error per bin
  lower_pad->cd();
  lower_pad->SetTickx();
  lower_pad->SetTicky();
 
  // Set axis titles
  error_bands->SetFillColor(38);
  error_bands->SetLineColor(38);
  error_bands->GetYaxis()->SetTitle("#sigma_{stat} (%)");
  error_bands->GetXaxis()->SetTitle(titles.names[i]+" ["+titles.units[i]+"]");

  double size_ratio = upper_pad->GetAbsHNDC()/lower_pad->GetAbsHNDC();
  // x axis config
  error_bands->GetXaxis()->SetTitleSize(1.1*size_ratio*error_bands->GetXaxis()->GetTitleSize());
  error_bands->GetXaxis()->SetLabelSize(size_ratio*error_bands->GetXaxis()->GetLabelSize());
  error_bands->GetXaxis()->SetLabelOffset(0.04);
  error_bands->GetXaxis()->SetTickLength(size_ratio*0.04);
  error_bands->SetTitleOffset(1.0, "x");
  // y axis config
  error_bands->GetYaxis()->SetTitleSize(1.1*size_ratio*error_bands->GetYaxis()->GetTitleSize());
  error_bands->GetYaxis()->SetLabelSize(size_ratio*error_bands->GetYaxis()->GetLabelSize());
  error_bands->GetYaxis()->CenterTitle();
  error_bands->GetYaxis()->SetTickLength(0.015);
  error_bands->SetNdivisions(105, "y");
  error_bands->SetTitleOffset(0.3, "y");

  // Draw the error bars
  if(error_bands->GetNbinsX() < 40) error_bands->Draw("B");
  else error_bands->Draw("C");
  
  TString output_file = fOutputFile;
  output_file.ReplaceAll(".","_"+name+".");
  canvas->SaveAs(output_file);
}

// Plot a 1D stacked hist
void Plot1D(THStack* hstack, TLegend* legend, Titles titles, TH1D* total_hist, size_t i, size_t j = -1, size_t k = -1){

  // Create the canvas
  TString name = hstack->GetName();
  TCanvas *canvas = new TCanvas("canvas"+name,"canvas");

  // Split the pad for histogram and error plot
  canvas->SetTopMargin(0.15);
  canvas->SetBottomMargin(0.15);
  canvas->SetLeftMargin(0.12);
  canvas->SetRightMargin(0.04);

  // Draw the stacked histogram and legend
  hstack->Draw("HIST");

  if(fShowErrorBars){
    total_hist->SetLineWidth(2);
    total_hist->SetMarkerStyle(1);
    total_hist->Draw("E1 X0 SAME");
  }

  if(fPlotStacked){
    legend->SetNColumns(legend->GetNRows());
    legend->SetFillStyle(0);
    legend->Draw();
    canvas->Update();
    legend->SetX1NDC(0.24);
    legend->SetY1NDC(0.85);
    legend->SetX2NDC(0.96);
    legend->SetY2NDC(0.91);
    canvas->Modified();
  }
  // Set the titles
  if(fPlotXSec){
    hstack->GetYaxis()->SetTitle(GetXSecTitle(titles, i, j, k));
  }
  else if(fMaxError > 0){
    hstack->GetYaxis()->SetTitle("Events (/Bin width)");
  }
  else{
    hstack->GetYaxis()->SetTitle("Events");
  }
  hstack->GetXaxis()->SetTitle(titles.names[i]+" ["+titles.units[i]+"]");
  // X axis config
  hstack->GetXaxis()->SetTitleOffset(1.);
  hstack->GetXaxis()->SetTickLength(0.02);
  hstack->GetXaxis()->SetTitleSize(1.1*hstack->GetXaxis()->GetTitleSize());
  // Y axis config
  hstack->GetYaxis()->SetTitleOffset(0.95);
  hstack->GetYaxis()->SetTickLength(0.015);
  double title_size = 1.1*hstack->GetYaxis()->GetTitleSize();
  if(fPlotXSec && fPlotVariables.size()==1){ 
    title_size = 1.0*hstack->GetYaxis()->GetTitleSize();
    hstack->GetYaxis()->SetTitleOffset(1.05);
  }
  if(fPlotXSec && fPlotVariables.size()==2){ 
    title_size = 0.8*hstack->GetYaxis()->GetTitleSize();
    hstack->GetYaxis()->SetTitleOffset(1.15);
  }
  if(fPlotXSec && fPlotVariables.size()==3){ 
    title_size = 0.6*hstack->GetYaxis()->GetTitleSize();
    hstack->GetYaxis()->SetTitleOffset(1.25);
  }
    
  hstack->GetYaxis()->SetTitleSize(title_size);
  hstack->GetYaxis()->SetNdivisions(110);
  if(fPlotXSec && fPlotVariables.size()==1)
  canvas->Modified();

  // Text position and content
  double width = 0.65*(hstack->GetXaxis()->GetXmax()-hstack->GetXaxis()->GetXmin())+hstack->GetXaxis()->GetXmin();
  double height = hstack->GetMaximum();
  double upper_text_size = 0.6*hstack->GetYaxis()->GetTitleSize();
  if(fShowInfo) DrawInfo(titles, width, height, upper_text_size);
  
  TString output_file = fOutputFile;
  output_file.ReplaceAll(".","_"+name+".");
  canvas->SaveAs(output_file);
}



// Draw efficiency/purity as function of some variable
void PlotEfficiency(TH1D* select, TH1D* total, TString name, TString xaxis, TString yaxis){

  TCanvas *canvas = new TCanvas(name+yaxis, "");
  canvas->SetTopMargin(0.1);
  canvas->SetBottomMargin(0.16);
  canvas->SetLeftMargin(0.14);
  canvas->SetRightMargin(0.04);

  TGraphAsymmErrors *graph = new TGraphAsymmErrors();

  graph->SetMarkerColor(46);
  graph->SetLineColor(46);
  graph->GetXaxis()->SetTitle(xaxis);
  graph->GetYaxis()->SetTitle(yaxis);
  // X axis config
  graph->GetXaxis()->SetTitleOffset(1.1);
  graph->GetXaxis()->SetTickLength(0.04);
  graph->GetXaxis()->SetTitleSize(1.1*graph->GetXaxis()->GetTitleSize());
  // Y axis config
  graph->GetYaxis()->SetTitleOffset(0.95);
  graph->GetYaxis()->SetTickLength(0.015);
  graph->GetYaxis()->SetTitleSize(1.1*graph->GetYaxis()->GetTitleSize());
  graph->GetYaxis()->SetNdivisions(108);

  graph->BayesDivide(select, total);
  graph->Draw("ap");
  graph->GetYaxis()->SetRangeUser(0, 1); 
  canvas->Modified();

  TString output_file = fOutputFile;
  output_file.ReplaceAll(".","_"+yaxis+".");
  output_file.ReplaceAll(".","_"+name+".");
  canvas->SaveAs(output_file);

}

// Plot a 2D histogram
void Plot2D(TH2D* hist, TString name, TString xaxis, TString yaxis){

  TCanvas *canvas = new TCanvas(name, "", 900, 600);
  canvas->SetFrameLineWidth(4.);
  canvas->SetLineWidth(4.);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->SetBottomMargin(0.16);
  canvas->SetLeftMargin(0.16);
  canvas->SetRightMargin(0.16);

  hist->GetXaxis()->SetTitle(xaxis);
  hist->GetYaxis()->SetTitle(yaxis);
  // X axis config
  hist->GetXaxis()->SetTitleOffset(1.1);
  hist->GetXaxis()->SetTickLength(0.04);
  hist->GetXaxis()->SetTitleSize(1.1*hist->GetXaxis()->GetTitleSize());
  hist->GetXaxis()->SetNdivisions(108);
  // Y axis config
  hist->GetYaxis()->SetTitleOffset(0.95);
  hist->GetYaxis()->SetTickLength(0.015);
  hist->GetYaxis()->SetTitleSize(1.1*hist->GetYaxis()->GetTitleSize());
  hist->GetYaxis()->SetNdivisions(108);

  hist->Draw("colz");

  TString output_file = fOutputFile;
  output_file.ReplaceAll(".","_"+name+".");
  canvas->SaveAs(output_file);
}

// Create stacked histogram and legend from data
std::pair<THStack*, TLegend*> StackHist1D(std::map<std::string, std::vector<std::vector<double>>> data, TString name, TString title, std::vector<std::vector<double>> bin_edges, int i, int j = -1, int bin_j = -1, int k = -1, int bin_k = -1){

  double edges_array[bin_edges[i].size()];
  std::copy(bin_edges[i].begin(), bin_edges[i].end(), edges_array);
  THStack *hstack = new THStack(name, title);
  TLegend *legend = new TLegend(0.14, 0., 0.94, 0.06);

  int index = 0;
  for(auto const& dat: data){
    TH1D* hist = new TH1D(name+dat.first.c_str(), title, bin_edges[i].size()-1, edges_array);
    for(size_t n = 0; n < dat.second.size(); n++){
      if(j==-1 && k == -1){
        hist->Fill(dat.second[n][i]);
      }
      else if(dat.second[n][j] >= bin_edges[j][bin_j] && dat.second[n][j] < bin_edges[j][bin_j+1]){
        if(k == -1){
          hist->Fill(dat.second[n][i]);
        }
        else{
          if(dat.second[n][k] >= bin_edges[k][bin_k] && dat.second[n][k] < bin_edges[k][bin_k+1]){
            hist->Fill(dat.second[n][i]);
          }
        }
      }
    }
    hist->Scale(fPotScaleFac);
    // If plotting cross section convert from rate
    if(fPlotXSec){
      double width = 1;
      if(j != -1) width = width * (bin_edges[j][bin_j+1] - bin_edges[j][bin_j]);
      if(k != -1) width = width * (bin_edges[k][bin_k+1] - bin_edges[k][bin_k]);
      double xsec_scale = 1e38/(width * fFlux * fTargets);
      hist->Scale(xsec_scale, "width");
    }
    // Else if max error used divide each bin by width
    else if (fMaxError > 0 || fBinEdges[i].size()>1){
      hist->Scale(1, "width");
    }
    hist->SetFillColor(fCols[index]);
    hist->SetLineColor(fCols[index]);
    if(!fPlotFilled){
      hist->SetFillColor(0);
      hist->SetLineWidth(3);
    }
    hstack->Add(hist);
    legend->AddEntry(hist, dat.first.c_str(), "lf");
    index++;
  }

  return std::make_pair(hstack, legend);
}

// Get the total (unstacked) histogram
TH1D* GetTotalHist(std::vector<std::vector<double>> data, TString name, std::vector<std::vector<double>> bin_edges, int i, int j = -1, int bin_j = -1, int k = -1, int bin_k = -1){
    
  double edges_array[bin_edges[i].size()];
  std::copy(bin_edges[i].begin(), bin_edges[i].end(), edges_array);
  TH1D *total_hist = new TH1D("total"+name, "hist", bin_edges[i].size()-1, edges_array);

  for (int n = 0; n < data.size(); n++){
    if(j == -1 && k == -1){
      total_hist->Fill(data[n][i]);
    }
    else if(data[n][j] >= bin_edges[j][bin_j] && data[n][j] < bin_edges[j][bin_j+1]){
      if(k == -1)
        total_hist->Fill(data[n][i]);
      else if(data[n][k] >= bin_edges[k][bin_k] && data[n][k] < bin_edges[k][bin_k+1])
        total_hist->Fill(data[n][i]);
    }
  }
  // Include scale factor bin by bin as ->Scale() won't change errors
  for(size_t n = 0; n <= total_hist->GetNbinsX(); n++){
    total_hist->SetBinContent(n, total_hist->GetBinContent(n)*fPotScaleFac);
  }
  // If plotting cross section convert from rate
  if(fPlotXSec){
    double width = 1;
    if(j != -1) width = width * (bin_edges[j][bin_j+1] - bin_edges[j][bin_j]);
    if(k != -1) width = width * (bin_edges[k][bin_k+1] - bin_edges[k][bin_k]);
    double xsec_scale = 1e38/(width * fFlux * fTargets);
    total_hist->Scale(xsec_scale, "width");
  }
  // Else if max error used divide each bin by width
  else if (fMaxError > 0 || fBinEdges[i].size()>1){
    total_hist->Scale(1, "width");
  }
  return total_hist;
}


// Get the percentage statistical error per bin
TH1D* GetErrorBand(TH1D* total_hist, TString name, std::vector<std::vector<double>> bin_edges, int i){

  double edges_array[bin_edges[i].size()];
  std::copy(bin_edges[i].begin(), bin_edges[i].end(), edges_array);
  TH1D *error_band = new TH1D("error_band"+name, "", bin_edges[i].size()-1, edges_array);

  // Set the bin errors on seperate plot
  for (int n = 1; n <= total_hist->GetNbinsX(); n++){
   error_band->SetBinContent(n, 0);
   if (total_hist->GetBinContent(n) > 0)
     error_band->SetBinContent(n, 100*total_hist->GetBinError(n)/total_hist->GetBinContent(n));
  }
  return error_band;

}

// Rebin to the maximum bin error for 2D histograms
std::vector<double> ChangeBinning2D(std::vector<std::vector<double>> data, std::vector<std::vector<double>> bin_edges, int i, int j, int bin_j, int k = -1, int bin_k = -1){

  double edges_array[bin_edges[i].size()];
  std::copy(bin_edges[i].begin(), bin_edges[i].end(), edges_array);
  TH1D *temp_hist = new TH1D("temp_hist", "", bin_edges[i].size()-1, edges_array);
  // Loop over data
  for(auto const& dat : data){
    // Fill temporary histogram
    if(dat[j] >= bin_edges[j][bin_j] && dat[j] < bin_edges[j][bin_j+1]){
      if(k == -1){
        temp_hist->Fill(dat[i]);
      }
      else if(dat[k] >= bin_edges[k][bin_k] && dat[k] < bin_edges[k][bin_k+1]){
        temp_hist->Fill(dat[i]);
      }
    }
  }
  // Include scale factor bin by bin as Scale() won't change errors
  for(size_t n = 1; n <= temp_hist->GetNbinsX(); n++){
    temp_hist->SetBinContent(n, temp_hist->GetBinContent(n)*fPotScaleFac);
  }
  // Change the binning so that all bin errors below maximum
  std::vector<double> bin_edges_new = ChangeBinning(temp_hist, bin_edges[i][bin_edges[i].size()-1]);
  delete temp_hist;
  return bin_edges_new;

}

// Make the efficiency and purity plots for up to 3 variables
void PlotEffPur(std::vector<Interaction> interactions, TString name, Titles titles, std::vector<std::vector<double>> bin_edges, int i, int j = -1, int bin_j = -1, int k = -1, int bin_k = -1){

  double edges_array[bin_edges[i].size()];
  std::copy(bin_edges[i].begin(), bin_edges[i].end(), edges_array);

  TH1D *eff_numerator = new TH1D("eff_numerator", "", bin_edges[i].size()-1, edges_array);
  TH1D *eff_denom = new TH1D("eff_denom", "", bin_edges[i].size()-1, edges_array);
  TH1D *pur_numerator = new TH1D("pur_numerator", "", bin_edges[i].size()-1, edges_array);
  TH1D *pur_denom = new TH1D("pur_denom", "", bin_edges[i].size()-1, edges_array);

  for (auto const& in : interactions){
    // Denominator of efficiency plot is all interactions that are selected in truth
    if(in.true_selected){ 
      if(j == -1 && k == -1){ 
        eff_denom->Fill(in.true_variables[i]);
      }
      else if(in.true_variables[j] >= bin_edges[j][bin_j] && in.true_variables[j] < bin_edges[j][bin_j+1]){
        if(k == -1){
          eff_denom->Fill(in.true_variables[i]);
        }
        else if(in.true_variables[k] >= bin_edges[k][bin_k] && in.true_variables[k] < bin_edges[k][bin_k+1]){
          eff_denom->Fill(in.true_variables[i]);
        }
      }
    }
    // Denominator of purity plot is all interaction that are selected after reconstruction
    if(in.selected){ 
      if(j == -1 && k == -1){ 
        pur_denom->Fill(in.variables[i]);
      }
      else if(in.variables[j] >= bin_edges[j][bin_j] && in.variables[j] < bin_edges[j][bin_j+1]){
        if(k == -1){
          pur_denom->Fill(in.variables[i]);
        }
        else if(in.variables[k] >= bin_edges[k][bin_k] && in.variables[k] < bin_edges[k][bin_k+1]){
          pur_denom->Fill(in.variables[i]);
        }
      }
    }
    // Numerator of efficiency and purity plots is all interactions that are selected in both truth and reco
    if(in.selected && in.true_selected){ 
      if(j == -1 && k == -1){
        eff_numerator->Fill(in.true_variables[i]);
        pur_numerator->Fill(in.variables[i]);
      }
      else{
        if(in.true_variables[j] >= bin_edges[j][bin_j] && in.true_variables[j] < bin_edges[j][bin_j+1]){
          if(k == -1){
            eff_numerator->Fill(in.true_variables[i]);
          }
          else if(in.true_variables[k] >= bin_edges[k][bin_k] && in.true_variables[k] < bin_edges[k][bin_k+1]){
            eff_numerator->Fill(in.true_variables[i]);
          }
        }
        if(in.variables[j] >= bin_edges[j][bin_j] && in.variables[j] < bin_edges[j][bin_j+1]){
          if(k == -1){
            pur_numerator->Fill(in.variables[i]);
          }
          else if(in.variables[k] >= bin_edges[k][bin_k] && in.variables[k] < bin_edges[k][bin_k+1]){
            pur_numerator->Fill(in.variables[i]);
          }
        }
      }
    }
  }

  // Efficiency: selected/total in true
  PlotEfficiency(eff_numerator, eff_denom, name, titles.names[i]+"^{true} ["+titles.units[i]+"]", "Efficiency");
  // Purity: correct selected/total selected
  PlotEfficiency(pur_numerator, pur_denom, name, titles.names[i]+"^{reco} ["+titles.units[i]+"]", "Purity");

  if(eff_numerator) delete eff_numerator; 
  if(eff_denom) delete eff_denom;
  if(pur_numerator) delete pur_numerator;
  if(pur_denom) delete pur_denom;

}

// Main
void PhysicsBookPlots(){

  // Get the configuration file
  std::cout<<"Reading the config file...\n";
  std::string input_file = "config.txt";
  Configure(input_file);
  GetMetaData();
  std::cout<<"...Finished.\n";

  // Get from configuration
  std::cout<<"Getting labels...\n";
  Titles titles = GetTitles();
  SetStyle();
  std::cout<<"...Finished.\n";

  // Get the plotting variable
  std::cout<<"Reading from the tree...\n";
  std::vector<Interaction> interactions = ReadData();
  std::map<std::string, std::vector<std::vector<double>>> stack_data;
  std::vector<std::vector<double>> total_data;
  for(auto const& in : interactions){
    if(!in.selected) continue;
    total_data.push_back(in.variables);
    if(!fPlotStacked){
      stack_data["all"].push_back(in.variables);
    }
    else if(fStackBy == "fsi"){
      stack_data[in.fsi].push_back(in.variables);
    }
    else if(fStackBy == "int"){
      stack_data[in.int_type].push_back(in.variables);
    }
    else if(fStackBy == "nu"){
      stack_data[in.nu_type].push_back(in.variables);
    }
  }
  std::cout<<"...Finished.\n";

  // Get the binning
  std::cout<<"Getting the correct binning...\n";
  std::vector<std::vector<double>> bin_edges = GetBinning(total_data);
  std::cout<<"...Finished.\n";

  // Calculate how many 1D histograms the choice of variables and binning would produce
  int n_hists = 0;
  if(fPlotVariables.size() == 1) n_hists = 1;
  else if(fPlotVariables.size() == 2){
    n_hists = bin_edges[0].size() + bin_edges[1].size();
  }
  else if(fPlotVariables.size() == 3){
    n_hists = 3 + (bin_edges[0].size()-1)*(bin_edges[1].size()-1) + (bin_edges[0].size()-1)*(bin_edges[2].size()-1) + (bin_edges[1].size()-1)*(bin_edges[2].size()-1);
  }

  // Ask user if they want to make that many histograms
  std::string response = "y";
  if(n_hists>10){
    std::cout<<"This will produce "<<n_hists<<" histograms, continue (y/n)? ";
    std::cin>>response;
  }
  if(response=="n") exit(1);

  // Loop over the number of variables - this is how many sets of 1D hists we will have
  std::cout<<"Making the plots...\n";
  for(size_t d_i = 0; d_i < fPlotVariables.size(); d_i++){

    // Don't do anything if we don't want to see the plots for this variable
    if(!fShowPlots[d_i]) continue;

    double edges_array[bin_edges[d_i].size()];
    std::copy(bin_edges[d_i].begin(), bin_edges[d_i].end(), edges_array);

    // Get the file name and title of the histogram
    TString name_1D = fPlotVariables[d_i];
    TString title_1D = titles.hist_titles[d_i];

    // Get the statistical errors per bin
    TH1D* total_hist = GetTotalHist(total_data, name_1D, bin_edges, d_i);
    //std::cout<<"Total cross section = "<<total_data.size()*fPotScaleFac*1e38/(fFlux*fTargets)<<" 10^{-38}\n";
    TH1D* error_band = GetErrorBand(total_hist, name_1D, bin_edges, d_i);
    // Create a total 1D stacked histogram for each of the variables
    std::pair<THStack*, TLegend*> stack = StackHist1D(stack_data, name_1D, title_1D, bin_edges, d_i);

    // Draw the plots
    if(fShowStatError) Plot1DWithErrors(stack.first, stack.second, error_band, titles, total_hist, d_i);
    else Plot1D(stack.first, stack.second, titles, total_hist, d_i);

    // Plot efficiency and purity if option selected and reconstruction selected
    if(fPlotEffPur && fStage == "reco"){
      PlotEffPur(interactions, name_1D, titles, bin_edges, d_i);
    }

    // Loop over the other variables - this is how many sets of 2D hists we will have
    for(size_t d_j = 0; d_j < fPlotVariables.size(); d_j++){
      if(d_j == d_i) continue;

      // FIXME bother with 2D slices of 3D hist?
      double edges_array2[bin_edges[d_j].size()];
      std::copy(bin_edges[d_j].begin(), bin_edges[d_j].end(), edges_array2);

      // Create a total 2D histogram for each combination of variables
      TH2D* hist_2D = new TH2D(name_1D+fPlotVariables[d_j], "", bin_edges[d_i].size()-1, edges_array, bin_edges[d_j].size()-1, edges_array2);
      for (int n = 0; n < total_data.size(); n++){
        hist_2D->Fill(total_data[n][d_i], total_data[n][d_j]);
      }
      Plot2D(hist_2D, fPlotVariables[d_i]+"_"+fPlotVariables[d_j], titles.names[d_i]+" ["+titles.units[d_i]+"]", titles.names[d_j]+" ["+titles.units[d_j]+"]");

      // Loop over the bins for variable 2
      for(size_t bin_j = 0; bin_j < bin_edges[d_j].size()-1; bin_j++){
        // Only make these plots for 2 variables
        if(fPlotVariables.size() != 2) continue;

        // Rebin so every bin below maximum error if set
        std::vector<std::vector<double>> bin_edges_copy = bin_edges;
        if(fMaxError>0){
          std::vector<double> bin_edges_new = ChangeBinning2D(total_data, bin_edges, d_i, d_j, bin_j);
          bin_edges_copy[d_i] = bin_edges_new;
        }

        // Get the file name and title of the histogram
        TString name_2D = fPlotVariables[d_i] +"_"
                          + fPlotVariables[d_j] +"_"+ Form("%.1f", bin_edges_copy[d_j][bin_j]) +"_"+ Form("%.1f", bin_edges_copy[d_j][bin_j+1]);
        TString title_2D = titles.hist_titles[d_j] 
                           +": ["+ Form("%.2f", bin_edges_copy[d_j][bin_j]) +", "+ Form("%.2f", bin_edges_copy[d_j][bin_j+1]) +"]";

        // Get the statistical errors per bin
        TH1D *total_hist_2D = GetTotalHist(total_data, name_2D, bin_edges_copy, d_i, d_j, bin_j);
        TH1D *error_band_2D = GetErrorBand(total_hist_2D, name_2D, bin_edges_copy, d_i);
        if(error_band_2D->Integral(0, error_band_2D->GetNbinsX()) == 0) continue;

        // Create a 1D stacked histogram for each of the bins
        std::pair<THStack*, TLegend*> stack_2D = StackHist1D(stack_data, name_2D, title_2D, bin_edges_copy, d_i, d_j, bin_j);

        // Draw the plots
        if(fShowStatError) Plot1DWithErrors(stack_2D.first, stack_2D.second, error_band_2D, titles, total_hist_2D, d_i, d_j);
        else Plot1D(stack_2D.first, stack_2D.second, titles, total_hist_2D, d_i, d_j);

        // Plot efficiency and purity if option selected and reconstruction selected
        if(fPlotEffPur && fStage == "reco"){
          PlotEffPur(interactions, name_2D, titles, bin_edges_copy, d_i, d_j, bin_j);
        }
      }

      // Loop over the other variables - this is how many sets of 3D hists we will have
      for(size_t d_k = 0; d_k < fPlotVariables.size(); d_k++){
        if(d_k == d_i || d_k <= d_j) continue;
        // Only produce these plots for 3 variables
        if(fPlotVariables.size() != 3) continue;

        // Loop over the bins for variable 3
        for(size_t bin_j = 0; bin_j < bin_edges[d_j].size()-1; bin_j++){

          // If maximum error set, get new k bins based on the j binning
          std::vector<std::vector<double>> bin_edges_copy = bin_edges;
          if(fMaxError>0){
            std::vector<double> bin_edges_new = ChangeBinning2D(total_data, bin_edges, d_k, d_j, bin_j);
            bin_edges_copy[d_k] = bin_edges_new;
          }

          // Loop over the bins for variable 2
          for(size_t bin_k = 0; bin_k < bin_edges_copy[d_k].size()-1; bin_k++){

            // If maximum error set, get new i bins based on the j and k binning
            if(fMaxError>0){
              std::vector<double> bin_edges_new = ChangeBinning2D(total_data, bin_edges_copy, d_i, d_j, bin_j, d_k, bin_k);
              bin_edges_copy[d_i] = bin_edges_new;
            }

            // Get the file name and title of the histogram
            TString name_3D = fPlotVariables[d_i] + "_"
                              + fPlotVariables[d_j] +"_"+ Form("%.1f", bin_edges_copy[d_j][bin_j]) +"_"+ Form("%.1f", bin_edges_copy[d_j][bin_j+1]) +"_"
                              + fPlotVariables[d_k] +"_"+ Form("%.1f", bin_edges_copy[d_k][bin_k]) +"_"+ Form("%.1f", bin_edges_copy[d_k][bin_k+1]);
            TString title_3D = titles.hist_titles[d_j] 
                               +": ["+ Form("%.2f", bin_edges_copy[d_j][bin_j]) +", "+ Form("%.2f", bin_edges_copy[d_j][bin_j+1]) +"] "
                               + titles.hist_titles[d_k] 
                               +": ["+ Form("%.2f", bin_edges_copy[d_k][bin_k]) +", "+ Form("%.2f", bin_edges_copy[d_k][bin_k+1]) +"]";

            // Get the statistical errors per bin
            TH1D *total_hist_3D = GetTotalHist(total_data, name_3D, bin_edges_copy, d_i, d_j, bin_j, d_k, bin_k);
            TH1D *error_band_3D = GetErrorBand(total_hist_3D, name_3D, bin_edges_copy, d_i);
            // Don't draw empty histograms
            if(error_band_3D->Integral(0, error_band_3D->GetNbinsX()) == 0) continue;

            // Create a 1D stacked histogram for each of the bins
            std::pair<THStack*, TLegend*> stack_3D = StackHist1D(stack_data, name_3D, title_3D, bin_edges_copy, d_i, d_j, bin_j, d_k, bin_k);

            // Draw the plots
            if(fShowStatError) Plot1DWithErrors(stack_3D.first, stack_3D.second, error_band_3D, titles, total_hist_3D, d_i, d_j, d_k);
            else Plot1D(stack_3D.first, stack_3D.second, titles, total_hist_3D, d_i, d_j, d_k);

            // Plot efficiency and purity if option selected and reconstruction selected
            if(fPlotEffPur && fStage == "reco"){
              PlotEffPur(interactions, name_3D, titles, bin_edges_copy, d_i, d_j, bin_j, d_k, bin_k);
            }
          }
        }
      }

    }
  }
  std::cout<<"...Finished.\n";

/*
  if(fPlotResponse){
    // Nj = number of events generated in bin j
    // Nij = number of events observed in bin i, generated in bin j
    // N0j = number of events not observed, generated in bin j
    // Response matrix: Aij = Nij/Nj
    // FIXME I don't think this is how you do it... units need to be bins (all same size)
    TH2D *response = new TH2D("response", "", bin_edges.size()-1, edges_array, bin_edges.size()-1, edges_array);
    TH1D *Nj = new TH1D("Nj", "", bin_edges.size()-1, edges_array);
    for(auto const& in : interactions){
      if(!in.true_selected) continue;
      Nj->Fill(in.true_variables[0]);
      if(in.selected) response->Fill(in.true_variables[0], in.variables[0]);
      else response->Fill(in.true_variables[0], -99999);
    }
    for(size_t i = 0; i <= response->GetNbinsX(); i++){
      for(size_t j = 1; j <= response->GetNbinsY(); j++){
        if(Nj->GetBinContent(j) != 0)
          response->SetBinContent(i, j, response->GetBinContent(i, j)/Nj->GetBinContent(j));
        else
          response->SetBinContent(i, j, 0);
      }
    }
    Plot2D(response, "Response", "True "+titles.units[0], "Reco "+titles.units[0]);
  }
  */

}

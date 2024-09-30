#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"

#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

void BasicVariables()
{
  const TString saveDir = "/exp/sbnd/data/users/hlay/ncpizero/plots/pizeroreco/basicvariables/tmp";
  const bool save = true;

  using namespace std;
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("pizeroana/pizeros");
  tree->Add("/exp/sbnd/data/users/hlay/ncpizero/pizeroana_sbnd.root");

  std::vector<Cut> cuts = { {"all", "", "", kMagenta+2},
                            {"neutrino", "mct_origin==1", "", kGreen+2},
                            {"cosmic", "mct_origin==2", "", kRed-2},
                            {"neutrino_av", "mct_origin==1 && mc_av", "", kGreen+2},
                            {"neutrino_fv", "mct_origin==1 && mc_fv", "", kGreen+2},
                            {"cosmic_av", "mct_origin==2 && mc_av", "", kRed-2},
                            {"cosmic_fv", "mct_origin==2 && mc_fv", "", kRed-2},
  };

  std::vector<Plot> plots = { {"mctruth_origin", "mct_origin", ";Origin;#pi^{0}s",
                               3, 0, 3, kBlack, false, "", true,
                               {"Unknown", "Neutrino Beam", "Cosmic Ray"}},
                              {"pizero_vx", "mc_vx", ";x position [cm];#pi^{0}s",
                               50, -400, 400},
                              {"pizero_vx_detector", "mc_vx", ";x position [cm];#pi^{0}s",
                               55, -220, 220},
                              {"pizero_vy", "mc_vy", ";y position [cm];#pi^{0}s",
                               50, -400, 400},
                              {"pizero_vy_detector", "mc_vy", ";y position [cm];#pi^{0}s",
                               55, -220, 220},
                              {"pizero_vz", "mc_vz", ";z position [cm];#pi^{0}s",
                               50, -250, 750},
                              {"pizero_vz_detector", "mc_vz", ";z position [cm];#pi^{0}s",
                               52, -10, 510},
                              {"pizero_vt", "mc_vt", ";time [ns];#pi^{0}s",
                               180, -1.8e6, 1.6e6},
                              {"pizero_vt_beam", "mc_vt", ";time [ns];#pi^{0}s",
                               60, -100, 1700},
                              {"pizero_theta", "mc_vtheta", ";#theta(#circ);#pi^{0}s", 45, 0, 180},
                              {"pizero_phi", "mc_vphi", ";#phi(#circ);#pi^{0}s", 45, -180, 180},
                              {"pizero_cos_theta", "cos(TMath::DegToRad() * mc_vtheta)",
                               ";cos(#theta);#pi^{0}s",
                               50, -1, 1},
                              {"pizero_cos_phi", "cos(TMath::DegToRad() * mc_vphi)",
                               ";cos(#phi);#pi^{0}s",
                               50, -1, 1},
                              {"pizero_energy", "mc_ve", ";E (GeV);#pi^{0}s",
                               100, 0, 2.5},
                              {"pizero_decay_open_angle", "mc_open_angle",
                               ";#theta_{#gamma#gamma} (#circ);#pi^{0}s",
                               45, 0, 180},
                              {"pizero_decay_cos_open_angle",
                               "cos(TMath::DegToRad() * mc_open_angle)",
                               ";cos(#theta_{#gamma#gamma});#pi^{0}s",
                               50, -1, 1},
                              {"leading_gamma_energy", "mc_gamma1_ve", ";E_{#gamma1} (GeV);#pi^{0}s",
                               100, 0, 2.5},
                              {"leading_gamma_theta", "mc_gamma1_vtheta",
                               ";#theta_{#gamma1};#pi^{0}s", 45, 0, 180,
                               kBlack, false, "mc_gamma1_pdg!=-999999"},
                              {"leading_gamma_phi", "mc_gamma1_vphi", ";#phi_{#gamma1};#pi^{0}s",
                               45, -180, 180, kBlack, false, "mc_gamma1_pdg!=-999999"},
                              {"leading_gamma_cos_theta",
                               "cos(TMath::DegToRad() * mc_gamma1_vtheta)",
                               ";cos(#theta_{#gamma1});#pi^{0}s",
                               50, -1, 1, kBlack, false, "mc_gamma1_pdg!=-999999"},
                              {"leading_gamma_cos_phi", "cos(TMath::DegToRad() * mc_gamma1_vphi)",
                               ";cos(#phi_{#gamma1});#pi^{0}s",
                               50, -1, 1, kBlack, false, "mc_gamma1_pdg!=-999999"},
                              {"subleading_gamma_energy", "mc_gamma2_ve",
                               ";E_{#gamma2} (GeV);#pi^{0}s",
                               100, 0, 2.5},
                              {"subleading_gamma_theta", "mc_gamma2_vtheta",
                               ";#theta_{#gamma2};#pi^{0}s", 45, 0, 180,
                               kBlack, false, "mc_gamma2_pdg!=-999999"},
                              {"subleading_gamma_phi", "mc_gamma2_vphi", ";#phi_{#gamma2};#pi^{0}s",
                               45, -180, 180, kBlack, false, "mc_gamma2_pdg!=-999999"},
                              {"subleading_gamma_cos_theta",
                               "cos(TMath::DegToRad() * mc_gamma2_vtheta)",
                               ";cos(#theta_{#gamma2});#pi^{0}s",
                               50, -1, 1, kBlack, false, "mc_gamma2_pdg!=-999999"},
                              {"subleading_gamma_cos_phi", "cos(TMath::DegToRad() * mc_gamma2_vphi)",
                               ";cos(#phi_{#gamma2});#pi^{0}s",
                               50, -1, 1, kBlack, false, "mc_gamma2_pdg!=-999999"},
  };
  
  std::vector<TwoDPlot> twodplots = {
    {"pizero_energy_versus_decay_open_angle", "mc_ve:mc_open_angle",
     ";#theta_{#gamma#gamma} (#circ);E_{#pi^{0}} (GeV);#pi^{0}s",
     45, 0, 180, 100, 0, 2.5, false, "mc_gamma1_pdg==22 && mc_gamma2_pdg==22" },
    {"pizero_energy_versus_decay_open_angle_log", "mc_ve:mc_open_angle",
     ";#theta_{#gamma#gamma} (#circ);E_{#pi^{0}} (GeV);#pi^{0}s",
     45, 0, 180, 100, 0, 2.5, true, "mc_gamma1_pdg==22 && mc_gamma2_pdg==22" },
  };
  
  for(auto const & cut : cuts)
    {
      gSystem->Exec("mkdir -p " + saveDir + "/" + cut.name);

      for(auto &plot : plots)
        {
          TCanvas *canvas = new TCanvas("c_" + plot.name + "_" + cut.name,
                                        "c_" + plot.name + "_" + cut.name);
          canvas->cd();

          plot.colour = cut.colour;

          MakePlot(canvas, tree, plot, cut);

          if(save)
            {
              canvas->SaveAs(saveDir + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".png");
              canvas->SaveAs(saveDir + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".pdf");
            }
          delete canvas;
        }

      for(auto const &plot : twodplots)
        {
          gStyle->SetNdivisions(505, "x");
          gStyle->SetPalette(kBlueRedYellow);

          TCanvas *canvas = new TCanvas("c_" + plot.name + "_" + cut.name,
                                        "c_" + plot.name + "_" + cut.name);
          canvas->cd();

          MakeTwoDPlot(canvas, tree, plot, cut);

          if(save)
            {
              canvas->SaveAs(saveDir + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".png");
              canvas->SaveAs(saveDir + "/" + cut.name + "/" + plot.name + "_" + cut.name + ".pdf");
            }
          delete canvas;
        }
    }
}

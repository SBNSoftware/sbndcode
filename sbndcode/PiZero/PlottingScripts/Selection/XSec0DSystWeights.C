#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Plots.h"
#include "Selections.h"
#include "Common.C"
#include "WeightNames.h"
#include "Enumerate.h"
#include "LatexHeaders.h"
#include "/exp/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv8/integrated_flux/FluxMap.h"

std::vector<int> *nu_event_type_incl = 0, *nu_event_type_0p0pi = 0, *nu_event_type_Np0pi = 0,
  *slc_true_event_type_incl = 0, *slc_true_event_type_0p0pi = 0, *slc_true_event_type_Np0pi = 0,
  *slc_n_primary_razzled_muons = 0, *slc_n_primary_razzled_photons = 0, *slc_n_primary_razzled_pions_thresh = 0,
  *slc_n_primary_razzled_protons_thresh = 0;
std::vector<bool> *slc_is_clear_cosmic = 0, *slc_is_fv = 0, *slc_best_pzc_good_kinematics = 0, *slc_all_other_trks_contained = 0;
std::vector<float> *slc_crumbs_score = 0;
std::vector<double> *slc_opt0_fracPE = 0, *slc_opt0_score = 0;

std::vector<std::vector<std::vector<double>>*> nu_parameter_weights, slc_parameter_weights;

void InitialiseTree(TChain *tree, std::vector<std::string> &weight_names);

void EvaluateTree(TChain *tree, const double &scaling, std::vector<TH1D*> &nominalSelectedHists, std::vector<TH1D*> &nominalTrueSignalHists,
                  std::vector<TH1D*> &nominalSelectedBackgroundHists, std::vector<std::vector<std::vector<TH1D*>>> &univTrueSignalHists,
                  std::vector<std::vector<std::vector<TH1D*>>> &univSelectedBackgroundHists, std::vector<std::vector<std::vector<TH1D*>>> &univSelectedSignalHists,
                  const std::vector<Cut> &signals, std::vector<std::string> weight_names, const uint n_univs, const bool combine);

void XSec0DSystWeights(const TString productionVersion, const TString saveDirExt, const std::vector<Cut> &signals,
                       std::vector<std::string> weight_names, const uint n_univs, const bool combine)
{
  nu_parameter_weights.resize(weight_names.size(), 0);
  slc_parameter_weights.resize(weight_names.size(), 0);

  const TString saveDir = baseSaveDir + "/" + productionVersion + "/xsec_zero_d_syst_weights/" + saveDirExt;
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_rockbox.root";
  const TString intimeFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_intime.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxEvents = new TChain("ncpizeroana/events");
  rockboxEvents->Add(rockboxFile);
  TChain *rockboxSubruns = new TChain("ncpizeroana/subruns");
  rockboxSubruns->Add(rockboxFile);

  TChain *intimeEvents = new TChain("ncpizeroana/events");
  intimeEvents->Add(intimeFile);
  TChain *intimeSubruns = new TChain("ncpizeroana/subruns");
  intimeSubruns->Add(intimeFile);

  double rockboxScaling, intimeScaling;
  GetScaling(rockboxSubruns, intimeSubruns, rockboxScaling, intimeScaling);

  std::vector<TH1D*> nominalSelectedHists, nominalTrueSignalHists, nominalSelectedBackgroundHists, nominalXSecHists;
  std::vector<std::vector<std::vector<TH1D*>>> univTrueSignalHists, univSelectedBackgroundHists, univSelectedSignalHists, univXSecHists;

  for(auto&& [signal_i, signal] : enumerate(signals))
    {
      nominalSelectedHists.push_back(new TH1D(Form("nominalSelectedHist_%lu", signal_i), "", 1, 0, 1));
      nominalTrueSignalHists.push_back(new TH1D(Form("nominalTrueSignalHist_%lu", signal_i), "", 1, 0, 1));
      nominalSelectedBackgroundHists.push_back(new TH1D(Form("nominalSelectedBackgroundHist_%lu", signal_i), "", 1, 0, 1));
      nominalXSecHists.push_back(new TH1D(Form("nominalXSecHist_%lu", signal_i), Form("%s;;#sigma (cm^{2}/nucleon)", signal.printed_name.Data()), 1, 0, 1));
      
      univTrueSignalHists.push_back(std::vector<std::vector<TH1D*>>());
      univSelectedBackgroundHists.push_back(std::vector<std::vector<TH1D*>>());
      univSelectedSignalHists.push_back(std::vector<std::vector<TH1D*>>());
      univXSecHists.push_back(std::vector<std::vector<TH1D*>>());

      for(auto&& [weight_i, weight] : enumerate(weight_names))
        {
          univTrueSignalHists[signal_i].push_back(std::vector<TH1D*>());
          univSelectedBackgroundHists[signal_i].push_back(std::vector<TH1D*>());
          univSelectedSignalHists[signal_i].push_back(std::vector<TH1D*>());
          univXSecHists[signal_i].push_back(std::vector<TH1D*>());

          for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
            {
              univTrueSignalHists[signal_i][weight_i].push_back(new TH1D(Form("univTrueSignalHist_%lu_%lu_%u", signal_i, weight_i, univ_i), "", 1, 0, 1));
              univSelectedBackgroundHists[signal_i][weight_i].push_back(new TH1D(Form("univSelectedBackgroundHist_%lu_%lu_%u", signal_i, weight_i, univ_i), "", 1, 0, 1));
              univSelectedSignalHists[signal_i][weight_i].push_back(new TH1D(Form("univSelectedSignalHist_%lu_%lu_%u", signal_i, weight_i, univ_i), "", 1, 0, 1));
              univXSecHists[signal_i][weight_i].push_back(new TH1D(Form("univXSecHist_%lu_%lu_%u", signal_i, weight_i, univ_i), "", 1, 0, 1));
            }
        }

      if(combine)
        {
          univTrueSignalHists[signal_i].push_back(std::vector<TH1D*>());
          univSelectedBackgroundHists[signal_i].push_back(std::vector<TH1D*>());
          univSelectedSignalHists[signal_i].push_back(std::vector<TH1D*>());
          univXSecHists[signal_i].push_back(std::vector<TH1D*>());

          for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
            {
              univTrueSignalHists[signal_i][weight_names.size()].push_back(new TH1D(Form("univTrueSignalHist_%lu_%lu_%u", signal_i, weight_names.size(), univ_i), "", 1, 0, 1));
              univSelectedBackgroundHists[signal_i][weight_names.size()].push_back(new TH1D(Form("univSelectedBackgroundHist_%lu_%lu_%u", signal_i, weight_names.size(), univ_i), "", 1, 0, 1));
              univSelectedSignalHists[signal_i][weight_names.size()].push_back(new TH1D(Form("univSelectedSignalHist_%lu_%lu_%u", signal_i, weight_names.size(), univ_i), "", 1, 0, 1));
              univXSecHists[signal_i][weight_names.size()].push_back(new TH1D(Form("univXSecHist_%lu_%lu_%u", signal_i, weight_names.size(), univ_i), "", 1, 0, 1));
            }
        }

    }

  InitialiseTree(rockboxEvents, weight_names);
  EvaluateTree(rockboxEvents, rockboxScaling, nominalSelectedHists, nominalTrueSignalHists, nominalSelectedBackgroundHists,
               univTrueSignalHists, univSelectedBackgroundHists, univSelectedSignalHists, signals, weight_names, n_univs, combine);

  InitialiseTree(intimeEvents, weight_names);
  EvaluateTree(intimeEvents, intimeScaling, nominalSelectedHists, nominalTrueSignalHists, nominalSelectedBackgroundHists,
               univTrueSignalHists, univSelectedBackgroundHists, univSelectedSignalHists, signals, weight_names, n_univs, combine);


  if(combine)
    weight_names.push_back("all");

  ofstream texFile;
  texFile.open(saveDir + "/systematic_fractional_errors.tex");

  texFile << docStart 
          << tableStart
          << "\\hline\n"
          << "Systematic & \\multicolumn{2}{|c|}{NC1$\\pi^{0}$} & \\multicolumn{2}{|c|}{NC1$\\pi^{0}$0p0$\\pi^{\\pm}$} & \\multicolumn{2}{|c|}{NC1$\\pi^{0}$Np0$\\pi^{\\pm}$} \\\\ \\hline"
          << "& Resolution (\\%) & Bias (\\%) & Resolution (\\%) & Bias (\\%) & Resolution (\\%) & Bias (\\%)"
          << "\\\\ \\hline" << std::endl;
  
  for(auto&& [ weight_i, weight ] : enumerate(weight_names))
    {
      texFile << "\\url{" << weight << "}";
      gSystem->Exec("mkdir -p " + saveDir + "/" + weight.c_str());

      TCanvas *canvas = new TCanvas(Form("canvas_%lu", weight_i), Form("canvas_%lu", weight_i));
      canvas->cd();
      canvas->Divide(3, 1);

      for(auto&& [ signal_i, signal ] : enumerate(signals))
        {
          canvas->cd(signal_i + 1);
          gPad->SetBottomMargin(0.05);
          gPad->SetTopMargin(0.1);
          gPad->SetLeftMargin(0.25);
          gPad->SetRightMargin(0.05);

          const double nomSig  = nominalSelectedHists[signal_i]->GetBinContent(1) - nominalSelectedBackgroundHists[signal_i]->GetBinContent(1);
          const double nomEff  = nomSig / nominalTrueSignalHists[signal_i]->GetBinContent(1);
          const double nomXSec = nomSig / (nomEff * nTargets * intFlux);

          const double selEntries  = nominalSelectedHists[signal_i]->GetEntries();
          const double backEntries = nominalSelectedBackgroundHists[signal_i]->GetEntries();
          const double statError   = ((std::sqrt(selEntries) + std::sqrt(backEntries)) / (selEntries - backEntries)) * nomXSec;

          nominalXSecHists[signal_i]->SetBinContent(1, nomXSec);
          nominalXSecHists[signal_i]->SetBinError(1, statError);

          nominalXSecHists[signal_i]->GetYaxis()->SetTitleOffset(1.9);
          nominalXSecHists[signal_i]->GetXaxis()->SetLabelSize(0);
          nominalXSecHists[signal_i]->GetXaxis()->SetLabelOffset(999);
          nominalXSecHists[signal_i]->Draw("histe][");
          gPad->Update();
          TPaveText* title = (TPaveText*)gPad->FindObject("title");
          title->SetY1NDC(0.92);
          title->SetY2NDC(1);
          title->SetX1NDC(0.45);
          title->SetX2NDC(.8);
          gPad->Modified();
          gPad->Update();
          nominalXSecHists[signal_i]->SetMinimum(0);
          nominalXSecHists[signal_i]->SetMaximum(5e-40);
        }

      canvas->SaveAs(saveDir + "/" + weight.c_str() + "/nominal.png");
      canvas->SaveAs(saveDir + "/" + weight.c_str() + "/nominal.pdf");

      for(auto&& [ signal_i, signal ] : enumerate(signals))
        {
          canvas->cd(signal_i + 1);
          nominalXSecHists[signal_i]->Draw("hist][");
          gPad->Update();
          TPaveText* title = (TPaveText*)gPad->FindObject("title");
          title->SetY1NDC(0.92);
          title->SetY2NDC(1);
          title->SetX1NDC(0.45);
          title->SetX2NDC(.8);
          gPad->Modified();
          gPad->Update();

          for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
            {
              // Option to use either background subtraction (commented out) or purity correction (currently in use) to calculate the xsec
              //              const double univSig  = nominalSelectedHists[signal_i]->GetBinContent(1) - univSelectedBackgroundHists[signal_i][weight_i][univ_i]->GetBinContent(1);
              //              const double univEff  = univSig / univTrueSignalHists[signal_i][weight_i][univ_i]->GetBinContent(1);

              const double univPur  = univSelectedSignalHists[signal_i][weight_i][univ_i]->GetBinContent(1) /
                (univSelectedSignalHists[signal_i][weight_i][univ_i]->GetBinContent(1) + univSelectedBackgroundHists[signal_i][weight_i][univ_i]->GetBinContent(1));
              const double univSig  = nominalSelectedHists[signal_i]->GetBinContent(1) * univPur;
              const double univEff  = univSelectedSignalHists[signal_i][weight_i][univ_i]->GetBinContent(1) / univTrueSignalHists[signal_i][weight_i][univ_i]->GetBinContent(1);

              const double flux     = saveDirExt == "flux" ? univsIntegratedFluxMap.at(univ_i) : intFlux;
              const double univXSec = univSig / (univEff * nTargets * flux);

              univXSecHists[signal_i][weight_i][univ_i]->SetBinContent(1, univXSec);

              univXSecHists[signal_i][weight_i][univ_i]->SetLineWidth(1);
              univXSecHists[signal_i][weight_i][univ_i]->SetLineColor(kMagenta-10);
              univXSecHists[signal_i][weight_i][univ_i]->Draw("hist][same");
            }

          nominalXSecHists[signal_i]->Draw("hist][same");
        }

      canvas->cd(1);

      TLegend *univsLeg = new TLegend(0.3, 0.4, 0.7, 0.5);
      univsLeg->AddEntry(nominalXSecHists[0], "Nominal", "l");
      univsLeg->AddEntry(univXSecHists[0][weight_i][0], "Universes", "l");
      univsLeg->SetTextSize(.06);
      univsLeg->Draw();

      canvas->SaveAs(saveDir + "/" + weight.c_str() + "/univs.png");
      canvas->SaveAs(saveDir + "/" + weight.c_str() + "/univs.pdf");

      for(auto&& [ signal_i, signal ] : enumerate(signals))
        {
          canvas->cd(signal_i + 1);
          nominalXSecHists[signal_i]->Draw("hist][");
          gPad->Update();
          TPaveText* title = (TPaveText*)gPad->FindObject("title");
          title->SetY1NDC(0.92);
          title->SetY2NDC(1);
          title->SetX1NDC(0.45);
          title->SetX2NDC(.8);
          gPad->Modified();
          gPad->Update();

          TH1D* cvErrHist = new TH1D(Form("cvErrHist_%lu_%lu", signal_i, weight_i), "", 1, 0, 1);

          TCanvas *fitCanvas = new TCanvas(Form("fitCanvas_%lu_%lu", signal_i, weight_i), Form("fitCanvas_%lu_%lu", signal_i, weight_i));
          fitCanvas->cd();

          double minVal = std::numeric_limits<double>::max(), maxVal = std::numeric_limits<double>::lowest();

          for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
            {
              const double value = univXSecHists[signal_i][weight_i][univ_i]->GetBinContent(1);

              if(value < minVal)
                minVal = value;

              if(value > maxVal)
                maxVal = value;
            }

          double diff = maxVal - minVal;

          if(diff > 1e-44)
            {
              TH1D* fitHist = new TH1D(Form("fitHist_%lu_%lu", signal_i, weight_i), ";#sigma (cm^{2}/nucleon);Universes", 10, minVal - 0.05 * diff, maxVal + 0.05 * diff);

              for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
                fitHist->Fill(univXSecHists[signal_i][weight_i][univ_i]->GetBinContent(1));

              TF1* fGaus = new TF1("fGaus", "gaus", minVal - 0.05 * diff, maxVal + 0.05 * diff);
              fGaus->SetLineColor(kRed+2);

              TFitResultPtr fitResult = fitHist->Fit(fGaus, "S");

              cvErrHist->SetBinContent(1, fGaus->GetParameter("Mean"));
              cvErrHist->SetBinError(1, fGaus->GetParameter("Sigma"));

              if(!fitResult->IsValid())
                {
                  gSystem->Exec("mkdir -p " + saveDir + "/" + weight.c_str() + "/bad_fits");
                  fitCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/bad_fits/" + Form("%s.png", signal.name.Data()));
                  fitCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/bad_fits/" + Form("%s.pdf", signal.name.Data()));
                }
              else
                {
                  gSystem->Exec("mkdir -p " + saveDir + "/" + weight.c_str() + "/fits");
                  fitCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/fits/" + Form("%s.png", signal.name.Data()));
                  fitCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/fits/" + Form("%s.pdf", signal.name.Data()));
                }
            }
          else
            {
              cvErrHist->SetBinContent(1, maxVal);
              cvErrHist->SetBinError(1, 0);
            }

          canvas->cd(signal_i + 1);

          cvErrHist->SetMarkerStyle(1);
          cvErrHist->SetLineColor(kBlue+2);
          cvErrHist->Draw("hist][esame");

          texFile << Form(" & %.2f & %.2f",
                          cvErrHist->GetBinError(1) / nominalXSecHists[signal_i]->GetBinContent(1) * 100.,
                          100. * (cvErrHist->GetBinContent(1) - nominalXSecHists[signal_i]->GetBinContent(1)) / nominalXSecHists[signal_i]->GetBinContent(1));
        }

      texFile << "\\\\ \\hline";
      
      canvas->cd(1);

      TH1D* tmpHist = new TH1D(Form("tmpHist_%lu", weight_i), "", 1, 0, 1);
      tmpHist->SetLineColor(kBlue+2);

      TLegend *cvLeg = new TLegend(0.3, 0.4, 0.7, 0.5);
      cvLeg->AddEntry(nominalXSecHists[0], "Nominal", "l");
      cvLeg->AddEntry(tmpHist, "CV #pm 1#sigma", "l");
      cvLeg->SetTextSize(.06);
      cvLeg->Draw();

      canvas->SaveAs(saveDir + "/" + weight.c_str() + "/cverr.png");
      canvas->SaveAs(saveDir + "/" + weight.c_str() + "/cverr.pdf");
    }

  texFile << tableEnd << docEnd;
  texFile.close();
  gSystem->Exec("pdflatex -output-directory " + saveDir + " " + saveDir + "/systematic_fractional_errors.tex");
}

void InitialiseTree(TChain *tree, std::vector<std::string> &weight_names)
{
  tree->SetBranchStatus("*", 0);

  tree->SetBranchAddress("nu_event_type_incl", &nu_event_type_incl);
  tree->SetBranchAddress("nu_event_type_0p0pi", &nu_event_type_0p0pi);
  tree->SetBranchAddress("nu_event_type_Np0pi", &nu_event_type_Np0pi);
  
  tree->SetBranchAddress("slc_true_event_type_incl", &slc_true_event_type_incl);
  tree->SetBranchAddress("slc_true_event_type_0p0pi", &slc_true_event_type_0p0pi);
  tree->SetBranchAddress("slc_true_event_type_Np0pi", &slc_true_event_type_Np0pi);
  
  tree->SetBranchAddress("slc_n_primary_razzled_muons", &slc_n_primary_razzled_muons);
  tree->SetBranchAddress("slc_n_primary_razzled_photons", &slc_n_primary_razzled_photons);
  tree->SetBranchAddress("slc_n_primary_razzled_pions_thresh", &slc_n_primary_razzled_pions_thresh);
  
  tree->SetBranchAddress("slc_n_primary_razzled_protons_thresh", &slc_n_primary_razzled_protons_thresh);
  tree->SetBranchAddress("slc_is_clear_cosmic", &slc_is_clear_cosmic);
  tree->SetBranchAddress("slc_is_fv", &slc_is_fv);
  tree->SetBranchAddress("slc_best_pzc_good_kinematics", &slc_best_pzc_good_kinematics);
  tree->SetBranchAddress("slc_all_other_trks_contained", &slc_all_other_trks_contained);
  tree->SetBranchAddress("slc_crumbs_score", &slc_crumbs_score);
  tree->SetBranchAddress("slc_opt0_fracPE", &slc_opt0_fracPE);
  tree->SetBranchAddress("slc_opt0_score", &slc_opt0_score);

  for(auto&& [i, name] : enumerate(weight_names))
    tree->SetBranchAddress(("nu_weight_" + name).c_str(), &nu_parameter_weights[i]);

  for(auto&& [i, name] : enumerate(weight_names))
    tree->SetBranchAddress(("slc_true_weight_" + name).c_str(), &slc_parameter_weights[i]);
}

void EvaluateTree(TChain *tree, const double &scaling, std::vector<TH1D*> &nominalSelectedHists, std::vector<TH1D*> &nominalTrueSignalHists,
                  std::vector<TH1D*> &nominalSelectedBackgroundHists, std::vector<std::vector<std::vector<TH1D*>>> &univTrueSignalHists,
                  std::vector<std::vector<std::vector<TH1D*>>> &univSelectedBackgroundHists, std::vector<std::vector<std::vector<TH1D*>>> &univSelectedSignalHists,
                  const std::vector<Cut> &signals, std::vector<std::string> weight_names, const uint n_univs, const bool combine)
{
  const uint N = tree->GetEntries();

  for(uint i = 0; i < N; ++i)
    {
      if(!(i%10000))
        std::cout << "Event: " << i << " / " << N << std::endl;
      tree->GetEntry(i);

      for(uint j = 0; j < nu_event_type_incl->size(); ++j)
        {
          for(auto&& [ signal_i, signal ] : enumerate(signals))
            {
              if((signal.name == "ncpizero_incl" && nu_event_type_incl->at(j) == 0) ||
                 (signal.name == "ncpizero_0p0pi" && nu_event_type_0p0pi->at(j) == 0) ||
                 (signal.name == "ncpizero_Np0pi" && nu_event_type_Np0pi->at(j) == 0))
                {
                  nominalTrueSignalHists[signal_i]->Fill(0.5, scaling);
                  
                  std::vector<float> all_parameter_weights(n_univs, 1.);

                  for(auto&& [ weight_i, weight ] : enumerate(weight_names))
                    {
                      for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
                        {
                          univTrueSignalHists[signal_i][weight_i][univ_i]->Fill(0.5, scaling * nu_parameter_weights[weight_i]->at(j)[univ_i]);
                          all_parameter_weights[univ_i] *= nu_parameter_weights[weight_i]->at(j)[univ_i];
                        }
                    }

                  if(combine)
                    {
                      for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
                        univTrueSignalHists[signal_i][weight_names.size()][univ_i]->Fill(0.5, scaling * all_parameter_weights[univ_i]);
                    }
                }
            }
        }

      for(uint j = 0; j < slc_true_event_type_incl->size(); ++j)
        {
          const bool sel_incl = !slc_is_clear_cosmic->at(j) && slc_is_fv->at(j) && slc_crumbs_score->at(j)>-0.195
            && slc_n_primary_razzled_muons->at(j) == 0 && slc_n_primary_razzled_photons->at(j)>1 && slc_best_pzc_good_kinematics->at(j)
            && slc_opt0_fracPE->at(j)<0.756 && slc_opt0_fracPE->at(j)>-0.7 && slc_opt0_score->at(j)>150 && slc_all_other_trks_contained->at(j);

          const bool sel_0p0pi = !slc_is_clear_cosmic->at(j) && slc_is_fv->at(j) && slc_crumbs_score->at(j)>-0.195
            && slc_n_primary_razzled_muons->at(j) == 0 && slc_n_primary_razzled_photons->at(j)>1 && slc_best_pzc_good_kinematics->at(j)
            && slc_n_primary_razzled_pions_thresh->at(j) == 0 && slc_n_primary_razzled_protons_thresh->at(j) == 0
            && slc_opt0_fracPE->at(j)<0.408 && slc_opt0_fracPE->at(j)>-0.7 && slc_opt0_score->at(j)>150;

          const bool sel_Np0pi = !slc_is_clear_cosmic->at(j) && slc_is_fv->at(j) && slc_crumbs_score->at(j)>-0.16
            && slc_n_primary_razzled_muons->at(j) == 0 && slc_n_primary_razzled_photons->at(j)>1 && slc_best_pzc_good_kinematics->at(j)
            && slc_n_primary_razzled_pions_thresh->at(j) == 0 && slc_n_primary_razzled_protons_thresh->at(j) > 0
            && slc_opt0_fracPE->at(j)<0.836 && slc_opt0_fracPE->at(j)>-0.376 && slc_opt0_score->at(j)>210;

          for(auto&& [ signal_i, signal ] : enumerate(signals))
            {
              if((signal.name == "ncpizero_incl" && sel_incl) ||
                 (signal.name == "ncpizero_0p0pi" && sel_0p0pi) ||
                 (signal.name == "ncpizero_Np0pi" && sel_Np0pi))
                {
                  nominalSelectedHists[signal_i]->Fill(0.5, scaling);

                  if((signal.name == "ncpizero_incl" && slc_true_event_type_incl->at(j) != 0) ||
                     (signal.name == "ncpizero_0p0pi" && slc_true_event_type_0p0pi->at(j) != 0) ||
                     (signal.name == "ncpizero_Np0pi" && slc_true_event_type_Np0pi->at(j) != 0))
                    {
                      nominalSelectedBackgroundHists[signal_i]->Fill(0.5, scaling);

                      std::vector<float> all_parameter_weights(n_univs, 1.);

                      for(auto&& [ weight_i, weight ] : enumerate(weight_names))
                        {
                          for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
                            {
                              const double w = slc_true_event_type_incl->at(j) > 5 ? 1. : slc_parameter_weights[weight_i]->at(j)[univ_i];
  
                              univSelectedBackgroundHists[signal_i][weight_i][univ_i]->Fill(0.5, scaling * w);
                              all_parameter_weights[univ_i] *= w;
                            }
                        }

                      if(combine)
                        {
                          for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
                            univSelectedBackgroundHists[signal_i][weight_names.size()][univ_i]->Fill(0.5, scaling * all_parameter_weights[univ_i]);
                        }
                    }
                  else
                    {
                      std::vector<float> all_parameter_weights(n_univs, 1.);

                      for(auto&& [ weight_i, weight ] : enumerate(weight_names))
                        {
                          for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
                            {
                              const double w = slc_true_event_type_incl->at(j) > 5 ? 1. : slc_parameter_weights[weight_i]->at(j)[univ_i];

                              univSelectedSignalHists[signal_i][weight_i][univ_i]->Fill(0.5, scaling * w);
                              all_parameter_weights[univ_i] *= w;
                            }
                        }

                      if(combine)
                        {
                          for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
                            univSelectedSignalHists[signal_i][weight_names.size()][univ_i]->Fill(0.5, scaling * all_parameter_weights[univ_i]);
                        }
                    }
                }
            }
        }
    }
}

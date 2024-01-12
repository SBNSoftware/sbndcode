#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Plots.h"
#include "Selections.h"
#include "Common.C"
#include "WeightNames.h"
#include "Enumerate.h"

void TrueEventMode2DPlotsSystWeights2(const TString productionVersion, const TString saveDirExt, const TwoDPlotSet &plotSet, const std::vector<Cut> &signals,
                                      std::vector<std::string> weight_names, const int n_univs, const bool combine)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/true_event_modes_two_d_syst_weights_2/" + saveDirExt;
  gSystem->Exec("mkdir -p " + saveDir);

  const TString ncpizeroFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_ncpizero.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *ncpizeroEvents = new TChain("ncpizeroana/events");
  ncpizeroEvents->Add(ncpizeroFile);

  TChain *ncpizeroSubruns = new TChain("ncpizeroana/subruns");
  ncpizeroSubruns->Add(ncpizeroFile);

  const double ncpizeroPOT     = GetPOT(ncpizeroSubruns);
  const double ncpizeroScaling = goalPOT / ncpizeroPOT;

  int run, subrun, event;
  std::vector<int> *nu_event_type_incl = 0, *nu_event_type_0p0pi = 0, *nu_event_type_Np0pi = 0, *nu_mode = 0;
  std::vector<std::vector<double>> *nu_pz_pizero_mom = 0, *nu_pz_cos_theta_pizero = 0;
  std::vector<std::vector<std::vector<double>>*> parameter_weights(weight_names.size(), 0);

  ncpizeroEvents->SetBranchStatus("*", 0);

  ncpizeroEvents->SetBranchAddress("run", &run);
  ncpizeroEvents->SetBranchAddress("subrun", &subrun);
  ncpizeroEvents->SetBranchAddress("event", &event);

  ncpizeroEvents->SetBranchAddress("nu_event_type_incl", &nu_event_type_incl);
  ncpizeroEvents->SetBranchAddress("nu_event_type_0p0pi", &nu_event_type_0p0pi);
  ncpizeroEvents->SetBranchAddress("nu_event_type_Np0pi", &nu_event_type_Np0pi);
  ncpizeroEvents->SetBranchAddress("nu_mode", &nu_mode);

  ncpizeroEvents->SetBranchAddress("nu_pz_pizero_mom", &nu_pz_pizero_mom);
  ncpizeroEvents->SetBranchAddress("nu_pz_cos_theta_pizero", &nu_pz_cos_theta_pizero);

  double bins1[plotSet.nbins1+1], bins2[plotSet.nbins2+1];

  for(int i = 0; i <= plotSet.nbins1; ++i)
    bins1[i] = plotSet.bins1.at(i);

  for(int i = 0; i <= plotSet.nbins2; ++i)
    bins2[i] = plotSet.bins2.at(i);

  for(auto&& [i, name] : enumerate(weight_names))
    ncpizeroEvents->SetBranchAddress(("nu_weight_" + name).c_str(), &parameter_weights[i]);

  const int N = ncpizeroEvents->GetEntries();

  std::vector<TH2D*> nominalHists;
  std::vector<std::vector<std::vector<TH2D*>>> univHists;

  for(auto&& [ signal_i, signal ] : enumerate(signals))
    {
      nominalHists.push_back(new TH2D(Form("nominalHist_%s", signal.name.Data()),
                                      Form(";%s;%s;Events", plotSet.axis1.Data(), plotSet.axis2.Data()),
                                      plotSet.nbins1, bins1, plotSet.nbins2, bins2));

      univHists.push_back(std::vector<std::vector<TH2D*>>());

      for(auto&& [ weight_i, weight ] : enumerate(weight_names))
        {
          univHists[signal_i].push_back(std::vector<TH2D*>());

          for(int univ_i = 0; univ_i < n_univs; ++univ_i)
            {
              univHists[signal_i][weight_i].push_back(new TH2D(Form("univHist_%s_%s_%i", signal.name.Data(), weight.c_str(), univ_i),
                                                               Form(";%s;%s;Events", plotSet.axis1.Data(), plotSet.axis2.Data()),
                                                               plotSet.nbins1, bins1, plotSet.nbins2, bins2));
            }
        }

      if(combine)
        {
          univHists[signal_i].push_back(std::vector<TH2D*>());

          for(int univ_i = 0; univ_i < n_univs; ++univ_i)
            {
              univHists[signal_i][weight_names.size()].push_back(new TH2D(Form("univHist_%s_%s_%i", signal.name.Data(), "all", univ_i),
                                                                          Form(";%s;%s;Events", plotSet.axis1.Data(), plotSet.axis2.Data()),
                                                                          plotSet.nbins1, bins1, plotSet.nbins2, bins2));
            }
	}
    }

  for(int i = 0; i < N; ++i)
    {
      if(!(i%1000))
        std::cout << "Event: " << i << " / " << N << std::endl;
      ncpizeroEvents->GetEntry(i);

      for(int j = 0; j < nu_event_type_incl->size(); ++j)
        {
          for(auto&& [ signal_i, signal ] : enumerate(signals))
            {
              if((signal.name == "ncpizero_incl" && nu_event_type_incl->at(j) == 0) ||
                 (signal.name == "ncpizero_0p0pi" && nu_event_type_0p0pi->at(j) == 0) ||
                 (signal.name == "ncpizero_Np0pi" && nu_event_type_Np0pi->at(j) == 0))
                {
                  nominalHists[signal_i]->Fill(1e3 * nu_pz_pizero_mom->at(j).at(0), nu_pz_cos_theta_pizero->at(j).at(0));

                  std::vector<float> all_parameter_weights(n_univs, 1.);

                  for(auto&& [ weight_i, weight ] : enumerate(weight_names))
                    {
                      for(int univ_i = 0; univ_i < n_univs; ++univ_i)
                        {
                          univHists[signal_i][weight_i][univ_i]->Fill(1e3 * nu_pz_pizero_mom->at(j).at(0), nu_pz_cos_theta_pizero->at(j).at(0), parameter_weights[weight_i]->at(j)[univ_i]);
                          all_parameter_weights[univ_i] *= parameter_weights[weight_i]->at(j)[univ_i];
                        }
                    }

                  if(combine)
                    {
                      for(int univ_i = 0; univ_i < n_univs; ++univ_i)
                        univHists[signal_i][weight_names.size()][univ_i]->Fill(1e3 * nu_pz_pizero_mom->at(j).at(0), nu_pz_cos_theta_pizero->at(j).at(0), all_parameter_weights[univ_i]);
                    }
                }
            }
        }
    }

  weight_names.push_back("all");

  for(auto&& [ signal_i, signal ] : enumerate(signals))
    {
      nominalHists[signal_i]->Scale(ncpizeroScaling);
      NormaliseEntriesByBinWidthX(nominalHists[signal_i], plotSet.scale);

      for(auto&& [ weight_i, weight ] : enumerate(weight_names))
        {
          gSystem->Exec("mkdir -p " + saveDir + "/" + weight.c_str());

          for(int univ_i = 0; univ_i < n_univs; ++univ_i)
            {
              univHists[signal_i][weight_i][univ_i]->Scale(ncpizeroScaling);
              NormaliseEntriesByBinWidthX(univHists[signal_i][weight_i][univ_i], plotSet.scale);
            }
        }
    }

  for(auto&& [ signal_i, signal ] : enumerate(signals))
    {
      for(auto&& [ weight_i, weight ] : enumerate(weight_names))
        {
          TCanvas *canvas = new TCanvas("canvas","canvas");
          canvas->cd();
          canvas->Draw();
          canvas->Divide(3, 3); // note this hardcodes for 9 bins or less for var2

          for(int frame_i = 0; frame_i < plotSet.nbins2; ++frame_i)
            {
              canvas->cd(frame_i+1);
              gPad->SetTopMargin(0.15);

              TH1D* nominalHist = nominalHists[signal_i]->ProjectionX(Form("nominalHist_%s_%s_%i", signal.name.Data(), weight.c_str(), frame_i), frame_i+1, frame_i+1);
              nominalHist->SetTitle(Form("%.2f < %s < %.2f;%s;Events / %s", plotSet.bins2[frame_i], plotSet.axis2.Data(), plotSet.bins2[frame_i+1], plotSet.axis1.Data(), plotSet.normalisationUnit.Data()));

              nominalHist->SetMaximum(1.2 * nominalHist->GetMaximum());
              nominalHist->Draw("histe");
              nominalHist->GetYaxis()->SetTitleOffset(1.3);
            }

          canvas->SaveAs(saveDir + "/" + weight.c_str() + "/" + Form("pizero_mom_cos_theta_%s_nominal.png", signal.name.Data()));
          canvas->SaveAs(saveDir + "/" + weight.c_str() + "/" + Form("pizero_mom_cos_theta_%s_nominal.pdf", signal.name.Data()));

          for(int frame_i = 0; frame_i < plotSet.nbins2; ++frame_i)
            {
              canvas->cd(frame_i+1);

              TH1D* nominalHist = nominalHists[signal_i]->ProjectionX(Form("nominalHist_%s_%s_%i", signal.name.Data(), weight.c_str(), frame_i), frame_i+1, frame_i+1);
              nominalHist->SetTitle(Form("%.2f < %s < %.2f;%s;Events / %s", plotSet.bins2[frame_i], plotSet.axis2.Data(), plotSet.bins2[frame_i+1], plotSet.axis1.Data(), plotSet.normalisationUnit.Data()));

              nominalHist->Draw("hist");
              nominalHist->GetYaxis()->SetTitleOffset(1.3);

              for(int univ_i = 0; univ_i < n_univs; ++univ_i)
                {
                  TH1D* univHist = univHists[signal_i][weight_i][univ_i]->ProjectionX(Form("univHist_%s_%s_%i_%i", signal.name.Data(), weight.c_str(), frame_i, univ_i), frame_i+1, frame_i+1);
                  univHist->SetTitle(Form("%.2f < %s < %.2f;%s;Events / %s", plotSet.bins2[frame_i], plotSet.axis2.Data(), plotSet.bins2[frame_i+1], plotSet.axis1.Data(), plotSet.normalisationUnit.Data()));

                  univHist->SetLineWidth(1);
                  univHist->SetLineColor(kMagenta-10);
                  univHist->Draw("histsame");
                }

              nominalHist->Draw("histsame");
            }

          canvas->SaveAs(saveDir + "/" + weight.c_str() + "/" + Form("pizero_mom_cos_theta_%s_%s_univs.png", signal.name.Data(), weight.c_str()));
          canvas->SaveAs(saveDir + "/" + weight.c_str() + "/" + Form("pizero_mom_cos_theta_%s_%s_univs.pdf", signal.name.Data(), weight.c_str()));

          for(int frame_i = 0; frame_i < plotSet.nbins2; ++frame_i)
            {
              canvas->cd(frame_i+1);

              TH1D* nominalHist = nominalHists[signal_i]->ProjectionX(Form("nominalHist_%s_%s_%i", signal.name.Data(), weight.c_str(), frame_i), frame_i+1, frame_i+1);
              nominalHist->SetTitle(Form("%.2f < %s < %.2f;%s;Events / %s", plotSet.bins2[frame_i], plotSet.axis2.Data(), plotSet.bins2[frame_i+1], plotSet.axis1.Data(), plotSet.normalisationUnit.Data()));

              nominalHist->Draw("hist");
              nominalHist->GetYaxis()->SetTitleOffset(1.3);

              TH1D* cvErrHist = new TH1D(Form("cvErrHist_%s_%s_%i", signal.name.Data(), weight.c_str(), frame_i), "",  plotSet.nbins1, bins1);

              for(int bin_i = 1; bin_i <= plotSet.nbins1; ++bin_i)
                {
                  TCanvas *binCanvas = new TCanvas(Form("binCanvas_%s_%s_%i_%i", signal.name.Data(), weight.c_str(), frame_i, bin_i),
                                                   Form("binCanvas_%s_%s_%i_%i", signal.name.Data(), weight.c_str(), frame_i, bin_i));
                  binCanvas->cd();

                  TH1D* binHist = new TH1D(Form("binHist_%s_%s_%i_%i", signal.name.Data(), weight.c_str(), frame_i, bin_i), "", 25,
                                           .75 * nominalHists[signal_i]->GetBinContent(bin_i, frame_i+1), 1.25 * nominalHists[signal_i]->GetBinContent(bin_i, frame_i+1));

                  for(int univ_i = 0; univ_i < n_univs; ++univ_i)
                    binHist->Fill(univHists[signal_i][weight_i][univ_i]->GetBinContent(bin_i, frame_i+1));

                  TF1* fGaus = new TF1("fGaus", "gaus", .75 * nominalHists[signal_i]->GetBinContent(bin_i, frame_i+1),
                                       1.25 * nominalHists[signal_i]->GetBinContent(bin_i, frame_i+1));
                  binHist->Fit(fGaus);

                  cvErrHist->SetBinContent(bin_i, fGaus->GetParameter("Mean"));
                  cvErrHist->SetBinError(bin_i, fGaus->GetParameter("Sigma"));
                }

              canvas->cd(frame_i+1);
              cvErrHist->SetMarkerStyle(1);
              cvErrHist->SetLineColor(kBlue+2);
              cvErrHist->Draw("histesame");
            }

          canvas->SaveAs(saveDir + "/" + weight.c_str() + "/" + Form("pizero_mom_cos_theta_%s_%s_cv_err.png", signal.name.Data(), weight.c_str()));
          canvas->SaveAs(saveDir + "/" + weight.c_str() + "/" + Form("pizero_mom_cos_theta_%s_%s_cv_err.pdf", signal.name.Data(), weight.c_str()));

          delete canvas;
        }
    }
}

#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Plots.h"
#include "Selections.h"
#include "Common.C"
#include "WeightNames.h"
#include "Enumerate.h"

void TrueEventModePlotsSystWeights(const TString productionVersion, const TString saveDirExt, const std::vector<Cut> &signals,
                                   std::vector<std::string> weight_names, const int n_univs, const bool combine)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/true_event_modes_syst_weights/" + saveDirExt;
  gSystem->Exec("mkdir -p " + saveDir);

  const TString ncpizeroFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_ncpizero.root";

  const TwoDPlotSet &plotSet = true_observables_twod_sets[0];

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

  std::vector<TH1D*> nominalMomHists, nominalCosThetaHists;
  std::vector<std::vector<std::vector<TH1D*>>> univMomHists, univCosThetaHists;

  for(auto&& [ signal_i, signal ] : enumerate(signals))
    {
      nominalMomHists.push_back(new TH1D(Form("nominalMomHist_%s", signal.name.Data()),
                                         Form(";%s;Events", plotSet.axis1.Data()), plotSet.nbins1, bins1));

      nominalCosThetaHists.push_back(new TH1D(Form("nominalCosThetaHist_%s", signal.name.Data()),
                                              Form(";%s;Events", plotSet.axis2.Data()), plotSet.nbins2, bins2));

      univMomHists.push_back(std::vector<std::vector<TH1D*>>());
      univCosThetaHists.push_back(std::vector<std::vector<TH1D*>>());

      for(auto&& [ weight_i, weight ] : enumerate(weight_names))
        {
          univMomHists[signal_i].push_back(std::vector<TH1D*>());
          univCosThetaHists[signal_i].push_back(std::vector<TH1D*>());

          for(int univ_i = 0; univ_i < n_univs; ++univ_i)
            {
              univMomHists[signal_i][weight_i].push_back(new TH1D(Form("univMomHist_%s_%s_%i", signal.name.Data(), weight.c_str(), univ_i),
                                                                  Form(";%s;Events", plotSet.axis1.Data()), plotSet.nbins1, bins1));

              univCosThetaHists[signal_i][weight_i].push_back(new TH1D(Form("univCosThetaHist_%s_%s_%i", signal.name.Data(), weight.c_str(), univ_i),
                                                                       Form(";%s;Events", plotSet.axis2.Data()), plotSet.nbins2, bins2));
            }
        }

      if(combine)
        {
          univMomHists[signal_i].push_back(std::vector<TH1D*>());
          univCosThetaHists[signal_i].push_back(std::vector<TH1D*>());

          for(int univ_i = 0; univ_i < n_univs; ++univ_i)
            {
              univMomHists[signal_i][weight_names.size()].push_back(new TH1D(Form("univMomHist_%s_%s_%i", signal.name.Data(), "all", univ_i),
                                                                             Form(";%s;Events", plotSet.axis1.Data()), plotSet.nbins1, bins1));

              univCosThetaHists[signal_i][weight_names.size()].push_back(new TH1D(Form("univCosThetaHist_%s_%s_%i", signal.name.Data(), "all", univ_i),
                                                                                  Form(";%s;Events", plotSet.axis2.Data()), plotSet.nbins2, bins2));
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
                  nominalMomHists[signal_i]->Fill(1e3 * nu_pz_pizero_mom->at(j).at(0));
                  nominalCosThetaHists[signal_i]->Fill(nu_pz_cos_theta_pizero->at(j).at(0));

                  std::vector<float> all_parameter_weights(n_univs, 1.);

                  for(auto&& [ weight_i, weight ] : enumerate(weight_names))
                    {
                      for(int univ_i = 0; univ_i < n_univs; ++univ_i)
                        {
                          if(parameter_weights[weight_i]->at(j)[univ_i] > 1e2)
                            parameter_weights[weight_i]->at(j)[univ_i] = 1.;

                          univMomHists[signal_i][weight_i][univ_i]->Fill(1e3 * nu_pz_pizero_mom->at(j).at(0), parameter_weights[weight_i]->at(j)[univ_i]);
                          univCosThetaHists[signal_i][weight_i][univ_i]->Fill(nu_pz_cos_theta_pizero->at(j).at(0), parameter_weights[weight_i]->at(j)[univ_i]);
                          all_parameter_weights[univ_i] *= parameter_weights[weight_i]->at(j)[univ_i];
                        }
                    }

                  if(combine)
                    {
                      for(int univ_i = 0; univ_i < n_univs; ++univ_i)
                        {
                          univMomHists[signal_i][weight_names.size()][univ_i]->Fill(1e3 * nu_pz_pizero_mom->at(j).at(0), all_parameter_weights[univ_i]);
                          univCosThetaHists[signal_i][weight_names.size()][univ_i]->Fill(nu_pz_cos_theta_pizero->at(j).at(0), all_parameter_weights[univ_i]);
                        }
                    }
                }
            }
        }
    }

  if(combine)
    weight_names.push_back("all");

  for(auto&& [ signal_i, signal ] : enumerate(signals))
    {
      nominalMomHists[signal_i]->Scale(ncpizeroScaling);
      NormaliseEntriesByBinWidth(nominalMomHists[signal_i], plotSet.scale1);

      nominalCosThetaHists[signal_i]->Scale(ncpizeroScaling);
      NormaliseEntriesByBinWidth(nominalCosThetaHists[signal_i], plotSet.scale2);

      for(auto&& [ weight_i, weight ] : enumerate(weight_names))
        {
          gSystem->Exec("mkdir -p " + saveDir + "/" + weight.c_str());

          for(int univ_i = 0; univ_i < n_univs; ++univ_i)
            {
              univMomHists[signal_i][weight_i][univ_i]->Scale(ncpizeroScaling);
              NormaliseEntriesByBinWidth(univMomHists[signal_i][weight_i][univ_i], plotSet.scale1);

              univCosThetaHists[signal_i][weight_i][univ_i]->Scale(ncpizeroScaling);
              NormaliseEntriesByBinWidth(univCosThetaHists[signal_i][weight_i][univ_i], plotSet.scale2);
            }
        }
    }

  for(auto&& [ signal_i, signal ] : enumerate(signals))
    {
      nominalMomHists[signal_i]->SetTitle(Form(";%s;Events / %s", plotSet.axis1.Data(), plotSet.normalisationUnit1.Data()));
      nominalMomHists[signal_i]->SetMaximum(1.2 * nominalMomHists[signal_i]->GetMaximum());
      nominalMomHists[signal_i]->SetMinimum(0.);

      nominalCosThetaHists[signal_i]->SetTitle(Form(";%s;Events / %s", plotSet.axis2.Data(), plotSet.normalisationUnit2.Data()));
      nominalCosThetaHists[signal_i]->SetMaximum(1.2 * nominalCosThetaHists[signal_i]->GetMaximum());
      nominalCosThetaHists[signal_i]->SetMinimum(0.);

      for(auto&& [ weight_i, weight ] : enumerate(weight_names))
        {
          TCanvas *momCanvas = new TCanvas("momCanvas","momCanvas");
          momCanvas->cd();

          nominalMomHists[signal_i]->Draw("histe");
          nominalMomHists[signal_i]->GetYaxis()->SetTitleOffset(1.3);

          momCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/" + Form("pizero_mom_%s_nominal.png", signal.name.Data()));
          momCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/" + Form("pizero_mom_%s_nominal.pdf", signal.name.Data()));

          nominalMomHists[signal_i]->Draw("hist");

          for(int univ_i = 0; univ_i < n_univs; ++univ_i)
            {
              univMomHists[signal_i][weight_i][univ_i]->SetTitle(Form(";%s;Events / %s", plotSet.axis1.Data(), plotSet.normalisationUnit1.Data()));
              univMomHists[signal_i][weight_i][univ_i]->SetLineWidth(1);
              univMomHists[signal_i][weight_i][univ_i]->SetLineColor(kMagenta-10);
              univMomHists[signal_i][weight_i][univ_i]->Draw("histsame");
            }

          nominalMomHists[signal_i]->Draw("histsame");

          momCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/" + Form("pizero_mom_%s_%s_univs.png", signal.name.Data(), weight.c_str()));
          momCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/" + Form("pizero_mom_%s_%s_univs.pdf", signal.name.Data(), weight.c_str()));

          ofstream logFile;
          logFile.open(saveDir + "/" + weight.c_str() + "/" + signal.name.Data() + "_log.txt");
          logFile << weight << "\n\n";

          nominalMomHists[signal_i]->Draw("hist");
          nominalMomHists[signal_i]->GetYaxis()->SetTitleOffset(1.3);

          TH1D* cvErrMomHist = new TH1D(Form("cvErrMomHist_%s_%s", signal.name.Data(), weight.c_str()), "",  plotSet.nbins1, bins1);

          for(int bin_i = 1; bin_i <= plotSet.nbins1; ++bin_i)
            {
              TCanvas *binMomCanvas = new TCanvas(Form("binMomCanvas_%s_%s_%i", signal.name.Data(), weight.c_str(), bin_i),
                                                  Form("binMomCanvas_%s_%s_%i", signal.name.Data(), weight.c_str(), bin_i));
              binMomCanvas->cd();

              double minVal = std::numeric_limits<double>::max(), maxVal = std::numeric_limits<double>::lowest();

              for(int univ_i = 0; univ_i < n_univs; ++univ_i)
                {
                  const double value = univMomHists[signal_i][weight_i][univ_i]->GetBinContent(bin_i);

                  if(value < minVal)
                    minVal = value;

                  if(value > maxVal)
                    maxVal = value;
                }

              double diff = maxVal - minVal;

              if(diff > std::numeric_limits<double>::epsilon())
                {
                  TH1D* binHist = new TH1D(Form("binMomHist_%s_%s_%i", signal.name.Data(), weight.c_str(), bin_i), "", 25, minVal - 0.05 * diff, maxVal + 0.05 * diff);

                  for(int univ_i = 0; univ_i < n_univs; ++univ_i)
                    binHist->Fill(univMomHists[signal_i][weight_i][univ_i]->GetBinContent(bin_i));

                  TF1* fGaus = new TF1("fGaus", "gaus", minVal - 0.05 * diff, maxVal + 0.05 * diff);

                  fGaus->SetLineColor(kRed+2);

                  TFitResultPtr fitResult = binHist->Fit(fGaus, "S");

                  cvErrMomHist->SetBinContent(bin_i, fGaus->GetParameter("Mean"));
                  cvErrMomHist->SetBinError(bin_i, fGaus->GetParameter("Sigma"));

                  logFile << Form("Mom, Bin %i, Sucess? ", bin_i);
                  if(fitResult->IsValid())
                    logFile << Form("Yes, Mean %f, Sigma %f\n", fGaus->GetParameter("Mean"), fGaus->GetParameter("Sigma"));
                  else
                    logFile << "No\n";

                  if(!fitResult->IsValid() || fGaus->GetParameter("Sigma") > 5000.)
                    {
                      gSystem->Exec("mkdir -p " + saveDir + "/" + weight.c_str() + "/bad_fits");
                      binMomCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/bad_fits/" + Form("mom_bin%i_%s_%s_cv_err.png", bin_i, signal.name.Data(), weight.c_str()));
                      binMomCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/bad_fits/" + Form("mom_bin%i_%s_%s_cv_err.pdf", bin_i, signal.name.Data(), weight.c_str()));
                    }
                }
              else
                {
                  cvErrMomHist->SetBinContent(bin_i, maxVal);
                  cvErrMomHist->SetBinError(bin_i, 0);

                  logFile << "Mom, Bin %i, No Variation - Error Set to 0";
                }
            }

          momCanvas->cd();

          cvErrMomHist->SetMarkerStyle(1);
          cvErrMomHist->SetLineColor(kBlue+2);
          cvErrMomHist->Draw("histesame");

          momCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/" + Form("pizero_mom_%s_%s_cv_err.png", signal.name.Data(), weight.c_str()));
          momCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/" + Form("pizero_mom_%s_%s_cv_err.pdf", signal.name.Data(), weight.c_str()));

          delete momCanvas;

          TCanvas *cosThetaCanvas = new TCanvas("cosThetaCanvas","cosThetaCanvas");
          cosThetaCanvas->cd();

          nominalCosThetaHists[signal_i]->Draw("histe");
          nominalCosThetaHists[signal_i]->GetYaxis()->SetTitleOffset(1.3);

          cosThetaCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/" + Form("pizero_cos_theta_%s_nominal.png", signal.name.Data()));
          cosThetaCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/" + Form("pizero_cos_theta_%s_nominal.pdf", signal.name.Data()));

          nominalCosThetaHists[signal_i]->Draw("hist");

          for(int univ_i = 0; univ_i < n_univs; ++univ_i)
            {
              univCosThetaHists[signal_i][weight_i][univ_i]->SetTitle(Form(";%s;Events / %s", plotSet.axis2.Data(), plotSet.normalisationUnit2.Data()));
              univCosThetaHists[signal_i][weight_i][univ_i]->SetLineWidth(1);
              univCosThetaHists[signal_i][weight_i][univ_i]->SetLineColor(kMagenta-10);
              univCosThetaHists[signal_i][weight_i][univ_i]->Draw("histsame");
            }

          nominalCosThetaHists[signal_i]->Draw("histsame");

          cosThetaCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/" + Form("pizero_cos_theta_%s_%s_univs.png", signal.name.Data(), weight.c_str()));
          cosThetaCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/" + Form("pizero_cos_theta_%s_%s_univs.pdf", signal.name.Data(), weight.c_str()));

          logFile << weight << "\n\n";

          nominalCosThetaHists[signal_i]->Draw("hist");
          nominalCosThetaHists[signal_i]->GetYaxis()->SetTitleOffset(1.3);

          TH1D* cvErrCosThetaHist = new TH1D(Form("cvErrCosThetaHist_%s_%s", signal.name.Data(), weight.c_str()), "",  plotSet.nbins2, bins2);

          for(int bin_i = 1; bin_i <= plotSet.nbins2; ++bin_i)
            {
              TCanvas *binCosThetaCanvas = new TCanvas(Form("binCosThetaCanvas_%s_%s_%i", signal.name.Data(), weight.c_str(), bin_i),
                                                       Form("binCosThetaCanvas_%s_%s_%i", signal.name.Data(), weight.c_str(), bin_i));
              binCosThetaCanvas->cd();

              double minVal = std::numeric_limits<double>::max(), maxVal = std::numeric_limits<double>::lowest();

              for(int univ_i = 0; univ_i < n_univs; ++univ_i)
                {
                  const double value = univCosThetaHists[signal_i][weight_i][univ_i]->GetBinContent(bin_i);

                  if(value < minVal)
                    minVal = value;

                  if(value > maxVal)
                    maxVal = value;
                }

              double diff = maxVal - minVal;

              if(diff > std::numeric_limits<double>::epsilon())
                {
                  TH1D* binHist = new TH1D(Form("binCosThetaHist_%s_%s_%i", signal.name.Data(), weight.c_str(), bin_i), "", 25, minVal - 0.05 * diff, maxVal + 0.05 * diff);

                  for(int univ_i = 0; univ_i < n_univs; ++univ_i)
                    binHist->Fill(univCosThetaHists[signal_i][weight_i][univ_i]->GetBinContent(bin_i));

                  TF1* fGaus = new TF1("fGaus", "gaus", minVal - 0.05 * diff, maxVal + 0.05 * diff);

                  fGaus->SetLineColor(kRed+2);

                  TFitResultPtr fitResult = binHist->Fit(fGaus, "S");

                  cvErrCosThetaHist->SetBinContent(bin_i, fGaus->GetParameter("Mean"));
                  cvErrCosThetaHist->SetBinError(bin_i, fGaus->GetParameter("Sigma"));

                  logFile << Form("CosTheta, Bin %i, Sucess? ", bin_i);
                  if(fitResult->IsValid())
                    logFile << Form("Yes, Mean %f, Sigma %f\n", fGaus->GetParameter("Mean"), fGaus->GetParameter("Sigma"));
                  else
                    logFile << "No\n";

                  if(!fitResult->IsValid() || fGaus->GetParameter("Sigma") > 500.)
                    {
                      gSystem->Exec("mkdir -p " + saveDir + "/" + weight.c_str() + "/bad_fits");
                      binCosThetaCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/bad_fits/" + Form("cos_theta_bin%i_%s_%s_cv_err.png", bin_i, signal.name.Data(), weight.c_str()));
                      binCosThetaCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/bad_fits/" + Form("cos_theta_bin%i_%s_%s_cv_err.pdf", bin_i, signal.name.Data(), weight.c_str()));
                    }
                }
              else
                {
                  cvErrCosThetaHist->SetBinContent(bin_i, maxVal);
                  cvErrCosThetaHist->SetBinError(bin_i, 0.);

                  logFile << "CosTheta, Bin %i, No Variation - Error Set to 0";
                }
            }

          cosThetaCanvas->cd();

          cvErrCosThetaHist->SetMarkerStyle(1);
          cvErrCosThetaHist->SetLineColor(kBlue+2);
          cvErrCosThetaHist->Draw("histesame");

          cosThetaCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/" + Form("pizero_cos_theta_%s_%s_cv_err.png", signal.name.Data(), weight.c_str()));
          cosThetaCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/" + Form("pizero_cos_theta_%s_%s_cv_err.pdf", signal.name.Data(), weight.c_str()));

          delete cosThetaCanvas;
        }
    }
}

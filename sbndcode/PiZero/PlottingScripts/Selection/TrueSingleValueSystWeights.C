#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Plots.h"
#include "Selections.h"
#include "Common.C"
#include "WeightNames.h"
#include "Enumerate.h"
#include "LatexHeaders.h"

void TrueSingleValueSystWeights(const TString productionVersion, const TString saveDirExt, const std::vector<Cut> &signals,
                                std::vector<std::string> weight_names, const uint n_univs, const bool combine)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/true_single_value_syst_weights/" + saveDirExt;
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
  const TString potString      = POTString();

  int run, subrun, event;
  std::vector<int> *nu_event_type_incl = 0, *nu_event_type_0p0pi = 0, *nu_event_type_Np0pi = 0;
  std::vector<std::vector<std::vector<double>>*> parameter_weights(weight_names.size(), 0);

  ncpizeroEvents->SetBranchStatus("*", 0);

  ncpizeroEvents->SetBranchAddress("run", &run);
  ncpizeroEvents->SetBranchAddress("subrun", &subrun);
  ncpizeroEvents->SetBranchAddress("event", &event);

  ncpizeroEvents->SetBranchAddress("nu_event_type_incl", &nu_event_type_incl);
  ncpizeroEvents->SetBranchAddress("nu_event_type_0p0pi", &nu_event_type_0p0pi);
  ncpizeroEvents->SetBranchAddress("nu_event_type_Np0pi", &nu_event_type_Np0pi);

  for(auto&& [i, name] : enumerate(weight_names))
    ncpizeroEvents->SetBranchAddress(("nu_weight_" + name).c_str(), &parameter_weights[i]);

  const uint N = ncpizeroEvents->GetEntries();

  std::vector<TH1D*> nominalHists;
  std::vector<std::vector<std::vector<TH1D*>>> univHists;

  for(auto&& [signal_i, signal] : enumerate(signals))
    {
      nominalHists.push_back(new TH1D(Form("nominalHist_%lu", signal_i), Form("%s;;Events%s", signal.printed_name.Data(), potString.Data()), 1, 0, 1));
      
      univHists.push_back(std::vector<std::vector<TH1D*>>());

      for(auto&& [weight_i, weight] : enumerate(weight_names))
        {
          univHists[signal_i].push_back(std::vector<TH1D*>());

          for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
            univHists[signal_i][weight_i].push_back(new TH1D(Form("univHist_%lu_%lu_%u", signal_i, weight_i, univ_i), "", 1, 0, 1));
        }

      if(combine)
        {
          univHists[signal_i].push_back(std::vector<TH1D*>());
          for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
            univHists[signal_i][weight_names.size()].push_back(new TH1D(Form("univHist_%lu_%lu_%u", signal_i, weight_names.size(), univ_i), "", 1, 0, 1));
        }

    }

  for(uint i = 0; i < N; ++i)
    {
      if(!(i%1000))
        std::cout << "Event: " << i << " / " << N << std::endl;
      ncpizeroEvents->GetEntry(i);

      for(uint j = 0; j < nu_event_type_incl->size(); ++j)
        {
          for(auto&& [ signal_i, signal ] : enumerate(signals))
            {
              if((signal.name == "ncpizero_incl" && nu_event_type_incl->at(j) == 0) ||
                 (signal.name == "ncpizero_0p0pi" && nu_event_type_0p0pi->at(j) == 0) ||
                 (signal.name == "ncpizero_Np0pi" && nu_event_type_Np0pi->at(j) == 0))
                {
                  nominalHists[signal_i]->Fill(0.5);
                  
                  std::vector<float> all_parameter_weights(n_univs, 1.);

                  for(auto&& [ weight_i, weight ] : enumerate(weight_names))
                    {
                      for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
                        {
                          univHists[signal_i][weight_i][univ_i]->Fill(0.5,parameter_weights[weight_i]->at(j)[univ_i]);
                          all_parameter_weights[univ_i] *= parameter_weights[weight_i]->at(j)[univ_i];
                        }
                    }

                  if(combine)
                    {
                      for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
                        univHists[signal_i][weight_names.size()][univ_i]->Fill(0.5, all_parameter_weights[univ_i]);
                    }
                }
            }
        }
    }

  if(combine)
    weight_names.push_back("all");

  for(auto&& [ signal_i, signal ] : enumerate(signals))
    {
      nominalHists[signal_i]->Scale(ncpizeroScaling);

      for(auto&& [ weight_i, weight ] : enumerate(weight_names))
        {
          gSystem->Exec("mkdir -p " + saveDir + "/" + weight.c_str());

          for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
            univHists[signal_i][weight_i][univ_i]->Scale(ncpizeroScaling);
        }
    }

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

          nominalHists[signal_i]->GetYaxis()->SetTitleOffset(1.9);
          nominalHists[signal_i]->GetXaxis()->SetLabelSize(0);
          nominalHists[signal_i]->GetXaxis()->SetLabelOffset(999);
          nominalHists[signal_i]->Draw("histe][");
          gPad->Update();
          TPaveText* title = (TPaveText*)gPad->FindObject("title");
          title->SetY1NDC(0.92);
          title->SetY2NDC(1);
          title->SetX1NDC(0.4);
          title->SetX2NDC(.75);
          gPad->Modified();
          gPad->Update();
          nominalHists[signal_i]->SetMinimum(0);
          nominalHists[signal_i]->SetMaximum(3.4e5);
        }

      canvas->SaveAs(saveDir + "/" + weight.c_str() + "/nominal.png");
      canvas->SaveAs(saveDir + "/" + weight.c_str() + "/nominal.pdf");

      for(auto&& [ signal_i, signal ] : enumerate(signals))
        {
          canvas->cd(signal_i + 1);
          nominalHists[signal_i]->Draw("hist][");
          gPad->Update();
          TPaveText* title = (TPaveText*)gPad->FindObject("title");
          title->SetY1NDC(0.92);
          title->SetY2NDC(1);
          title->SetX1NDC(0.4);
          title->SetX2NDC(.75);
          gPad->Modified();
          gPad->Update();

          for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
            {
              univHists[signal_i][weight_i][univ_i]->SetLineWidth(1);
              univHists[signal_i][weight_i][univ_i]->SetLineColor(kMagenta-10);
              univHists[signal_i][weight_i][univ_i]->Draw("hist][same");
            }

          nominalHists[signal_i]->Draw("hist][same");
        }

      canvas->cd(1);

      TLegend *univsLeg = new TLegend(0.3, 0.4, 0.7, 0.5);
      univsLeg->AddEntry(nominalHists[0], "Nominal", "l");
      univsLeg->AddEntry(univHists[0][weight_i][0], "Universes", "l");
      univsLeg->SetTextSize(.06);
      univsLeg->Draw();

      canvas->SaveAs(saveDir + "/" + weight.c_str() + "/univs.png");
      canvas->SaveAs(saveDir + "/" + weight.c_str() + "/univs.pdf");

      for(auto&& [ signal_i, signal ] : enumerate(signals))
        {
          canvas->cd(signal_i + 1);
          nominalHists[signal_i]->Draw("hist][");
          gPad->Update();
          TPaveText* title = (TPaveText*)gPad->FindObject("title");
          title->SetY1NDC(0.92);
          title->SetY2NDC(1);
          title->SetX1NDC(0.4);
          title->SetX2NDC(.75);
          gPad->Modified();
          gPad->Update();

          TH1D* cvErrHist = new TH1D(Form("cvErrHist_%lu_%lu", signal_i, weight_i), "", 1, 0, 1);

          TCanvas *fitCanvas = new TCanvas(Form("fitCanvas_%lu_%lu", signal_i, weight_i), Form("fitCanvas_%lu_%lu", signal_i, weight_i));
          fitCanvas->cd();

          double minVal = std::numeric_limits<double>::max(), maxVal = std::numeric_limits<double>::lowest();

          for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
            {
              const double value = univHists[signal_i][weight_i][univ_i]->GetBinContent(1);

              if(value < minVal)
                minVal = value;

              if(value > maxVal)
                maxVal = value;
            }

          double diff = maxVal - minVal;

          if(diff > std::numeric_limits<double>::epsilon())
            {
              TH1D* fitHist = new TH1D(Form("fitHist_%lu_%lu", signal_i, weight_i), "", 25, minVal - 0.05 * diff, maxVal + 0.05 * diff);

              for(uint univ_i = 0; univ_i < n_univs; ++univ_i)
                fitHist->Fill(univHists[signal_i][weight_i][univ_i]->GetBinContent(1));

              TF1* fGaus = new TF1("fGaus", "gaus", minVal - 0.05 * diff, maxVal + 0.05 * diff);
              fGaus->SetLineColor(kRed+2);

              TFitResultPtr fitResult = fitHist->Fit(fGaus, "S");

              cvErrHist->SetBinContent(1, fGaus->GetParameter("Mean"));
              cvErrHist->SetBinError(1, fGaus->GetParameter("Sigma"));

              if(!fitResult->IsValid() || fGaus->GetParameter("Sigma") > 20000.)
                {
                  gSystem->Exec("mkdir -p " + saveDir + "/" + weight.c_str() + "/bad_fits");
                  fitCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/bad_fits/" + Form("%s.png", signal.name.Data()));
                  fitCanvas->SaveAs(saveDir + "/" + weight.c_str() + "/bad_fits/" + Form("%s.pdf", signal.name.Data()));
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
                          cvErrHist->GetBinError(1) / nominalHists[signal_i]->GetBinContent(1) * 100.,
                          100. * (cvErrHist->GetBinContent(1) - nominalHists[signal_i]->GetBinContent(1)) / nominalHists[signal_i]->GetBinContent(1));
        }

      texFile << "\\\\ \\hline";
      
      canvas->cd(1);

      TH1D* tmpHist = new TH1D(Form("tmpHist_%lu", weight_i), "", 1, 0, 1);
      tmpHist->SetLineColor(kBlue+2);

      TLegend *cvLeg = new TLegend(0.3, 0.4, 0.7, 0.5);
      cvLeg->AddEntry(nominalHists[0], "Nominal", "l");
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

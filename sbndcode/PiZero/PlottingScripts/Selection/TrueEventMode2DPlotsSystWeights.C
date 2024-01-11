#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Plots.h"
#include "Selections.h"
#include "Common.C"
#include "WeightNames.h"

void TrueEventMode2DPlotsSystWeights(const TString productionVersion, const TwoDPlotSet &plotSet, const std::vector<Cut> &signals,
				     const std::vector<std::string> weight_names, const int n_univs)
{
  const TString saveDir = baseSaveDir + "/" + productionVersion + "/true_event_modes_two_d_syst_weights";
  gSystem->Exec("mkdir -p " + saveDir);

  const TString ncpizeroFile = baseFileDir + "/" + productionVersion + "/true_event_trees/" + productionVersion + "_ncpizero_true_event_trees.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *ncpizeroEvents = new TChain("nus");
  ncpizeroEvents->Add(ncpizeroFile);

  TChain *ncpizeroSubruns = new TChain("subruns");
  ncpizeroSubruns->Add(ncpizeroFile);

  const double ncpizeroPOT     = GetPOT(ncpizeroSubruns);
  const double ncpizeroScaling = goalPOT / ncpizeroPOT;

  for(const TString& weight_name : weight_names)
    {
      gSystem->Exec("mkdir -p " + saveDir + "/" + weight_name);

      for(auto const& signal : signals)
	{
	  double bins1[plotSet.nbins1 + 1];
	  std::copy(plotSet.bins1.begin(), plotSet.bins1.end(), bins1);

	  std::vector<THStack*> nominalStacks, weightedStacks;
	  for(int i = 0; i < plotSet.nbins2; ++i)
	    {
	      nominalStacks.push_back(new THStack(Form("nominalStack%i", i),
						  Form("%.2f < %s < %.2f;%s;Events / %s", plotSet.bins2[i],
						       plotSet.axis2.Data(), plotSet.bins2[i+1], plotSet.axis1.Data(),
						       plotSet.normalisationUnit.Data())));

	      weightedStacks.push_back(new THStack(Form("weightedStack%i", i),
						   Form("%.2f < %s < %.2f;%s;Events / %s", plotSet.bins2[i],
							plotSet.axis2.Data(), plotSet.bins2[i+1], plotSet.axis1.Data(),
							plotSet.normalisationUnit.Data())));
	    }

	  TLegend *lEventModes = new TLegend(.55, .28, .9, .8);

	  for(auto const& category : event_modes)
	    {
	      for(int i = 0; i < plotSet.nbins2; ++i)
		{
		  TH1F *hNominal = new TH1F(Form("hNominal%s%s%i%s", plotSet.name.Data(), signal.name.Data(), i, category.name.Data()), "", plotSet.nbins1, bins1);
		  TH1F *hWeighted = new TH1F(Form("hWeighted%s%s%i%s", plotSet.name.Data(), signal.name.Data(), i, category.name.Data()), "", plotSet.nbins1, bins1);
		      
		  ncpizeroEvents->Draw(Form("%s>>hNominal%s%s%i%s", plotSet.var1.Data(), plotSet.name.Data(), signal.name.Data(), i, category.name.Data()),
				       Form("((%s>%f && %s<%f) && %s && %s)", plotSet.var2.Data(), plotSet.bins2[i], plotSet.var2.Data(), plotSet.bins2[i+1], category.cut.GetTitle(), signal.cut.GetTitle()));
		      
		  hNominal->Scale(ncpizeroScaling);
		  NormaliseEntriesByBinWidth(hNominal, plotSet.scale);

		  std::vector<TH1F*> tempHists;
		  for(int univ = 0; univ < n_univs; ++univ)
		    {
		      TH1F *hTemp = new TH1F(Form("h%s%s%i%s%i", plotSet.name.Data(), signal.name.Data(), i, category.name.Data(), univ), "", plotSet.nbins1, bins1);

		      ncpizeroEvents->Draw(Form("%s>>h%s%s%i%s%i", plotSet.var1.Data(), plotSet.name.Data(), signal.name.Data(), i, category.name.Data(), univ),
					   Form("((%s>%f && %s<%f) && %s && %s) * %s", plotSet.var2.Data(), plotSet.bins2[i], plotSet.var2.Data(), plotSet.bins2[i+1],
						category.cut.GetTitle(), signal.cut.GetTitle(),
						Form("nu_weight_" + weight_name + "[%i]", univ)));

		      hTemp->Scale(ncpizeroScaling);
		      NormaliseEntriesByBinWidth(hTemp, plotSet.scale);
		      tempHists.push_back(hTemp);
		    }

		  std::vector<TH1F*> universesPerBin;

		  for(int bin = 0; bin < plotSet.nbins1; ++bin)
		    {
		      universesPerBin.push_back(new TH1F(Form("universesBin%i", bin+1), "", 25, .75 * hNominal->GetBinContent(bin+1), 1.25 * hNominal->GetBinContent(bin+1)));

		      for(int univ = 0; univ < n_univs; ++univ)
			universesPerBin[bin]->Fill(tempHists[univ]->GetBinContent(bin+1));

		      TCanvas *binCanvas = new TCanvas(Form("binCanvas%i", bin+1), Form("binCanvas%i", bin+1));
		      binCanvas->cd();

		      TF1 *fGaus = new TF1("fGaus", "gaus", .75 * hNominal->GetBinContent(bin+1), 1.25 * hNominal->GetBinContent(bin+1));
		      universesPerBin[bin]->Fit(fGaus);
		      fGaus->SetLineColor(kRed+2);

		      universesPerBin[bin]->Draw("hist");
		      fGaus->Draw("same");

		      binCanvas->SaveAs(saveDir + "/" + weight_name + Form("/%s_%s_bin_%i_%i.png", signal.name.Data(), category.name.Data(), i, bin+1));
		      binCanvas->SaveAs(saveDir + "/" + weight_name + Form("/%s_%s_bin_%i_%i.pdf", signal.name.Data(), category.name.Data(), i, bin+1));

		      hWeighted->SetBinContent(bin+1, fGaus->GetParameter("Mean"));
		      hWeighted->SetBinError(bin+1, fGaus->GetParameter("Sigma"));
		    }

		  hNominal->SetFillColorAlpha(category.colour, 0.4);
		  hNominal->SetLineColor(category.colour);
		  hNominal->SetLineWidth(2);
		  hNominal->SetMarkerStyle(1);

		  hWeighted->SetFillStyle(0);
		  hWeighted->SetLineColor(category.colour);
		  hWeighted->SetLineWidth(3);
		  hWeighted->SetLineStyle(7);
		  hWeighted->SetMarkerStyle(1);
		      
		  if(i == 0)
		    lEventModes->AddEntry(hNominal, category.name, "f");
		      
		  nominalStacks[i]->Add(hNominal);
		  weightedStacks[i]->Add(hWeighted);
		}
	    }

	  gStyle->SetPaperSize(45,80);

	  const int N = events->GetEntries();

	  for(int i = 0; i < N; ++i)
	    {
	      if(!(i%10000))
		std::cout << "Event: " << i << " / " << N << std::endl;
	      events->GetEntry(i);

	  TCanvas *canvas = new TCanvas("canvas","canvas");
	  canvas->cd();
	  canvas->Draw();
	  canvas->Divide(3, 3); // note this hardcodes for 9 bins or less for var2

	  for(int i = 0; i < plotSet.nbins2; ++i)
	    {
	      canvas->cd(i+1);
	      gPad->SetTopMargin(0.15);
	      nominalStacks[i]->SetMaximum(1.2 * nominalStacks[i]->GetMaximum());
	      nominalStacks[i]->Draw("histe");
	      nominalStacks[i]->GetYaxis()->SetTitleOffset(1.3);

	      weightedStacks[i]->Draw("histesame");

	      if(i==0)
		lEventModes->Draw();
	    }

	  canvas->SaveAs(saveDir + "/" + weight_name + Form("/%s_%s_by_event_mode_weight_%s.png", plotSet.name.Data(), signal.name.Data(), weight_name.Data()));
	  canvas->SaveAs(saveDir + "/" + weight_name + Form("/%s_%s_by_event_mode_weight_%s.pdf", plotSet.name.Data(), signal.name.Data(), weight_name.Data()));

	  delete canvas;
	}
    }
}

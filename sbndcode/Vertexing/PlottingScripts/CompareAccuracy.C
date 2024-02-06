void CompareAccuracy(const TString sample, const TString setupName, const bool cc_only, const int tpc = -1);

void CompareAccuracy()
{
  CompareAccuracy("rockbox", "all", true);
  CompareAccuracy("rockbox", "all", false);
  CompareAccuracy("rockbox", "simple", true);
  CompareAccuracy("rockbox", "simple", false);
  CompareAccuracy("rockbox", "centroid", true);
  CompareAccuracy("rockbox", "centroid", false);

  CompareAccuracy("rockbox", "all", true, 0);
  CompareAccuracy("rockbox", "all", false, 0);
  CompareAccuracy("rockbox", "simple", true, 0);
  CompareAccuracy("rockbox", "simple", false, 0);
  CompareAccuracy("rockbox", "centroid", true, 0);
  CompareAccuracy("rockbox", "centroid", false, 0);

  CompareAccuracy("rockbox", "all", true, 1);
  CompareAccuracy("rockbox", "all", false, 1);
  CompareAccuracy("rockbox", "simple", true, 1);
  CompareAccuracy("rockbox", "simple", false, 1);
  CompareAccuracy("rockbox", "centroid", true, 1);
  CompareAccuracy("rockbox", "centroid", false, 1);

  CompareAccuracy("intrnue", "all", true);
  CompareAccuracy("intrnue", "all", false);
  CompareAccuracy("intrnue", "simple", true);
  CompareAccuracy("intrnue", "simple", false);
  CompareAccuracy("intrnue", "centroid", true);
  CompareAccuracy("intrnue", "centroid", false);

  CompareAccuracy("intrnue", "all", true, 0);
  CompareAccuracy("intrnue", "all", false, 0);
  CompareAccuracy("intrnue", "simple", true, 0);
  CompareAccuracy("intrnue", "simple", false, 0);
  CompareAccuracy("intrnue", "centroid", true, 0);
  CompareAccuracy("intrnue", "centroid", false, 0);

  CompareAccuracy("intrnue", "all", true, 1);
  CompareAccuracy("intrnue", "all", false, 1);
  CompareAccuracy("intrnue", "simple", true, 1);
  CompareAccuracy("intrnue", "simple", false, 1);
  CompareAccuracy("intrnue", "centroid", true, 1);
  CompareAccuracy("intrnue", "centroid", false, 1);
}

void CompareAccuracy(const TString sample, const TString setupName, const bool cc_only, const int tpc)
{
  const TString saveDir = "/exp/sbnd/data/users/hlay/vertexing_fix/validation/" + sample;
  gSystem->Exec("mkdir -p " + saveDir);

  const TString file = "/pnfs/sbnd/persistent/users/hlay/vertexing_fix/vertexing_fix_" + sample + ".root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  struct Sample
  {
    TString name;
    TString treeName;
    TString legName;
    int     colour;
    TChain *tree = NULL;
  };

  std::vector<Sample> all_setups = { { "nominal", "vertexana/SliceTree", "Nominal", kBlack },
                                     { "correctedWeighting", "vertexanacorrectedrefinement/SliceTree", "Corrected Weighting", kMagenta+2 },
                                     { "noRefinement", "vertexananorefinement/SliceTree", "No Refinement", kSpring-1 },
                                     { "noWeighting", "vertexananoweights/SliceTree", "No Weighting", kRed+2 },
                                     { "smootherWeighting", "vertexanasmootherweights/SliceTree", "Smoother Weighting", kCyan+2 },
  };

  std::vector<Sample> centroid_setups = { { "nominal", "vertexana/SliceTree", "Nominal", kBlack },
                                          { "correctedWeighting", "vertexanacorrectedrefinement/SliceTree", "Corrected Weighting", kMagenta+2 },
                                          { "noRefinement", "vertexananorefinement/SliceTree", "No Refinement", kSpring-1 },
                                          { "noWeighting", "vertexananoweights/SliceTree", "No Weighting", kRed+2 },
                                          { "smootherWeighting", "vertexanasmootherweights/SliceTree", "Smoother Weighting", kCyan+2 },
                                          { "centroidPos", "vertexanacentroidpos/SliceTree", "Centroid Position", kViolet+3 },
  };

  std::vector<Sample> simple_setups = { { "nominal", "vertexana/SliceTree", "Nominal", kBlack },
                                        { "correctedWeighting", "vertexanacorrectedrefinement/SliceTree", "Corrected Weighting", kMagenta+2 },
  };

  std::vector<Sample> setups;

  if(setupName == "all")
    setups = all_setups;
  else if(setupName == "centroid")
    setups = centroid_setups;
  else if(setupName == "simple")
    setups = simple_setups;
  
  for(auto& setup : setups)
    {
      setup.tree = new TChain(setup.treeName);
      setup.tree->Add(file);
    }

  TCanvas *canvas = new TCanvas("canvas", "canvas");
  canvas->cd();

  TLegend *legend = new TLegend(.5, .8, .5, .8);
  
  std::vector<TH1F*> hists;
  float max = std::numeric_limits<float>::lowest();

  for(auto const& setup : setups)
    {
      TH1F *hist = new TH1F(Form("hist_%s", setup.name.Data()), ";Distance between true and reconstructed vertex (cm);Events",
                            40, 0, 10);

      TCut evSel = "is_fv";

      if(cc_only)
        {
          evSel += "cc";
          if(sample == "rockbox")
            evSel += "numu";
          else if(sample == "intrnue")
            evSel += "nue";
        }

      if(tpc == 0)
        evSel += "true_vtx_x<0";
      else if(tpc == 1)
        evSel += "true_vtx_x>0";

      setup.tree->Draw(Form("dr_x_corr>>hist_%s", setup.name.Data()), evSel);
      hist->SetLineColor(setup.colour);
      hists.push_back(hist);

      legend->AddEntry(hist, setup.legName, "le");

      if(hist->GetMaximum() > max)
        max = hist->GetMaximum();
    }

  max *= 1.3;

  for(auto const& hist : hists)
    {
      hist->SetMaximum(max);
      hist->GetYaxis()->SetTitleOffset(1.25);
      hist->GetXaxis()->SetTitleOffset(1.25);
      hist->GetXaxis()->SetTitleSize(.048);
      hist->Draw("histesame");
    }

  legend->Draw();

  TString addedName = "";

  if(cc_only)
    {
      addedName += "_cc";
      if(sample == "rockbox")
        addedName += "_numu";
      else if(sample == "intrnue")
        addedName += "_nue";
    }

  if(tpc != -1)
    addedName += Form("_tpc_%i", tpc);

  canvas->SaveAs(saveDir + "/vertex_accuracy_" + sample + "_" + setupName + "_variations" + addedName + ".png");
  canvas->SaveAs(saveDir + "/vertex_accuracy_" + sample + "_" + setupName + "_variations" + addedName + ".pdf");
}

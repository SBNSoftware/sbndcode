#include "Common.C"
#include "XSecPlot.h"
#include "XSecSample.h"
#include "Selection.h"
#include "WeightSet.h"
#include "Enumerate.h"
#include "SystNames.h"

XSecSamples SetupSamples(const TString productionVersion)
{
  const TString rockboxFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_rockbox.root";
  const TString intimeFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_intime.root";

  TChain *rockboxNus = new TChain("ncpizeroxsectrees/neutrinos");
  rockboxNus->Add(rockboxFile);
  TChain *rockboxSlices = new TChain("ncpizeroxsectrees/slices");
  rockboxSlices->Add(rockboxFile);
  TChain *rockboxSubruns = new TChain("ncpizeroxsectrees/subruns");
  rockboxSubruns->Add(rockboxFile);

  TChain *intimeNus = new TChain("ncpizeroxsectrees/neutrinos");
  intimeNus->Add(intimeFile);
  TChain *intimeSlices = new TChain("ncpizeroxsectrees/slices");
  intimeSlices->Add(intimeFile);
  TChain *intimeSubruns = new TChain("ncpizeroxsectrees/subruns");
  intimeSubruns->Add(intimeFile);

  double rockboxScaling, intimeScaling;
  GetScaling(rockboxSubruns, intimeSubruns, rockboxScaling, intimeScaling);

  XSecSamples samples = { { "Rockbox", rockboxNus, rockboxSlices, rockboxScaling },
                          { "Intime", intimeNus, intimeSlices, intimeScaling }
  };

  return samples;
}

void FillPlots(XSecSamples &samples, Selections &selections, WeightSets &weightSets)
{
  int nweights = 0;

  for(Selection &selection : selections)
    {
      XSecPlot *plot = selection.plot;

      for(Bin* bin : plot->GetBins())
        {
          for(WeightSet &weightSet : weightSets)
            {
              for(std::string &name : weightSet.list)
                {
                  bin->AddWeight(name, weightSet.nunivs);
                  ++nweights;
                }
            }
        }
    }

  nweights /= selections.size();

  for(XSecSample &sample : samples)
    {
      std::cout << '\n' << "Sample: " << sample.name << std::endl;

      const float scaling = sample.scaling / samples[0].scaling;

      std::deque<bool>                 sel(selections.size(), false);
      std::vector<int>                 event_type(selections.size(), -1);
      std::vector<std::vector<float>*> weights;
      float                            comp;
      double                           val0;
      double                           val1;

      for(int w = 0; w < nweights; ++w)
        weights.push_back(new std::vector<float>());

      sample.nutree->SetBranchStatus("*", 0);

      for(auto&& [ selection_i, selection ] : enumerate(selections))
        {
          XSecPlot *plot = selection.plot;

          const std::string var0 = plot->GetVar0();
          const std::string var1 = plot->GetVar1();

          sample.nutree->SetBranchAddress(selection.signal, &event_type[selection_i]);
          sample.nutree->SetBranchAddress(var0.c_str(), &val0);
          sample.nutree->SetBranchAddress(var1.c_str(), &val1);

          int weight_i = 0;

          for(auto&& [ weightSet_i, weightSet ] : enumerate(weightSets))
            {
              for(auto&& [ name_i, name ] : enumerate(weightSet.list))
                {
                  sample.nutree->SetBranchAddress(name.c_str(), &weights[weight_i]);
                  ++weight_i;
                }
            }
        }

      const int n_nu = sample.nutree->GetEntries();

      std::cout << "\n\tNeutrinos" << std::endl;

      for(int i = 0; i < n_nu; ++i)
        {
          if(!(i%100000))
            std::cout << '\t' << i << " / " << n_nu
                      << " (" << (100. * i) / n_nu << "%)" << std::endl;

          sample.nutree->GetEntry(i);

          for(auto&& [ selection_i, selection ] : enumerate(selections))
            {
              if(selection_i == 0 && event_type[selection_i] == 1)
                event_type[selection_i] = 0;

              if(event_type[selection_i] == 0)
                {
                  XSecPlot *plot = selection.plot;

                  for(Bin* bin : plot->GetBins())
                    {
                      if(!bin->InBin(val0, val1))
                        continue;

                      bin->IncrementNominalBinTrueSignal(scaling);

                      int weight_i = 0;

                      for(auto&& [ weightSet_i, weightSet ] : enumerate(weightSets))
                        {
                          for(auto&& [ name_i, name ] : enumerate(weightSet.list))
                            {
                              if(event_type[selection_i] == 7 || event_type[selection_i] == 8)
                                *weights[weight_i] = std::vector<float>(weightSet.nunivs, 1.);

                              for(int univ = 0; univ < weightSet.nunivs; ++univ)
                                bin->IncrementUniverseBinTrueSignal(name, univ, scaling * weights[weight_i]->at(univ));

                              ++weight_i;
                            }
                        }
                    }
                }
            }
        }  

      sample.slicetree->SetBranchStatus("*", 0);
      sample.slicetree->SetBranchAddress("comp", &comp);

      for(auto&& [ selection_i, selection ] : enumerate(selections))
        {
          XSecPlot *plot = selection.plot;

          const std::string var0 = plot->GetVar0();
          const std::string var1 = plot->GetVar1();

          sample.slicetree->SetBranchAddress(selection.signal, &event_type[selection_i]);
          sample.slicetree->SetBranchAddress(selection.cut, &sel[selection_i]);
          sample.slicetree->SetBranchAddress(("reco_" + var0).c_str(), &val0);
          sample.slicetree->SetBranchAddress(("reco_" + var1).c_str(), &val1);

          int weight_i = 0;

          for(auto&& [ weightSet_i, weightSet ] : enumerate(weightSets))
            {
              for(auto&& [ name_i, name ] : enumerate(weightSet.list))
                {
                  sample.slicetree->SetBranchAddress(name.c_str(), &weights[weight_i]);
                  ++weight_i;
                }
            }
        }

      const int n_slice = sample.slicetree->GetEntries();

      std::cout << "\n\tSlices" << std::endl;

      for(int i = 0; i < n_slice; ++i)
        {
          if(!(i%100000))
            std::cout << '\t' << i << " / " << n_slice
                      << " (" << (100. * i) / n_slice << "%)" << std::endl;

          sample.slicetree->GetEntry(i);

          for(auto&& [ selection_i, selection ] : enumerate(selections))
            {
              if(selection_i == 0 && event_type[selection_i] == 1)
                event_type[selection_i] = 0;

              if(sel[selection_i])
                {
                  XSecPlot *plot = selection.plot;

                  for(Bin* bin : plot->GetBins())
                    {
                      if(!bin->InBin(val0, val1))
                        continue;

                      bin->IncrementNominalBinCount(scaling);

                      if(!(event_type[selection_i] == 0 && comp > .5))
                        bin->IncrementNominalBinBkgdCount(scaling);

                      int weight_i = 0;

                      for(auto&& [ weightSet_i, weightSet ] : enumerate(weightSets))
                        {
                          for(auto&& [ name_i, name ] : enumerate(weightSet.list))
                            {
                              if(event_type[selection_i] == 7 || event_type[selection_i] == 8)
                                *weights[weight_i] = std::vector<float>(weightSet.nunivs, 1.);

                              for(int univ = 0; univ < weightSet.nunivs; ++univ)
                                {
                                  bin->IncrementUniverseBinCount(name, univ, scaling * weights[weight_i]->at(univ));

                                  if(!(event_type[selection_i] == 0 && comp > .5))
                                    bin->IncrementUniverseBinBkgdCount(name, univ, scaling * weights[weight_i]->at(univ));
                                }
                              ++weight_i;
                            }
                        }
                    }
                }
            }
        }      
    }

  for(Selection &selection : selections)
    {
      XSecPlot *plot = selection.plot;

      for(Bin* bin : plot->GetBins())
        {
          bin->SetScaleFactor(samples[0].scaling);
          bin->Update();
          bin->CalculateXSecPurity();
          //          bin->CalculateSystFracErrorsMedianPercentiles();
          bin->CalculateSystFracErrorsNominalSD();
        }

      plot->InsertSystFracFlatError("ntargets", 0.01);
      plot->InsertSystFracFlatError("pot", 0.02);
    }
}

std::vector<std::string> SystSetToWeightList(const std::vector<Syst> &systs)
{
  std::vector<std::string> list;

  for(auto const& syst : systs)
    list.push_back(syst.name);

  return list;
}

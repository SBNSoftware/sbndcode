#include "XSecPlot.h"
#include "XSecSample.h"
#include "Selection.h"
#include "WeightSet.h"
#include "Enumerate.h"

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


      for(int w = 0; w < nweights; ++w)
        weights.push_back(new std::vector<float>);

      sample.nutree->SetBranchStatus("*", 0);

      for(auto&& [ selection_i, selection ] : enumerate(selections))
        {
          XSecPlot *plot = selection.plot;

          sample.nutree->SetBranchAddress(selection.signal, &event_type[selection_i]);

          for(Bin* bin : plot->GetBins())
            {
              for(auto&& [ weightSet_i, weightSet ] : enumerate(weightSets))
                {
                  for(auto&& [ name_i, name ] : enumerate(weightSet.list))
                    {
                      int weight_i = weightSet_i * weightSets.size() + name_i;

                      sample.nutree->SetBranchAddress(name.c_str(), &weights[weight_i]);
                    }
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
              if(event_type[selection_i] == 0)
                {
                  XSecPlot *plot = selection.plot;

                  for(Bin* bin : plot->GetBins())
                    {
                      bin->IncrementNominalBinTrueSignal(scaling);

                      for(auto&& [ weightSet_i, weightSet ] : enumerate(weightSets))
                        {
                          for(auto&& [ name_i, name ] : enumerate(weightSet.list))
                            {
                              int weight_i = weightSet_i * weightSets.size() + name_i;

                              if(event_type[selection_i] == 7 || event_type[selection_i] == 8)
                                *weights[weight_i] = std::vector<float>(weightSet.nunivs, 1.);

                              for(int univ = 0; univ < weightSet.nunivs; ++univ)
                                bin->IncrementUniverseBinTrueSignal(name, univ, scaling * weights[weight_i]->at(univ));
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

          sample.slicetree->SetBranchAddress(selection.signal, &event_type[selection_i]);
          sample.slicetree->SetBranchAddress(selection.cut, &sel[selection_i]);

          for(Bin* bin : plot->GetBins())
            {
              for(auto&& [ weightSet_i, weightSet ] : enumerate(weightSets))
                {
                  for(auto&& [ name_i, name ] : enumerate(weightSet.list))
                    {
                      int weight_i = weightSet_i * weightSets.size() + name_i;

                      sample.slicetree->SetBranchAddress(name.c_str(), &weights[weight_i]);
                    }
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
              if(sel[selection_i])
                {
                  XSecPlot *plot = selection.plot;

                  for(Bin* bin : plot->GetBins())
                    {
                      bin->IncrementNominalBinCount(scaling);

                      if(!(event_type[selection_i] == 0 && comp > .5))
                        bin->IncrementNominalBinBkgdCount(scaling);

                      for(auto&& [ weightSet_i, weightSet ] : enumerate(weightSets))
                        {
                          for(auto&& [ name_i, name ] : enumerate(weightSet.list))
                            {
                              int weight_i = weightSet_i * weightSets.size() + name_i;

                              if(event_type[selection_i] == 7 || event_type[selection_i] == 8)
                                *weights[weight_i] = std::vector<float>(weightSet.nunivs, 1.);

                              for(int univ = 0; univ < weightSet.nunivs; ++univ)
                                {
                                  bin->IncrementUniverseBinCount(name, univ, scaling * weights[weight_i]->at(univ));

                                  if(!(event_type[selection_i] == 0 && comp > .5))
                                    bin->IncrementUniverseBinBkgdCount(name, univ, scaling * weights[weight_i]->at(univ));
                                }
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
          bin->CalculateSystFracErrors();
        }
    }
}

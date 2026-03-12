#include "Common.C"
#include "XSecPlot.h"
#include "XSecSample.h"
#include "Selection.h"
#include "WeightSet.h"
#include "Enumerate.h"
#include "SystNames.h"

XSecSamples SetupSamples(const TString productionVersion)
{
  const TString rockboxFile  = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_rockbox.root";
  const TString ncpizeroFile = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_ncpizero2.root";
  const TString intimeFile   = baseFileDir + "/" + productionVersion + "/" + productionVersion + "_intime.root";

  TChain *rockboxNus = new TChain("ncpizeroxsectrees/neutrinos");
  rockboxNus->Add(rockboxFile);
  TChain *rockboxSlices = new TChain("ncpizeroxsectrees/slices");
  rockboxSlices->Add(rockboxFile);
  TChain *rockboxSubruns = new TChain("ncpizeroxsectrees/subruns");
  rockboxSubruns->Add(rockboxFile);

  TChain *ncpizeroNus = new TChain("ncpizeroxsectrees/neutrinos");
  ncpizeroNus->Add(ncpizeroFile);
  TChain *ncpizeroSlices = new TChain("ncpizeroxsectrees/slices");
  ncpizeroSlices->Add(ncpizeroFile);
  TChain *ncpizeroSubruns = new TChain("ncpizeroxsectrees/subruns");
  ncpizeroSubruns->Add(ncpizeroFile);

  TChain *intimeNus = new TChain("ncpizeroxsectrees/neutrinos");
  intimeNus->Add(intimeFile);
  TChain *intimeSlices = new TChain("ncpizeroxsectrees/slices");
  intimeSlices->Add(intimeFile);
  TChain *intimeSubruns = new TChain("ncpizeroxsectrees/subruns");
  intimeSubruns->Add(intimeFile);

  double rockboxScaling, ncpizeroScaling, intimeScaling;
  GetScaling(rockboxSubruns, ncpizeroSubruns, intimeSubruns, rockboxScaling, ncpizeroScaling, intimeScaling);

  XSecSamples samples = { { "rockbox", rockboxNus, rockboxSlices, rockboxScaling },//, { 0, 1 } },
                          //              { "ncpizero", ncpizeroNus, ncpizeroSlices, ncpizeroScaling, { 2, 3, 4, 5, 6, 7, 8 } },
                          { "intime", intimeNus, intimeSlices, intimeScaling }
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
      double                           val;
      double                           trueVal;

      for(int w = 0; w < nweights; ++w)
        weights.push_back(new std::vector<float>());

      sample.nutree->SetBranchStatus("*", 0);

      for(auto&& [ selection_i, selection ] : enumerate(selections))
        {
          XSecPlot *plot = selection.plot;

          const std::string var = plot->GetVar();

          sample.nutree->SetBranchAddress(selection.signal, &event_type[selection_i]);
          sample.nutree->SetBranchAddress(var.c_str(), &val);

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
              if(sample.mask.count(event_type[selection_i]) != 0)
                continue;

              if(selection_i == 0 && event_type[selection_i] == 1)
                {
                  std::cout << "Why is this happening?" << std::endl;
                  event_type[selection_i] = 0;
                }

              if(event_type[selection_i] == 0)
                {
                  XSecPlot *plot = selection.plot;

                  for(Bin* bin : plot->GetBins())
                    {
                      if(!bin->InBin(val))
                        continue;

                      bin->IncrementNominalBinTrueSignal(scaling * extraSignalScaling.at(selection.name));

                      int weight_i = 0;

                      for(auto&& [ weightSet_i, weightSet ] : enumerate(weightSets))
                        {
                          for(auto&& [ name_i, name ] : enumerate(weightSet.list))
                            {
                              if(event_type[selection_i] == 7 || event_type[selection_i] == 8)
                                *weights[weight_i] = std::vector<float>(weightSet.nunivs, 1.);

                              for(int univ = 0; univ < weightSet.nunivs; ++univ)
                                bin->IncrementUniverseBinTrueSignal(name, univ, scaling * weights[weight_i]->at(univ) * extraSignalScaling.at(selection.name));

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

          std::string var = plot->GetVar();

          sample.slicetree->SetBranchAddress(var.c_str(), &trueVal);

          if(var == "pizero_mom")
            var += "_fit";

          sample.slicetree->SetBranchAddress(selection.signal, &event_type[selection_i]);
          sample.slicetree->SetBranchAddress(selection.cut, &sel[selection_i]);

          sample.slicetree->SetBranchAddress(("reco_" + var).c_str(), &val);

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
              if(sample.mask.count(event_type[selection_i]) != 0)
                continue;

              if(selection_i == 0 && event_type[selection_i] == 1)
                {
                  std::cout << "Why is this happening?" << std::endl;
                  event_type[selection_i] = 0;
                }

              if(sel[selection_i])
                {
                  XSecPlot *plot = selection.plot;

                  for(Bin* bin : plot->GetBins())
                    {
                      if(!bin->InBin(val) && !bin->InBin(trueVal))
                        continue;

                      if(bin->InBin(trueVal) && (event_type[selection_i] == 0 && comp > .5))
                        bin->IncrementNominalBinSelSignalTrueBin(scaling * extraSignalScaling.at(selection.name));

                      if(bin->InBin(val))
                        {
                          if(!(event_type[selection_i] == 0 && comp > .5))
                            {
                              bin->IncrementNominalBinCount(scaling);
                              bin->IncrementNominalBinBkgdCount(scaling);
                            }
                          else
                            {
                              bin->IncrementNominalBinCount(scaling * extraSignalScaling.at(selection.name));
                            }
                        }

                      int weight_i = 0;

                      for(auto&& [ weightSet_i, weightSet ] : enumerate(weightSets))
                        {
                          for(auto&& [ name_i, name ] : enumerate(weightSet.list))
                            {
                              if(event_type[selection_i] == 7 || event_type[selection_i] == 8)
                                *weights[weight_i] = std::vector<float>(weightSet.nunivs, 1.);

                              for(int univ = 0; univ < weightSet.nunivs; ++univ)
                                {
                                  if(bin->InBin(trueVal) && (event_type[selection_i] == 0 && comp > .5))
                                    bin->IncrementUniverseBinSelSignalTrueBin(name, univ, scaling * weights[weight_i]->at(univ) * extraSignalScaling.at(selection.name));

                                  if(bin->InBin(val))
                                    {
                                      if(!(event_type[selection_i] == 0 && comp > .5))
                                        {
                                          bin->IncrementUniverseBinCount(name, univ, scaling * weights[weight_i]->at(univ));
                                          bin->IncrementUniverseBinBkgdCount(name, univ, scaling * weights[weight_i]->at(univ));
                                        }
                                      else
                                        {
                                          bin->IncrementUniverseBinCount(name, univ, scaling * weights[weight_i]->at(univ) *  extraSignalScaling.at(selection.name));
                                        }
                                    }
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
          //bin->CalculateSystFracErrorsMedianPercentiles();
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

double Chi2Comp(TH1F* observed, TH1F* pred, TH2D* cov)
{
  double sum = 0.;

  TMatrixD covmatrix(observed->GetNbinsX(), observed->GetNbinsX());

  for(int i = 1; i < observed->GetNbinsX() + 1; ++i)
    {
      for(int j = 1; j < observed->GetNbinsX() + 1; ++j)
	{
	  covmatrix[i - 1][j - 1] = cov->GetBinContent(i, j);
	}
    }
  covmatrix *= 1e88;
  covmatrix.Invert();
  covmatrix *= 1e88;
  
  for(int i = 1; i < observed->GetNbinsX() + 1; ++i)
    {
      for(int j = 1; j < observed->GetNbinsX() + 1; ++j)
	{
	  sum += (observed->GetBinContent(i) - pred->GetBinContent(i)) * covmatrix[i - 1][j - 1] * (observed->GetBinContent(j) - pred->GetBinContent(j));
	}
    }

  return sum;
}

TH2D* NormMatrix(TH1F *p, TH2D *m)
{
  TH2D* n = new TH2D(Form("%s_norm", m->GetName()), ";Bin;Bin", m->GetNbinsX(), .5, m->GetNbinsX() + .5,
		     m->GetNbinsY(), .5, m->GetNbinsY() + .5);

  double p_sum = 0;

  for(int i = 1; i < p->GetNbinsX() + 1; ++i)
    {
      p_sum += p->GetBinContent(i);
    }

  double m_sum = 0;
  
  for(int i = 1; i < m->GetNbinsX() + 1; ++i)
    {
      for(int j = 1; j < m->GetNbinsY() + 1; ++j)
	{
	  m_sum += m->GetBinContent(i, j);
	}
    }

  for(int i = 1; i < m->GetNbinsX() + 1; ++i)
    {
      for(int j = 1; j < m->GetNbinsY() + 1; ++j)
	{
	  n->SetBinContent(i, j, p->GetBinContent(i) * p->GetBinContent(j) * m_sum / (p_sum * p_sum));
	}
    }

  return n;
}

TH2D* ShapeMatrix(TH1F *p, TH2D *m)
{
  TH2D* s = new TH2D(Form("%s_shape", m->GetName()), ";Bin;Bin", m->GetNbinsX(), .5, m->GetNbinsX() + .5,
		     m->GetNbinsY(), .5, m->GetNbinsY() + .5);

  double p_sum = 0;

  for(int i = 1; i < p->GetNbinsX() + 1; ++i)
    {
      p_sum += p->GetBinContent(i);
    }

  for(int i = 1; i < m->GetNbinsX() + 1; ++i)
    {
      for(int j = 1; j < m->GetNbinsY() + 1; ++j)
	{
	  double sum_m_ik = 0, sum_m_kj = 0, sum_m_kl = 0;

	  for(int k = 1; k < m->GetNbinsY() + 1; ++k)
	    sum_m_ik += m->GetBinContent(i, k);

	  for(int k = 1; k < m->GetNbinsX() + 1; ++k)
	    sum_m_kj += m->GetBinContent(k, j);

	  for(int k = 1; k < m->GetNbinsX() + 1; ++k)
	    {
	      for(int l = 1; l < m->GetNbinsY() + 1; ++l)
		{
		  sum_m_kl += m->GetBinContent(k, l);
		}
	    }

	  double term_1 = (p->GetBinContent(j) / p_sum) * sum_m_ik;
	  double term_2 = (p->GetBinContent(i) / p_sum) * sum_m_kj;
	  double term_3 = ((p->GetBinContent(i) * p->GetBinContent(j)) / (p_sum * p_sum)) * sum_m_kl;

	  s->SetBinContent(i, j, m->GetBinContent(i, j) - term_1 - term_2 + term_3);
	}
    }

  return s;
}

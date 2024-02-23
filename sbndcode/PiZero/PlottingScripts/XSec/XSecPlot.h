#pragma once

#include "Bin.h"

class XSecPlot {

 private:
  std::string         _name;
  std::string         _axes_labels;
  std::vector<Bin*>   _bins;
  std::vector<double> _var0Bins;
  std::vector<double> _var1Bins;
  std::string         _var0;
  std::string         _var1;

 public:

  XSecPlot(const std::string name, const std::string axes_labels, const int N,
           const std::vector<double> var0Bins, const std::vector<double> var1Bins,
           const std::string var0, const std::string var1,
           const double nTargets, const double intFlux)
    {
      _name        = name;
      _axes_labels = axes_labels;
      _var0Bins    = var0Bins;
      _var1Bins    = var1Bins;
      _var0        = var0;
      _var1        = var1;

      if((var0Bins.size() - 1) * (var1Bins.size() - 1) != N ||
         var0Bins.size() < 2 || var1Bins.size() < 2)
        throw std::runtime_error("Bin Sizes Messed Up");

      int n = 0;

      for(int i = 1; i < var0Bins.size(); ++i)
        {
          const double width0 = var0Bins.size() == 2 ? 1. : var0Bins[i] - var0Bins[i-1];

          for(int j = 1; j < var1Bins.size(); ++j)
            {
              const double width1 = var1Bins.size() == 2 ? 1. : var1Bins[j] - var1Bins[j-1];
              _bins.emplace_back(new Bin(n, width0 * width1, nTargets, intFlux,
                                         var0Bins[i-1], var0Bins[i], var1Bins[j-1], var1Bins[j]));

              ++n;
            }
        }
    }

  void SetName(const std::string name)
  {
    _name = name;
  }

  void SetAxesLabels(const std::string axes_labels)
  {
    _axes_labels = axes_labels;
  }

  std::string GetName() const
    {
      return _name;
    }

  std::string GetAxesLabels() const
    {
      return _axes_labels;
    }

  Bin* GetBin(const int i)
  {
    if(i >= _bins.size())
      throw std::runtime_error("Asking for non-existent bin");

    return _bins[i];
  }

  std::vector<Bin*> GetBins()
    {
      return _bins;
    }

  std::string GetVar0()
    {
      return _var0;
    }

  std::string GetVar1()
    {
      return _var1;
    }

  TH1F* GetNominalHist(const bool statErr = true)
  {
    TH1F* hist = new TH1F(_name.c_str(), _axes_labels.c_str(), 1, 0, 1);
    hist->SetBinContent(1, _bins[0]->GetNominalXSec());

    if(statErr)
      hist->SetBinError(1, _bins[0]->GetNominalFracStatErr() * _bins[0]->GetNominalXSec());
    else
      hist->SetBinError(1, 0);

    hist->SetLineColor(kBlack);

    return hist;
  }

  TH1F* GetUniverseHist(const std::string &weightName, const int univ)
  {
    TH1F* hist = new TH1F(Form("%s%s%i",_name.c_str(), weightName.c_str(), univ),
                          _axes_labels.c_str(), 1, 0, 1);
    hist->SetBinContent(1, _bins[0]->GetUniverseXSec(weightName, univ));
    hist->SetBinError(1, 0);

    hist->SetMarkerStyle(0);
    hist->SetLineColor(kMagenta-10);
    hist->SetLineWidth(1);

    return hist;
  }

  TGraphAsymmErrors* GetCVErrGraph(const std::string &weightName)
  {
    TGraphAsymmErrors* graph = new TGraphAsymmErrors(1);

    double low, cv, high;

    std::tie(low, cv, high) = _bins[0]->GetSystFracErrors(weightName);

    graph->SetPoint(0, 0.5, cv);
    graph->SetPointError(0, 0.48, 0.48, low, high);

    graph->SetMarkerStyle(1);
    graph->SetLineColor(kBlue+2);
    graph->SetLineWidth(5);

    return graph;
  }

  TGraphAsymmErrors* GetCVErrGraph(const std::vector<std::string> &weightNames)
  {
    TGraphAsymmErrors* graph = new TGraphAsymmErrors(1);

    const double nom  = _bins[0]->GetNominalXSec();
    const double bias = _bins[0]->GetFracSystBiasQuadSum(weightNames);
    const double res  = _bins[0]->GetFracSystResAveQuadSum(weightNames);

    graph->SetPoint(0, 0.5, nom * (1 - bias));
    graph->SetPointError(0, 0.48, 0.48, nom * res, nom * res);

    graph->SetMarkerStyle(1);
    graph->SetLineColor(kBlue+2);
    graph->SetLineWidth(5);

    return graph;
  }
};

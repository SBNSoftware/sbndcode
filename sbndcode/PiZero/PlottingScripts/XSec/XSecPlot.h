#pragma once

#include "Bin.h"

class XSecPlot {

 private:
  std::string       _name;
  std::string       _axes_labels;
  std::vector<Bin*> _bins;

 public:

  XSecPlot(const std::string name, const std::string axes_labels, const int n,
           const double nTargets, const double intFlux)
    {
      _name        = name;
      _axes_labels = axes_labels;

      for(int i = 0; i < n; ++i)
        _bins.emplace_back(new Bin(i, 1., nTargets, intFlux));
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
};

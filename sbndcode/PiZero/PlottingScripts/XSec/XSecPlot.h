#pragma once

#include "Bin.h"
#include "Enumerate.h"

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

  TH2F* GetNominalHist(const bool statErr = true)
  {
    double var0Bins[_var0Bins.size()];
    GetVar0BinsArray(var0Bins);
    double var1Bins[_var1Bins.size()];
    GetVar1BinsArray(var1Bins);

    TH2F* hist = new TH2F(_name.c_str(), _axes_labels.c_str(),
                          _var0Bins.size() - 1, var0Bins,
                          _var1Bins.size() - 1, var1Bins);

    for(int i = 0; i < _bins.size(); ++i)
      {
        const int binIndex0 = i / _var0Bins.size() + 1;
        const int binIndex1 = i % _var0Bins.size() + 1;

        hist->SetBinContent(binIndex0, binIndex1, _bins[i]->GetNominalXSec());

        if(statErr)
          hist->SetBinError(binIndex0, binIndex1, _bins[i]->GetNominalFracStatErr() * _bins[i]->GetNominalXSec());
        else
          hist->SetBinError(binIndex0, binIndex1, 0);
      }

    return hist;
  }

  TH1F* GetNominalHist0D(const bool statErr = true)
  {
    double var0Bins[_var0Bins.size()];
    GetVar0BinsArray(var0Bins);

    if(_bins.size() != 1)
      throw std::runtime_error("Asking for 0D hist but numbers of bins is not 1");

    TH2F* full_hist = GetNominalHist(statErr);

    TH1F* hist = new TH1F(_name.c_str(), _axes_labels.c_str(),
                          _var0Bins.size() - 1, var0Bins);

    hist->SetBinContent(1, full_hist->GetBinContent(1, 1));
    hist->SetBinError(1, full_hist->GetBinError(1, 1));
    hist->SetLineColor(kBlack);

    return hist;
  }

  TH1F* GetNominalHist1D(const bool statErr = true)
  {
    if(_var0Bins.size() != 2 && _var1Bins.size() != 2)
      throw std::runtime_error("Asking for 1D hist but neither set has 1 bin");

    TH2F* full_hist = GetNominalHist(statErr);

    if(_var0Bins.size() == 2)
      {
        double var1Bins[_var1Bins.size()];
        GetVar1BinsArray(var1Bins);

        TH1F* hist = new TH1F(_name.c_str(), _axes_labels.c_str(),
                              _var1Bins.size() - 1, var1Bins);

        for(auto&& [ binEdge_i, binEdge ] : enumerate(_var1Bins))
          {
            hist->SetBinContent(binEdge_i + 1, full_hist->GetBinContent(1, binEdge_i + 1));
            hist->SetBinError(binEdge_i + 1, full_hist->GetBinError(1, binEdge_i + 1));
            hist->SetLineColor(kBlack);
          }

        return hist;
      }

    if(_var1Bins.size() == 2)
      {
        double var0Bins[_var0Bins.size()];
        GetVar0BinsArray(var0Bins);

        TH1F* hist = new TH1F(_name.c_str(), _axes_labels.c_str(),
                              _var0Bins.size() - 1, var0Bins);

        for(auto&& [ binEdge_i, binEdge ] : enumerate(_var0Bins))
          {
            hist->SetBinContent(binEdge_i + 1, full_hist->GetBinContent(binEdge_i + 1, 1));
            hist->SetBinError(binEdge_i + 1, full_hist->GetBinError(binEdge_i + 1, 1));
            hist->SetLineColor(kBlack);
          }

        return hist;
      }

    return NULL;
  }

  TH2F* GetUniverseHist(const std::string &weightName, const int univ)
  {
    double var0Bins[_var0Bins.size()];
    GetVar0BinsArray(var0Bins);
    double var1Bins[_var1Bins.size()];
    GetVar1BinsArray(var1Bins);

    TH2F* hist = new TH2F(Form("General%s%s%i",_name.c_str(), weightName.c_str(), univ),
                          _axes_labels.c_str(),
                          _var0Bins.size() - 1, var0Bins,
                          _var1Bins.size() - 1, var1Bins);

    for(int i = 0; i < _bins.size(); ++i)
      {
        const int binIndex0 = i / _var0Bins.size() + 1;
        const int binIndex1 = i % _var0Bins.size() + 1;

        hist->SetBinContent(binIndex0, binIndex1, _bins[i]->GetUniverseXSec(weightName, univ));
        hist->SetBinError(binIndex0, binIndex1, 0);
      }

    return hist;
  }

  TH1F* GetUniverseHist0D(const std::string &weightName, const int univ)
  {
    double var0Bins[_var0Bins.size()];
    GetVar0BinsArray(var0Bins);

    if(_bins.size() != 1)
      throw std::runtime_error("Asking for 0D hist but numbers of bins is not 1");

    TH2F* full_hist = GetUniverseHist(weightName, univ);

    TH1F* hist = new TH1F(Form("OneD%s%s%i",_name.c_str(), weightName.c_str(), univ),
                          _axes_labels.c_str(),
                          _var0Bins.size() - 1, var0Bins);

    hist->SetBinContent(1, full_hist->GetBinContent(1, 1));
    hist->SetBinError(1, full_hist->GetBinError(1, 1));

    hist->SetMarkerStyle(0);
    hist->SetLineColor(kMagenta-10);
    hist->SetLineWidth(1);

    return hist;
  }

  TH1F* GetUniverseHist1D(const std::string &weightName, const int univ)
  {
    if(_var0Bins.size() != 2 && _var1Bins.size() != 2)
      throw std::runtime_error("Asking for 1D hist but neither set has 1 bin");

    TH2F* full_hist = GetUniverseHist(weightName, univ);

    if(_var0Bins.size() == 2)
      {
        double var1Bins[_var1Bins.size()];
        GetVar1BinsArray(var1Bins);

        TH1F* hist = new TH1F(_name.c_str(), _axes_labels.c_str(),
                              _var1Bins.size() - 1, var1Bins);

        for(auto&& [ binEdge_i, binEdge ] : enumerate(_var1Bins))
          {
            hist->SetBinContent(binEdge_i + 1, full_hist->GetBinContent(1, binEdge_i + 1));
            hist->SetBinError(binEdge_i + 1, full_hist->GetBinError(1, binEdge_i + 1));
            hist->SetLineColor(kBlack);
          }

        return hist;
      }

    if(_var1Bins.size() == 2)
      {
        double var0Bins[_var0Bins.size()];
        GetVar0BinsArray(var0Bins);

        TH1F* hist = new TH1F(_name.c_str(), _axes_labels.c_str(),
                              _var0Bins.size() - 1, var0Bins);

        for(auto&& [ binEdge_i, binEdge ] : enumerate(_var0Bins))
          {
            hist->SetBinContent(binEdge_i + 1, full_hist->GetBinContent(binEdge_i + 1, 1));
            hist->SetBinError(binEdge_i + 1, full_hist->GetBinError(binEdge_i + 1, 1));
            hist->SetLineColor(kBlack);
          }

        return hist;
      }

    return NULL;
  }

  TGraphAsymmErrors* GetCVErrGraph(const std::string &weightName)
  {
    TGraphAsymmErrors* graph = new TGraphAsymmErrors(1);

    for(auto&& [ bin_i, bin ] : enumerate(_bins))
      {
        double low, cv, high;
        std::tie(low, cv, high) = bin->GetSystFracErrors(weightName);

        graph->SetPoint(bin_i, bin->GetVar0Center(), cv);
        graph->SetPointError(0, 0.49 * bin->GetVar0Width(),
                             0.49 * bin->GetVar0Width(), low * cv, high * cv);

        graph->SetMarkerStyle(1);
        graph->SetLineColor(kBlue+2);
        graph->SetLineWidth(5);
      }

    return graph;
  }

  TGraphAsymmErrors* GetCVErrGraph(const std::vector<std::string> &weightNames)
  {
    TGraphAsymmErrors* graph = new TGraphAsymmErrors(1);

    for(auto&& [ bin_i, bin ] : enumerate(_bins))
      {
        const double nom  = bin->GetNominalXSec();
        const double bias = bin->GetFracSystBiasQuadSum(weightNames);
        const double res  = bin->GetFracSystResAveQuadSum(weightNames);

        graph->SetPoint(bin_i, bin->GetVar0Center(), nom * (1 - bias));
        graph->SetPointError(bin_i, 0.49 * bin->GetVar0Width(),
                             0.49 * bin->GetVar0Width(), nom * res, nom * res);

        graph->SetMarkerStyle(1);
        graph->SetLineColor(kBlue+2);
        graph->SetLineWidth(5);
      }

    return graph;
  }

  void GetVar0BinsArray(double var0Bins[])
  {
    for(auto&& [ binEdge_i, binEdge ] : enumerate(_var0Bins))
      var0Bins[binEdge_i] = binEdge;
  }

  void GetVar1BinsArray(double var1Bins[])
  {
    for(auto&& [ binEdge_i, binEdge ] : enumerate(_var1Bins))
      var1Bins[binEdge_i] = binEdge;
  }
};

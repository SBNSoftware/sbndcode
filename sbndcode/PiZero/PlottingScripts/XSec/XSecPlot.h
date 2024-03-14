#pragma once

#include "Bin.h"
#include "Enumerate.h"

// Flux Config I Front Face
/*
constexpr double numufluxscale  = 0.92341993;
constexpr double anumufluxscale = 0.070391454;
constexpr double nuefluxscale   = 0.0055868623;
constexpr double anuefluxscale  = 0.00060175334;
*/

// Flux Config I Effective Z (11227.8cm)
/*
constexpr double numufluxscale  = 0.92379188;
constexpr double anumufluxscale = 0.070015689;
constexpr double nuefluxscale   = 0.0055930597;
constexpr double anuefluxscale  = 0.00059937579;
*/

// Flux Config L Effective Z (11227.8cm)
constexpr double numufluxscale  = 0.92364301;
constexpr double anumufluxscale = 0.070125432;
constexpr double nuefluxscale   = 0.0056067340;
constexpr double anuefluxscale  = 0.00062482531;

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

      for(int i = 1; i < var1Bins.size(); ++i)
        {
          const double width1 = var1Bins.size() == 2 ? 1. : var1Bins[i] - var1Bins[i-1];

          for(int j = 1; j < var0Bins.size(); ++j)
            {
              const double width0 = var0Bins.size() == 2 ? 1. : var0Bins[j] - var0Bins[j-1];
              _bins.emplace_back(new Bin(n, width0 * width1, nTargets, intFlux,
                                         var0Bins[j-1], var0Bins[j], var1Bins[i-1], var1Bins[i]));

              ++n;
            }
        }
    }

  void Print()
  {
    for(Bin* bin : _bins)
      {
        bin->Print();
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

  void InsertSystFracFlatError(const std::string &name, const double &frac)
  {
    for(int i = 0; i < _bins.size(); ++i)
      {
        const double nom = _bins[i]->GetNominalXSec();
        _bins[i]->InsertSystFracError(name, nom, frac);
      }
  }

  TH2F* GetNominalHist(const bool statErr = true)
  {
    double var0Bins[_var0Bins.size()];
    GetVar0BinsArray(var0Bins);
    double var1Bins[_var1Bins.size()];
    GetVar1BinsArray(var1Bins);

    TH2F* hist = new TH2F(Form("General_%s", _name.c_str()), _axes_labels.c_str(),
                          _var0Bins.size() - 1, var0Bins,
                          _var1Bins.size() - 1, var1Bins);

    for(int i = 0; i < _bins.size(); ++i)
      {
        const int binIndex0 = i % ( _var0Bins.size() - 1 ) + 1;
        const int binIndex1 = i / ( _var0Bins.size() - 1 ) + 1;

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

    TH1F* hist = new TH1F(Form("ZeroD_%s", _name.c_str()), _axes_labels.c_str(),
                          _var0Bins.size() - 1, var0Bins);

    hist->SetMarkerStyle(0);
    hist->SetBinContent(1, full_hist->GetBinContent(1, 1));
    hist->SetBinError(1, full_hist->GetBinError(1, 1));
    hist->SetLineColor(kBlack);

    delete full_hist;
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

        TH1F* hist = new TH1F(Form("OneD_%s", _name.c_str()), _axes_labels.c_str(),
                              _var1Bins.size() - 1, var1Bins);

        for(auto&& [ binEdge_i, binEdge ] : enumerate(_var1Bins))
          {
            hist->SetBinContent(binEdge_i + 1, full_hist->GetBinContent(1, binEdge_i + 1));
            hist->SetBinError(binEdge_i + 1, full_hist->GetBinError(1, binEdge_i + 1));
          }

        hist->SetMarkerStyle(0);
        hist->SetLineColor(kBlack);

        delete full_hist;
        return hist;
      }

    if(_var1Bins.size() == 2)
      {
        double var0Bins[_var0Bins.size()];
        GetVar0BinsArray(var0Bins);

        TH1F* hist = new TH1F(Form("OneD_%s", _name.c_str()), _axes_labels.c_str(),
                              _var0Bins.size() - 1, var0Bins);

        for(auto&& [ binEdge_i, binEdge ] : enumerate(_var0Bins))
          {
            hist->SetBinContent(binEdge_i + 1, full_hist->GetBinContent(binEdge_i + 1, 1));
            hist->SetBinError(binEdge_i + 1, full_hist->GetBinError(binEdge_i + 1, 1));
          }

        hist->SetMarkerStyle(0);
        hist->SetLineColor(kBlack);

        delete full_hist;
        return hist;
      }

    return NULL;
  }

  TH1F* GetNominalHist2D(const bool statErr = true, const int var1Bin = 1)
  {
    TH2F* full_hist = GetNominalHist(statErr);

    double var0Bins[_var0Bins.size()];
    GetVar0BinsArray(var0Bins);

    TH1F* hist = new TH1F(Form("TwoD_%s_%i", _name.c_str(), var1Bin),
                          _axes_labels.c_str(),
                          _var0Bins.size() - 1, var0Bins);

    for(auto&& [ binEdge_i, binEdge ] : enumerate(_var0Bins))
      {
        hist->SetBinContent(binEdge_i + 1, full_hist->GetBinContent(binEdge_i + 1, var1Bin));
        hist->SetBinError(binEdge_i + 1, full_hist->GetBinError(binEdge_i + 1, var1Bin));
      }

    hist->SetMarkerStyle(0);
    hist->SetLineColor(kBlack);

    delete full_hist;
    return hist;
  }

  TH2F* GetUniverseHist(const std::string &weightName, const int univ)
  {
    double var0Bins[_var0Bins.size()];
    GetVar0BinsArray(var0Bins);
    double var1Bins[_var1Bins.size()];
    GetVar1BinsArray(var1Bins);

    TH2F* hist = new TH2F(Form("General_%s_%s_%i",_name.c_str(), weightName.c_str(), univ),
                          _axes_labels.c_str(),
                          _var0Bins.size() - 1, var0Bins,
                          _var1Bins.size() - 1, var1Bins);

    for(int i = 0; i < _bins.size(); ++i)
      {
        const int binIndex0 = i % ( _var0Bins.size() -1 ) + 1;
        const int binIndex1 = i / ( _var0Bins.size() -1 ) + 1;

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

    TH1F* hist = new TH1F(Form("ZeroD_%s_%s_%i",_name.c_str(), weightName.c_str(), univ),
                          _axes_labels.c_str(),
                          _var0Bins.size() - 1, var0Bins);

    hist->SetBinContent(1, full_hist->GetBinContent(1, 1));
    hist->SetBinError(1, full_hist->GetBinError(1, 1));

    hist->SetMarkerStyle(0);
    hist->SetLineColor(kMagenta-10);
    hist->SetLineWidth(1);

    delete full_hist;
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

        TH1F* hist = new TH1F(Form("OneD_%s_%s_%i",_name.c_str(), weightName.c_str(), univ),
                              _axes_labels.c_str(), _var1Bins.size() - 1, var1Bins);

        for(auto&& [ binEdge_i, binEdge ] : enumerate(_var1Bins))
          {
            hist->SetBinContent(binEdge_i + 1, full_hist->GetBinContent(1, binEdge_i + 1));
            hist->SetBinError(binEdge_i + 1, full_hist->GetBinError(1, binEdge_i + 1));
          }

        hist->SetMarkerStyle(0);
        hist->SetLineColor(kMagenta-10);
        hist->SetLineWidth(1);

        delete full_hist;
        return hist;
      }

    if(_var1Bins.size() == 2)
      {
        double var0Bins[_var0Bins.size()];
        GetVar0BinsArray(var0Bins);

        TH1F* hist = new TH1F(Form("OneD_%s_%s_%i",_name.c_str(), weightName.c_str(), univ),
                              _axes_labels.c_str(), _var0Bins.size() - 1, var0Bins);

        for(auto&& [ binEdge_i, binEdge ] : enumerate(_var0Bins))
          {
            hist->SetBinContent(binEdge_i + 1, full_hist->GetBinContent(binEdge_i + 1, 1));
            hist->SetBinError(binEdge_i + 1, full_hist->GetBinError(binEdge_i + 1, 1));
          }

        hist->SetMarkerStyle(0);
        hist->SetLineColor(kMagenta-10);
        hist->SetLineWidth(1);

        delete full_hist;
        return hist;
      }

    return NULL;
  }

  TH1F* GetUniverseHist2D(const std::string &weightName, const int univ, const int var1Bin)
  {
    TH2F* full_hist = GetUniverseHist(weightName, univ);

    double var0Bins[_var0Bins.size()];
    GetVar0BinsArray(var0Bins);

    TH1F* hist = new TH1F(Form("TwoD_%s_%s_%i_%i",_name.c_str(), weightName.c_str(), univ, var1Bin),
                          _axes_labels.c_str(), _var0Bins.size() - 1, var0Bins);

    for(auto&& [ binEdge_i, binEdge ] : enumerate(_var0Bins))
      {
        hist->SetBinContent(binEdge_i + 1, full_hist->GetBinContent(binEdge_i + 1, var1Bin));
        hist->SetBinError(binEdge_i + 1, full_hist->GetBinError(binEdge_i + 1, var1Bin));
      }

    hist->SetMarkerStyle(0);
    hist->SetLineColor(kMagenta-10);
    hist->SetLineWidth(1);

    delete full_hist;
    return hist;
  }

  TGraphAsymmErrors* GetCVErrGraph(const std::string &weightName, const int var1Bin = -1)
  {
    TGraphAsymmErrors* graph = new TGraphAsymmErrors();
    int bin_count = 0;

    for(auto&& [ bin_i, bin ] : enumerate(_bins))
      {
        const int binIndex1 = bin_i / (_var0Bins.size() - 1) + 1;

        if(var1Bin != -1 && var1Bin != binIndex1)
          continue;

        double low, cv, high;
        std::tie(low, cv, high) = bin->GetSystFracErrors(weightName);

        if(cv == 0.)
          continue;

        double centre, width;

        if(_var0Bins.size() == 2)
          {
            centre = bin->GetVar1Center();
            width  = bin->GetVar1Width();
          }
        else
          {
            centre = bin->GetVar0Center();
            width  = bin->GetVar0Width();
          }

        graph->SetPoint(bin_count, centre, cv);
        graph->SetPointError(bin_count, 0.49 * width, 0.49 * width, low * cv, high * cv);

        graph->SetMarkerStyle(1);
        graph->SetLineColor(kBlue+2);
        graph->SetLineWidth(5);

        ++bin_count;
      }

    return graph;
  }

  void CombineErrorsInQuaderature(const std::vector<std::string> &weightNames, const std::string &weightName)
  {
    for(auto const&  bin : _bins)
      {
        const double nom = bin->GetNominalXSec();

        const double bias = bin->GetFracSystBiasQuadSum(weightNames);
        const double res  = bin->GetFracSystResAveQuadSum(weightNames);

        bin->InsertSystFracError(weightName, (1 + bias) * nom, res);
      }
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

  std::string GetVar1BinString(const std::string var, const int var1Bin)
    {
      return Form("%.2f < %s < %.2f", _var1Bins[var1Bin-1], var.c_str(), _var1Bins[var1Bin]);
    }

  TH2D* CreateFractionalCovarianceMatrix(const std::string weightName)
  {
    TH2D *matrix = new TH2D(Form("FracCovMat%s", weightName.c_str()), ";Bin;Bin",
                            _bins.size(), 0.5, _bins.size() + 0.5,
                            _bins.size(), 0.5, _bins.size() + 0.5);

    for(auto&& [ binIndex_i, bin_i ] : enumerate(_bins))
      {
        const double xsec_nom_i = bin_i->GetNominalXSec();

        for(auto&& [ binIndex_j, bin_j ] : enumerate(_bins))
          {
            const int nunivs = bin_i->GetNUniverses(weightName);

            const double xsec_nom_j = bin_j->GetNominalXSec();

            double sum = 0.;

            for(int univ = 0; univ < nunivs; ++univ)
              {
                const double xsec_i = bin_i->GetUniverseXSec(weightName, univ);
                const double xsec_j = bin_j->GetUniverseXSec(weightName, univ);

                sum += (xsec_i - xsec_nom_i) * (xsec_j - xsec_nom_j);
              }

            matrix->SetBinContent(binIndex_i + 1, binIndex_j + 1, sum / (nunivs * xsec_nom_i * xsec_nom_j));
          }
      }

    return matrix;
  }

  TH2D* CreateCorrelationMatrix(const std::string weightName)
  {
    const double compfactor = 1;

    TH2D *cov_matrix = new TH2D(Form("CovMat%s", weightName.c_str()), ";Bin;Bin",
                                _bins.size(), 0.5, _bins.size() + 0.5,
                                _bins.size(), 0.5, _bins.size() + 0.5);

    for(auto&& [ binIndex_i, bin_i ] : enumerate(_bins))
      {
        const double xsec_nom_i = bin_i->GetNominalXSec();

        for(auto&& [ binIndex_j, bin_j ] : enumerate(_bins))
          {
            const int nunivs = bin_i->GetNUniverses(weightName);

            const double xsec_nom_j = bin_j->GetNominalXSec();

            double sum = 0.;

            for(int univ = 0; univ < nunivs; ++univ)
              {
                const double xsec_i = bin_i->GetUniverseXSec(weightName, univ);
                const double xsec_j = bin_j->GetUniverseXSec(weightName, univ);

                sum += (xsec_i - xsec_nom_i) * (xsec_j - xsec_nom_j);
              }

            cov_matrix->SetBinContent(binIndex_i + 1, binIndex_j + 1, (sum * compfactor) / nunivs);
          }
      }

    TH2D *corr_matrix = new TH2D(Form("CorrMat%s", weightName.c_str()), ";Bin;Bin",
                                 _bins.size(), 0.5, _bins.size() + 0.5,
                                 _bins.size(), 0.5, _bins.size() + 0.5);

    for(auto&& [ binIndex_i, bin_i ] : enumerate(_bins))
      {
        for(auto&& [ binIndex_j, bin_j ] : enumerate(_bins))
          {
            const double norm = TMath::Sqrt(cov_matrix->GetBinContent(binIndex_i + 1, binIndex_i + 1) * cov_matrix->GetBinContent(binIndex_j + 1, binIndex_j + 1));

            corr_matrix->SetBinContent(binIndex_i + 1, binIndex_j + 1, cov_matrix->GetBinContent(binIndex_i + 1, binIndex_j + 1) / norm);
          }
      }

    return corr_matrix;
  }

  TH1F* GetPredictedHist0D(const TString selName, const TString genName, /*const TString flavour, */const float scale = 1.)
  {
    //    TFile* xsecFile = new TFile("/exp/sbnd/data/users/hlay/ncpizero/generators/genie_xsec_" + flavour + ".root", "READ");
    TFile* foldFile = new TFile("/exp/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv17/forwardfoldingmatrices/forwardfoldingmatrices.root", "READ");

    TFile* xsecFileNuMu  = new TFile("/exp/sbnd/data/users/hlay/ncpizero/generators/genie_xsec_numu.root", "READ");
    TFile* xsecFileANuMu = new TFile("/exp/sbnd/data/users/hlay/ncpizero/generators/genie_xsec_anumu.root", "READ");
    TFile* xsecFileNuE   = new TFile("/exp/sbnd/data/users/hlay/ncpizero/generators/genie_xsec_nue.root", "READ");
    TFile* xsecFileANuE  = new TFile("/exp/sbnd/data/users/hlay/ncpizero/generators/genie_xsec_anue.root", "READ");

    //  TH1F* hist          = (TH1F*) xsecFile->Get(Form("%s_%s_%s", var.c_str(), selName.Data(), flavour.Data()));
    TH1F* histNuMu      = (TH1F*) xsecFileNuMu->Get(Form("%s_numu", selName.Data()));
    TH1F* histANuMu     = (TH1F*) xsecFileANuMu->Get(Form("%s_anumu", selName.Data()));
    TH1F* histNuE       = (TH1F*) xsecFileNuE->Get(Form("%s_nue", selName.Data()));
    TH1F* histANuE      = (TH1F*) xsecFileANuE->Get(Form("%s_anue", selName.Data()));

    histNuMu->Scale(numufluxscale);
    histANuMu->Scale(anumufluxscale);
    histNuE->Scale(nuefluxscale);
    histANuE->Scale(anuefluxscale);

    histNuMu->Add(histANuMu);
    histNuMu->Add(histNuE);
    histNuMu->Add(histANuE);

    histNuMu->Scale(scale);
    return histNuMu;
  }

  TH1F* GetPredictedHist1D(const TString selName, const TString genName, /*const TString flavour, */const float scale = 1.)
  {
    //    TFile* xsecFile = new TFile("/exp/sbnd/data/users/hlay/ncpizero/generators/genie_xsec_" + flavour + ".root", "READ");
    TFile* foldFile = new TFile("/exp/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv17/forwardfoldingmatrices/forwardfoldingmatrices.root", "READ");

    TFile* xsecFileNuMu  = new TFile("/exp/sbnd/data/users/hlay/ncpizero/generators/genie_xsec_numu.root", "READ");
    TFile* xsecFileANuMu = new TFile("/exp/sbnd/data/users/hlay/ncpizero/generators/genie_xsec_anumu.root", "READ");
    TFile* xsecFileNuE   = new TFile("/exp/sbnd/data/users/hlay/ncpizero/generators/genie_xsec_nue.root", "READ");
    TFile* xsecFileANuE  = new TFile("/exp/sbnd/data/users/hlay/ncpizero/generators/genie_xsec_anue.root", "READ");

    std::string var;

    if(_var0Bins.size() == 2)
      var = _var1;
    else if(_var1Bins.size() == 2)
      var = _var0;
    else
      throw std::runtime_error("This is not a 1D setup");

    //  TH1F* hist          = (TH1F*) xsecFile->Get(Form("%s_%s_%s", var.c_str(), selName.Data(), flavour.Data()));
    TH1F* histNuMu      = (TH1F*) xsecFileNuMu->Get(Form("%s_%s_numu", var.c_str(), selName.Data()));
    TH1F* histANuMu     = (TH1F*) xsecFileANuMu->Get(Form("%s_%s_anumu", var.c_str(), selName.Data()));
    TH1F* histNuE       = (TH1F*) xsecFileNuE->Get(Form("%s_%s_nue", var.c_str(), selName.Data()));
    TH1F* histANuE      = (TH1F*) xsecFileANuE->Get(Form("%s_%s_anue", var.c_str(), selName.Data()));
    TH2D* foldingmatrix = (TH2D*) foldFile->Get(Form("hForwardFold_%s_%s", var.c_str(), selName.Data()));

    histNuMu->Scale(numufluxscale);
    histANuMu->Scale(anumufluxscale);
    histNuE->Scale(nuefluxscale);
    histANuE->Scale(anuefluxscale);

    histNuMu->Add(histANuMu);
    histNuMu->Add(histNuE);
    histNuMu->Add(histANuE);

    //TH1F* foldedHist = Fold(histNuMu, foldingmatrix);

    TH1F* foldedHist = histNuMu;
    foldedHist->Scale(scale);
    return foldedHist;
  }

  TH1F* GetFracErrorHist0D(const std::string &name)
  {
    double var0Bins[_var0Bins.size()];
    GetVar0BinsArray(var0Bins);

    TH1F* hist = new TH1F(Form("ZeroD_%s_frac_error_%s", _name.c_str(), name.c_str()), _axes_labels.c_str(),
                          _var0Bins.size() - 1, var0Bins);

    hist->SetBinContent(1, _bins[0]->GetFracSystResAve(name));

    return hist;
  }

  TH1F* GetFracErrorHist1D(const std::string &name)
  {
    if(_var0Bins.size() != 2 && _var1Bins.size() != 2)
      throw std::runtime_error("Asking for 1D hist but neither set has 1 bin");

    if(_var0Bins.size() == 2)
      {
        double var1Bins[_var1Bins.size()];
        GetVar1BinsArray(var1Bins);

        TH1F* hist = new TH1F(Form("OneD_%s_%s", _name.c_str(), name.c_str()), _axes_labels.c_str(),
                              _var1Bins.size() - 1, var1Bins);

        for(auto&& [ bin_i, bin ] : enumerate(_bins))
          {
            hist->SetBinContent(bin_i + 1, bin->GetFracSystResAve(name));
            hist->SetBinError(bin_i + 1, 0);
          }

        hist->SetMarkerStyle(0);
        hist->SetLineColor(kBlack);

        return hist;
      }

    if(_var1Bins.size() == 2)
      {
        double var0Bins[_var0Bins.size()];
        GetVar0BinsArray(var0Bins);

        TH1F* hist = new TH1F(Form("OneD_%s_%s", _name.c_str(), name.c_str()), _axes_labels.c_str(),
                              _var0Bins.size() - 1, var0Bins);

        for(auto&& [ bin_i, bin ] : enumerate(_bins))
          {
            hist->SetBinContent(bin_i + 1, bin->GetFracSystResAve(name));
            hist->SetBinError(bin_i + 1, 0);
          }

        hist->SetMarkerStyle(0);
        hist->SetLineColor(kBlack);

        return hist;
      }

    return NULL;
  }
};

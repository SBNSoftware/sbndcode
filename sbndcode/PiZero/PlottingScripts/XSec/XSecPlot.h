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

// Flux Config L Effective Z (11227.8cm) - Subset
/*
  constexpr double numufluxscale  = 0.92364301;
  constexpr double anumufluxscale = 0.070125432;
  constexpr double nuefluxscale   = 0.0056067340;
  constexpr double anuefluxscale  = 0.00062482531;
*/

// Flux Config L Ray-traced
constexpr double numufluxscale  = 0.92536980;
constexpr double anumufluxscale = 0.068619804;
constexpr double nuefluxscale   = 0.0054254241;
constexpr double anuefluxscale  = 0.00058497585;

/*
  const std::map<TString, double> extraSignalScaling = { { "ncpizero_incl", 1.01763 },
  { "ncpizero_0p0pi", 1.01166 },
  { "ncpizero_Np0pi", 1.0251 }
  };
*/

const std::map<TString, double> extraSignalScaling = { { "ncpizero_incl", 1. },
                                                       { "ncpizero_0p0pi", 1. },
                                                       { "ncpizero_Np0pi", 1. }
};

const TString foldFileName   = "/exp/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroBv3/forwardfoldingmatrices/forwardfoldingmatrices.root";
const TString unfoldFileName = "/exp/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroBv3/unfoldingmatrices/unfoldingmatrices.root";

class XSecPlot {

private:
  std::string         _name;
  std::string         _axes_labels;
  std::vector<Bin*>   _bins;
  std::vector<double> _varBins;
  std::string         _var;

public:

  XSecPlot(const std::string name, const std::string axes_labels, const int N,
           const std::vector<double> varBins, const std::string var,
           const double nTargets, const double intFlux)
  {
    _name        = name;
    _axes_labels = axes_labels;
    _varBins     = varBins;
    _var         = var;

    if(varBins.size() - 1 != N)
      throw std::runtime_error("Bin Sizes Messed Up");

    for(int i = 1; i < varBins.size(); ++i)
      {
        const double width = varBins.size() == 2 ? 1. : varBins[i] - varBins[i-1];
        _bins.emplace_back(new Bin(i - 1, width, nTargets, intFlux, varBins[i-1], varBins[i]));
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

  std::string GetVar()
  {
    return _var;
  }

  std::vector<float> GetNominalXSecs()
  {
    std::vector<float> values;

    for(auto const& bin : _bins)
      values.push_back(bin->GetNominalXSec());

    return values;
  }

  std::vector<float> GetUniverseXSecs(const std::string &weightName, const int univ)
  {
    std::vector<float> values;

    for(auto const& bin : _bins)
      values.push_back(bin->GetUniverseXSec(weightName, univ));

    return values;
  }

  std::vector<float> GetFracErrorsXSec(const std::string &name)
  {
    std::vector<float> values;

    for(auto const& bin : _bins)
      values.push_back(bin->GetFracSystResAveXSec(name));

    return values;
  }

  std::vector<float> GetFracErrorsEfficiency(const std::string &name)
  {
    std::vector<float> values;

    for(auto const& bin : _bins)
      values.push_back(bin->GetFracSystResAveEfficiency(name));

    return values;
  }

  std::vector<float> GetFracErrorsPurity(const std::string &name)
  {
    std::vector<float> values;

    for(auto const& bin : _bins)
      values.push_back(bin->GetFracSystResAvePurity(name));

    return values;
  }

  std::vector<float> GetFracErrorsBkgdCount(const std::string &name, const bool scale = false)
  {
    std::vector<float> values;

    for(auto const& bin : _bins)
      values.push_back(bin->GetFracSystResAveBkgdCount(name, scale));

    return values;
  }

  std::vector<float> GetUnderlyingMCXSecs()
  {
    std::vector<float> values;

    for(auto const& bin : _bins)
      values.push_back(bin->GetUnderlyingMCXSec());

    return values;
  }

  std::vector<float> GetNominalEfficiencies()
  {
    std::vector<float> values;

    for(auto const& bin : _bins)
      values.push_back(bin->GetNominalEfficiency());

    return values;
  }

  std::vector<float> GetNominalStatErrors(const std::vector<float> &values)
  {
    std::vector<float> errors;

    for(int i = 0; i < _bins.size(); ++i)
      errors.push_back(_bins[i]->GetNominalFracStatErr() * values[i]);

    return errors;
  }

  std::vector<float> GetUniverseEfficiencies(const std::string &weightName, const int univ)
  {
    std::vector<float> values;

    for(auto const& bin : _bins)
      values.push_back(bin->GetUniverseEfficiency(weightName, univ));

    return values;
  }

  std::vector<float> GetNominalPurities()
  {
    std::vector<float> values;

    for(auto const& bin : _bins)
      values.push_back(bin->GetNominalPurity());

    return values;
  }

  std::vector<float> GetUniversePurities(const std::string &weightName, const int univ)
  {
    std::vector<float> values;

    for(auto const& bin : _bins)
      values.push_back(bin->GetUniversePurity(weightName, univ));

    return values;
  }

  std::vector<float> GetNominalBkgdCounts(const bool scale = false)
  {
    std::vector<float> values;

    for(auto const& bin : _bins)
      values.push_back(bin->GetNominalBkgdCount(scale));

    return values;
  }

  std::vector<float> GetUniverseBkgdCounts(const std::string &weightName, const int univ, const bool scale = false)
  {
    std::vector<float> values;

    for(auto const& bin : _bins)
      values.push_back(bin->GetUniverseBkgdCount(weightName, univ, scale));

    return values;
  }

  void InsertSystFracFlatError(const std::string &name, const double &frac)
  {
    for(int i = 0; i < _bins.size(); ++i)
      {
        const double nomxsec   = _bins[i]->GetNominalXSec();
        const double nomeff    = _bins[i]->GetNominalEfficiency();
        const double nompur    = _bins[i]->GetNominalPurity();
        const double nombcount = _bins[i]->GetNominalBkgdCount();

        _bins[i]->InsertSystFracErrorXSec(name, nomxsec, frac);
        _bins[i]->InsertSystFracErrorEfficiency(name, nomeff, frac);
        _bins[i]->InsertSystFracErrorPurity(name, nompur, frac);
        _bins[i]->InsertSystFracErrorBkgdCount(name, nombcount, frac);
      }
  }

  TH1F* MakeHist(const TString &name, const std::vector<float> &values, const std::vector<float> &errors)
  {
    double varBins[_varBins.size()];
    GetVarBinsArray(varBins);

    if(values.size() != _bins.size())
      std::runtime_error("Bin sizes in compatible");

    TH1F* hist = new TH1F(name, _axes_labels.c_str(), _varBins.size() - 1, varBins);

    for(int i = 0; i < _bins.size(); ++i)
      {
        hist->SetBinContent(i + 1, values[i]);
        hist->SetBinError(i + 1, errors[i]);
      }

    return hist;
  }

  TH1F* GetNominalXSecHist(const bool statErr = true, const bool unfold = false, const TString selName = "")
  {
    if(_varBins.size() == 2 && unfold)
      std::runtime_error("Can't unfold 0D result");

    const std::vector<float> values = GetNominalXSecs();
    const std::vector<float> errors = statErr ? GetNominalStatErrors(values) : std::vector<float>(_bins.size(), 0.);

    TH1F* hist = MakeHist(Form("Nominal_XSec_%s", _name.c_str()), values, errors);

    TFile* foldFile = new TFile(unfoldFileName, "READ");
    TH2D* foldingmatrix = (TH2D*) foldFile->Get(Form("hUnfold_%s_%s", _var.c_str(), selName.Data()));

    TH1F *effHist = GetNominalEfficiencyHist();
    TH1F* unfoldedHist = unfold ? Fold(hist, effHist, foldingmatrix) : hist;

    unfoldedHist->SetMarkerStyle(0);
    unfoldedHist->SetLineColor(kBlack);

    return unfoldedHist;
  }

  TH1F* GetNominalEfficiencyHist(const bool statErr = true)
  {
    const std::vector<float> values = GetNominalEfficiencies();
    const std::vector<float> errors = statErr ? GetNominalStatErrors(values) : std::vector<float>(_bins.size(), 0.);

    TH1F* hist = MakeHist(Form("Nominal_Efficiency_%s", _name.c_str()), values, errors);

    return hist;
  }

  TH1F* GetNominalPurityHist(const bool statErr = true)
  {
    const std::vector<float> values = GetNominalPurities();
    const std::vector<float> errors = statErr ? GetNominalStatErrors(values) : std::vector<float>(_bins.size(), 0.);

    TH1F* hist = MakeHist(Form("Nominal_Purity_%s", _name.c_str()), values, errors);

    return hist;
  }

  TH1F* GetNominalBkgdCountHist(const bool statErr = true, const bool scale = false)
  {
    const std::vector<float> values = GetNominalBkgdCounts(scale);
    const std::vector<float> errors = statErr ? GetNominalStatErrors(values) : std::vector<float>(_bins.size(), 0.);

    TH1F* hist = MakeHist(Form("Nominal_BkgdCount_%s", _name.c_str()), values, errors);

    return hist;
  }

  TH1F* GetUniverseXSecHist(const std::string &weightName, const int univ)
  {
    const std::vector<float> values = GetUniverseXSecs(weightName, univ);
    const std::vector<float> errors = std::vector<float>(_bins.size(), 0.);

    TH1F* hist = MakeHist(Form("Universe_XSec_%s_%s_%d", _name.c_str(), weightName.c_str(), univ), values, errors);

    hist->SetMarkerStyle(0);
    hist->SetLineColor(kMagenta-10);
    hist->SetLineWidth(1);

    return hist;
  }

  TH1F* GetUniverseEfficiencyHist(const std::string &weightName, const int univ)
  {
    const std::vector<float> values = GetUniverseEfficiencies(weightName, univ);
    const std::vector<float> errors = std::vector<float>(_bins.size(), 0.);

    TH1F* hist = MakeHist(Form("Universe_Efficiency_%s_%s_%d", _name.c_str(), weightName.c_str(), univ), values, errors);

    hist->SetMarkerStyle(0);
    hist->SetLineColor(kMagenta-10);
    hist->SetLineWidth(1);

    return hist;
  }

  TH1F* GetUniversePurityHist(const std::string &weightName, const int univ)
  {
    const std::vector<float> values = GetUniversePurities(weightName, univ);
    const std::vector<float> errors = std::vector<float>(_bins.size(), 0.);

    TH1F* hist = MakeHist(Form("Universe_Purity_%s_%s_%d", _name.c_str(), weightName.c_str(), univ), values, errors);

    hist->SetMarkerStyle(0);
    hist->SetLineColor(kMagenta-10);
    hist->SetLineWidth(1);

    return hist;
  }

  TH1F* GetUniverseBkgdCountHist(const std::string &weightName, const int univ, const bool scale = false)
  {
    const std::vector<float> values = GetUniverseBkgdCounts(weightName, univ, scale);
    const std::vector<float> errors = std::vector<float>(_bins.size(), 0.);

    TH1F* hist = MakeHist(Form("Universe_BkgdCount_%s_%s_%d", _name.c_str(), weightName.c_str(), univ), values, errors);

    hist->SetMarkerStyle(0);
    hist->SetLineColor(kMagenta-10);
    hist->SetLineWidth(1);

    return hist;
  }

  TH1F* GetFracErrorXSecHist(const std::string &name)
  {
    const std::vector<float> values = GetFracErrorsXSec(name);
    const std::vector<float> errors = std::vector<float>(_bins.size(), 0.);

    TH1F* hist = MakeHist(Form("Frac_Error_XSec_%s_%s", _name.c_str(), name.c_str()), values, errors);

    hist->SetMarkerStyle(0);
    hist->SetLineColor(kBlack);

    return hist;
  }

  TH1F* GetUnderlyingMCXSecHist(const TString selName = "", const bool fold = false)
  {
    if(_varBins.size() == 2 && fold)
      std::runtime_error("Can't fold 0D result");

    const std::vector<float> values = GetUnderlyingMCXSecs();
    const std::vector<float> errors = std::vector<float>(_bins.size(), 0.);

    TH1F* hist = MakeHist(Form("Underlying_MC_XSec_%s", _name.c_str()), values, errors);

    TFile* foldFile = new TFile(foldFileName, "READ");
    TH2D* foldingmatrix = (TH2D*) foldFile->Get(Form("hForwardFold_%s_%s", _var.c_str(), selName.Data()));

    TH1F *effHist = GetNominalEfficiencyHist();
    TH1F* foldedHist = fold ? Fold(hist, effHist, foldingmatrix) : hist;

    foldedHist->SetMarkerStyle(0);
    foldedHist->SetLineColor(kBlack);

    return foldedHist;
  }

  TGraphAsymmErrors* GetCVErrXSecGraph(const std::string &weightName)
  {
    TGraphAsymmErrors* graph = new TGraphAsymmErrors();

    for(auto&& [ bin_i, bin ] : enumerate(_bins))
      {
        double low, cv, high;
        std::tie(low, cv, high) = bin->GetSystFracErrorsXSec(weightName);

        if(cv == 0.)
          continue;

        const double centre = bin->GetVarCenter();
        const double width  = bin->GetVarWidth();

        graph->SetPoint(bin_i, centre, cv);
        graph->SetPointError(bin_i, 0.49 * width, 0.49 * width, low * cv, high * cv);

      }

    graph->SetMarkerStyle(1);
    graph->SetLineColor(kBlue+2);
    graph->SetLineWidth(5);

    return graph;
  }

  TGraphAsymmErrors* GetCVErrEfficiencyGraph(const std::string &weightName)
  {
    TGraphAsymmErrors* graph = new TGraphAsymmErrors();

    for(auto&& [ bin_i, bin ] : enumerate(_bins))
      {
        double low, cv, high;
        std::tie(low, cv, high) = bin->GetSystFracErrorsEfficiency(weightName);

        if(cv == 0.)
          continue;

        const double centre = bin->GetVarCenter();
        const double width  = bin->GetVarWidth();

        graph->SetPoint(bin_i, centre, cv);
        graph->SetPointError(bin_i, 0.49 * width, 0.49 * width, low * cv, high * cv);

      }

    graph->SetMarkerStyle(1);
    graph->SetLineColor(kBlue+2);
    graph->SetLineWidth(5);

    return graph;
  }

  TGraphAsymmErrors* GetCVErrPurityGraph(const std::string &weightName)
  {
    TGraphAsymmErrors* graph = new TGraphAsymmErrors();

    for(auto&& [ bin_i, bin ] : enumerate(_bins))
      {
        double low, cv, high;
        std::tie(low, cv, high) = bin->GetSystFracErrorsPurity(weightName);

        if(cv == 0.)
          continue;

        const double centre = bin->GetVarCenter();
        const double width  = bin->GetVarWidth();

        graph->SetPoint(bin_i, centre, cv);
        graph->SetPointError(bin_i, 0.49 * width, 0.49 * width, low * cv, high * cv);

      }

    graph->SetMarkerStyle(1);
    graph->SetLineColor(kBlue+2);
    graph->SetLineWidth(5);

    return graph;
  }

  TGraphAsymmErrors* GetCVErrBkgdCountGraph(const std::string &weightName, const bool scale = false)
  {
    TGraphAsymmErrors* graph = new TGraphAsymmErrors();

    for(auto&& [ bin_i, bin ] : enumerate(_bins))
      {
        double low, cv, high;
        std::tie(low, cv, high) = bin->GetSystFracErrorsBkgdCount(weightName, scale);

        if(cv == 0.)
          continue;

        const double centre = bin->GetVarCenter();
        const double width  = bin->GetVarWidth();

        graph->SetPoint(bin_i, centre, cv);
        graph->SetPointError(bin_i, 0.49 * width, 0.49 * width, low * cv, high * cv);

      }

    graph->SetMarkerStyle(1);
    graph->SetLineColor(kBlue+2);
    graph->SetLineWidth(5);

    return graph;
  }

  void CombineErrorsInQuaderature(const std::vector<std::string> &weightNames, const std::string &weightName)
  {
    for(auto const&  bin : _bins)
      {
        const double nomxsec   = bin->GetNominalXSec();
        const double nomeff    = bin->GetNominalEfficiency();
        const double nompur    = bin->GetNominalPurity();
        const double nombcount = bin->GetNominalBkgdCount();

        const double biasxsec = bin->GetFracSystBiasQuadSumXSec(weightNames);
        const double resxsec  = bin->GetFracSystResAveQuadSumXSec(weightNames);

        const double biaseff = bin->GetFracSystBiasQuadSumEfficiency(weightNames);
        const double reseff  = bin->GetFracSystResAveQuadSumEfficiency(weightNames);

        const double biaspur = bin->GetFracSystBiasQuadSumPurity(weightNames);
        const double respur  = bin->GetFracSystResAveQuadSumPurity(weightNames);

        const double biasbcount = bin->GetFracSystBiasQuadSumBkgdCount(weightNames);
        const double resbcount  = bin->GetFracSystResAveQuadSumBkgdCount(weightNames);

        bin->InsertSystFracErrorXSec(weightName, (1 + biasxsec) * nomxsec, resxsec);
        bin->InsertSystFracErrorEfficiency(weightName, (1 + biaseff) * nomeff, reseff);
        bin->InsertSystFracErrorPurity(weightName, (1 + biaspur) * nompur, respur);
        bin->InsertSystFracErrorBkgdCount(weightName, (1 + biasbcount) * nombcount, resbcount);
      }
  }

  void GetVarBinsArray(double varBins[])
  {
    for(auto&& [ binEdge_i, binEdge ] : enumerate(_varBins))
      varBins[binEdge_i] = binEdge;
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

            cov_matrix->SetBinContent(binIndex_i + 1, binIndex_j + 1, sum / nunivs);
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

  TH1F* GetPredictedXSecHist(const TString selName, const TString genName, const float scale = 1., const bool fold = false)
  {
    if(_varBins.size() == 2 && fold)
      std::runtime_error("Can't fold 0D result");

    TFile* xsecFileNuMu  = new TFile("/exp/sbnd/data/users/hlay/ncpizero/generators/genie_xsec_numu.root", "READ");
    TFile* xsecFileANuMu = new TFile("/exp/sbnd/data/users/hlay/ncpizero/generators/genie_xsec_anumu.root", "READ");
    TFile* xsecFileNuE   = new TFile("/exp/sbnd/data/users/hlay/ncpizero/generators/genie_xsec_nue.root", "READ");
    TFile* xsecFileANuE  = new TFile("/exp/sbnd/data/users/hlay/ncpizero/generators/genie_xsec_anue.root", "READ");

    TH1F* histNuMu;
    TH1F* histANuMu;
    TH1F* histNuE;
    TH1F* histANuE;

    if(_varBins.size() == 2)
      {
        histNuMu  = (TH1F*) xsecFileNuMu->Get(Form("%s_numu", selName.Data()));
        histANuMu = (TH1F*) xsecFileANuMu->Get(Form("%s_anumu", selName.Data()));
        histNuE   = (TH1F*) xsecFileNuE->Get(Form("%s_nue", selName.Data()));
        histANuE  = (TH1F*) xsecFileANuE->Get(Form("%s_anue", selName.Data()));
      }
    else
      {
        histNuMu  = (TH1F*) xsecFileNuMu->Get(Form("%s_%s_numu", _var.c_str(), selName.Data()));
        histANuMu = (TH1F*) xsecFileANuMu->Get(Form("%s_%s_anumu", _var.c_str(), selName.Data()));
        histNuE   = (TH1F*) xsecFileNuE->Get(Form("%s_%s_nue", _var.c_str(), selName.Data()));
        histANuE  = (TH1F*) xsecFileANuE->Get(Form("%s_%s_anue", _var.c_str(), selName.Data()));
      }

    histNuMu->Scale(numufluxscale);
    histANuMu->Scale(anumufluxscale);
    histNuE->Scale(nuefluxscale);
    histANuE->Scale(anuefluxscale);

    histNuMu->Add(histANuMu);
    histNuMu->Add(histNuE);
    histNuMu->Add(histANuE);

    TFile* foldFile = new TFile(foldFileName, "READ");
    TH2D* foldingmatrix = (TH2D*) foldFile->Get(Form("hForwardFold_%s_%s", _var.c_str(), selName.Data()));
    TH1F *effHist = GetNominalEfficiencyHist();
    TH1F* foldedHist = fold ? Fold(histNuMu, effHist, foldingmatrix) : histNuMu;

    foldedHist->Scale(scale);
    return foldedHist;
  }
};

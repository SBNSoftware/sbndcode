#pragma once

#include "/exp/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroBv1/integrated_flux_eff_z/FluxMap.h"
#include "UniverseBin.h"

class Bin {

private:
  int                                                        _index;
  UniverseBin                                               *_nominalBin;
  std::map<std::string, std::vector<UniverseBin*>>          *_universeBins;
  std::map<std::string, std::tuple<double, double, double>> *_systFracErrors;
  double                                                     _varLow;
  double                                                     _varHigh;

public:
  
  Bin(const int index, const double binWidth, const double nTargets, const double intFlux,
      const double varLow, const double varHigh)
  {
    _index          = index;
    _nominalBin     = new UniverseBin(binWidth, nTargets, intFlux);
    _universeBins   = new std::map<std::string, std::vector<UniverseBin*>>();
    _systFracErrors = new std::map<std::string, std::tuple<double, double, double>>();
    _varLow         = varLow;
    _varHigh        = varHigh;
  }

  void Print()
  {
    std::cout << "\n======================================\n"
              << "Bin Index: " << _index << '\n'
              << "Var: " << _varLow << " --> " << _varHigh  << '\n';

    _nominalBin->Print();

    std::cout <<"======================================\n" << std::endl;
  }

  void Print(const std::string weightName, const int univ)
  {
    std::cout << "\n======================================\n"
              << "Bin Index: " << _index << " (" << weightName << ", " << univ << ")" << '\n'
              << "Var: " << _varLow << " --> " << _varHigh  << '\n';

    _nominalBin->Print();
    _universeBins->at(weightName).at(univ)->Print();

    std::cout <<"======================================\n" << std::endl;
  }

  void Update()
  {
    _nominalBin->Update();
    
    for(auto&& [ name, univBins ] : *_universeBins)
      {
        for(UniverseBin* bin : univBins)
          bin->Update();
      }
  }

  void CalculateXSecPurity()
  {
    _nominalBin->CalculateXSecPurity();
    
    for(auto&& [ name, univBins ] : *_universeBins)
      {
        for(UniverseBin* bin : univBins)
          bin->CalculateXSecPurity(_nominalBin->GetCount());
      }
  }

  void SetScaleFactor(const double &scaleFactor)
  {
    _nominalBin->SetScaleFactor(scaleFactor);

    for(auto&& [ name, univBins ] : *_universeBins)
      {
        for(UniverseBin* bin : univBins)
          bin->SetScaleFactor(scaleFactor);
      }
  }

  double GetNominalXSec()
  {
    return _nominalBin->GetXSec();
  }

  double GetUnderlyingMCXSec()
  {
    return _nominalBin->GetUnderlyingMCXSec();
  }

  double GetNominalFracStatErr()
  {
    return _nominalBin->GetFracStatErr();
  }

  double GetNominalEfficiency()
  {
    return _nominalBin->GetEfficiency();
  }

  void AddWeight(const std::string weightName, const int nunivs)
  {
    const double binWidth = _nominalBin->GetBinWidth();
    const double nTargets = _nominalBin->GetNTargets();
    double intFlux  = _nominalBin->GetIntFlux();

    _universeBins->insert({weightName, std::vector<UniverseBin*>()});

    for(int univ = 0; univ < nunivs; ++univ)
      {
        if(univsIntegratedFluxMap.find(weightName) !=univsIntegratedFluxMap.end() &&
           univ < univsIntegratedFluxMap.at(weightName).size())
          intFlux = univsIntegratedFluxMap.at(weightName).at(univ);
      
        _universeBins->at(weightName).push_back(new UniverseBin(binWidth, nTargets, intFlux));
      }
  }

  double GetUniverseXSec(const std::string &weightName, const int univ)
  {
    return _universeBins->at(weightName).at(univ)->GetXSec();
  }

  double GetNUniverses(const std::string &weightName)
  {
    return _universeBins->at(weightName).size();
  }

  void IncrementNominalBinCount(const double &increment)
  {
    _nominalBin->IncrementCount(increment);
  }

  void IncrementNominalBinBkgdCount(const double &increment)
  {
    _nominalBin->IncrementBkgdCount(increment);
  }

  void IncrementNominalBinTrueSignal(const double &increment)
  {
    _nominalBin->IncrementTrueSignal(increment);
  }

  void IncrementNominalBinSelSignalTrueBin(const double &increment)
  {
    _nominalBin->IncrementSelSignalTrueBin(increment);
  }

  void IncrementUniverseBinCount(const std::string &weightName, const int univ, 
                                 const double &increment)
  {
    _universeBins->at(weightName).at(univ)->IncrementCount(increment);
  }

  void IncrementUniverseBinBkgdCount(const std::string &weightName, const int univ, 
                                     const double &increment)
  {
    _universeBins->at(weightName).at(univ)->IncrementBkgdCount(increment);
  }

  void IncrementUniverseBinTrueSignal(const std::string &weightName, const int univ,
                                      const double &increment)
  {
    _universeBins->at(weightName).at(univ)->IncrementTrueSignal(increment);
  }

  void IncrementUniverseBinSelSignalTrueBin(const std::string &weightName, const int univ,
                                            const double &increment)
  {
    _universeBins->at(weightName).at(univ)->IncrementSelSignalTrueBin(increment);
  }

  void CalculateSystFracErrorsMedianPercentiles()
  {
    for(auto&& [ name, univBins ] : *_universeBins)
      {
        std::vector<double> xsecs;

        for(UniverseBin* bin : univBins)
          xsecs.push_back(bin->GetXSec());

        std::sort(xsecs.begin(), xsecs.end());

        const int n         = univBins.size();
        const double cv_n   = n * (1/2.);
        const double low_n  = n * (1/6.);
        const double high_n = n * (5/6.);

        double cv   = (xsecs[std::floor(cv_n)] + xsecs[std::ceil(cv_n)]) / 2.;
        double low  = (xsecs[std::floor(low_n)] + xsecs[std::ceil(low_n)]) / 2.;
        double high = (xsecs[std::floor(high_n)] + xsecs[std::ceil(high_n)]) / 2.;

        low  = std::abs(low - cv) / cv;
        high = std::abs(high - cv) / cv;

        _systFracErrors->insert({ name, { low, cv, high }});
      }
  }

  void CalculateSystFracErrorsNominalSD()
  {
    for(auto&& [ name, univBins ] : *_universeBins)
      {
        const double nom = GetNominalXSec();

        double sum = 0.;

        for(UniverseBin* bin : univBins)
          sum += TMath::Power((bin->GetXSec() - nom), 2);

        sum /= univBins.size();

        const double sd     = TMath::Sqrt(sum);
        const double fracsd = sd / nom;

        _systFracErrors->insert({ name, { fracsd, nom, fracsd }});
      }
  }

  std::tuple<double, double, double> GetSystFracErrors(const std::string &weightName)
  {
    return _systFracErrors->at(weightName);
  }

  void InsertSystFracError(const std::string &name, const double &centre, const double &frac)
  {
    _systFracErrors->insert({ name, { frac, centre, frac }});
  }

  bool InBin(const double &val)
  {
    return val > _varLow && val < _varHigh;
  }

  double GetFracSystResAve(const std::string &weightName)
  {
    double low, cv, high;

    std::tie(low, cv, high) = GetSystFracErrors(weightName);

    const double nom = _nominalBin->GetXSec();

    return (low + high) / 2.;
  }

  double GetFracSystBias(const std::string &weightName)
  {
    double low, cv, high;

    std::tie(low, cv, high) = GetSystFracErrors(weightName);

    const double nom = _nominalBin->GetXSec();
    return (cv - nom) / nom;
  }

  double GetFracSystResAveQuadSum(const std::vector<std::string> &weightNames)
  {
    double quadSum = 0.;
    for(auto const& weightName : weightNames)
      quadSum += TMath::Power(GetFracSystResAve(weightName), 2);

    return TMath::Sqrt(quadSum);
  }

  double GetFracSystBiasQuadSum(const std::vector<std::string> &weightNames)
  {
    double quadSum = 0.;
    for(auto const& weightName : weightNames)
      quadSum += TMath::Power(GetFracSystBias(weightName), 2);

    return TMath::Sqrt(quadSum);
  }

  double GetVarCenter()
  {
    return (_varLow + _varHigh) / 2.;
  }

  double GetVarWidth()
  {
    return _varHigh - _varLow;
  }
};

#pragma once

#include "/exp/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroBv1/integrated_flux_eff_z/FluxMap.h"
#include "UniverseBin.h"

class Bin {

private:
  int                                                        _index;
  UniverseBin                                               *_nominalBin;
  std::map<std::string, std::vector<UniverseBin*>>          *_universeBins;
  std::map<std::string, std::tuple<double, double, double>> *_systFracErrorsXSec;
  std::map<std::string, std::tuple<double, double, double>> *_systFracErrorsEfficiency;
  std::map<std::string, std::tuple<double, double, double>> *_systFracErrorsPurity;
  std::map<std::string, std::tuple<double, double, double>> *_systFracErrorsBkgdCount;
  double                                                     _varLow;
  double                                                     _varHigh;

public:
  
  Bin(const int index, const double binWidth, const double nTargets, const double intFlux,
      const double varLow, const double varHigh)
  {
    _index                    = index;
    _nominalBin               = new UniverseBin(binWidth, nTargets, intFlux);
    _universeBins             = new std::map<std::string, std::vector<UniverseBin*>>();
    _systFracErrorsXSec       = new std::map<std::string, std::tuple<double, double, double>>();
    _systFracErrorsEfficiency = new std::map<std::string, std::tuple<double, double, double>>();
    _systFracErrorsPurity     = new std::map<std::string, std::tuple<double, double, double>>();
    _systFracErrorsBkgdCount  = new std::map<std::string, std::tuple<double, double, double>>();
    _varLow                   = varLow;
    _varHigh                  = varHigh;
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

  double GetNominalPurity()
  {
    return _nominalBin->GetPurity();
  }

  double GetNominalBkgdCount()
  {
    return _nominalBin->GetBkgdCount();
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

  double GetUniverseEfficiency(const std::string &weightName, const int univ)
  {
    return _universeBins->at(weightName).at(univ)->GetEfficiency();
  }

  double GetUniversePurity(const std::string &weightName, const int univ)
  {
    return _universeBins->at(weightName).at(univ)->GetPurity();
  }

  double GetUniverseBkgdCount(const std::string &weightName, const int univ)
  {
    return _universeBins->at(weightName).at(univ)->GetBkgdCount();
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


  std::tuple<double, double, double> MedianPercentiles(std::vector<double> &values)
  {
    std::sort(values.begin(), values.end());

    const int n         = values.size();
    const double cv_n   = n * (1/2.);
    const double low_n  = n * (1/6.);
    const double high_n = n * (5/6.);

    double cv   = (values[std::floor(cv_n)] + values[std::ceil(cv_n)]) / 2.;
    double low  = (values[std::floor(low_n)] + values[std::ceil(low_n)]) / 2.;
    double high = (values[std::floor(high_n)] + values[std::ceil(high_n)]) / 2.;

    low  = std::abs(low - cv) / cv;
    high = std::abs(high - cv) / cv;

    return { low, cv, high };
  }

  void CalculateSystFracErrorsMedianPercentiles()
  {
    for(auto&& [ name, univBins ] : *_universeBins)
      {
        std::vector<double> xsecs;
        std::vector<double> effs;
        std::vector<double> purs;
        std::vector<double> bcounts;

        for(UniverseBin* bin : univBins)
          {
            xsecs.push_back(bin->GetXSec());
            effs.push_back(bin->GetEfficiency());
            purs.push_back(bin->GetPurity());
            bcounts.push_back(bin->GetBkgdCount());
          }

        _systFracErrorsXSec->insert({ name, MedianPercentiles(xsecs)});
        _systFracErrorsEfficiency->insert({ name, MedianPercentiles(effs)});
        _systFracErrorsPurity->insert({ name, MedianPercentiles(purs)});
        _systFracErrorsBkgdCount->insert({ name, MedianPercentiles(bcounts)});
      }
  }

  void CalculateSystFracErrorsNominalSD()
  {
    for(auto&& [ name, univBins ] : *_universeBins)
      {
        const double nomxsec   = GetNominalXSec();
        const double nomeff    = GetNominalEfficiency();
        const double nompur    = GetNominalPurity();
        const double nombcount = GetNominalBkgdCount();

        double sumxsec   = 0.;
        double sumeff    = 0.;
        double sumpur    = 0.;
        double sumbcount = 0.;

        for(UniverseBin* bin : univBins)
          {
            sumxsec   += TMath::Power((bin->GetXSec() - nomxsec), 2);
            sumeff    += TMath::Power((bin->GetEfficiency() - nomeff), 2);
            sumpur    += TMath::Power((bin->GetPurity() - nompur), 2);
            sumbcount += TMath::Power((bin->GetBkgdCount() - nombcount), 2);
          }

        sumxsec   /= univBins.size();
        sumeff    /= univBins.size();
        sumpur    /= univBins.size();
        sumbcount /= univBins.size();

        const double sdxsec     = TMath::Sqrt(sumxsec);
        const double fracsdxsec = sdxsec / nomxsec;

        const double sdeff     = TMath::Sqrt(sumeff);
        const double fracsdeff = sdeff / nomeff;

        const double sdpur     = TMath::Sqrt(sumpur);
        const double fracsdpur = sdpur / nompur;

        const double sdbcount     = TMath::Sqrt(sumbcount);
        const double fracsdbcount = sdbcount / nombcount;

        _systFracErrorsXSec->insert({ name, { fracsdxsec, nomxsec, fracsdxsec }});
        _systFracErrorsEfficiency->insert({ name, { fracsdeff, nomeff, fracsdeff }});
        _systFracErrorsPurity->insert({ name, { fracsdpur, nompur, fracsdpur }});
        _systFracErrorsBkgdCount->insert({ name, { fracsdbcount, nombcount, fracsdbcount }});
      }
  }

  std::tuple<double, double, double> GetSystFracErrorsXSec(const std::string &weightName)
  {
    return _systFracErrorsXSec->at(weightName);
  }

  std::tuple<double, double, double> GetSystFracErrorsEfficiency(const std::string &weightName)
  {
    return _systFracErrorsEfficiency->at(weightName);
  }

  std::tuple<double, double, double> GetSystFracErrorsPurity(const std::string &weightName)
  {
    return _systFracErrorsPurity->at(weightName);
  }

  std::tuple<double, double, double> GetSystFracErrorsBkgdCount(const std::string &weightName)
  {
    return _systFracErrorsBkgdCount->at(weightName);
  }

  void InsertSystFracErrorXSec(const std::string &name, const double &centre, const double &frac)
  {
    _systFracErrorsXSec->insert({ name, { frac, centre, frac }});
  }

  void InsertSystFracErrorEfficiency(const std::string &name, const double &centre, const double &frac)
  {
    _systFracErrorsEfficiency->insert({ name, { frac, centre, frac }});
  }

  void InsertSystFracErrorPurity(const std::string &name, const double &centre, const double &frac)
  {
    _systFracErrorsPurity->insert({ name, { frac, centre, frac }});
  }

  void InsertSystFracErrorBkgdCount(const std::string &name, const double &centre, const double &frac)
  {
    _systFracErrorsBkgdCount->insert({ name, { frac, centre, frac }});
  }

  bool InBin(const double &val)
  {
    return val > _varLow && val < _varHigh;
  }

  double GetFracSystResAveXSec(const std::string &weightName)
  {
    double low, cv, high;

    std::tie(low, cv, high) = GetSystFracErrorsXSec(weightName);

    const double nom = _nominalBin->GetXSec();

    return (low + high) / 2.;
  }

  double GetFracSystBiasXSec(const std::string &weightName)
  {
    double low, cv, high;

    std::tie(low, cv, high) = GetSystFracErrorsXSec(weightName);

    const double nom = _nominalBin->GetXSec();
    return (cv - nom) / nom;
  }

  double GetFracSystResAveQuadSumXSec(const std::vector<std::string> &weightNames)
  {
    double quadSum = 0.;
    for(auto const& weightName : weightNames)
      quadSum += TMath::Power(GetFracSystResAveXSec(weightName), 2);

    return TMath::Sqrt(quadSum);
  }

  double GetFracSystBiasQuadSumXSec(const std::vector<std::string> &weightNames)
  {
    double quadSum = 0.;
    for(auto const& weightName : weightNames)
      quadSum += TMath::Power(GetFracSystBiasXSec(weightName), 2);

    return TMath::Sqrt(quadSum);
  }

  double GetFracSystResAveEfficiency(const std::string &weightName)
  {
    double low, cv, high;

    std::tie(low, cv, high) = GetSystFracErrorsEfficiency(weightName);

    const double nom = _nominalBin->GetEfficiency();

    return (low + high) / 2.;
  }

  double GetFracSystBiasEfficiency(const std::string &weightName)
  {
    double low, cv, high;

    std::tie(low, cv, high) = GetSystFracErrorsEfficiency(weightName);

    const double nom = _nominalBin->GetEfficiency();
    return (cv - nom) / nom;
  }

  double GetFracSystResAveQuadSumEfficiency(const std::vector<std::string> &weightNames)
  {
    double quadSum = 0.;
    for(auto const& weightName : weightNames)
      quadSum += TMath::Power(GetFracSystResAveEfficiency(weightName), 2);

    return TMath::Sqrt(quadSum);
  }

  double GetFracSystBiasQuadSumEfficiency(const std::vector<std::string> &weightNames)
  {
    double quadSum = 0.;
    for(auto const& weightName : weightNames)
      quadSum += TMath::Power(GetFracSystBiasEfficiency(weightName), 2);

    return TMath::Sqrt(quadSum);
  }

  double GetFracSystResAvePurity(const std::string &weightName)
  {
    double low, cv, high;

    std::tie(low, cv, high) = GetSystFracErrorsPurity(weightName);

    const double nom = _nominalBin->GetPurity();

    return (low + high) / 2.;
  }

  double GetFracSystBiasPurity(const std::string &weightName)
  {
    double low, cv, high;

    std::tie(low, cv, high) = GetSystFracErrorsPurity(weightName);

    const double nom = _nominalBin->GetPurity();
    return (cv - nom) / nom;
  }

  double GetFracSystResAveQuadSumPurity(const std::vector<std::string> &weightNames)
  {
    double quadSum = 0.;
    for(auto const& weightName : weightNames)
      quadSum += TMath::Power(GetFracSystResAvePurity(weightName), 2);

    return TMath::Sqrt(quadSum);
  }

  double GetFracSystBiasQuadSumPurity(const std::vector<std::string> &weightNames)
  {
    double quadSum = 0.;
    for(auto const& weightName : weightNames)
      quadSum += TMath::Power(GetFracSystBiasPurity(weightName), 2);

    return TMath::Sqrt(quadSum);
  }

  double GetFracSystResAveBkgdCount(const std::string &weightName)
  {
    double low, cv, high;

    std::tie(low, cv, high) = GetSystFracErrorsBkgdCount(weightName);

    const double nom = _nominalBin->GetBkgdCount();

    return (low + high) / 2.;
  }

  double GetFracSystBiasBkgdCount(const std::string &weightName)
  {
    double low, cv, high;

    std::tie(low, cv, high) = GetSystFracErrorsBkgdCount(weightName);

    const double nom = _nominalBin->GetBkgdCount();
    return (cv - nom) / nom;
  }

  double GetFracSystResAveQuadSumBkgdCount(const std::vector<std::string> &weightNames)
  {
    double quadSum = 0.;
    for(auto const& weightName : weightNames)
      quadSum += TMath::Power(GetFracSystResAveBkgdCount(weightName), 2);

    return TMath::Sqrt(quadSum);
  }

  double GetFracSystBiasQuadSumBkgdCount(const std::vector<std::string> &weightNames)
  {
    double quadSum = 0.;
    for(auto const& weightName : weightNames)
      quadSum += TMath::Power(GetFracSystBiasBkgdCount(weightName), 2);

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

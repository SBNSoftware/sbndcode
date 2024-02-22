#pragma once

#include "/exp/sbnd/data/users/hlay/ncpizero/plots/NCPiZeroAv16/integrated_flux/FluxMapTmp.h"
#include "UniverseBin.h"

class Bin {

 private:
  int                                                        _index;
  UniverseBin                                               *_nominalBin;
  std::map<std::string, std::vector<UniverseBin*>>          *_universeBins;
  std::map<std::string, std::tuple<double, double, double>> *_systFracErrors;
  double                                                     _var0Low;
  double                                                     _var0High;
  double                                                     _var1Low;
  double                                                     _var1High;

 public:
  
  Bin(const int index, const double binWidth, const double nTargets, const double intFlux,
      const double var0Low, const double var0High, const double var1Low, const double var1High)
    {
      _index          = index;
      _nominalBin     = new UniverseBin(binWidth, nTargets, intFlux);
      _universeBins   = new std::map<std::string, std::vector<UniverseBin*>>();
      _systFracErrors = new std::map<std::string, std::tuple<double, double, double>>();
      _var0Low        = var0Low;
      _var0High       = var0High;
      _var1Low        = var1Low;
      _var1High       = var1High;
    }

  void Print()
  {
    std::cout << "\n======================================\n"
              << "Bin Index: " << _index << '\n';

    _nominalBin->Print();

    std::cout <<"======================================\n" << std::endl;
  }

  void Print(const std::string weightName, const int univ)
  {
    std::cout << "\n======================================\n"
              << "Bin Index: " << _index << '\n';

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
          bin->CalculateXSecPurity();
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

  double GetNominalFracStatErr()
  {
    return _nominalBin->GetFracStatErr();
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

  void CalculateSystFracErrors()
  {
    for(auto&& [ name, univBins ] : *_universeBins)
      {
        std::vector<double> xsecs;

        for(UniverseBin* bin : univBins)
          xsecs.push_back(bin->GetXSec());

        std::sort(xsecs.begin(), xsecs.end());

        const int n         = univBins.size();
        const double cv_n   = n * (1/2.);
        const double low_n  = n * (1/3.);
        const double high_n = n * (2/3.);

        double cv   = (xsecs[std::floor(cv_n)] + xsecs[std::ceil(cv_n)]) / 2.;
        double low  = (xsecs[std::floor(low_n)] + xsecs[std::ceil(low_n)]) / 2.;
        double high = (xsecs[std::floor(high_n)] + xsecs[std::ceil(high_n)]) / 2.;

        low  = std::abs(low - cv);
        high = std::abs(high - cv);

        _systFracErrors->insert({ name, { low, cv, high }});
      }
  }

  std::tuple<double, double, double> GetSystFracErrors(const std::string &weightName)
  {
    return _systFracErrors->at(weightName);
  }

  bool InBin(const double &val0, const double &val1)
  {
    return val0 > _var0Low && val0 < _var0High && val1 > _var1Low && val1 < _var1High;
  }
};

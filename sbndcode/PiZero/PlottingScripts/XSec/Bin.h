#pragma once

constexpr int def_int       = std::numeric_limits<int>::min();
constexpr double def_double = std::numeric_limits<double>::lowest();

class Bin {

 private:
  int                           _index;
  double                        _xsec;
  double                        _rawCount;
  double                        _scaledCount;
  double                        _bkgdCount;
  double                        _scaledBkgdCount;
  double                        _purity;
  double                        _efficiency;
  double                        _binWidth;
  double                        _fracStatErr;
  std::map<std::string, double> _fracSystErr;

 public:
  
  Bin()
    {
      _index           = def_int;
      _xsec            = def_double;
      _rawCount        = def_double;
      _scaledCount     = def_double;
      _bkgdCount       = def_double;
      _scaledBkgdCount = def_double;
      _purity          = def_double;
      _efficiency      = def_double;
      _binWidth        = 1.;
      _fracStatErr     = def_double;
      _fracSystErr     = std::map<std::string, double>();
    }

 Bin(const int index)
   : Bin()
    {
      _index = index;
    }

  void Print()
  {
    std::cout << "\n======================================\n"
              << "Bin Index: " << _index << '\n'
              << "Count: " << _rawCount << " (" << _scaledCount << ")\n"
              << "Background Count: " << _bkgdCount << " (" << _scaledBkgdCount << ")\n"
              << "Purity: " << _purity << '\n'
              << "Efficiency: " << _efficiency << '\n'
              << "Bin Width: " << _binWidth << '\n'
              << "Frac Stat Err: " << _fracStatErr << '\n'
              << "--------------------------------------\n"
              << "XSec: " << _xsec << '\n'
              <<"======================================\n" << std::endl;
  }

  void SetIndex(const int index)
  {
    _index = index;
  }

  void SetRawCount(const double rawCount)
  {
    _rawCount = rawCount;
  }

  void SetScaledCount(const double scaledCount)
  {
    _scaledCount = scaledCount;
  }

  void SetBackgroundCount(const double bkgdCount)
  {
    _bkgdCount = bkgdCount;
  }

  void SetScaledBackgroundCount(const double scaledBkgdCount)
  {
    _scaledBkgdCount = scaledBkgdCount;
  }

  void SetPurity(const double purity)
  {
    _purity = purity;
  }

  void SetEfficiency(const double efficiency)
  {
    _efficiency = efficiency;
  }

  void SetBinWidth(const double binWidth)
  {
    _binWidth = binWidth;
  }

  void SetFracStatErr(const double fracStatErr)
  {
    _fracStatErr = fracStatErr;
  }

  void InsertFracSystErr(const std::string name, const double error)
  {
    _fracSystErr[name] = error;
  }

  int GetIndex() const
  {
    return _index;
  }

  double GetXSec() const
  {
    return _xsec;
  }

  double GetRawCount() const
  {
    return _rawCount;
  }

  double GetScaledCount() const
  {
    return _scaledCount;
  }

  double GetBackgroundCount() const
  {
    return _bkgdCount;
  }

  double GetScaledBackgroundCount() const
  {
    return _scaledBkgdCount;
  }

  double GetPurity() const
  {
    return _purity;
  }

  double GetEfficiency() const
  {
    return _efficiency;
  }

  double GetBinWidth() const
  {
    return _binWidth;
  }

  double GetFracStatErr() const
  {
    return _fracStatErr;
  }

  double GetFracSystErr(const std::string &name) const
  {
    return _fracSystErr.at(name);
  }

  void CalculateXSecPurity(const double &nTargets, const double &intFlux)
  {
    _xsec = (_scaledCount * _purity) / (_efficiency * nTargets * intFlux * _binWidth);
  }

  void CalculateXSecBackgroundSubtraction(const double &nTargets, const double &intFlux)
  {
    _xsec = (_scaledCount - _scaledBkgdCount) / (_efficiency * nTargets * intFlux * _binWidth);
  }
};

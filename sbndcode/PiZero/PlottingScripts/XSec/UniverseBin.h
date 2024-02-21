#pragma once

constexpr int def_int       = std::numeric_limits<int>::min();
constexpr double def_double = std::numeric_limits<double>::lowest();

class UniverseBin {

 private:
  double _xsec;
  double _count;
  double _bkgdCount;
  double _trueSignal;
  double _scaleFactor;
  double _efficiency;
  double _purity;
  double _binWidth;
  double _nTargets;
  double _intFlux;
  double _fracStatErr;

 public:
  
  UniverseBin(const double binWidth, const double nTargets, const double intFlux)
    {
      _xsec        = def_double;
      _count       = 0.;
      _bkgdCount   = 0.;
      _trueSignal  = 0.;
      _scaleFactor = def_double;
      _purity      = def_double;
      _efficiency  = def_double;
      _fracStatErr = def_double;
      _binWidth    = binWidth;
      _nTargets    = nTargets;
      _intFlux     = intFlux;
    }

  void Print()
  {
    std::cout << "Count: " << _count << " (" << _count * _scaleFactor << ")\n"
              << "Background Count: " << _bkgdCount << " (" << _bkgdCount * _scaleFactor << ")\n"
              << "Purity: " << _purity << '\n'
              << "Efficiency: " << _efficiency << '\n'
              << "Bin Width: " << _binWidth << '\n'
              << "N Targets: " << _nTargets << '\n'
              << "Int Flux: " << _intFlux << '\n'
              << "--------------------------------------\n"
              << "XSec: " << _xsec << '\n'
              << "--------------------------------------\n";
  }

  void SetScaleFactor(const double &scaleFactor)
  {
    _scaleFactor = scaleFactor;
  }

  void Update()
  {
    _purity      = (_count - _bkgdCount) / _count;
    _efficiency  = (_count - _bkgdCount) / _trueSignal;
    _fracStatErr = std::sqrt(_count) / _count;
  }

  void CalculateXSecPurity()
  {
    _xsec = (_count * _scaleFactor * _purity) / (_efficiency * _nTargets * _intFlux * _binWidth);
  }

  void CalculateXSecBackgroundSubtraction()
  {
    _xsec = ((_count - _bkgdCount) * _scaleFactor) / (_efficiency * _nTargets * _intFlux * _binWidth);
  }

  double GetXSec()
  {
    return _xsec;
  }

  double GetFracStatErr()
  {
    return _fracStatErr;
  }

  double GetBinWidth()
  {
    return _binWidth;
  }

  double GetNTargets()
  {
    return _nTargets;
  }

  double GetIntFlux()
  {
    return _intFlux;
  }

  void IncrementCount(const double &increment)
  {
    _count += increment;
  }

  void IncrementBkgdCount(const double &increment)
  {
    _bkgdCount += increment;
  }

  void IncrementTrueSignal(const double &increment)
  {
    _trueSignal += increment;
  }
};

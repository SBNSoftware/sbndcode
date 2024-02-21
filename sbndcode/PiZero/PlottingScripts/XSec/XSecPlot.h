#pragma once

#include "Bin.h"

class XSecPlot {

 private:
  std::string      _name;
  std::string      _axes_labels;
  std::vector<Bin> _bins;

 public:

  XSecPlot(const std::string name, const std::string axes_labels, const int n)
    {
      _name        = name;
      _axes_labels = axes_labels;

      for(int i = 0; i < n; ++i)
        _bins.emplace_back(i);
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

  Bin& GetBin(const int i)
    {
      if(i >= _bins.size())
        throw std::runtime_error("Asking for non-existent bin");

      return _bins[i];
    }
};

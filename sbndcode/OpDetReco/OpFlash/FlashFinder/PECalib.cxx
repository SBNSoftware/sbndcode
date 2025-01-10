#ifndef PECALIB_CXX
#define PECALIB_CXX

#include "PECalib.h"

namespace lightana {

  PECalib::PECalib()
  {}

  void PECalib::Configure(const Config_t &pset)
  {

    _spe_area_gain_v.clear();
    _spe_area_gain_v = pset.get<std::vector<double> >("SPEAreaGainList",_spe_area_gain_v);

    if(_spe_area_gain_v.empty()) {
      double spe_area_gain = pset.get<double>("SPEAreaGain");
      _spe_area_gain_v.resize(NOpDets(),spe_area_gain);
    }

    if(_spe_area_gain_v.size() != NOpDets()) {
      std::cerr << "SPEAreaGain array size (" << _spe_area_gain_v.size()
                << ") != NOpDets (" << NOpDets() << ")..." << std::endl;
      throw std::exception();
    }

    _relative_qe_v.clear();
    _relative_qe_v = pset.get<std::vector<double> >("RelativeQEList",_relative_qe_v);

    if(_relative_qe_v.empty())
      _relative_qe_v.resize(NOpDets(),1.0);

    if(_relative_qe_v.size() != NOpDets()) {
      std::cerr << "RelativeQE array size (" << _relative_qe_v.size()
                << ") != NOpDets (" << NOpDets() << ")..." << std::endl;
      throw std::exception();
    }
  }

  double PECalib::Calibrate(const size_t opdet, const double area) const
  {
    if( opdet > NOpDets() ) {
      std::cerr << "OpDet ID " << opdet << " exceeding max # of OpDet (" << NOpDets() << ")" << std::endl;
      throw std::exception();
    }

    double area_pe = area / _spe_area_gain_v[opdet] * _relative_qe_v[opdet];

    return area_pe;
  }

}
#endif


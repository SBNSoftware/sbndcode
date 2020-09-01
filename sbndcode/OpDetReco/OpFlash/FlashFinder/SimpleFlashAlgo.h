#ifndef SIMPLEFLASHALGO_H
#define SIMPLEFLASHALGO_H

#include "FlashAlgoBase.h"
#include "FlashAlgoFactory.h"
#include <map>

namespace lightana
{

  class SimpleFlashAlgo : public FlashAlgoBase {

  public:

    SimpleFlashAlgo(const std::string name);

    void Configure(const Config_t &p);
    
    virtual ~SimpleFlashAlgo();

    LiteOpFlashArray_t RecoFlash(const LiteOpHitArray_t ophits);

    bool Veto(double t) const;

    const std::vector<double>& PESumArray() const { return _pesum_v; }

    const double TimeRes() const { return _time_res; }

  private:

    double TotalCharge(const std::vector<double>& PEs);
    double _min_pe_flash;   // minimum PE to make a flash
    double _min_pe_coinc;   // minimum PE to make a flash candidate
    double _min_mult_coinc; // minimum PE to use an OpHit
    double _integral_time;  // integral period
    double _veto_time;      // veto time
    double _time_res;       // time resolution of pe sum
    double _pre_sample;     // time pre-sample

    std::vector<double> _pesum_v;        // pw aum array
    std::vector<double> _pe_baseline_v;  // calibration: PEs to be subtracted from each opdet

    std::map<double,double> _flash_veto_range_m;  // veto window start

    bool _debug;            // debug mode flag

    // list of opchannel to use
    std::vector<int> _opch_to_index_v;
    std::vector<int> _index_to_opch_v;
                             
  };

  /**
    \class lightana::SimpleFlashAlgoFactory
    \brief A concrete factory class for lightana::SimpleFlashAlgo
  */
  class SimpleFlashAlgoFactory : public FlashAlgoFactoryBase {

  public:
    /// ctor
    SimpleFlashAlgoFactory() { FlashAlgoFactory::get().add_factory("SimpleFlashAlgo",this); }
    /// dtor
    ~SimpleFlashAlgoFactory() {}
    /// creation method
    FlashAlgoBase* create(const std::string instance_name) { return new SimpleFlashAlgo(instance_name); }
  };

}
#endif

/** @} */ // end of doxygen group

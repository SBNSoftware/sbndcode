////////////////////////////////////////////////////////////////////////
// File:        sbndPDMapAlg.hh
// Authors: Laura Paulucci, Franciole Marinho, and Iker de Icaza
//
// Updates: 2020-03, v08_45_00. Iker de Icaza icaza@fnal.gov
//          Added properties to PDS MAP and functions to access these.
//
// This class stores the SBND PDS Map and channel's properties;
// also implements functions to access these.
//
// As of version v09_08_00 each entry of the PDS Map
// has the following characteristics:
//
// channel: 0 to 303
// pd_type: pmt_coated, pmt_uncoated, xarapuca_vuv, xarapuca_vis, arapuca_vuv, arapuca_vis
// pds_box: 0 to 23
// sensible_to_vis: true or false
// sensible_to_vuv: true or false
// tpc: 0, 1
////////////////////////////////////////////////////////////////////////

#ifndef SBND_OPDETSIM_SBNDPDMAPALG_HH
#define SBND_OPDETSIM_SBNDPDMAPALG_HH

#include "sbncode/OpDet/PDMapAlg.h"
//#include "art/Utilities/ToolMacros.h"
//#include "art/Utilities/make_tool.h"

#include <algorithm>
#include <fstream>
#include <map>
#include <string>

#include "art_root_io/TFileService.h"

#include "json.hpp"

namespace opdet {

  class sbndPDMapAlg : PDMapAlg{

  public:
    //Default constructor
    explicit sbndPDMapAlg(const fhicl::ParameterSet& pset);
    sbndPDMapAlg() : sbndPDMapAlg(fhicl::ParameterSet()) {}
    //Default destructor
    ~sbndPDMapAlg();

    nlohmann::json getCollectionWithProperty(std::string property);
    template<typename T> nlohmann::json getCollectionWithProperty(std::string property, T property_value);
    template<typename Func> nlohmann::json getCollectionFromCondition(Func condition);

    // struct Config {};

    //  sbndPDMapAlg(Config const&) {}

    // void setup() {}

    bool isPDType(size_t ch, std::string pdname) const override;
    std::string pdType(size_t ch) const override;
    std::vector<int> getChannelsOfType(std::string pdname) const;
    size_t size() const;
    auto getChannelEntry(size_t ch) const;

  private:
    nlohmann::json PDmap;

  }; // class sbndPDMapAlg

  template<typename T>
  nlohmann::json sbndPDMapAlg::getCollectionWithProperty(std::string property, T property_value)
  {
    nlohmann::json subSetPDmap;
    std::copy_if (PDmap.begin(), PDmap.end(), std::back_inserter(subSetPDmap),
                  [property, property_value](auto const& e)->bool
                    {return e[property] == property_value;} );
    return subSetPDmap;
  }

  template<typename Func>
  nlohmann::json sbndPDMapAlg::getCollectionFromCondition(Func condition)
  {
    nlohmann::json subSetPDmap;
    std::copy_if (PDmap.begin(), PDmap.end(), std::back_inserter(subSetPDmap),
                  condition);
    return subSetPDmap;
  }

} // namespace

#endif // SBND_OPDETSIM_SBNDPDMAPALG_HH

////////////////////////////////////////////////////////////////////////
// File:        sbndPDMapAlg.h
// Authors: Laura Paulucci, Franciole Marinho, and Iker de Icaza
//
// Updates: 2020-03, v08_45_00. Iker de Icaza icaza@fnal.gov
//          Added properties to PDS MAP and functions to access these.
//
// This class stores the SBND PDS Map and channel's properties;
// also implements functions to access these.
//
// As of version v08_45_00 the PDS Map has:
// channel: 0 to 503
// pd_type: pmt_coated, pmt_uncoated, xarapuca, xarapuca_vuv, xarapuca_vis, arapuca_vuv, arapuca_vis
// pds_box: -12, to 12; skipping 0
// sensible_to: VUV or VIS
// tpc: 0, 1
// xarapuca_pos: top, bottom, null
////////////////////////////////////////////////////////////////////////

#ifndef SBND_OPDETSIM_SBNDPDMAPALG_H
#define SBND_OPDETSIM_SBNDPDMAPALG_H

#include <algorithm>
#include <fstream>
#include <map>
#include <string>

#include "art_root_io/TFileService.h"

#include "json.hpp"

namespace opdet {

  class sbndPDMapAlg {

  public:
    //Default constructor
    sbndPDMapAlg();
    //Default destructor
    ~sbndPDMapAlg();

    template<typename T> nlohmann::json getCollectionWithProperty(std::string property, T property_value);

    // struct Config {};

    //  sbndPDMapAlg(Config const&) {}

    // void setup() {}

    bool isPDType(size_t ch, std::string pdname) const;
    std::string pdType(size_t ch) const;
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
                  [property, property_value](const nlohmann::json e)->bool
                    {return e[property] == property_value;} );
    return subSetPDmap;
  }

} // namespace

#endif // SBND_OPDETSIM_SBNDPDMAPALG_H

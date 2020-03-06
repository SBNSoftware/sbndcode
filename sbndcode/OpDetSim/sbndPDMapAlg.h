////////////////////////////////////////////////////////////////////////
// File:        sbndPDMapAlg.h
// Authors: Laura Paulucci and Franciole Marinho
//
// This class stores the SBND PDS type map and implements a few functions
//
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

    nlohmann::json getCollectionWithProperty(std::string property, std::string property_value);
    nlohmann::json getCollectionWithProperty(std::string property, int property_value);
    // template<typename T> nlohmann::json getCollectionWithProperty(std::string property, T property_value);

    // struct Config {};

    //  sbndPDMapAlg(Config const&) {}

    // void setup() {}

    bool isPDType(size_t ch, std::string pdname) const;
    std::string pdType(size_t ch) const;
    size_t size() const;

  private:
    nlohmann::json PDmap;
    nlohmann::json subSetPDmap;

  }; // class sbndPDMapAlg

} // namespace

#endif // SBND_OPDETSIM_SBNDPDMAPALG_H

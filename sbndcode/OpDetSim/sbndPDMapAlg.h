////////////////////////////////////////////////////////////////////////
// File:        sbndPDMapAlg.h
// Authors: Laura Paulucci and Franciole Marinho
//
// This class stores the SBND PDS type map and implements a few functions
//
////////////////////////////////////////////////////////////////////////

#ifndef SBND_OPDETSIM_SBNDPDMAPALG_H
#define SBND_OPDETSIM_SBNDPDMAPALG_H

// LArSoft libraries

// framework libraries
#include <string>
#include <map>

namespace opdet {

  class sbndPDMapAlg {

  public:
    //Default constructor
    sbndPDMapAlg();
    //Default destructor
    ~sbndPDMapAlg();

    // struct Config {};

    //  sbndPDMapAlg(Config const&) {}

    // void setup() {}

    bool pdType(int ch, std::string pdname) const;
    std::string pdName(int ch) const;
    int size() const;

  private:
    std::map<int, std::string> PDmap;

  }; // class sbndPDMapAlg

} // namespace

#endif // SBND_OPDETSIM_SBNDPDMAPALG_H

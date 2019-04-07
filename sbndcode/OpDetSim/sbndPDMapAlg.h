////////////////////////////////////////////////////////////////////////
// File:        sbndPDMapAlg.h
// Authors: Laura Paulucci and Franciole Marinho
//
// This class stores the SBND PDS type map and implements a few functions
// 
////////////////////////////////////////////////////////////////////////

#ifndef SBNDPDMAPALG_H
#define SBNDPDMAPALG_H

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

  bool pdType(int ch, std::string pdname);
  std::string pdName(int ch);
  int size();
        
  private:	 
  std::map<int, std::string> PDmap;
          
  }; // class sbndPDMapAlg
             
} // namespace 
    
#endif // SBNDPDMAPALG_H

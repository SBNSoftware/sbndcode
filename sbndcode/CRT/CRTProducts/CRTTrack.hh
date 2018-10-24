/**
 * \class CRTData
 *
 * \ingroup crt
 *
 * \brief CRT Track Info
 *
 * \author $Author: David Lorca $
 *
 */

#ifndef CRTTrack_hh_
#define CRTTrack_hh_
#include <cstdint>

#include <vector>
#include <map>

namespace sbnd {
namespace crt {
  
  struct CRTTrack{

    std::vector<uint8_t> feb_id;
    std::map< uint8_t, std::vector<std::pair<int,float> > > pesmap;
    float peshit;
    uint32_t ts0_s;
    uint16_t ts0_s_err;
    uint32_t ts0_ns;
    uint16_t ts0_ns_err;
    int32_t ts1_ns; 
    uint16_t ts1_ns_err;                                                                                              
    int plane1;
    int plane2;
                           
    float x1_pos;
    float x1_err;
    float y1_pos;
    float y1_err;
    float z1_pos;
    float z1_err;
    float x2_pos;
    float x2_err;
    float y2_pos;
    float y2_err;
    float z2_pos;
    float z2_err;
    float length;
    float thetaxy;
    float phizy;
    uint32_t ts0_ns_h1;
    uint16_t ts0_ns_err_h1;
    uint32_t ts0_ns_h2;
    uint16_t ts0_ns_err_h2;

    bool complete;
       
    CRTTrack() {}
    
  };
  
}
}
#endif

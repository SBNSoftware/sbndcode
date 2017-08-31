#include "sbndcode/CRT/CRTData.hh"

namespace sbnd{
namespace crt{

  CRTData::CRTData(): fChannel(0), fT0(0), fT1(0){
  }
  CRTData::CRTData(uint32_t channel, uint32_t t0, 
    uint32_t t1, uint32_t adc):
    fChannel(channel),
    fT0(t0),
    fT1(t1),
    fADC(adc) {
    }
  CRTData::~CRTData(){
  }
  uint32_t CRTData::Channel() const {
    return fChannel;
  }
  uint32_t CRTData::T0() const {
    return fT0;
  }
  uint32_t CRTData::T1() const {
    return fT1;
  }
  uint32_t CRTData::ADC() const {
    return fADC;
  }
} // namespace crt
} // namespace sbnd

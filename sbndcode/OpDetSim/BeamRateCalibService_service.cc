/**
 * @file    BeamRateCalibService.cc
 * @brief   Service implementation associated with BeamRateCalibService.h
 * @author  Nikki Pallat (palla110@umn.edu)
 * @date    July 7, 2025
 * @see     BeamRateCalibService.h
 * @ingroup BeamRateCalibService
 *
 */

#include "BeamRateCalibService.h"

//namespace lar {
namespace calib {

  BeamRateCalibService::BeamRateCalibService(fhicl::ParameterSet const& pset, art::ActivityRegistry& amp)
    : fMonWidth(pset.get<int>("MonWidth", 3)),
      fTotalCAENBoards(pset.get<int>("TotalCAENBoards", 10)),
      PMTPerBoard(pset.get<int>("PMTPerBoard", 16)),
      Baseline(pset.get<int>("Baseline", 14250)),
      fMC(pset.get<bool>("MC", false))
  {}
  BeamRateCalibService::~BeamRateCalibService()
  {}

  void BeamRateCalibService::ConstructMonPulse(
    art::Handle< std::vector< raw::OpDetWaveform > > &waveHandle, 
    int MonThreshold, 
    std::vector<int> *MonPulse, 
    bool Saving, 
    int FlashCounter
  ) 
  {
    MonPulse->resize(20);
    for(int i=0; i<20; i++) (*MonPulse)[i]=12;

  } // ConstructMonPulse


  std::vector<bool> BeamRateCalibService::ConstructBinaryResponse(const raw::OpDetWaveform &wvf, int MonThreshold) const
  {
    std::vector<bool> BinaryResponse(wvf.size());

    return BinaryResponse; 
  } // ConstructBinaryResponse

} // namespace calib



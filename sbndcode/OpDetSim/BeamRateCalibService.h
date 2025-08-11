/**
 * @file    BeamRateCalibService.h
 * @brief   Provider modifying given waveform
 * @author  Nikki Pallat (palla110@umn.edu)
 * @date    July 7, 2025
 * @see     BeamRateCalibService.cxx BeamRateCalibService.h
 * @ingroup BeamRateCalibService
 *
 */

#ifndef BEAMRATECALIB_SERVICE_H
#define BEAMRATECALIB_SERVICE_H


// potential support libraries 
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Principal/Handle.h"
//#include "canvas/Persistency/Common/Handle.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "larcore/CoreUtils/ServiceUtil.h"

// C/C++ standard libraries
#include <string>
#include <vector>

//namespace lar {
   namespace calib {

       class BeamRateCalibService { 
         public:
           BeamRateCalibService(fhicl::ParameterSet const&pset, art::ActivityRegistry&amp);
           ~BeamRateCalibService();
           void ConstructMonPulse(
             art::Handle< std::vector< raw::OpDetWaveform > > &waveHandle, 
             int MonThreshold, 
             std::vector<int> *MonPulse, 
             bool Saving, 
             int FlashCounter
           ) const;

         //---------------------------------------------------------------------
         private:
           std::vector<bool> ConstructBinaryResponse(const raw::OpDetWaveform &wvf, int MonThreshold) const;

         int fMonWidth;
         int fTotalCAENBoards;
         int PMTPerBoard;
         int Baseline;
         bool fMC;
       }; // BeamRateCalibService

   } // namespace calib

DECLARE_ART_SERVICE(calib::BeamRateCalibService, LEGACY)

#endif // BEAMRATECALIB_SERVICE_H



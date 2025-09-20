/**
 * @file    TriggerEmulationService.h
 * @brief   Provider modifying given waveform
 * @author  Nikki Pallat (palla110@umn.edu)
 * @date    July 7, 2025
 * @see     TriggerEmulationService.cc TriggerEmulationService_service.cc
 * @ingroup TriggerEmulationService
 *
 */

#ifndef TRIGGEREMULATION_SERVICE_H
#define TRIGGEREMULATION_SERVICE_H


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

       class TriggerEmulationService { 
         public:
           TriggerEmulationService(fhicl::ParameterSet const&pset, art::ActivityRegistry&amp);
           ~TriggerEmulationService();
           void ConstructMonPulse(
             std::vector<raw::OpDetWaveform> fWaveforms,
             int MonThreshold, 
             std::vector<int> *MonPulse, 
             bool Saving, 
             int FlashCounter,
             int *numPairsOverThreshold = nullptr
           );

           int getTotalCAENBoards() const { return fTotalCAENBoards; } 
           int getPMTPerBoard() const { return PMTPerBoard; } 

         //---------------------------------------------------------------------
         private:
           std::vector<bool> ConstructBinaryResponse(const raw::OpDetWaveform &wvf, int MonThreshold);

         int fMonWidth;
         int fTotalCAENBoards;
         int PMTPerBoard;
         int Baseline;
         bool fMC;
       }; // TriggerEmulationService

   } // namespace calib

DECLARE_ART_SERVICE(calib::TriggerEmulationService, LEGACY)
#endif // TRIGGEREMULATION_SERVICE_H



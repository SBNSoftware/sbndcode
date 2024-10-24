//////////////////////////////////////////////////////////////////////////////////////////////////
// Class:       SBND::CRTCalibService
// Module type: service
// File:        SBND::CRTCalibService_service.cc
// Author:      Henry Lay, May 2024.
//
// Implementation of service for inputing CRT calibration values to reconstruction
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>

#include "CRTCalibService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
  
SBND::CRTCalibService::CRTCalibService(fhicl::ParameterSet const& pset)
{
  const std::string timingOffsetFile = pset.get<std::string>("TimingOffsetFileName");
  const std::string pedestalFile     = pset.get<std::string>("PedestalFileName");

  std::string timingOffsetPath, pedestalPath;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(timingOffsetFile, timingOffsetPath);
  sp.find_file(pedestalFile, pedestalPath);

  if(timingOffsetPath.empty())
    {
      std::cout << "SBND::CRTCalibService Input file " << timingOffsetFile << " not found" << std::endl;
      throw cet::exception("File not found");
    }

  if(pedestalPath.empty())
    {
      std::cout << "SBND::CRTCalibService Input file " << pedestalFile << " not found" << std::endl;
      throw cet::exception("File not found");
    }

  std::cout << "SBND CRT Channel Map: Building map from files...\n"
            << "\t Timing Offsets: " << timingOffsetFile << '\n'
            << "\t Pedestals: " << pedestalFile << std::endl;

  std::ifstream timingOffsetStream(timingOffsetPath, std::ios::in);
  std::string line;

  while(std::getline(timingOffsetStream, line))
    {
      std::stringstream linestream(line);
    
      unsigned int mac5, offset;
      linestream 
        >> mac5
        >> offset;

      fTimingOffsetFromFEBMAC5[mac5] = offset;
    }
  
  timingOffsetStream.close();

  std::ifstream pedestalStream(pedestalPath, std::ios::in);

  while(std::getline(pedestalStream, line))
    {
      std::stringstream linestream(line);
    
      unsigned int mac5, ch, pedestal;
      linestream 
        >> mac5
        >> ch
        >> pedestal;

      fPedestalFromFEBMAC5AndChannel[mac5][ch] = pedestal;
    }
  
  pedestalStream.close();
}

SBND::CRTCalibService::CRTCalibService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
  : SBND::CRTCalibService(pset)
{
}

double SBND::CRTCalibService::GetTimingOffsetFromFEBMAC5(unsigned int feb_mac5) const
{
  auto iter = fTimingOffsetFromFEBMAC5.find(feb_mac5);

  if(iter == fTimingOffsetFromFEBMAC5.end())
    {
      mf::LogInfo("SBND CRT Calibration Service") << "Asked for FEB with MAC5: " << feb_mac5 << '\n'
                                                  << "This FEB does not appear in the channel map."
                                                  << std::endl;

      return 0.;
    }

  return iter->second;
}

double SBND::CRTCalibService::GetPedestalFromFEBMAC5AndChannel(unsigned int feb_mac5, unsigned int ch) const
{
  auto iter = fPedestalFromFEBMAC5AndChannel.find(feb_mac5);

  if(iter == fPedestalFromFEBMAC5AndChannel.end())
    {
      mf::LogInfo("SBND CRT Calibration Service") << "Asked for FEB with MAC5: " << feb_mac5 << '\n'
                                                  << "This FEB does not appear in the calibration service."
                                                  << std::endl;

      return 0.;
    }

  auto subIter = iter->second.find(ch);

  if(subIter == iter->second.end())
    {
      mf::LogInfo("SBND CRT Calibration Service") << "Asked for channel: " << ch << '\n'
                                                  << "in FEB with MAC5: " << feb_mac5 << '\n'
                                                  << "This channel does not appear in the calibration service."
                                                  << std::endl;

      return 0.;
    }

  return subIter->second;
}

DEFINE_ART_SERVICE(SBND::CRTCalibService)

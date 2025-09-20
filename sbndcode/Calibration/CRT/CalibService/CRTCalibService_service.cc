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
  const std::string t0CableOffsetFile = pset.get<std::string>("T0CableOffsetFileName");
  const std::string t1CableOffsetFile = pset.get<std::string>("T1CableOffsetFileName");
  const std::string t0CalibOffsetFile = pset.get<std::string>("T0CalibOffsetFileName");
  const std::string t1CalibOffsetFile = pset.get<std::string>("T1CalibOffsetFileName");
  const std::string pedestalFile      = pset.get<std::string>("PedestalFileName");
  const std::string badChannelsFile   = pset.get<std::string>("BadChannelsFileName");

  std::string t0CableOffsetPath, t1CableOffsetPath, t0CalibOffsetPath, t1CalibOffsetPath,
    pedestalPath, badChannelsPath;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(t0CableOffsetFile, t0CableOffsetPath);
  sp.find_file(t1CableOffsetFile, t1CableOffsetPath);
  sp.find_file(t0CalibOffsetFile, t0CalibOffsetPath);
  sp.find_file(t1CalibOffsetFile, t1CalibOffsetPath);
  sp.find_file(pedestalFile, pedestalPath);
  sp.find_file(badChannelsFile, badChannelsPath);

  if(t0CableOffsetPath.empty())
    {
      std::cout << "SBND::CRTCalibService Input file " << t0CableOffsetFile << " not found" << std::endl;
      throw cet::exception("File not found");
    }

  if(t1CableOffsetPath.empty())
    {
      std::cout << "SBND::CRTCalibService Input file " << t1CableOffsetFile << " not found" << std::endl;
      throw cet::exception("File not found");
    }

  if(t0CalibOffsetPath.empty())
    {
      std::cout << "SBND::CRTCalibService Input file " << t0CalibOffsetFile << " not found" << std::endl;
      throw cet::exception("File not found");
    }

  if(t1CalibOffsetPath.empty())
    {
      std::cout << "SBND::CRTCalibService Input file " << t1CalibOffsetFile << " not found" << std::endl;
      throw cet::exception("File not found");
    }

  if(pedestalPath.empty())
    {
      std::cout << "SBND::CRTCalibService Input file " << pedestalFile << " not found" << std::endl;
      throw cet::exception("File not found");
    }

  if(badChannelsPath.empty())
    {
      std::cout << "SBND::CRTCalibService Input file " << badChannelsFile << " not found" << std::endl;
      throw cet::exception("File not found");
    }

  std::cout << "SBND CRT Channel Map: Building map from files...\n"
            << "\t T0 Offsets Cable: " << t0CableOffsetFile << '\n'
            << "\t T1 Offsets Cable: " << t1CableOffsetFile << '\n'
            << "\t T0 Offsets Calib: " << t0CalibOffsetFile << '\n'
            << "\t T1 Offsets Calib: " << t1CalibOffsetFile << '\n'
            << "\t Pedestals: " << pedestalFile << '\n'
            << "\t Bad Channels: " << badChannelsFile << std::endl;

  std::ifstream t0CableOffsetStream(t0CableOffsetPath, std::ios::in);
  std::string line;

  while(std::getline(t0CableOffsetStream, line))
    {
      std::stringstream linestream(line);

      unsigned int mac5;
      double offset;
      linestream
        >> mac5
        >> offset;

      fT0CableOffsetFromFEBMAC5[mac5] = offset;
    }

  t0CableOffsetStream.close();

  std::ifstream t1CableOffsetStream(t1CableOffsetPath, std::ios::in);

  while(std::getline(t1CableOffsetStream, line))
    {
      std::stringstream linestream(line);

      unsigned int mac5;
      double offset;
      linestream
        >> mac5
        >> offset;

      fT1CableOffsetFromFEBMAC5[mac5] = offset;
    }

  t1CableOffsetStream.close();

  std::ifstream t0CalibOffsetStream(t0CalibOffsetPath, std::ios::in);

  while(std::getline(t0CalibOffsetStream, line))
    {
      std::stringstream linestream(line);

      unsigned int mac5;
      double offset;
      linestream
        >> mac5
        >> offset;

      fT0CalibOffsetFromFEBMAC5[mac5] = offset;
    }

  t0CalibOffsetStream.close();

  std::ifstream t1CalibOffsetStream(t1CalibOffsetPath, std::ios::in);

  while(std::getline(t1CalibOffsetStream, line))
    {
      std::stringstream linestream(line);

      unsigned int mac5;
      double offset;
      linestream
        >> mac5
        >> offset;

      fT1CalibOffsetFromFEBMAC5[mac5] = offset;
    }

  t1CalibOffsetStream.close();

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

  std::ifstream badChannelsStream(badChannelsPath, std::ios::in);

  while(std::getline(badChannelsStream, line))
    {
      std::stringstream linestream(line);

      unsigned int mac5, ch;
      std::string status_string;

      linestream
        >> mac5
        >> ch
        >> status_string;

      sbnd::crt::CRTChannelStatus status      = sbnd::crt::CRTChannelStatus::kGoodChannel;
      sbnd::crt::CRTChannelStatus status_pair = sbnd::crt::CRTChannelStatus::kGoodChannel;

      if(status_string.compare("kDeadChannel") == 0)
        {
          status      = sbnd::crt::CRTChannelStatus::kDeadChannel;
          status_pair = sbnd::crt::CRTChannelStatus::kDeadNeighbourChannel;
        }
      else if(status_string.compare("kQuietChannel") == 0)
        {
          status      = sbnd::crt::CRTChannelStatus::kQuietChannel;
          status_pair = sbnd::crt::CRTChannelStatus::kQuietNeighbourChannel;
        }
      else
        {
          std::cout << "SBND::CRTCalibService unknown channel status " << status_string << std::endl;
          throw cet::exception("Unknown status");
        }

      fChannelStatusFromFEBMAC5AndChannel[mac5][ch] = status;

      unsigned int ch_pair = ch % 2 ? ch - 1 : ch + 1;
      fChannelStatusFromFEBMAC5AndChannel[mac5][ch_pair] = status_pair;
    }

  badChannelsStream.close();
}

SBND::CRTCalibService::CRTCalibService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
  : SBND::CRTCalibService(pset)
{
}

double SBND::CRTCalibService::GetT0CableOffsetFromFEBMAC5(unsigned int feb_mac5) const
{
  auto iter = fT0CableOffsetFromFEBMAC5.find(feb_mac5);

  if(iter == fT0CableOffsetFromFEBMAC5.end())
    {
      mf::LogInfo("SBND CRT Calibration Service") << "Asked for FEB with MAC5: " << feb_mac5 << '\n'
                                                  << "This FEB does not appear in the channel map."
                                                  << std::endl;

      return 0.;
    }

  return iter->second;
}

double SBND::CRTCalibService::GetT1CableOffsetFromFEBMAC5(unsigned int feb_mac5) const
{
  auto iter = fT1CableOffsetFromFEBMAC5.find(feb_mac5);

  if(iter == fT1CableOffsetFromFEBMAC5.end())
    {
      mf::LogInfo("SBND CRT Calibration Service") << "Asked for FEB with MAC5: " << feb_mac5 << '\n'
                                                  << "This FEB does not appear in the channel map."
                                                  << std::endl;

      return 0.;
    }

  return iter->second;
}

double SBND::CRTCalibService::GetT0CalibOffsetFromFEBMAC5(unsigned int feb_mac5) const
{
  auto iter = fT0CalibOffsetFromFEBMAC5.find(feb_mac5);

  if(iter == fT0CalibOffsetFromFEBMAC5.end())
    {
      mf::LogInfo("SBND CRT Calibration Service") << "Asked for FEB with MAC5: " << feb_mac5 << '\n'
                                                  << "This FEB does not appear in the channel map."
                                                  << std::endl;

      return 0.;
    }

  return iter->second;
}

double SBND::CRTCalibService::GetT1CalibOffsetFromFEBMAC5(unsigned int feb_mac5) const
{
  auto iter = fT1CalibOffsetFromFEBMAC5.find(feb_mac5);

  if(iter == fT1CalibOffsetFromFEBMAC5.end())
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

enum sbnd::crt::CRTChannelStatus SBND::CRTCalibService::GetChannelStatusFromFEBMAC5AndChannel(unsigned int feb_mac5,
                                                                                              unsigned int ch) const
{
  auto iter = fChannelStatusFromFEBMAC5AndChannel.find(feb_mac5);

  if(iter == fChannelStatusFromFEBMAC5AndChannel.end())
    return sbnd::crt::CRTChannelStatus::kGoodChannel;

  auto subIter = iter->second.find(ch);

  if(subIter == iter->second.end())
    return sbnd::crt::CRTChannelStatus::kGoodChannel;

  return subIter->second;
}

DEFINE_ART_SERVICE(SBND::CRTCalibService)

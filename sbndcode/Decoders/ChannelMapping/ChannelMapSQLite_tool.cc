/*
 *  @file   ChannelMapSQLite_tool.cc
 *
 *  @brief  This tool converts from daq to LArSoft format with noise filtering
 *
 */

// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "wda.h"

// LArSoft includes
#include "sbndcode/Decoders/ChannelMapping/IChannelMapping.h"
#include "larevt/CalibrationDBI/Providers/DBFolder.h"

// std includes
#include <string>
#include <iostream>
#include <memory>
#include <sqlite3.h> 
#include <time.h>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace sbndDB {
/**
 *  @brief  ChannelMapSQLite class definiton
 */
class ChannelMapSQLite : virtual public IChannelMapping
{
public:
  /**
   *  @brief  Constructor
   *
   *  @param  pset
   */
  explicit ChannelMapSQLite(fhicl::ParameterSet const &pset);
  
  /**
   *  @brief  Destructor
   */
  ~ChannelMapSQLite();
  
  /**
   *  @brief Define the returned data structures for a mapping between TPC Fragment IDs
   *         and the related crate and readout information. 
   *         Then define the function interface to fill these data structures 
   */
  virtual int BuildTPCFragmentIDToReadoutIDMap(TPCFragmentIDToReadoutIDMap&) const override;
  
  /**
   *  @brief Define the returned data structures for a mapping between TPC readout boards
   *         and the channel information 
   *         Then define the function interface to fill these data structures 
   */
  virtual int BuildTPCReadoutBoardToChannelMap(TPCReadoutBoardToChannelMap&) const override;
  
  /**
   *  @brief Define the returned data structures for a mapping between PMT Fragment IDs
   *         and the related crate and readout information. 
   *         Then define the function interface to fill these data structures 
   */
  virtual int BuildFragmentToDigitizerChannelMap(FragmentToDigitizerChannelMap&) const override;
  
  
  /**
   *  @brief Define the returned data structures for a mapping between CRT hardware mac_address
   *         to the simulated mac_address.
   *         Then define the function interface to fill these data structures
   */
  virtual int BuildCRTChannelIDToHWtoSimMacAddressPairMap(CRTChannelIDToHWtoSimMacAddressPairMap&) const override;
  virtual int BuildTopCRTHWtoSimMacAddressPairMap(TopCRTHWtoSimMacAddressPairMap&) const override;

  virtual int BuildSideCRTCalibrationMap(SideCRTChannelToCalibrationMap&) const override;  
  
private:
  // Recover data from postgres database
  int GetDataset(const std::string&, int func(void*,int,char**,char**), void*) const;
  

  std::string fDBFileName;          //< File name of our sqlite database
  std::string fCalibDBFileName;     //< File name of our side crt calibration sqlite database
  std::string fTag;                 //< Tag for conditioned database
};

ChannelMapSQLite::ChannelMapSQLite(fhicl::ParameterSet const &pset)
{
    fDBFileName      = pset.get<std::string>("DBFileName");
    fCalibDBFileName = pset.get<std::string>("CalibDBFileName");
    fTag             = pset.get<std::string>("Tag");

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

ChannelMapSQLite::~ChannelMapSQLite()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

int callback(void *data, int argc, char **argv, char **azColName)
{
   int i;
   
   for(i = 0; i<argc; i++){
      std::cout << "column: " << azColName[i] << "- value: " << (argv[i] ? argv[i] : "NULL") << " ";
   }
   std::cout << std::endl;
   
   std::cout << "\n" << std::endl;
   return 0;
}

// -----------------------------------------------------
// This Function does the basic information retrieval 
// One passes in the type of data requested and a reference to the data holder
// The function connects via the libwda function to recover the
// data and checks to insure there were not connection errors. 
// The function returns the error status, if null then success
//-----------------------------------------------------
  int ChannelMapSQLite::GetDataset(const std::string& table, int func(void*,int,char**,char**), void* data) const
  {
    std::string fullFileName;
    cet::search_path searchPath("FW_SEARCH_PATH");
    
    if (!searchPath.find_file(fDBFileName, fullFileName))
      throw cet::exception("ChannelMapSQLite::GetDataset") << "Can't find input file: '" << fDBFileName << "'\n";
    
    // Set up to open the database
    sqlite3* database;
    
    int rc = sqlite3_open(fullFileName.c_str(), &database);
    
    if (rc)
      throw cet::exception("ChannelMapSQLite::GetDataset") << "Can't open the database, return code:" << sqlite3_errmsg(database) << "'\n";
    
    // Test reading the database
    std::string select = "SELECT * FROM " + table;
    
    char*       zErrMsg = 0;
    
    rc = sqlite3_exec(database, select.c_str(), func, data, &zErrMsg);
    
    if( rc != SQLITE_OK ) 
      {
	std::cout << "ChannelMapSQLite::GetDataset: SQL error: " << zErrMsg << std::endl;
	sqlite3_free(zErrMsg);
      } 
    else 
      {
        std::cout << "ChannelMapSQLite::GetDataset: Successfully read database" << std::endl;
      }
    
    sqlite3_close(database);
    
    return 0;
  }
  
// -----------------------------------------------------
// The aim of this function is to build a map between the
// TPC Fragment IDs and the readout board IDs. Here we 
// expect there will be some number of boards per Fragment
//-----------------------------------------------------
int buildTPCFragmentIDToReadoutIDMap_callback(void* dataIn, int argc, char**argv, char** azColName)
{
    const unsigned int tpcIdentifier(0x00001000);

    IChannelMapping::TPCFragmentIDToReadoutIDMap& fragmentBoardMap = *(IChannelMapping::TPCFragmentIDToReadoutIDMap*)dataIn;

    // Include a by hand mapping of fragement ID to crate
    // Note that we now know we can get this from the "flanges" table... so an upgrade coming soon...
    using FlangeIDToCrateMap = std::map<size_t,std::string>;
    FlangeIDToCrateMap flangeIDToCrateMap;
    flangeIDToCrateMap[19]  = "WW01T";
    flangeIDToCrateMap[68]  = "WW01M";
    flangeIDToCrateMap[41]  = "WW01B";
    flangeIDToCrateMap[11]  = "WW02";
    flangeIDToCrateMap[17]  = "WW03";
    flangeIDToCrateMap[36]  = "WW04";
    flangeIDToCrateMap[18]  = "WW05";
    flangeIDToCrateMap[58]  = "WW06";
    flangeIDToCrateMap[71]  = "WW07";
    flangeIDToCrateMap[14]  = "WW08";
    flangeIDToCrateMap[25]  = "WW09";
    flangeIDToCrateMap[34]  = "WW10";
    flangeIDToCrateMap[67]  = "WW11";
    flangeIDToCrateMap[33]  = "WW12";
    flangeIDToCrateMap[87]  = "WW13";
    flangeIDToCrateMap[10]  = "WW14";
    flangeIDToCrateMap[59]  = "WW15";
    flangeIDToCrateMap[95]  = "WW16";
    flangeIDToCrateMap[22]  = "WW17";
    flangeIDToCrateMap[91]  = "WW18";
    flangeIDToCrateMap[61]  = "WW19";
    flangeIDToCrateMap[55]  = "WW20T";
    flangeIDToCrateMap[97]  = "WW20M";
    flangeIDToCrateMap[100] = "WW20B";
    flangeIDToCrateMap[83]  = "WE01T";
    flangeIDToCrateMap[85]  = "WE01M";
    flangeIDToCrateMap[7]   = "WE01B";
    flangeIDToCrateMap[80]  = "WE02";
    flangeIDToCrateMap[52]  = "WE03";
    flangeIDToCrateMap[32]  = "WE04";
    flangeIDToCrateMap[70]  = "WE05";
    flangeIDToCrateMap[74]  = "WE06";
    flangeIDToCrateMap[46]  = "WE07";
    flangeIDToCrateMap[81]  = "WE08";
    flangeIDToCrateMap[63]  = "WE09";
    flangeIDToCrateMap[30]  = "WE10";
    flangeIDToCrateMap[51]  = "WE11";
    flangeIDToCrateMap[90]  = "WE12";
    flangeIDToCrateMap[23]  = "WE13";
    flangeIDToCrateMap[93]  = "WE14";
    flangeIDToCrateMap[92]  = "WE15";
    flangeIDToCrateMap[88]  = "WE16";
    flangeIDToCrateMap[73]  = "WE17";
    flangeIDToCrateMap[1]   = "WE18";
    flangeIDToCrateMap[66]  = "WE19";
    flangeIDToCrateMap[48]  = "WE20T";
    flangeIDToCrateMap[13]  = "WE20M";
    flangeIDToCrateMap[56]  = "WE20B";
    flangeIDToCrateMap[94]  = "EW01T";
    flangeIDToCrateMap[77]  = "EW01M";
    flangeIDToCrateMap[72]  = "EW01B";
    flangeIDToCrateMap[65]  = "EW02";
    flangeIDToCrateMap[4]   = "EW03";
    flangeIDToCrateMap[89]  = "EW04";
    flangeIDToCrateMap[37]  = "EW05";
    flangeIDToCrateMap[76]  = "EW06";
    flangeIDToCrateMap[49]  = "EW07";
    flangeIDToCrateMap[60]  = "EW08";
    flangeIDToCrateMap[21]  = "EW09";
    flangeIDToCrateMap[6]   = "EW10";
    flangeIDToCrateMap[62]  = "EW11";
    flangeIDToCrateMap[2]   = "EW12";
    flangeIDToCrateMap[29]  = "EW13";
    flangeIDToCrateMap[44]  = "EW14";
    flangeIDToCrateMap[9]   = "EW15";
    flangeIDToCrateMap[31]  = "EW16";
    flangeIDToCrateMap[98]  = "EW17";
    flangeIDToCrateMap[38]  = "EW18";
    flangeIDToCrateMap[99]  = "EW19";
    flangeIDToCrateMap[53]  = "EW20T";
    flangeIDToCrateMap[82]  = "EW20M";
    flangeIDToCrateMap[35]  = "EW20B";
    flangeIDToCrateMap[96]  = "EE01T";
    flangeIDToCrateMap[28]  = "EE01M";
    flangeIDToCrateMap[16]  = "EE01T";
    flangeIDToCrateMap[69]  = "EE02";
    flangeIDToCrateMap[20]  = "EE02";
    flangeIDToCrateMap[79]  = "EE02";
    flangeIDToCrateMap[50]  = "EE02";
    flangeIDToCrateMap[45]  = "EE02";
    flangeIDToCrateMap[84]  = "EE02";
    flangeIDToCrateMap[42]  = "EE02";
    flangeIDToCrateMap[39]  = "EE02";
    flangeIDToCrateMap[26]  = "EE02";
    flangeIDToCrateMap[64]  = "EE02";
    flangeIDToCrateMap[43]  = "EE02";
    flangeIDToCrateMap[47]  = "EE02";
    flangeIDToCrateMap[15]  = "EE02";
    flangeIDToCrateMap[3]   = "EE02";
    flangeIDToCrateMap[27]  = "EE02";
    flangeIDToCrateMap[24]  = "EE02";
    flangeIDToCrateMap[40]  = "EE02";
    flangeIDToCrateMap[75]  = "EE02";
    flangeIDToCrateMap[86]  = "EE20T";
    flangeIDToCrateMap[54]  = "EE20M";
    flangeIDToCrateMap[8]   = "EE20B";

    unsigned int fragmentID = std::stol(argv[8],nullptr,16);
    
    if (fragmentID & tpcIdentifier)
    {
        if (fragmentBoardMap.find(fragmentID) == fragmentBoardMap.end())
        {
            unsigned int flangeID = std::stol(argv[1]);
            fragmentBoardMap[fragmentID].first = flangeIDToCrateMap[flangeID];
        }

        unsigned int readoutID = std::stol(argv[0]);
        fragmentBoardMap[fragmentID].second.emplace_back(readoutID);

    }

    return 0;
}

  // -----------------------------------------------------
  // The aim of this function is to build a map between the
  // TPC Fragment IDs and the readout board IDs. Here we 
  // expect there will be some number of boards per Fragment
  //-----------------------------------------------------
  int ChannelMapSQLite::BuildTPCFragmentIDToReadoutIDMap(TPCFragmentIDToReadoutIDMap& fragmentBoardMap) const
  {
    const std::string dataType("readout_boards");
    
    // Recover the data from the database
    int error = GetDataset(dataType,buildTPCFragmentIDToReadoutIDMap_callback,&fragmentBoardMap);
    
    // If there was an error the function above would have printed a message so bail out
    if (error)
      throw cet::exception("ChannelMapSQLite::BuildTPCFragmentIDToReadoutIDMap") << "Encountered error in reading the database: '" << error << "'\n";
    
    return error;
  }
  
  // -----------------------------------------------------
  // The aim of this function is to build a map between the
  // TPC readout board IDs and the associated channels. So for
  // each readout board ID we expect a number of channels back
  // from the data base. So the returned data structure will 
  // be a map of readout ID to a vector of channels
  //-----------------------------------------------------
  const unsigned int CHANNELSPERBOARD = 64;
  int buildTPCReadoutBoardToChannelMap_callback(void* dataIn, int argc, char**argv, char** azColName)
  {
    IChannelMapping::TPCReadoutBoardToChannelMap& rbChanMap = *(IChannelMapping::TPCReadoutBoardToChannelMap*)dataIn;
    
    unsigned int readoutBoardID = std::stol(argv[2]);
    
    if (rbChanMap.find(readoutBoardID) == rbChanMap.end())
      {
        unsigned int readoutBoardSlot = std::stol(argv[4]);
	
        rbChanMap[readoutBoardID].first = readoutBoardSlot;
        rbChanMap[readoutBoardID].second.resize(CHANNELSPERBOARD);
      }
    
    unsigned int channelNum = std::stol(argv[5]);
    unsigned int channelID  = std::stol(argv[0]);
    
    std::string fragmentBuffer = argv[10];
    
    // Make sure lower case... (sigh...)
    std::transform(fragmentBuffer.begin(),fragmentBuffer.end(),fragmentBuffer.begin(),[](char c){return std::tolower(c);});
    
    unsigned int plane(3);
    
    if      (fragmentBuffer.find("collection")  != std::string::npos) plane = 2;
    else if (fragmentBuffer.find("induction 2") != std::string::npos) plane = 1;
    else if (fragmentBuffer.find("induction 1") != std::string::npos) plane = 0;

    if (plane > 2) std::cout << "YIKES!!! Plane is " << plane << " for channel " << channelID << " with type " << std::string(fragmentBuffer) << std::endl;
    
    rbChanMap[readoutBoardID].second[channelNum] = IChannelMapping::ChannelPlanePair(channelID,plane);
    
    return 0;
  }
  
  int ChannelMapSQLite::BuildTPCReadoutBoardToChannelMap(TPCReadoutBoardToChannelMap& rbChanMap) const
  {
    const std::string  dataType("daq_channels");
    
    // Recover the data from the database
    int error = GetDataset(dataType,buildTPCReadoutBoardToChannelMap_callback,&rbChanMap);
    
    // If there was an error the function above would have printed a message so bail out
    if (error)
      throw cet::exception("ChannelMapSQLite::BuildTPCReadoutBoardToChannelMap") << "Encountered error in reading the database: '" << error << "'\n";
    
    return error;
  }
  
  //******************* PMT Channel Mapping ***********************
  int buildFragmentToDigitizerChannelMap_callback(void* dataIn, int argc, char**argv, char** azColName)
  {
    IChannelMapping::FragmentToDigitizerChannelMap& fragmentToDigitizerChannelMap = *(IChannelMapping::FragmentToDigitizerChannelMap*)dataIn;
    
    // Start extracting info
    //    std::string  digitizerBuffer    = argv[8];
    unsigned int fragmentID         = std::stol(argv[18]);
    unsigned int digitizerChannelNo = std::stol(argv[9]);
    unsigned int channelID          = std::stol(argv[17]);

    // Read the laser channel 
    std::string laserChannelLabel  = argv[7]; // format is L-<number> . <number> is int from [1-41]
    unsigned int laserChannel      = std::stol(laserChannelLabel.substr(2,4));  // try-catch syntax for stol or not necessary ? 


    // Fill the map
    fragmentToDigitizerChannelMap[fragmentID].emplace_back(digitizerChannelNo,channelID,laserChannel);
    
    return 0;
  }
  
  int ChannelMapSQLite::BuildFragmentToDigitizerChannelMap(FragmentToDigitizerChannelMap& fragmentToDigitizerChannelMap) const
  {
    // clearing is cleansing
    fragmentToDigitizerChannelMap.clear();
    // Recover the information from the database on the mapping 
    const std::string  dataType("pmt_placements");
    
    // Recover the data from the database
    int error = GetDataset(dataType,buildFragmentToDigitizerChannelMap_callback,&fragmentToDigitizerChannelMap);
    
    // If there was an error the function above would have printed a message so bail out
    if (error)
      throw cet::exception("ChannelMapSQLite::BuildFragmentToDigitizerChannelMap") << "Encountered error in reading the database: '" << error << "'\n";
    
    return error;
  }

  
  //******************* CRT Channel Mapping ***********************
  
  int buildCRTChannelIDToHWtoSimMacAddressPairMap_callback(void* dataIn, int argc, char**argv, char** azColName)
  {
    IChannelMapping::CRTChannelIDToHWtoSimMacAddressPairMap& crtChannelIDToHWtoSimMacAddressPairMap = *(IChannelMapping::CRTChannelIDToHWtoSimMacAddressPairMap*)dataIn;

    /*    
    // Start extracting info
    unsigned int channelID         = std::stol(argv[10]);
    unsigned int simmacaddress     = std::stol(argv[11]);
    unsigned int hwmacaddress      = std::stol(argv[12]);
    */
     
    unsigned int channelID         = strcmp(argv[10],"None")==0 ? 0 : std::stol(argv[10]);
    unsigned int simmacaddress     = strcmp(argv[11],"None")==0 ? 0 : std::stol(argv[11]);
    unsigned int hwmacaddress      = strcmp(argv[12],"None")==0 ? 0 : std::stol(argv[12]); 
   
    // Fill the map
    crtChannelIDToHWtoSimMacAddressPairMap[channelID]=std::make_pair(hwmacaddress, simmacaddress);

    return 0;
  }
  
  int ChannelMapSQLite::BuildCRTChannelIDToHWtoSimMacAddressPairMap(CRTChannelIDToHWtoSimMacAddressPairMap& crtChannelIDToHWtoSimMacAddressPairMap) const
  {
    // clearing is cleansing
    crtChannelIDToHWtoSimMacAddressPairMap.clear();
    // Recover the information from the database on the mapping
    const std::string  dataType("feb_channels");
    
    // Recover the data from the database
    int error = GetDataset(dataType,buildCRTChannelIDToHWtoSimMacAddressPairMap_callback,&crtChannelIDToHWtoSimMacAddressPairMap);
    
    // If there was an error the function above would have printed a message so bail out
    if (error)
      throw cet::exception("ChannelMapSQLite::BuildCRTChannelIDToHWtoSimMacAddressPairMap") << "Encountered error in reading the database: '" << error << "'\n";
    
    return error;
  }

  // Top CRT
  //---------------------------------- 
  int buildTopCRTHWtoSimMacAddressPairMap_callback(void* dataIn, int argc, char**argv, char** azColName)
  {
    IChannelMapping::TopCRTHWtoSimMacAddressPairMap& topcrtHWtoSimMacAddressPairMap = *(IChannelMapping::TopCRTHWtoSimMacAddressPairMap*)dataIn;
    
    // Start extracting info
    unsigned int simmacaddress     =  strcmp(argv[41],"None")==0 ? 0 : std::stol(argv[41]);
    unsigned int hwmacaddress      =  strcmp(argv[3], "None")==0 ? 0 : std::stol(argv[3]);
        
    // Fill the map
    topcrtHWtoSimMacAddressPairMap[hwmacaddress] = simmacaddress;

    return 0;
  }

  int ChannelMapSQLite::BuildTopCRTHWtoSimMacAddressPairMap(TopCRTHWtoSimMacAddressPairMap& topcrtHWtoSimMacAddressPairMap) const
  {
    // clearing is cleansing
    topcrtHWtoSimMacAddressPairMap.clear();
    // Recover the information from the database on the mapping
    const std::string  dataType("crtfeb");
    
    // Recover the data from the database
    int error = GetDataset(dataType,buildTopCRTHWtoSimMacAddressPairMap_callback,&topcrtHWtoSimMacAddressPairMap);
    
    // If there was an error the function above would have printed a message so bail out
    if (error)
      throw cet::exception("ChannelMapSQLite::BuildTopCRTHWtoSimMacAddressPairMap") << "Encountered error in reading the database: '" << error << "'\n";
    
    return error;
  }
    

  //******************* Accessing Side CRT Calibration Database ***********************

  int ChannelMapSQLite::BuildSideCRTCalibrationMap(SideCRTChannelToCalibrationMap& sideCRTChannelToCalibrationMap) const
  {
    // clearing is cleansing
    sideCRTChannelToCalibrationMap.clear();
        
    std::string fullFileName;
    cet::search_path searchPath("FW_SEARCH_PATH");

    
    if (!searchPath.find_file(fCalibDBFileName+".db", fullFileName)){
      std::cout << "******* Succesfully found calibration input file:  " << fCalibDBFileName <<  std::endl;
      throw cet::exception("ChannelMapSQLite::GetDataset") << "Can't find calibration input file: '" << fCalibDBFileName << "'\n";
    }
    
    lariov::DBFolder database(fCalibDBFileName, "", "", fTag, true, false);
    
    database.UpdateData(1638918271*1e9);
    
    std::vector<unsigned int> channels;
    database.GetChannelList(channels);
    
    for (auto it = channels.begin(); it != channels.end(); ++it) {
      
      long mac5, chan;
      double  gain, ped;
      
      // Start extracting info
      database.GetNamedChannelData(*it, "mac5", mac5);
      database.GetNamedChannelData(*it, "localchannel", chan);
      database.GetNamedChannelData(*it, "gain", gain);
      database.GetNamedChannelData(*it, "pedestal", ped);
      
      // Fill the map
      sideCRTChannelToCalibrationMap.insert(std::make_pair(std::make_pair((int)mac5,(int)chan), std::make_pair(gain,ped)));
    }

    return 0;
  }

  
  DEFINE_ART_CLASS_TOOL(ChannelMapSQLite)
} // namespace lar_cluster3d

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
   *  @brief Define the returned data structures for a mapping between PMT Fragment IDs
   *         and the related crate and readout information. 
   *         Then define the function interface to fill these data structures 
   */
  virtual int BuildFragmentToDigitizerChannelMap(FragmentToDigitizerChannelMap&) const override;
  
private:
  // Recover data from postgres database
  int GetDataset(const std::string&, int func(void*,int,char**,char**), void*) const;
  

  std::string fDBFileName;          //< File name of our sqlite database
  std::string fTag;                 //< Tag for conditioned database
};

ChannelMapSQLite::ChannelMapSQLite(fhicl::ParameterSet const &pset)
{
    fDBFileName      = pset.get<std::string>("DBFileName");
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

  DEFINE_ART_CLASS_TOOL(ChannelMapSQLite)
} // namespace lar_cluster3d

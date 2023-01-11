 /*
 * @file   ChannelMapPostGres_tool.cc
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

// std includes
#include <string>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace sbndDB {
/**
 *  @brief  ChannelMapPostGres class definiton
 */
class ChannelMapPostGres : virtual public IChannelMapping
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit ChannelMapPostGres(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~ChannelMapPostGres();

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

  /**
   *  @brief Define the returned data structures for a mapping between Side CRT Channels and
   *         their calibration values.
   *         Then define the function interface to fill these data structures
   */
  virtual int BuildSideCRTCalibrationMap(SideCRTChannelToCalibrationMap&) const override;

private:
    // Recover data from postgres database
    int GetDataset(const std::string&, const std::string&, const std::string&, Dataset&) const;
    int GetCRTCaldata(const std::string&, const std::string&, Dataset&) const;
    uint32_t fNothing;     //< Nothing

};

ChannelMapPostGres::ChannelMapPostGres(fhicl::ParameterSet const &pset)
{
    fNothing = pset.get<uint32_t>("Nothing");

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

ChannelMapPostGres::~ChannelMapPostGres()
{
}

// -----------------------------------------------------
// This Function does the basic information retrieval 
// One passes in the type of data requested and a reference to the data holder
// The function connects via the libwda function to recover the
// data and checks to insure there were not connection errors. 
// The function returns the error status, if null then success
//-----------------------------------------------------
int ChannelMapPostGres::GetDataset(const std::string& name, const std::string& url, const std::string& dataType, Dataset& dataSet) const
{
    const int   timeout(200);
    std::string dburl = url + "&t=" + dataType;
    int error(0);
    dataSet = getDataWithTimeout(dburl.data(),name.data(),timeout,&error);
    if (error)
    {
        std::string errorMsg = "Database access GetDataset failed with error " + std::to_string(error) + "\nDB url: " 
                             + dburl + ", name: " + name + ", type: " + dataType;
        std::cout << "****> Database retrieval error, code: " << error << std::endl;
        throw std::runtime_error(errorMsg);
    }

    return error;
}

  //Adding this so I can control how this workflow goes -TB
  int ChannelMapPostGres::GetCRTCaldata(const std::string& name, const std::string& url, Dataset& dataSet) const
  {
    const int   timeout(200);
    int error(0);
    dataSet = getDataWithTimeout(url.c_str(),name.c_str(),timeout,&error);
    if (error)
      {
	std::string errorMsg = "Database access GetDataset failed with error " + std::to_string(error) + "\nurl: "
	  + url + ", name: " + name;
	std::cout << "****> Database retrieval error, code: " << error << std::endl;
        throw std::runtime_error(errorMsg);
      }

    return error;
  }


// -----------------------------------------------------
// The aim of this function is to build a map between the
// TPC Fragment IDs and the readout board IDs. Here we 
// expect there will be some number of boards per Fragment
//-----------------------------------------------------
int ChannelMapPostGres::BuildTPCFragmentIDToReadoutIDMap(TPCFragmentIDToReadoutIDMap& fragmentBoardMap) const
{
    const unsigned int tpcIdentifier(0x00001000);
    const std::string  name("sbnd_hw_readoutboard");
    const std::string  dburl("https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query?dbname=sbnd_hardware_prd");
    const std::string  dataType("readout_boards");
    Dataset            dataset;
    // Recover the data from the database
    int error = GetDataset(name,dburl,dataType,dataset);
    // Include a by hand mapping of fragement ID to crate
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
    // If there was an error the function above would have printed a message so bail out
    if (error) throw(std::exception());
    // Loop through the data to recover the channels
    // NOTE that we skip the first row because that is just the labels
    for(int row = 1; row < getNtuples(dataset); row++)
    {
        // Recover the row
        Tuple tuple = getTuple(dataset, row);
        if (tuple != NULL)
        {
            // Note that the fragment ID is stored in the database as a string which reads as a hex number
            // Meaning we have to read back as a string and decode to get the numerical value. 
            char fragmentBuffer[16];
            getStringValue(tuple, 8, fragmentBuffer, sizeof(fragmentBuffer), &error);
            if (error) throw std::runtime_error("Encountered error in trying to recover FragmentID from database");
            std::string fragmentIDString(fragmentBuffer,4);
            unsigned int fragmentID = std::stol(fragmentIDString,nullptr,16);
            if (!(fragmentID & tpcIdentifier)) continue;
            if (fragmentBoardMap.find(fragmentID) == fragmentBoardMap.end())
            {
                unsigned int flangeID = getLongValue(tuple, 1, &error);
                if (error) throw std::runtime_error("Encountered error in trying to recover Board Flange ID from database");
                fragmentBoardMap[fragmentID].first = flangeIDToCrateMap[flangeID];
            }
            unsigned int readoutID = getLongValue(tuple, 0, &error);
            if (error) throw std::runtime_error("Encountered error in trying to recover Board ReadoutID from database");
            fragmentBoardMap[fragmentID].second.emplace_back(readoutID);
            releaseTuple(tuple);
        }
    }
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
int ChannelMapPostGres::BuildTPCReadoutBoardToChannelMap(TPCReadoutBoardToChannelMap& rbChanMap) const
{
    const std::string  name("sbnd_hardware_prd");
    const std::string  dburl("https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query?dbname=sbnd_hardware_prd");
    const std::string  dataType("daq_channels");
    Dataset            dataset;
    // Recover the data from the database
    int error = GetDataset(name,dburl,dataType,dataset);
    // If there was an error the function above would have printed a message so bail out
    if (error) throw std::runtime_error("Encountered error accessing the database with GetDataset");
    // Loop through the data to recover the channels, making sure to skip the first (header) row
    for(int row = 1; row < getNtuples(dataset); row++)
    {
        // Recover the row
        Tuple tuple = getTuple(dataset, row);
        if (tuple != NULL)
        {
            unsigned int readoutBoardID   = getLongValue(tuple, 2, &error);
            if (error) throw std::runtime_error("Encountered error when trying to read Board ReadoutID");
            if (rbChanMap.find(readoutBoardID) == rbChanMap.end())
            {
                unsigned int readoutBoardSlot = getLongValue(tuple, 4, &error);
                if (error) throw std::runtime_error("Encountered error when trying to read Board Readout slot");
                rbChanMap[readoutBoardID].first = readoutBoardSlot;
                rbChanMap[readoutBoardID].second.resize(CHANNELSPERBOARD);
            }
            unsigned int channelNum = getLongValue(tuple, 5, &error);
            if (error) throw std::runtime_error("Encountered error when trying to read channel number");
            unsigned int channelID = getLongValue(tuple, 0, &error);
            if (error) throw std::runtime_error("Encountered error when recovering the channel ID list");
            // Recover the plane identifier 
            char fragmentBuffer[16];
            getStringValue(tuple, 10, fragmentBuffer, sizeof(fragmentBuffer), &error);
            if (error) throw std::runtime_error("Encountered error when trying to read plane type");
            // Make sure lower case... (sigh...)
            for(size_t charIdx = 0; charIdx < sizeof(fragmentBuffer); charIdx++) fragmentBuffer[charIdx] = tolower(fragmentBuffer[charIdx]);
            unsigned int plane(3);
            if      (strstr(fragmentBuffer,"collection"))  plane = 2;
            else if (strstr(fragmentBuffer,"induction 2")) plane = 1;
            else if (strstr(fragmentBuffer,"induction 1")) plane = 0;
            if (plane > 2) std::cout << "YIKES!!! Plane is " << plane << " for channel " << channelID << " with type " << std::string(fragmentBuffer) << std::endl;
            rbChanMap[readoutBoardID].second[channelNum] = ChannelPlanePair(channelID,plane);
        }
    }

    return error;
}

//******************* PMT Channel Mapping ***********************
int ChannelMapPostGres::BuildFragmentToDigitizerChannelMap(FragmentToDigitizerChannelMap& fragmentToDigitizerChannelMap) const
{
    // clearing is cleansing
    fragmentToDigitizerChannelMap.clear();
    // Recover the information from the database on the mapping 
    const std::string  name("Pmt_placement");
    const std::string  dburl("https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query?dbname=sbnd_hardware_prd");
    const std::string  dataType("pmt_placements");
    Dataset            dataset;
    // Recover the data from the database
    int error = GetDataset(name,dburl,dataType,dataset);
    // If there was an error the function above would have printed a message so bail out
    if (error) throw(std::exception());
    // Ok, now we can start extracting the information
    // We do this by looping through the database and building the map from that
    for(int row = 1; row < getNtuples(dataset); row++)
    {
        // Recover the row
        Tuple tuple = getTuple(dataset, row);
        if (tuple != NULL)
        {
            char digitizerBuffer[10];
            // Recover the digitizer label first 
            getStringValue(tuple, 8, digitizerBuffer, sizeof(digitizerBuffer), &error);
            if (error) throw std::runtime_error("Encountered error when trying to recover the PMT digitizer label");
            std::string digitizerLabel(digitizerBuffer, 8); //sizeof(digitizerBuffer));
            // Recover the fragment id
            unsigned fragmentID = getLongValue(tuple, 18, &error);
            if (error) throw std::runtime_error("Encountered error when trying to recover the PMT fragment id");
            // Now recover the digitizer channel number
            unsigned int digitizerChannelNo = getLongValue(tuple, 9, &error);
            if (error) throw std::runtime_error("Encountered error when trying to recover the PMT digitizer channel number");
            // Get the LArsoft channel ID
            unsigned int channelID = getLongValue(tuple, 17, &error);
            if (error) throw std::runtime_error("Encountered error when trying to recover the PMT channel ID");
            // Get the laserChannelNumber 
            char laserChannelBuffer[10];
            getStringValue(tuple, 7, laserChannelBuffer, sizeof(laserChannelBuffer), &error);
            if (error) throw std::runtime_error("Encountered error when trying to recover the PMT laser channel label");
            std::string laserChannelLabel(laserChannelBuffer, 2, sizeof(laserChannelBuffer)); //sizeof(digitizerBuffer));
            unsigned int laserChannel = std::stol(laserChannelLabel);  // try-catch syntax for stol or not necessary ? 

            // Fill the map
            fragmentToDigitizerChannelMap[fragmentID].emplace_back(digitizerChannelNo,channelID,laserChannel);
            releaseTuple(tuple);
        }
    }

    return error;
}

  
//******************* CRT Channel Mapping ***********************
  
  int ChannelMapPostGres::BuildCRTChannelIDToHWtoSimMacAddressPairMap(CRTChannelIDToHWtoSimMacAddressPairMap& crtChannelIDToHWtoSimMacAddressPairMap) const
  {
    // clearing is cleansing
    crtChannelIDToHWtoSimMacAddressPairMap.clear();
    // Recover the information from the database on the mapping 
    const std::string  name("Feb_channels");
    const std::string  dburl("https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query?dbname=sbnd_hardware_prd");
    const std::string  dataType("feb_channels");
    Dataset            dataset;
    // Recover the data from the database
    int error = GetDataset(name,dburl,dataType,dataset);
    // If there was an error the function above would have printed a message so bail out
    if (error) throw(std::exception());
    // Ok, now we can start extracting the information
    // We do this by looping through the database and building the map from that
    for(int row = 1; row < getNtuples(dataset); row++)
    {
        // Recover the row
        Tuple tuple = getTuple(dataset, row);
        if (tuple != NULL)
        {
	  // Recover the simmacaddress
            unsigned int simmacaddress = getLongValue(tuple, 11, &error);
            if (error) throw std::runtime_error("Encountered error when trying to recover the CRT simmacaddress");
            // Now recover the hwmacaddress
            unsigned int hwmacaddress = getLongValue(tuple, 12, &error);
            if (error) throw std::runtime_error("Encountered error when trying to recover the CRT hwmacaddress");
            // Finally, get the LArsoft channel ID
            unsigned int channelID = getLongValue(tuple, 10, &error);
            if (error) throw std::runtime_error("Encountered error when trying to recover the CRT channel ID");
            // Fill the map
            crtChannelIDToHWtoSimMacAddressPairMap[channelID]=std::make_pair(hwmacaddress,simmacaddress);
            releaseTuple(tuple);
        }
    }

    return error;
  }

   
  /// Top CRT harware mac5 to software mac5 relation
  //----------------------------------------------------
  int ChannelMapPostGres::BuildTopCRTHWtoSimMacAddressPairMap(TopCRTHWtoSimMacAddressPairMap& topcrtHWtoSimMacAddressPairMap) const
  {
    // clearing is cleansing
    topcrtHWtoSimMacAddressPairMap.clear();
    // Recover the information from the database on the mapping 
    const std::string  name("topcrt_febs");
    const std::string  dburl("https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query?dbname=sbnd_hardware_prd");
    const std::string  dataType("crtfeb");
    Dataset            dataset;
    // Recover the data from the database
    int error = GetDataset(name,dburl,dataType,dataset);
    // If there was an error the function above would have printed a message so bail out
    if (error) throw(std::exception());
    // Ok, now we can start extracting the information
    // We do this by looping through the database and building the map from that
    for(int row = 1; row < getNtuples(dataset); row++)
    {
        // Recover the row
        Tuple tuple = getTuple(dataset, row);
        if (tuple != NULL)
        {
	  // Recover the simmacaddress
            unsigned int simmacaddress = getLongValue(tuple, 41, &error);
            if (error) throw std::runtime_error("Encountered error when trying to recover the CRT simmacaddress");
            // Now recover the hwmacaddress
            unsigned int hwmacaddress = getLongValue(tuple, 3, &error);
            if (error) throw std::runtime_error("Encountered error when trying to recover the CRT hwmacaddress");

            // Fill the map
            topcrtHWtoSimMacAddressPairMap[hwmacaddress] = simmacaddress;
            releaseTuple(tuple);
        }
    }

    return error;
  }
  

  /// Accesing the Side CRT charge calibration from the postgresql database
  //------------------------------------------------------------------------

  int ChannelMapPostGres::BuildSideCRTCalibrationMap(SideCRTChannelToCalibrationMap& sideCRTChannelToCalibrationMap) const {
    //    sideCRTChannelToCalibrationMap.clear();
    const std::string  name("SideCRT_calibration_data");
    const std::string dburl("https://dbdata0vm.fnal.gov:9443/sbnd_con_prod/app/data?f=crt_gain_reco_data&t=1638918270");
    Dataset            ds;
    // Recover the data from the database
    int error = GetCRTCaldata(name,dburl,ds);
    if (error)
      { fprintf(stderr, "error code=%d\n", error);    perror("error message"); }
    if (getHTTPstatus(ds) != 200)
      { fprintf(stderr, "HTTP code=%ld, message: '%s'\n", getHTTPstatus(ds), getHTTPmessage(ds)); }
    int mac5, chan;
    double gain, ped;
    int err;
    int nrows =  getNtuples(ds);
    int ncols;
    Tuple tup;
    for (int rows = 0; rows < nrows; rows++ ){
      tup = getTuple(ds, rows);                                           // Get the row with double array
      ncols = getNfields(tup);//check number of columns
      if(ncols <5) continue;//first few rows aren't actual data and have ncols==1, this excludes those
      if (tup != NULL)/*check that the tup is not empty before proceeding*/ {
        //assign values from the database to variables so we can use them
        mac5 = (int)getDoubleValue(tup,1,&err);
        chan = (int)getDoubleValue(tup,2,&err);
        gain = getDoubleValue(tup,3,&err);
        ped = getDoubleValue(tup,4,&err);
        //This line adds the association between the two pairs to the map object
        sideCRTChannelToCalibrationMap.insert(std::make_pair(std::make_pair(mac5,chan), std::make_pair(gain,ped)));

        releaseTuple(tup);
      }//end if(tup != NULL)
      else releaseTuple(tup);
    }//end loop over rows

    return error;
  }

 
DEFINE_ART_CLASS_TOOL(ChannelMapPostGres)
} // namespace sbndDB

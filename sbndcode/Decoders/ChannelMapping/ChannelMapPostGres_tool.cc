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
     *  @brief Define the returned data structures for a mapping between PMT Fragment IDs
     *         and the related crate and readout information. 
     *         Then define the function interface to fill these data structures 
     */
    virtual int BuildFragmentToDigitizerChannelMap(FragmentToDigitizerChannelMap&) const override;

private:
    // Recover data from postgres database
    int GetDataset(const std::string&, const std::string&, const std::string&, Dataset&) const;
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
//    const int   timeout(200);
    const int   timeout(10);
    std::string dburl = url + "&t=" + dataType;
    int error(0);

    std::cout << "url = " << dburl.data() << std::endl; 

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


//******************* PMT Channel Mapping ***********************
int ChannelMapPostGres::BuildFragmentToDigitizerChannelMap(FragmentToDigitizerChannelMap& fragmentToDigitizerChannelMap) const
{

    

    // clearing is cleansing
    fragmentToDigitizerChannelMap.clear();
    // Recover the information from the database on the mapping 
    const std::string  name("Pmt_placement");
//    const std::string  dburl("https://dbdata0vm.fnal.gov:9443/QE/hw/app/SQ/query?dbname=sbnd_hardware_prd");
//    const std::string  dataType("pmt_placements");
    
    //LAN: No online database so use google drive for now
    const std::string  dburl("https://docs.google.com/spreadsheets/d/183Ri3wMdBIsXpFP-jGajZZXQ-9_YCt-aifPQ5MB_VSQ/gviz/tq?tqx=out:csv&sheet=pmt_placements");
   const std::string dataType("");

    Dataset            dataset;
    // Recover the data from the database
    int error = GetDataset(name,dburl,dataType,dataset);
    // If there was an error the function above would have printed a message so bail out
    if (error) throw(std::exception());
    // Ok, now we can start extracting the information
    // We do this by looping through the database and building the map from that

//    for(int row = 1; row < getNtuples(dataset); row++)
    for(int row = 1; row < getNtuples(dataset) + 2; row++)
    {
        // Recover the row
        Tuple tuple = getTuple(dataset, row);
        if (tuple != NULL)
        {
            // Recover the digitizer label first 
            char digitizerBuffer[10];
            getStringValue(tuple, 1, digitizerBuffer, sizeof(digitizerBuffer), &error);
            if (error) throw std::runtime_error("Encountered error when trying to recover the PMT digitizer label");
            std::string digitizerLabel(digitizerBuffer, 8); //sizeof(digitizerBuffer));
          
            // Recover the fragment id
            unsigned fragmentID = getLongValue(tuple, 2, &error);
            if (error) throw std::runtime_error("Encountered error when trying to recover the PMT fragment id");
 
            // Now recover the digitizer channel number
            unsigned int digitizerChannelNo = getLongValue(tuple, 3, &error);
            if (error) throw std::runtime_error("Encountered error when trying to recover the PMT digitizer channel number");

            // Get the LArsoft channel ID
            unsigned int channelID = getLongValue(tuple, 4, &error);
            if (error) throw std::runtime_error("Encountered error when trying to recover the PMT channel ID");

            // Get the laserChannelNumber 
            char laserChannelBuffer[10];
            getStringValue(tuple, 5, laserChannelBuffer, sizeof(laserChannelBuffer), &error);
            if (error) throw std::runtime_error("Encountered error when trying to recover the PMT laser channel label");

            std::string laserChannelLabel(laserChannelBuffer, 2, sizeof(laserChannelBuffer)); //sizeof(digitizerBuffer));
            unsigned int laserChannel = std::stol(laserChannelLabel);  // try-catch syntax for stol or not necessary ? 

            // Fill the map
            fragmentToDigitizerChannelMap[fragmentID].emplace_back(digitizerChannelNo,channelID,laserChannel);
            releaseTuple(tuple);

            std::cout << "row #: " << row << ", fragmentID: " << fragmentID << std::endl; 
       }
    }

    return error;
}

DEFINE_ART_CLASS_TOOL(ChannelMapPostGres)
} // namespace sbndDB

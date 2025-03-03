/*
 *  Service for the PMT Calibration Database.
 *  Andrea Scarpelli (ascarpell@bnl.gov), Matteo Vicenzi (mvicenzi@bnl.gov)
 */
// Ported from icaruscode to SBND by Alejandro Sanchez-Castillo, Jan. 2025


// Framework includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Principal/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Local
#include "sbndcode/DatabaseInterface/PMTCalibrationDatabase.h"
#include "sbndcode/DatabaseInterface/PMTCalibrationDatabaseProvider.h"

// Database interface helpers
#include "larevt/CalibrationDBI/Providers/DBFolder.h"
#include "larevt/CalibrationDBI/IOVData/TimeStampDecoder.h"

// C/C++ standard libraries
#include <string>
#include <vector>

//--------------------------------------------------------------------------------

sbndDB::PMTCalibrationDatabaseProvider::PMTCalibrationDatabaseProvider
    (const fhicl::ParameterSet& pset) 
    : fVerbose{ pset.get<bool>("Verbose", false) }
    , fLogCategory{ pset.get<std::string>("LogCategory", "PMTTimingCorrection") }
    { 
        fhicl::ParameterSet const tags{ pset.get<fhicl::ParameterSet>("CorrectionTags") };
        fCablesTag  = tags.get<std::string>("CablesTag");
        fDatabaseTimeStamp = tags.get<long>("DatabaseTimeStamp");
        fTableName = tags.get<std::string>("TableName");
        fVariabletoread = tags.get<std::string>("VariableToRead");
        fSERLength = tags.get<size_t>("SERLength");
        if( fVerbose ) mf::LogInfo(fLogCategory) << "Database tags for timing corrections:\n"
						 << "Cables corrections  " << fCablesTag << "\n"  ;             
    }

// -------------------------------------------------------------------------------

uint64_t sbndDB::PMTCalibrationDatabaseProvider::RunToDatabaseTimestamp( uint32_t run ) const{

   // Run number to timestamp used in the db
   // DBFolder.h only takes 19 digit (= timestamp in nano second),
   // but ICARUS tables are currently using run numbers
   // Step 1) Add 1000000000 to the run number; e.g., run XXXXX -> 10000XXXXX
   // Step 2) Multiply 1000000000
   uint64_t runNum = uint64_t(run);
   uint64_t timestamp = runNum+1000000000;
   timestamp *= 1000000000;

   if( fVerbose ) mf::LogInfo(fLogCategory) << "Run " << runNum << " corrections from DB timestamp " << timestamp;
   
   return timestamp;
}

// -------------------------------------------------------------------------------

/// Function to look up the calibration database at the table holding the pmt hardware cables corrections
void sbndDB::PMTCalibrationDatabaseProvider::ReadPMTCalibration( uint32_t run ) { 
    
    // pmt_cables_delay: delays of the cables relative to trigger 
    // and reset distribution

    const std::string dbname(fTableName);
    lariov::DBFolder db(dbname, "", "", fCablesTag, true, false);

    bool ret = db.UpdateData( fDatabaseTimeStamp ); // select table based on run number  
    mf::LogDebug(fLogCategory) << dbname + " corrections" << (ret? "": " not") << " updated for run " << run;
    mf::LogTrace(fLogCategory)
           << "Fetched IoV [ " << db.CachedStart().DBStamp() << " ; " << db.CachedEnd().DBStamp()
           << " ] to cover t=" << RunToDatabaseTimestamp(run)
           << " [=" << lariov::TimeStampDecoder::DecodeTimeStamp(RunToDatabaseTimestamp(run)).DBStamp() << "]";
    
    std::vector<unsigned int> channelList;
    if (int res = db.GetChannelList(channelList); res != 0) {
        throw cet::exception
        ( "PMTTimingCorrectionsProvider" ) << "GetChannelList() returned " << res << " on run " << run << " query in " << dbname << "\n";
    }
    
    if (channelList.empty()) {
        throw cet::exception("PMTTimingCorrectionsProvider") << "Got an empty channel list for run " << run << " in " << dbname << "\n";
    }
    
    for( auto channel : channelList ) {
        // Read breakout box
        long _breakoutbox = 0;
        int error  = db.GetNamedChannelData( channel, "breakout_box", _breakoutbox );
        if( error ) throw cet::exception( "PMTTimingCorrectionsProvider" ) << "Encountered error (code " << error << ") while trying to access 'breakout_box' on table " << dbname << "\n";
        fPMTCalibrationData[channel].breakoutBox = static_cast<int>(_breakoutbox);

        // Read caen digitizer
        long _caen_digitizer=0;
        error  = db.GetNamedChannelData( channel, "caen_digitizer", _caen_digitizer );
        if( error ) throw cet::exception( "PMTTimingCorrectionsProvider" ) << "Encountered error (code " << error << ") while trying to access 'caen_digitizer' on table " << dbname << "\n";
        fPMTCalibrationData[channel].caenDigitizer = static_cast<int>(_caen_digitizer);

        // Read caen digitizer channel
        long _caen_digitizer_channel=0;
        error  = db.GetNamedChannelData( channel, "caen_digitizer_channel", _caen_digitizer_channel );
        if( error ) throw cet::exception( "PMTTimingCorrectionsProvider" ) << "Encountered error (code " << error << ") while trying to access 'caen_digitizer_channel' on table " << dbname << "\n";
        fPMTCalibrationData[channel].caenDigitizerChannel = static_cast<int>(_caen_digitizer_channel);

        // Read total transit time
        double _total_transit_time =0.;
        error  = db.GetNamedChannelData( channel, "total_transit_time", _total_transit_time );
        if( error ) throw cet::exception( "PMTTimingCorrectionsProvider" ) << "Encountered error (code " << error << ") while trying to access 'total_transit_time' on table " << dbname << "\n";
        fPMTCalibrationData[channel].totalTransitTime = _total_transit_time;
       
        // Read spe amplitude
        double _spe_amplitude=0.;
        error  = db.GetNamedChannelData( channel, "spe_amp", _spe_amplitude );
        if( error ) throw cet::exception( "PMTTimingCorrectionsProvider" ) << "Encountered error (code " << error << ") while trying to access 'spe_amplitude' on table " << dbname << "\n";
        fPMTCalibrationData[channel].spe_amplitude = _spe_amplitude;
   
        // Read gauss filter power
        double _gauss_w_wc_power=0;
        error  = db.GetNamedChannelData( channel, "gauss_w_wc_power", _gauss_w_wc_power );
        if( error ) throw cet::exception( "PMTTimingCorrectionsProvider" ) << "Encountered error (code " << error << ") while trying to access 'gauss_w_wc_power' on table " << dbname << "\n";
        fPMTCalibrationData[channel].gauss_wc_power = _gauss_w_wc_power;

        // Read gauss filter wc
        double _gauss_wc = 0.;
        error  = db.GetNamedChannelData( channel, "gauss_wc", _gauss_wc );
        if( error ) throw cet::exception( "PMTTimingCorrectionsProvider" ) << "Encountered error (code " << error << ") while trying to access 'gauss_wc' on table " << dbname << "\n";
        fPMTCalibrationData[channel].gauss_wc = _gauss_wc;

        // Read SER
        std::vector<double> _ser;
        std::string name_base = "ser_vec_";
        for(size_t i=0; i<fSERLength; i++)
        {
            std::string entry_num = std::to_string(i);
            std::string entry_name = name_base+entry_num;
            double _ser_component = 0.;
            error  = db.GetNamedChannelData( channel, entry_name, _ser_component );
            if( error ) throw cet::exception( "PMTTimingCorrectionsProvider" ) << "Encountered error (code " << error << ") while trying to access 'ser_vec' on table " << dbname << "\n";
            _ser.push_back(_ser_component);
        }
        fPMTCalibrationData[channel].ser = _ser;
    }
}

// -----------------------------------------------------------------------------

/// Read all the corrections from the database and save them inside a map, whose index 
/// is the PMT channel number
void sbndDB::PMTCalibrationDatabaseProvider::readPMTCalibrationDatabase(const art::Run& run){

    // Clear before the run
    fPMTCalibrationData.clear();

    ReadPMTCalibration(run.id().run());

    if( fVerbose ) {
        mf::LogInfo(fLogCategory) << "Dump information from database " << std::endl;
        mf::LogVerbatim(fLogCategory) << "channel, trigger cable delay, reset cable delay, laser corrections, muons corrections" << std::endl;
        for( auto const & [key, value] : fPMTCalibrationData ){
            mf::LogVerbatim(fLogCategory) 
               << key << " " 
               << value.breakoutBox << "," 
               << std::endl; 
        }
    }
}
/**
 * @file   icaruscode/Timing/PMTTimingCorrectionsProvider.cxx
 * @brief  Service for the PMT timing corrections.
 * @author Andrea Scarpelli (ascarpell@bnl.gov), Matteo Vicenzi (mvicenzi@bnl.gov)
 */

// Framework includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Principal/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Local
#include "sbndcode/DatabaseInterface/PMTTimingCorrections.h"
#include "sbndcode/DatabaseInterface/PMTTimingCorrectionsProvider.h"

// Database interface helpers
#include "larevt/CalibrationDBI/Providers/DBFolder.h"
#include "larevt/CalibrationDBI/IOVData/TimeStampDecoder.h"

// C/C++ standard libraries
#include <string>
#include <vector>

//--------------------------------------------------------------------------------

sbndDB::PMTTimingCorrectionsProvider::PMTTimingCorrectionsProvider
    (const fhicl::ParameterSet& pset) 
    : fVerbose{ pset.get<bool>("Verbose", false) }
    , fLogCategory{ pset.get<std::string>("LogCategory", "PMTTimingCorrection") }
    { 
        fhicl::ParameterSet const tags{ pset.get<fhicl::ParameterSet>("CorrectionTags") };
        fCablesTag  = tags.get<std::string>("CablesTag");
        fLaserTag   = tags.get<std::string>("LaserTag");
        fCosmicsTag = tags.get<std::string>("CosmicsTag");
        if( fVerbose ) mf::LogInfo(fLogCategory) << "Database tags for timing corrections:\n"
						 << "Cables corrections  " << fCablesTag << "\n"  
						 << "Laser corrections   " << fLaserTag  << "\n"
						 << "Cosmics corrections " << fCosmicsTag;
    }

// -------------------------------------------------------------------------------

uint64_t sbndDB::PMTTimingCorrectionsProvider::RunToDatabaseTimestamp( uint32_t run ) const{

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
void sbndDB::PMTTimingCorrectionsProvider::ReadPMTCablesCorrections( uint32_t run ) { 

    // pmt_cables_delay: delays of the cables relative to trigger 
    // and reset distribution
    const std::string dbname("pmt_cables_delays_data");
    lariov::DBFolder db(dbname, "", "", fCablesTag, true, false);

    bool ret = db.UpdateData( RunToDatabaseTimestamp(run) ); // select table based on run number   
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
        
        // PPS reset correction
        double reset_distribution_delay = 0;
        int error  = db.GetNamedChannelData( channel, "reset_distribution_delay", reset_distribution_delay );
        if( error ) throw cet::exception("PMTTimingCorrectionsProvider") << "Encountered error (code " << error << ") while trying to access 'reset_distribution_delay' on table " << dbname << "\n";

        // Trigger cable delay
        double trigger_reference_delay = 0;
        error  = db.GetNamedChannelData( channel, "trigger_reference_delay", trigger_reference_delay );
        if( error ) throw cet::exception( "PMTTimingCorrectionsProvider" ) << "Encountered error (code " << error << ") while trying to access 'trigger_reference_delay' on table " << dbname << "\n";

        // Phase correction
        double phase_correction = 0;
	      error = db.GetNamedChannelData( channel, "phase_correction", phase_correction );
        if( error ) throw cet::exception( "PMTTimingCorrectionsProvider" ) << "Encountered error (code " << error <<  ") while trying to access 'phase_correction' on table " << dbname << "\n";
   
        /// This is the delay due to the cables connecting the 'global' trigger crate FPGA to the spare channel of the first digitizer in each VME crates. 
        /// The phase correction is an additional fudge factor 
        /// holding local variation due to constant clock offsets (it can be up to 8 ns)
        /// It can be absorbed within other corrections if necessary
        /// Corrections are saved in ns, but icaruscode wants us
        /// Correction are saved with the sing correspoding to their time direction 
        fDatabaseTimingCorrections[channel].triggerCableDelay = -(trigger_reference_delay-phase_correction)/1000.;

        /// This is the delay along the distribution line of the TTT reset
        /// The phase correction is an additional fudge factor 
        /// holding local variation due to constant clock offsets (it can be up to 8 ns)
        /// It can be absorbed within other corrections if necessary
        /// Corrections are saved in ns, but icaruscode wants us
        /// Corrections are additive! 
        fDatabaseTimingCorrections[channel].resetCableDelay = (reset_distribution_delay-phase_correction)/1000.; 
    }

}


// -----------------------------------------------------------------------------

/// Function to look up the calibration database at the table holding the pmt timing corrections measured using the laser
void sbndDB::PMTTimingCorrectionsProvider::ReadLaserCorrections( uint32_t run ) { 

    const std::string dbname("pmt_laser_timing_data");
    lariov::DBFolder db(dbname, "", "", fLaserTag, true, false);

    bool ret = db.UpdateData( RunToDatabaseTimestamp(run) ); // select table based on run number   
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
      throw cet::exception("PMTTimingCorrectionsProvider") << "got an empty channel list for run " << run << " in " << dbname << "\n";
    }

    for( auto channel : channelList ) {
        
        // Laser correction
        double t_signal = 0;
        int error  = db.GetNamedChannelData( channel, "t_signal", t_signal );
        if( error ) throw cet::exception( "PMTTimingCorrectionsProvider" ) << "Encountered error (code " << error << ") while trying to access 't_signal' on table " << dbname << "\n";

        /// pmt_laser_delay: delay from the Electron Transit time inside the PMT 
        /// and the PMT signal cable 
        /// corrections are saved in ns, but icaruscode wants us
        /// Corrections are additive! 
        fDatabaseTimingCorrections[channel].laserCableDelay = -t_signal/1000.; 
    }
}

// -----------------------------------------------------------------------------

/// Function to look up the calibration database at the table holding the pmt timing corrections measured using cosmic muons
void sbndDB::PMTTimingCorrectionsProvider::ReadCosmicsCorrections( uint32_t run ) { 

    const std::string dbname("pmt_cosmics_timing_data");
    lariov::DBFolder db(dbname, "", "", fCosmicsTag, true, false);

    bool ret = db.UpdateData( RunToDatabaseTimestamp(run) ); // select table based on run number   
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
        
        // Cosmics correction
 	double mean_residual_ns = 0;
	int error = db.GetNamedChannelData( channel, "mean_residual_ns", mean_residual_ns );
        if( error ) throw cet::exception( "PMTTimingCorrectionsProvider" ) << "Encountered error (code " << error << ") while trying to access 'mean_residual_ns' on table " << dbname << "\n";

        /// pmt_cosmics_residual: time residuals from downward going cosmics tracks 
        /// correcting for point-like laser emission and pmts that do not see laser light
        /// corrections are saved in ns, but icaruscode wants us
        /// Corrections are additive! 
        fDatabaseTimingCorrections[channel].cosmicsCorrections = -mean_residual_ns/1000.; 
    }

}


// -----------------------------------------------------------------------------


/// Read all the corrections from the database and save them inside a map, whose index 
/// is the PMT channel number
void sbndDB::PMTTimingCorrectionsProvider::readTimeCorrectionDatabase(const art::Run& run){

    // Clear before the run
    fDatabaseTimingCorrections.clear();

    ReadPMTCablesCorrections(run.id().run());
    ReadLaserCorrections(run.id().run());
    ReadCosmicsCorrections(run.id().run());

    if( fVerbose ) {

        mf::LogInfo(fLogCategory) << "Dump information from database " << std::endl;
        mf::LogVerbatim(fLogCategory) << "channel, trigger cable delay, reset cable delay, laser corrections, muons corrections" << std::endl;
        for( auto const & [key, value] : fDatabaseTimingCorrections ){
            mf::LogVerbatim(fLogCategory) 
               << key << " " 
               << value.triggerCableDelay << "," 
               << value.resetCableDelay << ", " 
               << value.laserCableDelay << ", "
               << value.cosmicsCorrections << ","
               << std::endl; 
        }
    }

}
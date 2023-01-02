/**
 * @file   sbndcode/Decoders/DecoderTools/Dumpers/PMTTimingCorrectionsProvider.cxx
 * @brief  Service for the PMT timing corrections.
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 */

// Framework includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Principal/Run.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "sbndcode/Decoders/DecoderTools/Dumpers/PMTTimingCorrections.h"
#include "sbndcode/Decoders/DecoderTools/Dumpers/PMTTimingCorrectionsProvider.h"

// Database interface helpers
#include "wda.h"

// C/C++ standard libraries
#include <memory> // std::unique_ptr<>
#include <optional>
#include <string>
#include <utility> // std::move()
#include <cassert>
#include <tuple>

//--------------------------------------------------------------------------------

sbndDB::PMTTimingCorrectionsProvider::PMTTimingCorrectionsProvider
    (const fhicl::ParameterSet& pset) 
    : fUrl{ pset.get<std::string>("DatabaseURL", "https://dbdata0vm.fnal.gov:9443/sbnd_con_prod/app/data?") }
    , fTimeout{ pset.get<unsigned int>("Timeout", 1000) }
    , fVerbose{ pset.get<bool>("Verbose", false) }
    , fLogCategory{ pset.get<std::string>("LogCategory", "PMTTimingCorrection") }
    , fClocksData{ art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob() }
    { }


// -------------------------------------------------------------------------------

/// Access the PostgreSQL calibration database for sbnd
/// Query depends on the run 
int sbndDB::PMTTimingCorrectionsProvider::ConnectToDataset(
    const std::string& name, uint32_t run, Dataset& dataset ) const {

    int error = 0;
    std::string url = fUrl + "f="+name +"&t="+std::to_string(run);
    dataset = getDataWithTimeout( url.c_str(), "", fTimeout, &error );
    if ( error ){
        throw cet::exception(fLogCategory) 
          << "Calibration database access failed. URL (" << url 
          << ") Error code: " << error;
    }
    if ( getHTTPstatus(dataset) !=200 ){
        throw cet::exception(fLogCategory)
            << "Calibration database access failed. URL (" << url
            << "). HTTP error status: " << getHTTPstatus(dataset) 
            << ". HTTP error message: " << getHTTPmessage(dataset); 
    }

    return error;

}


// -----------------------------------------------------------------------------

/// Function to look up the calibration database at the table holding the pmt hardware cables corrections
void sbndDB::PMTTimingCorrectionsProvider::ReadPMTCablesCorrections( uint32_t run ) { 

    // pmt_cables_delay: delays of the cables relative to trigger 
    // and reset distribution
    const std::string name("pmt_cables_delays_data");
    Dataset dataset;
    int error = ConnectToDataset( name, run, dataset );

    if (error) throw(std::exception());

    for( int row=4; row < getNtuples(dataset); row++ ) {

        Tuple tuple = getTuple(dataset, row);

        if( tuple != NULL ) {

            unsigned int channel_id = getLongValue( tuple, 0, &error );
            if( error ) throw std::runtime_error( "Encountered error while trying to access 'channel_id' on table " + name );
            
            // PPS reset correction
            double reset_distribution_delay = getDoubleValue( tuple, 11, &error );
            if( error ) throw std::runtime_error( "Encountered error '" + std::to_string(error) + "' while trying to access 'reset_distribution_delay' on table " + name );

            // Trigger cable delay
            double trigger_reference_delay = getDoubleValue( tuple, 10, &error );
            if( error ) throw std::runtime_error( "Encountered error while trying to access 'trigger_reference_delay' on table " + name );

            // Phase correction
            double phase_correction = getDoubleValue( tuple, 13, &error );
            if( error ) throw std::runtime_error( "Encountered error while trying to access 'phase_correction' on table " + name );
   
            /// This is the delay due to the cables connecting the 'global' FPGA of the trigger crate to the spare channel of the first digitizer in each VME crates. 
            ////The phase correction is an additional fudge factor 
            /// holding local variation due to constant clock offsets (it can be up to 8 ns)
            /// It can be absorbed within other corrections if necessary
            /// Corrections are saved in ns, but sbndcode wants us
            /// Correction are saved with the sing correspoding to their time direction 
            fDatabaseTimingCorrections[channel_id].triggerCableDelay = -(trigger_reference_delay-phase_correction)/1000.;

            /// This is the delay along the distribution line of the TTT reset
            /// The phase correction is an additional fudge factor 
            /// holding local variation due to constant clock offsets (it can be up to 8 ns)
            /// It can be absorbed within other corrections if necessary
            /// Corrections are saved in ns, but sbndcode wants us
            /// Corrections are additive! 
            fDatabaseTimingCorrections[channel_id].resetCableDelay = (reset_distribution_delay-phase_correction)/1000.; 

            releaseTuple(tuple);
        }

    }

}


// -----------------------------------------------------------------------------

/// Function to look up the calibration database at the table holding the pmt timing corrections measured using the laser
void sbndDB::PMTTimingCorrectionsProvider::ReadLaserCorrections( uint32_t run ) { 

    const std::string name("pmt_laser_timing_data");
    Dataset dataset;
    int error = ConnectToDataset( name, run, dataset );

    if (error) throw(std::exception());

    for( int row=4; row < getNtuples(dataset); row++ ) {

        Tuple tuple = getTuple(dataset, row);
        if( tuple != NULL ) {

            unsigned int channel_id = getLongValue( tuple, 0, &error );
            if( error ) throw std::runtime_error( "Encountered error while trying to access channel_id on table " + name );

            double t_signal = getDoubleValue( tuple, 5, &error );
            if( error ) throw std::runtime_error( "Encountered error while trying to access 't_signal' on table " + name );

            /// pmt_laser_delay: delay from the Electron Transit time inside the PMT 
            /// and the PMT signal cable 
            /// corrections are saved in ns, but sbndcode wants us
            /// Corrections are additive! 
            fDatabaseTimingCorrections[channel_id].laserCableDelay = -t_signal/1000.; 

            releaseTuple(tuple);
        }

    }

}

// -----------------------------------------------------------------------------

/*
/// Function to look up the calibration database at the table holding the pmt timing corrections measured using cosmic muons
void sbnd::PMTTimingCorrectionsProvider::ReadCosmicsCorrections( uint32_t run ) { 

    const std::string name("pmt_muons_timing_data");
    Dataset dataset;
    int error = ConnectToDataset( name, run, dataset );

    if (error) throw(std::exception());

    for( int row=4; row < getNtuples(dataset); row++ ) {

        Tuple tuple = getTuple(dataset, row);
        if( tuple != NULL ) {

            unsigned int channel_id = getLongValue( tuple, 0, &error );
            if( error ) throw std::runtime_error( "Encountered error while trying to access channel_id on table " + name );

            double t_signal = getDoubleValue( tuple, 5, &error );
            if( error ) throw std::runtime_error( "Encountered error while trying to access 't_signal' on table " + name );

            /// pmt_laser_delay: delay from the Electron Transit time inside the PMT 
            /// and the PMT signal cable 
            /// corrections are saved in ns, but sbndcode wants us
            /// Corrections are additive! 
            fDatabaseTimingCorrections[channel_id].cosmicsCorrections = t_signal/1000.; 

            releaseTuple(tuple);
        }

    }

}
*/


// -----------------------------------------------------------------------------


/// Read all the corrections from the database and save them inside a map, whose index 
/// is the PMT channel number
void sbndDB::PMTTimingCorrectionsProvider::readTimeCorrectionDatabase(const art::Run& run){

    // Clear before the run
    fDatabaseTimingCorrections.clear();

    ReadPMTCablesCorrections(run.id().run());

    ReadLaserCorrections(run.id().run());

    //ReadCosmicsCorrections(run.id().run());

    if( fVerbose ) {

        mf::LogInfo(fLogCategory) << "Dump information from database " << std::endl;
        mf::LogInfo(fLogCategory) << "channel, trigger cable delay, reset cable delay, laser corrections, muons correctsions" << std::endl;

        for( auto const & [key, value] : fDatabaseTimingCorrections ){
            mf::LogInfo(fLogCategory) << key << " " 
                  << value.triggerCableDelay << "," 
                  << value.resetCableDelay << ", " 
                  << value.laserCableDelay << ", "
                  << value.cosmicsCorrections << ","
                  << std::endl; 
        }
    }

}

#ifndef SIMPLEFLASHALGO_CXX
#define SIMPLEFLASHALGO_CXX

#include "SimpleFlashAlgo.h"
#include <set>
#include <algorithm>

namespace lightana{

    static SimpleFlashAlgoFactory __SimpleFlashAlgoFactoryStaticObject__;

    SimpleFlashAlgo::SimpleFlashAlgo(const std::string name)
    : FlashAlgoBase(name)
    {}

    size_t SimpleFlashAlgo::_nopdet_maxsize = 0;

    /*
    CONFIGURE:
    - Reads some configuration parameters
    - Makes sure the configuration is valid (the parameters make sense)
    */
    void SimpleFlashAlgo::Configure(const Config_t &p)
    {
        Reset();
        _debug          = p.get<bool>("DebugMode",false);       // Print couts debug info
        _min_pe_flash   = p.get<double>("PEThreshold",10);      // Minimum PE to declare a flash (20 in config file)
        _min_pe_coinc   = p.get<double>("MinPECoinc",  5);      // Minimum PE in coincidence window (6 in config file)
        _min_mult_coinc = p.get<double>("MinMultCoinc", 2);     // Minimum multiplicity in coincidence window (3 in config file)
        _integral_time = p.get<double>("IntegralTime",8);       // Default integration time for flash
        _pre_sample    = p.get<double>("PreSample",0.1);        // Time before the peak time to start integration
        _veto_time     = p.get<double>("VetoSize",8.);          // Veto time after a flash is found
        _time_res      = p.get<double>("TimeResolution",0.03);  // Time resolution = bin size (0.01(us) in config file)
        _tpc           = p.get<int>("TPC");                     // TPC to associate the flash with (0 or 1)
        //_pe_baseline_v.clear();
        //_pe_baseline_v = p.get<std::vector<double> >("PEBaseline",_pe_baseline_v);

        _min_pe_repeated = p.get<double>("MinPECoincRepeated", 20);              // Minimum PE in one bin to declare a repeated flash during an existing OpFlash
        _min_time_before = p.get<double>("MinTimeBefore", 0.5);             // minimum time separation to declare a repeated flash before an existing OpFlash
        _time_dif_flash_before = p.get<double>("TimeDifferenceFlashBefore", 0.05);             // minimum time separation to declare a repeated flash before an existing OpFlash

        // Check that integral_time > veto_time (they are set equal)
        if(_integral_time > _veto_time) {
            std::cerr << "Integral time cannot exceed veto time!" << std::endl;
            throw std::exception();
        }

        // Check that _min_pe_repeated >= _min_pe_coinc
        if(_min_pe_coinc > _min_pe_repeated) {
            std::cerr << "MinPERepeated cannot be smaller than MinPECoinc!" << std::endl;
            throw std::exception();
        }

        /*
        Configuration of GENERAL hit veto ranges. These veto windows are applied to ALL flashes.
        They are DIFFERENT from the veto windows that are applied after a flash is found.
        Right now, in the configuration file they are left empty (no general veto ranges)
        */
        auto const range_start_v = p.get<std::vector<double> >("HitVetoRangeStart");    // Veto ranges for hits
        auto const range_end_v   = p.get<std::vector<double> >("HitVetoRangeEnd");      // Veto ranges for hits
        // Check veto range config validity: start has the same length as end, start < end
        if(range_start_v.size() != range_end_v.size()) {
            std::cerr << "OpHit veto range start and end config value array have different length!" << std::endl;
            throw std::exception();
        }
        for(size_t i=0; i<range_start_v.size(); ++i) {
            if(range_start_v[i] >= range_end_v[i]) {
                std::cerr << "OpHit veto range element " << i
                << " have start time @ " << range_start_v[i]
                << " being later than end time @ " << range_end_v[i]
                << std::endl;
                throw std::exception();
            }
            _flash_veto_range_m.emplace(range_end_v[i],range_start_v[i]);
        }

        // Check what PD to use: ["pmt_coated", "pmt_uncoated", "xarapuca", "xarapuca_vuv", "xarapuca_vis"]
        // In the configuration file, only PMTs are included
        std::vector<std::string> pd_to_use;
        pd_to_use = p.get<std::vector<std::string>>("PD", pd_to_use);
        // Transform the list of strings to the ID of the PMTs
        std::vector<int> opch_to_use = PDNamesToList(pd_to_use);

        // opch_to_use --> channels that are allowed to be used
        // _index_to_opch_v --> channels that will actually be used, in the order of index
        _index_to_opch_v.clear();
        // Try to read the OpChannel array first
        _index_to_opch_v =  p.get<std::vector<int> >("OpChannel",_index_to_opch_v);
        // If not given (_index_to_opch_v is empty), try to read OpChannelRange or TPC
        if(_index_to_opch_v.empty()) {
            // Get the channels from this TPC --> USED CURRENTLY
            if(_tpc>=0) {
                auto const opch_v = ListOpChannelsByTPC(_tpc);
                _index_to_opch_v.reserve(opch_v.size());
                for(auto const& v : opch_v) {
                    auto iter = std::find(opch_to_use.begin(), opch_to_use.end(), v);
                    if (iter != opch_to_use.end()) {
                        _index_to_opch_v.push_back(v);
                        // std::cout << "Going to use ch " << v << std::endl;
                    }
                }

            // Get the channels from OpChannelRange: all channels that are in the range and also in opch_to_use. NOT USED CURRENTLY
            }else{
                auto opch_range = p.get<std::vector<int> >("OpChannelRange");
                // Check that the range is valid
                if(opch_range.size()!=2) {
                    std::cerr << "OpChannelRange must be a vector of length two!" << std::endl;
                    throw std::exception();
                }else if(opch_range[0] > opch_range[1]) {
                    std::cerr << "OpChannelRange 0th element (" << opch_range[0]
                    << ") must be larger than the 1st element (" << opch_range[1]
                    << ")" << std::endl;
                    throw std::exception();
                }
                // Fill _index_to_opch_v: look for the channels in the range that are also in opch_to_use
                _index_to_opch_v.reserve(opch_range[1]-opch_range[0]+1);
                for(int i=opch_range[0]; i<=opch_range[1]; ++i) {
                    auto iter = std::find(opch_to_use.begin(), opch_to_use.end(), i);
                    if (iter != opch_to_use.end()) {
                        _index_to_opch_v.push_back(i);
                    }
                }
            }
        }

        /* Now build the _opch_to_index_v mapping:
        Creates an internal index of the channels
        Example:
        IDs of the channels to use: [2, 3, 7, 10, ...]
        Internal new index (_index_to_opch_v): [-1, -1, 0, 1, -1, -1, -1, 2, -1, -1, 3, ...]
        -1: unused
        0, 1, 2, 3, ...: used channels such that _index_to_opch_v[2] = 0, _index_to_opch_v[3] = 1, _index_to_opch_v[7] = 2, _index_to_opch_v[10] = 3, ...
        */
        size_t valid_id=0;
        std::set<size_t> duplicate;
        for(auto const& ch : _index_to_opch_v) {
            if(ch >= (int)(_opch_to_index_v.size())) _opch_to_index_v.resize(ch+1,-1);
            // Check there are no duplicated events
            if(duplicate.find(ch) != duplicate.end()) {
                std::cerr << "Channel number " << ch << " is duplicated!" << std::endl;
                throw std::exception();
            }
            _opch_to_index_v[ch] = valid_id;
            valid_id += 1;
            duplicate.insert(ch);
        }

        if(_opch_to_index_v.empty()) {
            std::cerr << "Length of OpChannel array parameter is 0..." << std::endl;
            throw std::exception();
        }
        // std::cout << "Number of opch used to construct flashes: " << _index_to_opch_v.size() << std::endl;
        /*
        if(_pe_baseline_v.size() != duplicate.size()) {
          std::cout << "PEBaseline array length (" << _pe_baseline_v.size()
          << ") is not same as OpDet ID count (" << duplicate.size() << ")!" << std::endl;
          throw std::exception();
        }
        */

        if(_index_to_opch_v.size()>_nopdet_maxsize) _nopdet_maxsize = _index_to_opch_v.size();

    }

    // Method to check whether a hit at time t is vetoed by the GENERAL VETO RANGES
    // Different from the veto windows applied after a flash is found
    bool SimpleFlashAlgo::Veto(double t) const
    {
        auto iter = _flash_veto_range_m.lower_bound(t);
        if(iter == _flash_veto_range_m.end()) return false;
        return (t >= (*iter).second);
    }

    /*
    We end up with:
    - _min_pe_flash: minimum PE to declare a flash: TOTAL OR PER PMT?
    - _min_pe_coinc: minimum PE in coincidence window to declare a flash candidate: TOTAL OR PER PMT?
    - _min_mult_coinc: minimum multiplicity in coincidence window to declare a flash candidate
    - _integral_time: integration time for flash
    - _pre_sample: time before the peak time to start integration
    - _veto_time: veto time after a flash is found
    - _time_res: time resolution = bin size
    - _tpc: TPC to associate the flash with
    - **_index_to_opch_v**: mapping from internal index to OpChannel ID: array with the allowed channels [2, 4, 7, ...]
    - **_opch_to_index_v**: mapping from OpChannel ID to internal index: position (index) in the array = OpChannel ID [-1, -1, 0, -1, 1, -1, -1, 2, ...]
    - _flash_veto_range_m: map of general veto ranges for hits. NOW EMPTY
    */

    SimpleFlashAlgo::~SimpleFlashAlgo()
    {}

    LiteOpFlashArray_t SimpleFlashAlgo::RecoFlash(const LiteOpHitArray_t ophits) {

        Reset();        // Reset internal variables of the SimpleFlashAlgo class
        size_t max_ch = _opch_to_index_v.size() - 1;    // Maximum channel number (ID) to consider
        size_t NOpDet = _index_to_opch_v.size();        // Number of OpDet being used

        //static std::vector<double> pesum_v;
        static std::vector<double> mult_v;                          //< this is not strictly a multiplicity of PMTs, but multiplicity of OpHits in each bin
        static std::vector<std::vector<double> > pespec_v;          // per-channel PE spectrum per bin
        static std::vector<std::vector<unsigned int> > hitidx_v;    // indices of OpHits contributing to each time-bin
        // Find the time span of the hits
        double min_time=1.1e20;
        double max_time=1.1e20;
        for(auto const& oph : ophits) {
            if(max_time > 1.e20 || oph.peak_time > max_time) max_time = oph.peak_time;
            if(min_time > 1.e20 || oph.peak_time < min_time) min_time = oph.peak_time;
        }
        // Add some margin to the left and right
        min_time -= 10* _time_res;
        max_time += 10* _time_res;
        if(_debug)
            std::cout << "T span: " << min_time << " => " << max_time << " ... " << (size_t)((max_time - min_time) / _time_res) << std::endl;

        // calculate the number of bins needed of size _time_res = 10ns
        size_t nbins_pesum_v = (size_t)((max_time - min_time) / _time_res) + 1;
        // resize static vectors to size = number of bins. Fill them with value 0
        if(_pesum_v.size() < nbins_pesum_v) _pesum_v.resize(nbins_pesum_v,0);
        if(mult_v.size()   < nbins_pesum_v) mult_v.resize(nbins_pesum_v,0);
        if(pespec_v.size() < nbins_pesum_v) pespec_v.resize(nbins_pesum_v,std::vector<double>(_nopdet_maxsize));
        if(hitidx_v.size() < nbins_pesum_v) hitidx_v.resize(nbins_pesum_v,std::vector<unsigned int>());
        // reset pe_sum_v, set to 0
        for(size_t i=0; i<_pesum_v.size(); ++i) {
            _pesum_v[i] = 0;
        }
        // reset static vectors, set to 0
        for(size_t i=0; i<mult_v.size(); ++i) {
            mult_v[i]  = 0;
            hitidx_v[i].clear();
            for(auto& v : pespec_v[i]) v=0;
        }

        // Fill _pesum_v: number of PE in each time bin
        // Loop over OpHits
        for(size_t hitidx = 0; hitidx < ophits.size(); ++hitidx) {
            // Get the OpHit
            auto const& oph = ophits[hitidx];
            // Check that the channel is valid
            if(oph.channel > max_ch || _opch_to_index_v[oph.channel] < 0) {
                if(_debug) std::cout << "Ignoring OpChannel " << oph.channel << std::endl;
                continue;
            }
            // Check if the hit is vetoed by the GENERAL veto ranges
            if(Veto(oph.peak_time)) {
                if(_debug) std::cout << "Ignoring hit @ time " << oph.peak_time << std::endl;
                continue;
            }
            // Find the time bin index
            size_t index = (size_t)((oph.peak_time - min_time) / _time_res);
            // std::cout << "Ophit from ch " << oph.channel << " at time " << oph.peak_time << " with PE " << oph.pe << ", index " << index << std::endl;
            // Fill _pesum_v, mult_v, pespec_v, hitidx_v
            _pesum_v[index] += oph.pe;
            mult_v[index] += 1;
            pespec_v[index][_opch_to_index_v[oph.channel]] += oph.pe;
            hitidx_v[index].push_back(hitidx);
        }


        // Order by pe (above threshold)
        // To do that, we create a MAP of (1/pe_sum) --> index of the bin
        // We do 1/pe_sum so that the map is ordered from highest to lowest pe_sum
        std::map<double,size_t> pesum_idx_map;  // map: 1/pe_sum --> index of then bin
        // Loop over all time bins (already filled)
        for(size_t idx=0; idx<nbins_pesum_v; ++idx) {
            // std::cout <<  "    _pesum_v at " << idx << " is " << _pesum_v[idx] << ", _min_pe_coinc is " << _min_pe_coinc << std::endl;
            // Check if the bin is above threshold
            if(_pesum_v[idx] < _min_pe_coinc   ) continue;  // at least _min_pe_coinc PE = 6PE in the bin
            // std::cout <<  "    mult_v at " << idx << " is " << mult_v[idx] << ", _min_mult_coinc is " << _min_mult_coinc << std::endl;
            if(mult_v[idx]  < _min_mult_coinc ) continue;   // at least _min_mult_coinc = 3 OpHits in the bin
            // Add to the map
            pesum_idx_map[1./(_pesum_v[idx])] = idx;
        }

        // Get candidate flash times
        std::vector<std::pair<size_t,size_t> > flash_period_v;  // stores the time of the already defined OpFlashes (in bins) = their integration period (in bins)
        std::vector<size_t> flash_time_v;                       // stores the 'peak' time (index) of the flash
        size_t veto_ctr = (size_t)(_veto_time / _time_res);     // number of bins of the veto window of the OpFlashes
        size_t default_integral_ctr = (size_t)(_integral_time / _time_res); // number of bins of the integration window of the OpFlashes
        size_t precount = (size_t)(_pre_sample / _time_res);    // number of bins to go back from the peak time to start the integration
        // Reserve space in the vectors
        flash_period_v.reserve(pesum_idx_map.size());          
        flash_time_v.reserve(pesum_idx_map.size());

        double sum_baseline = 0;    // Total baseline to add to the PE threshold (=0 --> NO baseline)
        //for(auto const& v : _pe_baseline_v) sum_baseline += v;

        // Loop over the candidate flash times (ordered by PE): first = PE; second = index of the bin
        for(auto const& pe_idx : pesum_idx_map) {

            auto const& pe  = 1./(pe_idx.first);

            // Get the index of the bin
            auto const& idx = pe_idx.second;

            // Determine the start time of the integration window
            size_t start_time = idx;
            // Go back 'precount' bins
            if(start_time < precount) start_time = 0;
            else start_time = idx - precount;

            // see if this idx can be used
            bool skip=false;         // flag to skip this candidate flash
            size_t integral_ctr = default_integral_ctr; // For most OpFlashes, integral time = default_integral_ctr
            // Loop over the already defined flashes (flash_period_v) to see if there is overlap
            for(auto const& used_period : flash_period_v) {
                // If a flash candidate is inside an existing flash, veto it
                if( used_period.first <= start_time && start_time < (used_period.first + veto_ctr) ) {
                    skip=true;
                    break;
                }

                // If there is an OpFlash candidate that starts before an existing one,
                // and ends after the start of the existing one 
                // create a shortened integration window and a new flash
                if( (start_time <= used_period.first && (start_time + veto_ctr) > used_period.first) && (pe >= _min_pe_repeated) ) {
                    // To avoid super small integration windows (<0.5us), check that the shortened window is at least _min_time_before long
                    if( (start_time + (size_t)(_min_time_before/_time_res)) <= used_period.first ) {
                        if(_debug) 
                            std::cout << "Possible shortened OpFlash before existing OpFlash at time " << min_time + start_time * _time_res
                            << " (previous flash @ " << min_time + used_period.first*_time_res << " until " << min_time + (used_period.first + used_period.second)*_time_res
                            <<") "<< std::endl;
                        integral_ctr = used_period.first - start_time - (size_t)(_time_dif_flash_before/_time_res);  // _time_dif_flash_before/_time_res = number of bins before the existing flash to stop the integration = 0.05us before the existing flash
                    // If the shortened window would be too small, skip this candidate
                    }else{
                        skip=true;
                        break;
                    }
                }else if( (start_time <= used_period.first && (start_time + veto_ctr) > used_period.first) && (pe < _min_pe_repeated) ){
                    skip=true;
                    break;
                }

                /*
                // Same condition as before, but using the integration window instead of the veto window
                // As they are currently set equal, this condition never happens
                if( used_period.first >= start_time && used_period.first < (start_time + integral_ctr) ) {
                    if(_debug) std::cout << "Truncating flash @ " << start_time
                        << " (previous flash @ " << used_period.first
                        << ") ... integral ctr change: " << integral_ctr
                        << " => " << used_period.first - start_time << std::endl;

                    integral_ctr = used_period.first - start_time;
                }
                */
                if(_debug) {
                    std::cout << "Flash @ " << min_time + used_period.first * _time_res
                    << " => " << min_time + (used_period.first + used_period.second) * _time_res
                    << " does not interfare with THIS flash @ " << min_time + start_time * _time_res
                    << " => " << min_time + (start_time + integral_ctr) * _time_res << std::endl;
                }
            }
            if(skip) {
                if(_debug){
                    std::cout << "Skipping a candidate @ " << min_time + start_time * _time_res << " as it is in a veto window!" <<std::endl;
                } 
                continue;
            }
            // Check the thresholds in PE yo see if this candidate can be an OpFlash:
            double pesum = 0;   // Total PE in the integration window
            // Integrate from start_time to start_time + integral_ctr and calulate the total PE in pesum
            for(size_t i=start_time; i<std::min(nbins_pesum_v,(start_time+integral_ctr)); ++i)
                pesum += _pesum_v[i];

            // Check that the total PE is larger than the threshold
            if(pesum < (_min_pe_flash + sum_baseline)) {
                if(_debug) std::cout << "Skipping a candidate @ " << start_time  << " => " << start_time + integral_ctr
                    << " as it got " << pesum
                    << " PE which is lower than threshold " << (_min_pe_flash + sum_baseline) << std::endl;
                continue;
            }

            if(_debug) std::cout << "Claiming a flash candidate @ " << min_time + start_time * _time_res
                << " => " << min_time + (start_time + integral_ctr) * _time_res
                << " with PE = " << pesum
                << std::endl;

            // Fill the vector with the existing flashes
            flash_period_v.push_back(std::pair<size_t,size_t>(start_time,integral_ctr));
            flash_time_v.push_back(idx);
        }

        // Construct flash objects
        LiteOpFlashArray_t res;
        // Loop over the found flashes (their time)
        for(size_t flash_idx=0; flash_idx<flash_period_v.size(); ++flash_idx) {

            auto const& start  = flash_period_v[flash_idx].first;   // start of the OpFlash
            auto const& period = flash_period_v[flash_idx].second;  // period (in bins) of the OpFlash
            auto const& time   = flash_time_v[flash_idx];           // 'peak' time (index) of the OpFlash

            // Build the PE spectrum (PE in each PMT) for this flash
            std::vector<double> pe_v(max_ch+1,0);   // initialize to 0
            // Sum the PE in the integration window for each PMT
            // loop over time bins inside the integration window
            for(size_t index=start; index<(start+period) && index<pespec_v.size(); ++index) {
                // Loop over the PMTs
                for(size_t pmt_index=0; pmt_index<NOpDet; ++pmt_index)
                    // Sum the PE in each PMT
                    pe_v[_index_to_opch_v[pmt_index]] += pespec_v[index][pmt_index];
            }
            // Set negative PE to 0
            // Loop over OpChannels
            for(size_t opch=0; opch<max_ch; ++opch) {
                // Skip unused channels
                if(_opch_to_index_v[opch]<0) continue;

                //pe_v[opch] -= _pe_baseline_v[_opch_to_index_v[opch]];

                if(pe_v[opch]<0) pe_v[opch]=0;

            }

            // Get the associated OpHits for this flash
            std::vector<unsigned int> asshit_v;
            // Loop over time bins inside the integration window
            for(size_t index=start; index<(start+period) && index<pespec_v.size(); ++index) {
                // Loop over the OpHits in this time bin
                for(auto const& idx : hitidx_v[index])
                    // Add the OpHit index to the associated OpHits vector
                    asshit_v.push_back(idx);
            }

            if(_debug) {
                std::cout << "Claiming a flash @ " << min_time + time * _time_res << " => " << min_time + (time + period) * _time_res
                << " : " << std::flush;
                double tmpsum=0;
                for(auto const& v : pe_v) { std::cout << v << " "; tmpsum +=v; }
                std::cout << " ... sum = " << tmpsum << std::endl;
            }

            // Define the flash:
            /*
            - time: min_time + time * _time_res --> time of the start of the flash 
            - time width: period * _time_res / 2. = half of the integration window
            - TPC: _tpc --> TPC where the interactio occurs
            - PE spectrum: pe_v --> vector with the PE in each PMT
            - Associated OpHits: asshit_v --> vector with the indices of the associated OpHits
            */
            LiteOpFlash_t flash( min_time + time * _time_res,
                                period * _time_res / 2.,
                                _tpc,
                                std::move(pe_v),
                                std::move(asshit_v));
            res.emplace_back( std::move(flash) );

        }
        if(_debug) std::cout << std::endl;
        return res;
    }

}
#endif

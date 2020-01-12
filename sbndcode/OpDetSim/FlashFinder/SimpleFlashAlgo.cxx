#ifndef SIMPLEFLASHALGO_CXX
#define SIMPLEFLASHALGO_CXX

#include "SimpleFlashAlgo.h"
#include <set>

namespace lightana{
    
    static SimpleFlashAlgoFactory __SimpleFlashAlgoFactoryStaticObject__;
    
    SimpleFlashAlgo::SimpleFlashAlgo(const std::string name)
    : FlashAlgoBase(name)
    {}
    
    void SimpleFlashAlgo::Configure(const Config_t &p)
    {
        Reset();
        _debug          = p.get<bool>("DebugMode",false);
        _min_pe_flash   = p.get<double>("PEThreshold",10);
        _min_pe_coinc   = p.get<double>("MinPECoinc",  5);
        _min_mult_coinc = p.get<double>("MinMultCoinc", 2);
        _integral_time = p.get<double>("IntegralTime",8);
        _pre_sample    = p.get<double>("PreSample",0.1);
        _veto_time     = p.get<double>("VetoSize",8.);
        _time_res      = p.get<double>("TimeResolution",0.03);
        //_pe_baseline_v.clear();
        //_pe_baseline_v = p.get<std::vector<double> >("PEBaseline",_pe_baseline_v);
        
        if(_integral_time > _veto_time) {
            std::cerr << "Integral time cannot exceed veto time!" << std::endl;
            throw std::exception();
        }
        auto const range_start_v = p.get<std::vector<double> >("HitVetoRangeStart");
        auto const range_end_v   = p.get<std::vector<double> >("HitVetoRangeEnd");
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
        
        _index_to_opch_v.clear();
        _index_to_opch_v =  p.get<std::vector<int> >("OpChannel",_index_to_opch_v);
        if(_index_to_opch_v.empty()) {
            int cryostat = p.get<int>("Cryostat",-1);
            if(cryostat>=0) {
                auto const opch_v = ListOpChannels(cryostat);
                _index_to_opch_v.reserve(opch_v.size());
                for(auto const& v : opch_v) _index_to_opch_v.push_back(v);
            }else{
                auto opch_range = p.get<std::vector<int> >("OpChannelRange");
                if(opch_range.size()!=2) {
                    std::cerr << "OpChannelRange must be a vector of length two!" << std::endl;
                    throw std::exception();
                }else if(opch_range[0] > opch_range[1]) {
                    std::cerr << "OpChannelRange 0th element (" << opch_range[0]
                    << ") must be larger than the 1st element (" << opch_range[1]
                    << ")" << std::endl;
                    throw std::exception();
                }
                _index_to_opch_v.reserve(opch_range[1]-opch_range[0]+1);
                for(int i=opch_range[0]; i<=opch_range[1]; ++i)
                    _index_to_opch_v.push_back(i);
            }
        }
        size_t valid_id=0;
        std::set<size_t> duplicate;
        for(auto const& ch : _index_to_opch_v) {
            if(ch >= (int)(_opch_to_index_v.size())) _opch_to_index_v.resize(ch+1,-1);
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
        /*
        if(_pe_baseline_v.size() != duplicate.size()) {
          std::cout << "PEBaseline array length (" << _pe_baseline_v.size()
          << ") is not same as OpDet ID count (" << duplicate.size() << ")!" << std::endl;
          throw std::exception();
        }
        */
    }
    
    bool SimpleFlashAlgo::Veto(double t) const
    {
        auto iter = _flash_veto_range_m.lower_bound(t);
        if(iter == _flash_veto_range_m.end()) return false;
        return (t >= (*iter).second);
    }
    
    SimpleFlashAlgo::~SimpleFlashAlgo()
    {}
    
    LiteOpFlashArray_t SimpleFlashAlgo::RecoFlash(const LiteOpHitArray_t ophits) {
        
        Reset();
        size_t max_ch = _opch_to_index_v.size() - 1;
        size_t NOpDet = _index_to_opch_v.size();
        
        //static std::vector<double> pesum_v;
        static std::vector<double> mult_v;  //< this is not strictly a multiplicity of PMTs, but multiplicity of hits
        static std::vector<std::vector<double> > pespec_v;
        static std::vector<std::vector<unsigned int> > hitidx_v;
        double min_time=1.1e20;
        double max_time=1.1e20;
        for(auto const& oph : ophits) {
            if(max_time > 1.e20 || oph.peak_time > max_time) max_time = oph.peak_time;
            if(min_time > 1.e20 || oph.peak_time < min_time) min_time = oph.peak_time;
        }
        min_time -= 10* _time_res;
        max_time += 10* _time_res;
        if(_debug)
            std::cout << "T span: " << min_time << " => " << max_time << " ... " << (size_t)((max_time - min_time) / _time_res) << std::endl;
        
        size_t nbins_pesum_v = (size_t)((max_time - min_time) / _time_res) + 1;
        if(_pesum_v.size() < nbins_pesum_v) _pesum_v.resize(nbins_pesum_v,0);
        if(mult_v.size()   < nbins_pesum_v) mult_v.resize(nbins_pesum_v,0);
        if(pespec_v.size() < nbins_pesum_v) pespec_v.resize(nbins_pesum_v,std::vector<double>(NOpDet));
        if(hitidx_v.size() < nbins_pesum_v) hitidx_v.resize(nbins_pesum_v,std::vector<unsigned int>());
        for(size_t i=0; i<_pesum_v.size(); ++i) {
            _pesum_v[i] = 0;
            mult_v[i]  = 0;
            hitidx_v[i].clear();
            for(auto& v : pespec_v[i]) v=0;
        }
        
        // Fill _pesum_v
        for(size_t hitidx = 0; hitidx < ophits.size(); ++hitidx) {
            auto const& oph = ophits[hitidx];
            if(oph.channel > max_ch || _opch_to_index_v[oph.channel] < 0) {
                if(_debug) std::cout << "Ignoring OpChannel " << oph.channel << std::endl;
                continue;
            }
            if(Veto(oph.peak_time)) {
                if(_debug) std::cout << "Ignoring hit @ time " << oph.peak_time << std::endl;
                continue;
            }
            size_t index = (size_t)((oph.peak_time - min_time) / _time_res);
            _pesum_v[index] += oph.pe;
            mult_v[index] += 1;
            pespec_v[index][_opch_to_index_v[oph.channel]] += oph.pe;
            hitidx_v[index].push_back(hitidx);
        }
        
        // Order by pe (above threshold)
        std::map<double,size_t> pesum_idx_map;
        for(size_t idx=0; idx<nbins_pesum_v; ++idx) {
            if(_pesum_v[idx] < _min_pe_coinc   ) continue;
            if(mult_v[idx]  < _min_mult_coinc ) continue;
            pesum_idx_map[1./(_pesum_v[idx])] = idx;
        }
        
        // Get candidate flash times
        std::vector<std::pair<size_t,size_t> > flash_period_v;
        std::vector<size_t> flash_time_v;
        size_t veto_ctr = (size_t)(_veto_time / _time_res);
        size_t default_integral_ctr = (size_t)(_integral_time / _time_res);
        size_t precount = (size_t)(_pre_sample / _time_res);
        flash_period_v.reserve(pesum_idx_map.size());
        flash_time_v.reserve(pesum_idx_map.size());
        
        double sum_baseline = 0;
        //for(auto const& v : _pe_baseline_v) sum_baseline += v;

        for(auto const& pe_idx : pesum_idx_map) {
                        
          //auto const& pe  = 1./(pe_idx.first);
            auto const& idx = pe_idx.second;
            
            size_t start_time = idx;
            if(start_time < precount) start_time = 0;
            else start_time = idx - precount;
            
            // see if this idx can be used
            bool skip=false;
            size_t integral_ctr = default_integral_ctr;
            for(auto const& used_period : flash_period_v) {
                if( start_time <= used_period.first && (start_time + veto_ctr) > used_period.first ) {
                    skip=true;
                    break;
                }
                if( used_period.first <= start_time && start_time < (used_period.first + veto_ctr) ) {
                    skip=true;
                    break;
                }
                if( used_period.first >= start_time && used_period.first < (start_time + integral_ctr) ) {
                    if(_debug) std::cout << "Truncating flash @ " << start_time
                        << " (previous flash @ " << used_period.first
                        << ") ... integral ctr change: " << integral_ctr
                        << " => " << used_period.first - start_time << std::endl;
                    
                    integral_ctr = used_period.first - start_time;
                }
                if(_debug) {
                    std::cout << "Flash @ " << min_time + used_period.first * _time_res
                    << " => " << min_time + (used_period.first + used_period.second) * _time_res
                    << " does not interfare with THIS flash @ " << min_time + start_time * _time_res
                    << " => " << min_time + (start_time + integral_ctr) * _time_res << std::endl;
                }
            }
            if(skip) {
                if(_debug) std::cout << "Skipping a candidate @ " << min_time + start_time * _time_res << " as it is in a veto window!" <<std::endl;
                continue;
            }
            
            // See if this flash is declarable
            double pesum = 0;
            for(size_t i=start_time; i<std::min(nbins_pesum_v,(start_time+integral_ctr)); ++i)
                
                pesum += _pesum_v[i];
            
            if(pesum < (_min_pe_flash + sum_baseline)) {
                if(_debug) std::cout << "Skipping a candidate @ " << start_time  << " => " << start_time + integral_ctr
                    << " as it got " << pesum
                    << " PE which is lower than threshold " << (_min_pe_flash + sum_baseline) << std::endl;
                continue;
            }
            
            flash_period_v.push_back(std::pair<size_t,size_t>(start_time,integral_ctr));
            flash_time_v.push_back(idx);
        }
        
        // Construct flash
        LiteOpFlashArray_t res;
        for(size_t flash_idx=0; flash_idx<flash_period_v.size(); ++flash_idx) {
            
            auto const& start  = flash_period_v[flash_idx].first;
            auto const& period = flash_period_v[flash_idx].second;
            auto const& time   = flash_time_v[flash_idx];
            
            std::vector<double> pe_v(max_ch+1,0);
            for(size_t index=start; index<(start+period) && index<pespec_v.size(); ++index) {
                
                for(size_t pmt_index=0; pmt_index<NOpDet; ++pmt_index)
                    
                    pe_v[_index_to_opch_v[pmt_index]] += pespec_v[index][pmt_index];
                
            }
            
            for(size_t opch=0; opch<max_ch; ++opch) {
                
                if(_opch_to_index_v[opch]<0) continue;
                
                //pe_v[opch] -= _pe_baseline_v[_opch_to_index_v[opch]];
                
                if(pe_v[opch]<0) pe_v[opch]=0;
                
            }
            
            std::vector<unsigned int> asshit_v;
            for(size_t index=start; index<(start+period) && index<pespec_v.size(); ++index) {
                for(auto const& idx : hitidx_v[index])
                    asshit_v.push_back(idx);
            }
            
            if(_debug) {
                std::cout << "Claiming a flash @ " << min_time + time * _time_res
                << " : " << std::flush;
                double tmpsum=0;
                for(auto const& v : pe_v) { std::cout << v << " "; tmpsum +=v; }
                std::cout << " ... sum = " << tmpsum << std::endl;
            }
            
            LiteOpFlash_t flash( min_time + time * _time_res,
                                period * _time_res / 2.,
                                std::move(pe_v),
                                std::move(asshit_v));
            res.emplace_back( std::move(flash) );
            
        }
        if(_debug) std::cout << std::endl;
        return res;
    }
    
}
#endif



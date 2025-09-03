#ifndef PIXELMAPVARS_H
#define PIXELMAPVARS_H

#include <vector>
#include <map>

struct PixelMapVars {
    std::vector<std::vector<float>> flash_ophit_pe;
    std::vector<std::vector<int>>   flash_ophit_ch;
    std::vector<std::vector<float>> flash_ophit_time;
    std::vector<double> nuvT;
    std::vector<double> dEpromx;
    std::vector<double> dEpromy;
    std::vector<double> dEpromz;
    std::vector<double> dEtpc;
    std::vector<double> nuvZ;
    
    // TensorFlow predictions for dEprom* (position reconstruction)
    std::vector<double> dEpromx_pred;
    std::vector<double> dEpromy_pred;
    std::vector<double> dEpromz_pred;
    
    // Differences: prediction - ground truth
    std::vector<double> dEpromx_diff;
    std::vector<double> dEpromy_diff;
    std::vector<double> dEpromz_diff;
    
    // Event information and performance metrics
    int run_id;
    int subrun_id; 
    int event_id;
    bool passed_filters;
    std::vector<double> error_3d; // 3D distance error for each reconstructed event
    
    // Channel dictionary mapping OpDetID to OpDetType
    std::map<int, int> channel_dict;
    
    // PE matrix [event][channel(312)]
    std::vector<std::vector<float>> pe_matrix;
    
    // PMT maps [ch_y][ch_z] 
    std::vector<std::vector<int>> coated_pmt_map;
    std::vector<std::vector<int>> uncoated_pmt_map;
    
    // Generated images [event][ch_y/2][ch_z][map_count]
    std::vector<std::vector<std::vector<std::vector<float>>>> pe_images;
    
    // Default constructor required by ROOT dictionary
    PixelMapVars() = default;
};

#endif // PIXELMAPVARS_H
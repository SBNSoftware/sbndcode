#ifndef PIXELMAPVARS_H
#define PIXELMAPVARS_H

#include <vector>

struct PixelMapVars {
    std::vector<std::vector<double>> flash_ophit_pe;
    std::vector<std::vector<int>>   flash_ophit_ch;
    std::vector<std::vector<double>> flash_ophit_time;
    std::vector<double> nuvT;
    std::vector<double> dEpromx;
    std::vector<double> dEpromy;
    std::vector<double> dEpromz;
    std::vector<double> dEtpc;
    std::vector<double> nuvZ;
    
    // Default constructor required by ROOT dictionary
    PixelMapVars() = default;
};

#endif // PIXELMAPVARS_H
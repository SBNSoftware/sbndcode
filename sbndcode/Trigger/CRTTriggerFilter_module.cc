////////////////////////////////////////////////////////////////////////
// Class:       MetricFilter
// Plugin Type: filter (Unknown Unknown)
// File:        MetricFilter_module.cc
//
// Generated at Fri Dec 13 06:19:04 2024 by Lynn Tung using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "sbndcode/Decoders/PTB/sbndptb.h"

#include <memory>

std::vector<int> bitsUP(int trigger_word)
{
    std::string binary = "";

    // Convierte el número a binario
    while (trigger_word > 0)
    {
        binary = std::to_string(trigger_word % 2) + binary;
        trigger_word /= 2;
    }

    // Si el número es 0, el binario es "0"
    if (binary.empty())
    {
        binary = "0";
    }

    // Encuentra las posiciones de los bits que están arriba (con valor 1)
    std::vector<int> bits_up;

    // Recorre la cadena de binario desde la derecha
    for (int i = binary.size() - 1; i >= 0; --i)
    {
        if (binary[i] == '1')
        {
            bits_up.push_back(binary.size() - 1 - i);
        }
    }

    return bits_up;
}

namespace sbndaq
{
    class CRTTriggerFilter;
}

class sbndaq::CRTTriggerFilter : public art::EDFilter
{
public:
    explicit CRTTriggerFilter(fhicl::ParameterSet const &p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    CRTTriggerFilter(CRTTriggerFilter const &) = delete;
    CRTTriggerFilter(CRTTriggerFilter &&) = delete;
    CRTTriggerFilter &operator=(CRTTriggerFilter const &) = delete;
    CRTTriggerFilter &operator=(CRTTriggerFilter &&) = delete;

    // Required functions.
    bool filter(art::Event &e) override;

private:
    std::string fPTBLabel;
    std::string fCRTTriggerType;
    int FilterHLT;
    // Declare member data here.
};

sbndaq::CRTTriggerFilter::CRTTriggerFilter(fhicl::ParameterSet const &p)
    : EDFilter{p} // ,
                  // More initializers here.
{
    fPTBLabel = p.get<std::string>("PTBLabel");
    fCRTTriggerType = p.get<std::string>("CRTTriggerType");
    if(fCRTTriggerType=="EW" || fCRTTriggerType=="WE") FilterHLT=15;
    else if(fCRTTriggerType=="NS" || fCRTTriggerType=="SN") FilterHLT=14;
    else
        throw std::runtime_error("CRT Filter not chosen correctly...");
}

bool sbndaq::CRTTriggerFilter::filter(art::Event &e)
{
    art::Handle<std::vector<raw::ptb::sbndptb>> ptbHandle_2;
    e.getByLabel(fPTBLabel, ptbHandle_2);
    for (int index = 0; index < int(ptbHandle_2->size()); index++)
    {
        auto ptb = (*ptbHandle_2)[index];
        auto hltrigs = ptb.GetHLTriggers();
        for (int HLT = 0; HLT < int(hltrigs.size()); HLT++)
        {
            std::vector<int> bits_up = bitsUP(hltrigs[HLT].trigger_word);
            for (size_t i = 0; i < bits_up.size(); i++)
            {
                //std::cout << " bit up " << bits_up[i] << std::endl;
                if(bits_up[i]==FilterHLT) return true; 
            }
        }
    }
    return false;
}

DEFINE_ART_MODULE(sbndaq::CRTTriggerFilter)
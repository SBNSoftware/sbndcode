#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"


//------------------------------------------------------------------------------
//--- opdet::sbndPDMapAlg implementation
//------------------------------------------------------------------------------

namespace opdet {

  sbndPDMapAlg::sbndPDMapAlg(const fhicl::ParameterSet&)
  {
    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file("sbnd_pds_mapping.json", fname);
    std::ifstream i(fname, std::ifstream::in);
    i >> PDmap;
    i.close();
  }

  sbndPDMapAlg::~sbndPDMapAlg()
  { }

  bool sbndPDMapAlg::isPDType(size_t ch, std::string pdname) const
  {
    if(PDmap.at(ch)["pd_type"] == std::string(pdname)) return true;
    return false;
  }

  std::string sbndPDMapAlg::pdType(size_t ch) const
  {
    if(ch < PDmap.size()) return PDmap.at(ch)["pd_type"];
    return "There is no such channel";
  }

  std::vector<int> sbndPDMapAlg::getChannelsOfType(std::string pdname) const
  {
    std::vector<int> out_ch_v;
    for (size_t ch = 0; ch < PDmap.size(); ch++) {
      if (PDmap.at(ch)["pd_type"] == pdname) out_ch_v.push_back(ch);
    }
    return out_ch_v;
  }

  size_t sbndPDMapAlg::size() const
  {
    return PDmap.size();
  }

  auto sbndPDMapAlg::getChannelEntry(size_t ch) const
  {
    return PDmap.at(ch);
  }

  // template<>
  nlohmann::json sbndPDMapAlg::getCollectionWithProperty(std::string property)
  {
    nlohmann::json subSetPDmap;
    std::copy_if (PDmap.begin(), PDmap.end(), std::back_inserter(subSetPDmap),
                  [property](const nlohmann::json e)->bool
                    {return e[property];} );
    return subSetPDmap;
  }

  // Look in the header for the implementation:
  // template<typename T>
  // nlohmann::json sbndPDMapAlg::getCollectionWithProperty(std::string property, T property_value)

DEFINE_ART_CLASS_TOOL(sbndPDMapAlg)
}

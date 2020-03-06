#include "sbndcode/OpDetSim/sbndPDMapAlg.h"


//------------------------------------------------------------------------------
//--- opdet::sbndPDMapAlg implementation
//------------------------------------------------------------------------------

namespace opdet {

  sbndPDMapAlg::sbndPDMapAlg()
  {
    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file("sbnd_pds_mapping.json", fname);
    std::ifstream i(fname);
    i >> PDmap;
  }

  sbndPDMapAlg::~sbndPDMapAlg()
  { }

  // TODO: I suspect this is not the most efficient way to do it, for
  // one it would be more memory efficient to have a container of
  // refs. Also this only allows to have one subcollection at a
  // time. ~icaza

  // TODO: instead of overloading, a better way would be to use a template. ~icaza
  nlohmann::json sbndPDMapAlg::getCollectionWithProperty(std::string property,
                                                         std::string property_value)
  {
    subSetPDmap.clear();
    std::copy_if (PDmap.begin(), PDmap.end(), std::back_inserter(subSetPDmap),
                  [property, property_value](const nlohmann::json e)->bool
                  {return e[property] == property_value;} );
    return subSetPDmap;
  }
  nlohmann::json sbndPDMapAlg::getCollectionWithProperty(std::string property,
                                                         int property_value)
  {
    subSetPDmap.clear();
    std::copy_if (PDmap.begin(), PDmap.end(), std::back_inserter(subSetPDmap),
                  [property, property_value](const nlohmann::json e)->bool
                  {return e[property] == property_value;} );
    return subSetPDmap;
  }
  // template<typename T>
  // nlohmann::json sbndPDMapAlg::getCollectionWithProperty(std::string property, T property_value)


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

  size_t sbndPDMapAlg::size() const
  {
    return PDmap.size();
  }

  auto sbndPDMapAlg::getChannelEntry(size_t ch) const
  {
    return PDmap.at(ch);
  }

}

#include "WeightNames.h"

struct WeightSet {
  std::string              name;
  std::vector<std::string> list;
  int                      nunivs;
};

typedef std::vector<WeightSet> WeightSets;

WeightSets weightSets = { { "flux", flux_weight_names_simple, 1000 },
                          { "genie", genie_weight_names, 500 },
                          { "geant4", geant4_weight_names, 1000 },
};

WeightSets weightSetsTest = { { "flux", tmp_flux_weight_names, 1000 },
};

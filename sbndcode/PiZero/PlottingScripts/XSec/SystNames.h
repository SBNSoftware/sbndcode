#pragma once

struct Syst {
  std::string name;
  int colour;
  std::string printName;
};


const std::vector<Syst> all_systs = { { "flux_weights_all", kMagenta+2, "Flux" },
                                      { "genie_all", kRed+2, "GENIE" },
                                      { "reinteractions_Geant4", kYellow+2, "Geant4" },
                                      { "ntargets", kGreen+2, "# Targets" },
                                      { "pot", kCyan+2, "POT" }
};

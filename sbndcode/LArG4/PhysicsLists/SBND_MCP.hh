#ifndef SBNDCODE_SBND_MCP_HH
#define SBNDCODE_SBND_MCP_HH 1

#include <Geant4/G4VModularPhysicsList.hh>

class SBND_MCP : public G4VModularPhysicsList {
public:
  SBND_MCP(int verbose = 1);
  virtual ~SBND_MCP() = default;

  SBND_MCP(const SBND_MCP &) = delete;
  SBND_MCP &operator=(const SBND_MCP &) = delete;
};


#endif //SBNDCODE_SBND_MCP_HH                                                                                                


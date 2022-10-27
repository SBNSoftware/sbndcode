/**
 * \brief Import CRT data from standalone DAQ into LArSoft
 *
 * \details To be written.
 *
 * \author Marco Del Tutto
 */

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Principal/EventPrincipal.h"

#include <fstream>
#include <vector>
#include <map>

#include "TFile.h"
#include "TTree.h"

class TTree;
class TFile;

namespace crt {
  class ImportCRTData {
  public:
    // Required constructor
    ImportCRTData(fhicl::ParameterSet const &pset,
                  art::ProductRegistryHelper &helper,
                  art::SourceHelper const &pm);

    // Required by FileReaderSource:
    void closeCurrentFile();
    void readFile(std::string const &name,
                  art::FileBlock* &fb);
    bool readNext(art::RunPrincipal* const &inR,
                  art::SubRunPrincipal* const &inSR,
                  art::RunPrincipal* &outR,
                  art::SubRunPrincipal* &outSR,
                  art::EventPrincipal* &outE);

  private:
    art::SourceHelper const      &fSourceHelper;
    art::SubRunID                 fSubRunID;

    bool                          fVerbose;
    uint32_t                      fEventCounter;
    int                           fMaxEvents;    // fhicl parameter.  Maximum number of events.
    float                         fPOT;
    float                         fCurrentPOT;
    float                         fSpills;
    float                         fCurrentSpills;
    uint32_t                      fTotalTreeEvents;

    std::map<unsigned, unsigned>               fMac5ToGeoID;
    std::vector<std::pair<unsigned, unsigned>> fMac5ToGeoIDVec;

    unsigned fT1Offset;

    TFile*                        fCRTInputFile;

    art::TypeLabel                fTLfebdata;

    TTree * fTree;
    std::vector<int> *fUnixS = nullptr; 
    std::vector<double> *fHit1Feb = nullptr;
    std::vector<UShort_t> *fHit1Flags = nullptr;
    std::vector<double> *fHit1T0 = nullptr;
    std::vector<double> *fHit1T1 = nullptr;
    std::vector<std::vector<uint16_t>> *fHit1Adc = nullptr;
    std::vector<double> *fHit2Feb = nullptr;
    std::vector<UShort_t> *fHit2Flags = nullptr;
    std::vector<double> *fHit2T0 = nullptr;
    std::vector<double> *fHit2T1 = nullptr;
    std::vector<std::vector<uint16_t>> *fHit2Adc = nullptr;
  };  // ImportCRTData
}



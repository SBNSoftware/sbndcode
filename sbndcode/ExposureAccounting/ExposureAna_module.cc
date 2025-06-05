////////////////////////////////////////////////////////////////////////
// Class:       ExposureAna
// Plugin Type: analyzer
// File:        ExposureAna_module.cc
// Author:      Henry Lay (h.lay@sheffield.ac.uk)
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Provenance/ProcessConfiguration.h"
#include "canvas/Persistency/Provenance/ProcessHistory.h"

#include "TTree.h"

#include "larcoreobj/SummaryData/POTSummary.h"

#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h"
#include "sbnobj/Common/POTAccounting/EXTCountInfo.h"

namespace sbnd {
  class ExposureAna;
}

class sbnd::ExposureAna : public art::EDAnalyzer {
public:
  explicit ExposureAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ExposureAna(ExposureAna const&) = delete;
  ExposureAna(ExposureAna&&) = delete;
  ExposureAna& operator=(ExposureAna const&) = delete;
  ExposureAna& operator=(ExposureAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void endSubRun(art::SubRun const& sr) override;

private:

  std::string fBNBSpillInfoModuleLabel, fEXTCountInfoModuleLabel, fMCPOTModuleLabel;
  TTree* fTree;

  // Tree variables

  int _run, _subrun;

  std::vector<double> _bnb_spill_pot;
  std::vector<float> _ext_count_gates_since_last_trigger;
  double _mc_pot;
  int _mc_spills, _mc_ngenevts;
};

sbnd::ExposureAna::ExposureAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  fBNBSpillInfoModuleLabel = p.get<std::string>("BNBSpillInfoModuleLabel");
  fEXTCountInfoModuleLabel = p.get<std::string>("EXTCountInfoModuleLabel");
  fMCPOTModuleLabel        = p.get<std::string>("MCPOTModuleLabel");

  art::ServiceHandle<art::TFileService> fs;

  fTree = fs->make<TTree>("tree","");
  fTree->Branch("run", &_run);
  fTree->Branch("subrun", &_subrun);

  fTree->Branch("bnb_spill_pot", &_bnb_spill_pot);
  fTree->Branch("ext_count_gates_since_last_trigger", &_ext_count_gates_since_last_trigger);
  fTree->Branch("mc_pot", &_mc_pot);
  fTree->Branch("mc_spills", &_mc_spills);
  fTree->Branch("mc_ngenevts", &_mc_ngenevts);
}

void sbnd::ExposureAna::analyze(art::Event const &e)
{
  _mc_ngenevts = 0;

  for(const art::ProcessConfiguration &process: e.processHistory())
    {
      std::optional<fhicl::ParameterSet> genConfig = e.getProcessParameterSet(process.processName());
      
      if(genConfig && genConfig->has_key("source") && genConfig->has_key("source.maxEvents") && genConfig->has_key("source.module_type"))
        {
          int maxEvents = genConfig->get<int>("source.maxEvents");
          std::string moduleType = genConfig->get<std::string>("source.module_type");

          if(moduleType == "EmptyEvent")
            {
              _mc_ngenevts += maxEvents;
            }
        }
    }
}

void sbnd::ExposureAna::endSubRun(art::SubRun const& sr)
{
  _run = sr.id().run();
  _subrun = sr.id().subRun();

  _bnb_spill_pot.clear();
  _ext_count_gates_since_last_trigger.clear();

  _mc_pot    = std::numeric_limits<double>::lowest();
  _mc_spills = std::numeric_limits<int>::lowest();

  // Get Data POT
  art::Handle<std::vector<sbn::BNBSpillInfo>> BNBSpillInfoHandle;
  sr.getByLabel(fBNBSpillInfoModuleLabel, BNBSpillInfoHandle);

  if(BNBSpillInfoHandle.isValid())
    {
      std::vector<art::Ptr<sbn::BNBSpillInfo>> BNBSpillInfoVec;
      art::fill_ptr_vector(BNBSpillInfoVec, BNBSpillInfoHandle);

      _bnb_spill_pot.resize(BNBSpillInfoVec.size());

      for(unsigned int i = 0; i < BNBSpillInfoVec.size(); ++i)
        _bnb_spill_pot[i] = BNBSpillInfoVec[i]->POT();
    }

  art::Handle<std::vector<sbn::EXTCountInfo>> EXTCountInfoHandle;
  sr.getByLabel(fEXTCountInfoModuleLabel, EXTCountInfoHandle);

  // Get Data Off-Beam Spills
  if(EXTCountInfoHandle.isValid())
    {
      std::vector<art::Ptr<sbn::EXTCountInfo>> EXTCountInfoVec;
      art::fill_ptr_vector(EXTCountInfoVec, EXTCountInfoHandle);

      _ext_count_gates_since_last_trigger.resize(EXTCountInfoVec.size());

      for(unsigned int i = 0; i < EXTCountInfoVec.size(); ++i)
        _ext_count_gates_since_last_trigger[i] = EXTCountInfoVec[i]->gates_since_last_trigger;
    }

  // Get MC POT
  art::Handle<sumdata::POTSummary> MCPOTHandle;
  sr.getByLabel(fMCPOTModuleLabel, MCPOTHandle);

  if(MCPOTHandle.isValid())
    {
      _mc_pot    = MCPOTHandle->totpot;
      _mc_spills = MCPOTHandle->totspills;
    }

  fTree->Fill();
}

DEFINE_ART_MODULE(sbnd::ExposureAna)

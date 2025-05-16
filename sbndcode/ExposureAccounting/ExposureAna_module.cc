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

#include "TTree.h"

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
  void analyze(art::Event const& e) override {};

  void endSubRun(art::SubRun const& sr) override;

private:

  std::string fBNBSpillInfoModuleLabel, fEXTCountInfoModuleLabel;
  TTree* fTree;

  // Tree variables

  int _run;
  int _subrun;

  std::vector<double> _bnb_spill_pot;
  std::vector<float> _ext_count_gates_since_last_trigger;
};

sbnd::ExposureAna::ExposureAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  fBNBSpillInfoModuleLabel = p.get<std::string>("BNBSpillInfoModuleLabel");
  fEXTCountInfoModuleLabel = p.get<std::string>("EXTCountInfoModuleLabel");

  art::ServiceHandle<art::TFileService> fs;

  fTree = fs->make<TTree>("tree","");
  fTree->Branch("run", &_run);
  fTree->Branch("subrun", &_subrun);

  fTree->Branch("bnb_spill_pot", &_bnb_spill_pot);
  fTree->Branch("ext_count_gates_since_last_trigger", &_ext_count_gates_since_last_trigger);
}

void sbnd::ExposureAna::endSubRun(art::SubRun const& sr)
{
  _run = sr.id().run();
  _subrun = sr.id().subRun();

  _bnb_spill_pot.clear();
  _ext_count_gates_since_last_trigger.clear();

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

  if(EXTCountInfoHandle.isValid())
    {
      std::vector<art::Ptr<sbn::EXTCountInfo>> EXTCountInfoVec;
      art::fill_ptr_vector(EXTCountInfoVec, EXTCountInfoHandle);

      _ext_count_gates_since_last_trigger.resize(EXTCountInfoVec.size());

      for(unsigned int i = 0; i < EXTCountInfoVec.size(); ++i)
        _ext_count_gates_since_last_trigger[i] = EXTCountInfoVec[i]->gates_since_last_trigger;
    }

  fTree->Fill();
}
DEFINE_ART_MODULE(sbnd::ExposureAna)

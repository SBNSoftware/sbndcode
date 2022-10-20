////////////////////////////////////////////////////////////////////////
// Class:       CRTSharpsAnalysis
// Plugin Type: analyzer
// File:        CRTSharpsAnalysis_module.cc
// Author:      Henry Lay (h.lay@lancaster.ac.uk
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
#include "canvas/Persistency/Common/FindManyP.h"

#include "TTree.h"

#include "sbnobj/SBND/CRT/FEBData.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

class CRTSharpsAnalysis;


class CRTSharpsAnalysis : public art::EDAnalyzer {
public:
  explicit CRTSharpsAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTSharpsAnalysis(CRTSharpsAnalysis const&) = delete;
  CRTSharpsAnalysis(CRTSharpsAnalysis&&) = delete;
  CRTSharpsAnalysis& operator=(CRTSharpsAnalysis const&) = delete;
  CRTSharpsAnalysis& operator=(CRTSharpsAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void AnalyseFEBDatas(std::vector<art::Ptr<sbnd::crt::FEBData>> &FEBDataVec);

private:

  sbnd::CRTGeoAlg fCRTGeoAlg;

  std::string fFEBDataModuleLabel;
  bool fDebug;

  TTree* fTree;

  // Tree variables

  int _run;
  int _subrun;
  int _event;

  std::vector<uint16_t>              _feb_mac5;
  std::vector<uint16_t>              _feb_flags;
  std::vector<uint32_t>              _feb_ts0;
  std::vector<uint32_t>              _feb_ts1;
  std::vector<uint32_t>              _feb_unixs;
  std::vector<std::vector<uint16_t>> _feb_adc;
  std::vector<uint32_t>              _feb_coinc;
};


CRTSharpsAnalysis::CRTSharpsAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  fFEBDataModuleLabel = p.get<std::string>("FEBDataLabel", "importdata");
  fDebug              = p.get<bool>("Debug", false);

  art::ServiceHandle<art::TFileService> fs;

  fTree = fs->make<TTree>("tree","");
  fTree->Branch("run", &_run);
  fTree->Branch("subrun", &_subrun);
  fTree->Branch("event", &_event);

  fTree->Branch("feb_mac5", "std::vector<uint16_t>", &_feb_mac5);
  fTree->Branch("feb_flags", "std::vector<uint16_t>", &_feb_flags);
  fTree->Branch("feb_ts0", "std::vector<uint32_t>", &_feb_ts0);
  fTree->Branch("feb_ts1", "std::vector<uint32_t>", &_feb_ts1);
  fTree->Branch("feb_unixs", "std::vector<uint32_t>", &_feb_unixs);
  fTree->Branch("feb_adc", "std::vector<std::vector<uint16_t>>", &_feb_adc);
  fTree->Branch("feb_coinc", "std::vector<uint32_t>", &_feb_coinc);

  if(fDebug)
    {
      for(auto const &[name, tagger] : fCRTGeoAlg.GetTaggers())
	{
	  std::cout << "Tagger:  " << tagger.name << '\n'
		    << "X - Min: " << tagger.minX << " Max: " << tagger.maxX << '\n'
		    << "Y - Min: " << tagger.minY << " Max: " << tagger.maxY << '\n'
		    << "Z - Min: " << tagger.minZ << " Max: " << tagger.maxZ << '\n' << std::endl;
	}
    }
}

void CRTSharpsAnalysis::analyze(art::Event const& e)
{
  _run = e.id().run();
  _subrun = e.id().subRun();
  _event =  e.id().event();

  if(fDebug) std::cout << "This is event " << _run << "-" << _subrun << "-" << _event << std::endl;

  // Get FEB Data
  art::Handle<std::vector<sbnd::crt::FEBData>> FEBDataHandle;
  e.getByLabel(fFEBDataModuleLabel, FEBDataHandle);
  if(!FEBDataHandle.isValid()){
    std::cout << "FEBData product " << fFEBDataModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<sbnd::crt::FEBData>> FEBDataVec;
  art::fill_ptr_vector(FEBDataVec, FEBDataHandle);

  // Fill FEBData variables
  AnalyseFEBDatas(FEBDataVec);

  // Fill the Tree
  fTree->Fill();
}

void CRTSharpsAnalysis::AnalyseFEBDatas(std::vector<art::Ptr<sbnd::crt::FEBData>> &FEBDataVec)
{
  unsigned nFEBData = FEBDataVec.size();

  _feb_mac5.resize(nFEBData);
  _feb_flags.resize(nFEBData);
  _feb_ts0.resize(nFEBData);
  _feb_ts1.resize(nFEBData);
  _feb_unixs.resize(nFEBData);
  _feb_adc.resize(nFEBData, std::vector<uint16_t>(32));
  _feb_coinc.resize(nFEBData);

  for (size_t i = 0; i < nFEBData; ++i)
    {
      auto FEBData = FEBDataVec[i];
      
      _feb_mac5[i]  = FEBData->Mac5();
      _feb_flags[i] = FEBData->Flags();
      _feb_ts0[i]   = FEBData->Ts0();
      _feb_ts1[i]   = FEBData->Ts1();
      _feb_unixs[i] = FEBData->UnixS();
      _feb_coinc[i] = FEBData->Coinc();

      for (size_t j = 0; j < 32; j++)
	_feb_adc[i][j] = FEBData->ADC(j);
    }
}

DEFINE_ART_MODULE(CRTSharpsAnalysis)

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
#include "sbnobj/Common/CRT/CRTHit.hh"

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

  void AnalyseCRTHits(std::vector<art::Ptr<sbn::crt::CRTHit>> &CRTHitVec, art::FindManyP<sbnd::crt::FEBData> &CRTHitToFEBData);

private:

  sbnd::CRTGeoAlg fCRTGeoAlg;

  std::string fFEBDataModuleLabel, fCRTHitModuleLabel;
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

  std::vector<double>                _chit_x;
  std::vector<double>                _chit_y;
  std::vector<double>                _chit_z;
  std::vector<double>                _chit_t0;
  std::vector<double>                _chit_t1;
  std::vector<double>                _chit_t1_diff;
  std::vector<uint64_t>              _chit_unix_s;
  std::vector<double>                _chit_h1_t0;
  std::vector<double>                _chit_h2_t0;
  std::vector<double>                _chit_h1_t1;
  std::vector<double>                _chit_h2_t1;
  std::vector<double>                _chit_pes;
  std::vector<int>                   _chit_plane;
  std::vector<std::vector<uint16_t>> _chit_sipm_adc;
  std::vector<std::vector<uint16_t>> _chit_sipm_corr_adc;
  std::vector<std::vector<uint16_t>> _chit_sipm_channel_id;
  std::vector<std::vector<uint16_t>> _chit_sipm_feb_mac5;
};


CRTSharpsAnalysis::CRTSharpsAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fCRTGeoAlg(p.get<fhicl::ParameterSet>("CRTGeoAlg", fhicl::ParameterSet()))
{
  fFEBDataModuleLabel = p.get<std::string>("FEBDataLabel", "importdata");
  fCRTHitModuleLabel  = p.get<std::string>("CRTHitLabel", "crthit");
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

  fTree->Branch("chit_x", "std::vector<double>", &_chit_x);
  fTree->Branch("chit_y", "std::vector<double>", &_chit_y);
  fTree->Branch("chit_z", "std::vector<double>", &_chit_z);
  fTree->Branch("chit_t0", "std::vector<double>", &_chit_t0);
  fTree->Branch("chit_t1", "std::vector<double>", &_chit_t1);
  fTree->Branch("chit_t1_diff", "std::vector<double>", &_chit_t1_diff);
  fTree->Branch("chit_unix_s", "std::vector<uint64_t>", &_chit_unix_s);
  fTree->Branch("chit_h1_t0", "std::vector<double>", &_chit_h1_t0);
  fTree->Branch("chit_h2_t0", "std::vector<double>", &_chit_h2_t0);
  fTree->Branch("chit_h1_t1", "std::vector<double>", &_chit_h1_t1);
  fTree->Branch("chit_h2_t1", "std::vector<double>", &_chit_h2_t1);
  fTree->Branch("chit_pes", "std::vector<double>", &_chit_pes);
  fTree->Branch("chit_plane", "std::vector<int>", &_chit_plane);
  fTree->Branch("chit_sipm_adc", "std::vector<std::vector<uint16_t> >", &_chit_sipm_adc);
  fTree->Branch("chit_sipm_corr_adc", "std::vector<std::vector<uint16_t> >", &_chit_sipm_corr_adc);
  fTree->Branch("chit_sipm_channel_id", "std::vector<std::vector<uint16_t> >", &_chit_sipm_channel_id);
  fTree->Branch("chit_sipm_feb_mac5", "std::vector<std::vector<uint16_t> >", &_chit_sipm_feb_mac5);

  if(fDebug)
    {
      for(auto const &[name, tagger] : fCRTGeoAlg.GetTaggers())
	{
	  std::cout << "Tagger:  " << tagger.name << '\n'
		    << "X - Min: " << tagger.minX << " Max: " << tagger.maxX << '\n'
		    << "Y - Min: " << tagger.minY << " Max: " << tagger.maxY << '\n'
		    << "Z - Min: " << tagger.minZ << " Max: " << tagger.maxZ << '\n' << std::endl;
	}

      std::cout << std::endl;

      for(auto const &[name, sipm] : fCRTGeoAlg.GetSiPMs())
	{
	  std::cout << "SiPM:  " << sipm.channel << " (" << sipm.channel/32 << " - " << sipm.channel%32 << ")" << '\n'
		    << "x: " << sipm.x << " y: " << sipm.y << " z: " << sipm.z << std::endl;
	}
    }
}

void CRTSharpsAnalysis::analyze(art::Event const& e)
{
  _run = e.id().run();
  _subrun = e.id().subRun();
  _event =  e.id().event();

  if(fDebug) std::cout << "This is event " << _run << "-" << _subrun << "-" << _event << std::endl;

  // Get FEBDatas
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

  // Get CRTHits
  art::Handle<std::vector<sbn::crt::CRTHit>> CRTHitHandle;
  e.getByLabel(fCRTHitModuleLabel, CRTHitHandle);
  if(!CRTHitHandle.isValid()){
    std::cout << "CRTHit product " << fCRTHitModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<sbn::crt::CRTHit>> CRTHitVec;
  art::fill_ptr_vector(CRTHitVec, CRTHitHandle);

  // Get CRTHit -> FEBData Assn
  art::FindManyP<sbnd::crt::FEBData> CRTHitToFEBData(CRTHitHandle, e, fCRTHitModuleLabel);

  // Fill CRTHit variables
  AnalyseCRTHits(CRTHitVec, CRTHitToFEBData);

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
      auto data = FEBDataVec[i];
      
      _feb_mac5[i]  = data->Mac5();
      _feb_flags[i] = data->Flags();
      _feb_ts0[i]   = data->Ts0();
      _feb_ts1[i]   = data->Ts1();
      _feb_unixs[i] = data->UnixS();
      _feb_coinc[i] = data->Coinc();

      for (size_t j = 0; j < 32; j++)
	_feb_adc[i][j] = data->ADC(j);
    }
}

void CRTSharpsAnalysis::AnalyseCRTHits(std::vector<art::Ptr<sbn::crt::CRTHit>> &CRTHitVec, art::FindManyP<sbnd::crt::FEBData> &CRTHitToFEBData)
{
  unsigned nCRTHits = CRTHitVec.size();

  _chit_x.resize(nCRTHits);
  _chit_y.resize(nCRTHits);
  _chit_z.resize(nCRTHits);
  _chit_t0.resize(nCRTHits);
  _chit_t1.resize(nCRTHits);
  _chit_t1_diff.resize(nCRTHits);
  _chit_unix_s.resize(nCRTHits);
  _chit_h1_t0.resize(nCRTHits);
  _chit_h2_t0.resize(nCRTHits);
  _chit_h1_t1.resize(nCRTHits);
  _chit_h2_t1.resize(nCRTHits);
  _chit_pes.resize(nCRTHits);
  _chit_plane.resize(nCRTHits);
  _chit_sipm_adc.resize(nCRTHits);
  _chit_sipm_corr_adc.resize(nCRTHits);
  _chit_sipm_channel_id.resize(nCRTHits);
  _chit_sipm_feb_mac5.resize(nCRTHits);

  for (size_t i = 0; i < nCRTHits; ++i)
    {
      auto hit = CRTHitVec[i];
      
      _chit_x[i]       = hit->x_pos;
      _chit_y[i]       = hit->y_pos;
      _chit_z[i]       = hit->z_pos;
      _chit_t0[i]      = hit->ts0_ns;
      _chit_t1[i]      = hit->ts1_ns;
      _chit_t1_diff[i] = hit->ts0_ns_corr; // the variable name in the object is old and is just a placeholder for diff, don't worry!
      _chit_unix_s[i]  = hit->ts0_s;
      _chit_pes[i]     = hit->peshit;

      if(hit->tagger == "volTESTTaggerSouth_0")
	_chit_plane[i] = 0; // upstream
      else
	_chit_plane[i] = 1; // downstream
      
      _chit_sipm_adc[i].resize(4);
      _chit_sipm_corr_adc[i].resize(4);
      const std::array<uint16_t, 4> adcs      = hit->raw_adcs;
      const std::array<uint16_t, 4> corr_adcs = hit->adcs;

      for(unsigned adc_i = 0; adc_i < 4; ++adc_i)
	{
	  _chit_sipm_adc[i][adc_i]      = adcs[adc_i];
	  _chit_sipm_corr_adc[i][adc_i] = corr_adcs[adc_i];
	}

      std::vector<art::Ptr<sbnd::crt::FEBData>> FEBDataVec = CRTHitToFEBData.at(hit.key());
      if(FEBDataVec.size() != 2)
	std::cout << "ERROR: CRTHit associated to " << FEBDataVec.size() << " FEBDatas" << std::endl;

      _chit_sipm_feb_mac5[i].resize(2);
      _chit_sipm_feb_mac5[i][0] = FEBDataVec[0]->Mac5();
      _chit_sipm_feb_mac5[i][1] = FEBDataVec[1]->Mac5();

      _chit_h1_t0[i] = FEBDataVec[0]->Ts0();
      _chit_h1_t1[i] = FEBDataVec[0]->Ts1();
      _chit_h2_t0[i] = FEBDataVec[1]->Ts0();
      _chit_h2_t1[i] = FEBDataVec[1]->Ts1();

      _chit_sipm_channel_id[i].resize(2);
      _chit_sipm_channel_id[i][0] = hit->channel0;
      _chit_sipm_channel_id[i][1] = hit->channel1;
    }
}

DEFINE_ART_MODULE(CRTSharpsAnalysis)

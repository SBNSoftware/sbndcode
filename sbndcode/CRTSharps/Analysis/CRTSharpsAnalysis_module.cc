////////////////////////////////////////////////////////////////////////
// Class:       CRTSharpsAnalysis
// Plugin Type: analyzer
// File:        CRTSharpsAnalysis_module.cc
// Author:      Henry Lay (h.lay@lancaster.ac.uk)
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
#include "sbnobj/SBND/Timing/DAQTimestamp.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"

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

  void AnalyseSPECTDCTimestamps(std::vector<art::Ptr<sbnd::timing::DAQTimestamp>> &DAQTimestampVec);

  void AnalyseCRTHits(std::vector<art::Ptr<sbn::crt::CRTHit>> &CRTHitVec, art::FindManyP<sbnd::crt::FEBData> &CRTHitToFEBData);

  void AnalyseCRTTracks(std::vector<art::Ptr<sbn::crt::CRTTrack>> &CRTTrackVec, art::FindManyP<sbn::crt::CRTHit> &CRTTrackToCRTHits);

private:

  sbnd::CRTGeoAlg fCRTGeoAlg;

  std::string fFEBDataModuleLabel, fSPECTDCModuleLabel, fCRTHitModuleLabel, fCRTTrackModuleLabel;
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

  std::vector<uint32_t>    _tdc_channel;
  std::vector<uint64_t>    _tdc_timestamp;
  std::vector<uint64_t>    _tdc_offset;
  std::vector<std::string> _tdc_name;

  std::vector<double>                _chit_x;
  std::vector<double>                _chit_y;
  std::vector<double>                _chit_z;
  std::vector<double>                _chit_ex;
  std::vector<double>                _chit_ey;
  std::vector<double>                _chit_ez;
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
  std::vector<std::vector<uint16_t>> _chit_sipm_raw_adc;
  std::vector<std::vector<uint16_t>> _chit_sipm_adc;
  std::vector<std::vector<uint16_t>> _chit_sipm_corr_adc;
  std::vector<std::vector<uint16_t>> _chit_sipm_channel_id;
  std::vector<std::vector<uint16_t>> _chit_sipm_feb_mac5;

  std::vector<double>                 _ct_time;
  std::vector<double>                 _ct_pes;
  std::vector<double>                 _ct_length;
  std::vector<double>                 _ct_tof;
  std::vector<double>                 _ct_hit1_x;
  std::vector<double>                 _ct_hit1_y;
  std::vector<double>                 _ct_hit1_z;
  std::vector<double>                 _ct_hit2_x;
  std::vector<double>                 _ct_hit2_y;
  std::vector<double>                 _ct_hit2_z;
  std::vector<double>                 _ct_hit1_ex;
  std::vector<double>                 _ct_hit1_ey;
  std::vector<double>                 _ct_hit1_ez;
  std::vector<double>                 _ct_hit2_ex;
  std::vector<double>                 _ct_hit2_ey;
  std::vector<double>                 _ct_hit2_ez;
  std::vector<double>                 _ct_hit1_t0;
  std::vector<double>                 _ct_hit1_t1;
  std::vector<double>                 _ct_hit2_t0;
  std::vector<double>                 _ct_hit2_t1;
  std::vector<uint16_t>               _ct_hit1_nhits;
  std::vector<uint16_t>               _ct_hit2_nhits;
  std::vector<std::vector<uint16_t> > _ct_hit1_sipm_raw_adc;
  std::vector<std::vector<uint16_t> > _ct_hit1_sipm_adc;
  std::vector<std::vector<uint16_t> > _ct_hit1_sipm_corr_adc;
  std::vector<std::vector<uint16_t> > _ct_hit2_sipm_raw_adc;
  std::vector<std::vector<uint16_t> > _ct_hit2_sipm_adc;
  std::vector<std::vector<uint16_t> > _ct_hit2_sipm_corr_adc;
};


CRTSharpsAnalysis::CRTSharpsAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fCRTGeoAlg(p.get<fhicl::ParameterSet>("CRTGeoAlg", fhicl::ParameterSet()))
  {
    fFEBDataModuleLabel  = p.get<std::string>("FEBDataModuleLabel", "importcrt");
    fSPECTDCModuleLabel  = p.get<std::string>("SPECTDCModuleLabel", "importspectdc");
    fCRTHitModuleLabel   = p.get<std::string>("CRTHitModuleLabel", "crthit");
    fCRTTrackModuleLabel = p.get<std::string>("CRTTrackModuleLabel", "crttrack");
    fDebug               = p.get<bool>("Debug", false);

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

    fTree->Branch("tdc_channel", "std::vector<uint32_t>", &_tdc_channel);
    fTree->Branch("tdc_timestamp", "std::vector<uint64_t>", &_tdc_timestamp);
    fTree->Branch("tdc_offset", "std::vector<uint64_t>", &_tdc_offset);
    fTree->Branch("tdc_name", "std::vector<std::string>", &_tdc_name);

    fTree->Branch("chit_x", "std::vector<double>", &_chit_x);
    fTree->Branch("chit_y", "std::vector<double>", &_chit_y);
    fTree->Branch("chit_z", "std::vector<double>", &_chit_z);
    fTree->Branch("chit_ex", "std::vector<double>", &_chit_ex);
    fTree->Branch("chit_ey", "std::vector<double>", &_chit_ey);
    fTree->Branch("chit_ez", "std::vector<double>", &_chit_ez);
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
    fTree->Branch("chit_sipm_raw_adc", "std::vector<std::vector<uint16_t> >", &_chit_sipm_raw_adc);
    fTree->Branch("chit_sipm_adc", "std::vector<std::vector<uint16_t> >", &_chit_sipm_adc);
    fTree->Branch("chit_sipm_corr_adc", "std::vector<std::vector<uint16_t> >", &_chit_sipm_corr_adc);
    fTree->Branch("chit_sipm_channel_id", "std::vector<std::vector<uint16_t> >", &_chit_sipm_channel_id);
    fTree->Branch("chit_sipm_feb_mac5", "std::vector<std::vector<uint16_t> >", &_chit_sipm_feb_mac5);

    fTree->Branch("ct_time", "std::vector<double>", &_ct_time);
    fTree->Branch("ct_pes", "std::vector<double>", &_ct_pes);
    fTree->Branch("ct_length", "std::vector<double>", &_ct_length);
    fTree->Branch("ct_tof", "std::vector<double>", &_ct_tof);
    fTree->Branch("ct_hit1_x", "std::vector<double>", &_ct_hit1_x);
    fTree->Branch("ct_hit1_y", "std::vector<double>", &_ct_hit1_y);
    fTree->Branch("ct_hit1_z", "std::vector<double>", &_ct_hit1_z);
    fTree->Branch("ct_hit2_x", "std::vector<double>", &_ct_hit2_x);
    fTree->Branch("ct_hit2_y", "std::vector<double>", &_ct_hit2_y);
    fTree->Branch("ct_hit2_z", "std::vector<double>", &_ct_hit2_z);
    fTree->Branch("ct_hit1_ex", "std::vector<double>", &_ct_hit1_ex);
    fTree->Branch("ct_hit1_ey", "std::vector<double>", &_ct_hit1_ey);
    fTree->Branch("ct_hit1_ez", "std::vector<double>", &_ct_hit1_ez);
    fTree->Branch("ct_hit2_ex", "std::vector<double>", &_ct_hit2_ex);
    fTree->Branch("ct_hit2_ey", "std::vector<double>", &_ct_hit2_ey);
    fTree->Branch("ct_hit2_ez", "std::vector<double>", &_ct_hit2_ez);
    fTree->Branch("ct_hit1_t0", "std::vector<double>", &_ct_hit1_t0);
    fTree->Branch("ct_hit1_t1", "std::vector<double>", &_ct_hit1_t1);
    fTree->Branch("ct_hit2_t0", "std::vector<double>", &_ct_hit2_t0);
    fTree->Branch("ct_hit2_t1", "std::vector<double>", &_ct_hit2_t1);
    fTree->Branch("ct_hit1_nhits", "std::vector<uint16_t>", &_ct_hit1_nhits);
    fTree->Branch("ct_hit2_nhits", "std::vector<uint16_t>", &_ct_hit2_nhits);
    fTree->Branch("ct_hit1_sipm_raw_adc", "std::vector<std::vector<uint16_t> >", &_ct_hit1_sipm_raw_adc);
    fTree->Branch("ct_hit1_sipm_adc", "std::vector<std::vector<uint16_t> >", &_ct_hit1_sipm_adc);
    fTree->Branch("ct_hit1_sipm_corr_adc", "std::vector<std::vector<uint16_t> >", &_ct_hit1_sipm_corr_adc);
    fTree->Branch("ct_hit2_sipm_raw_adc", "std::vector<std::vector<uint16_t> >", &_ct_hit2_sipm_raw_adc);
    fTree->Branch("ct_hit2_sipm_adc", "std::vector<std::vector<uint16_t> >", &_ct_hit2_sipm_adc);
    fTree->Branch("ct_hit2_sipm_corr_adc", "std::vector<std::vector<uint16_t> >", &_ct_hit2_sipm_corr_adc);

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

        for(auto const &[name, module] : fCRTGeoAlg.GetModules())
          {
            std::cout << "Module:  " << module.name << '\n'
                      << "X - Min: " << module.minX << " Max: " << module.maxX << '\n'
                      << "Y - Min: " << module.minY << " Max: " << module.maxY << '\n'
                      << "Z - Min: " << module.minZ << " Max: " << module.maxZ << '\n' << std::endl;
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

  // Get DAQTimestamps
  art::Handle<std::vector<sbnd::timing::DAQTimestamp>> DAQTimestampHandle;
  e.getByLabel(fSPECTDCModuleLabel, DAQTimestampHandle);
  if(!DAQTimestampHandle.isValid()){
    std::cout << "DAQTimestamp product " << fSPECTDCModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<sbnd::timing::DAQTimestamp>> DAQTimestampVec;
  art::fill_ptr_vector(DAQTimestampVec, DAQTimestampHandle);

  // Fill SPECTDC variables
  AnalyseSPECTDCTimestamps(DAQTimestampVec);

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

  // Get CRTTracks
  art::Handle<std::vector<sbn::crt::CRTTrack>> CRTTrackHandle;
  e.getByLabel(fCRTTrackModuleLabel, CRTTrackHandle);
  if(!CRTTrackHandle.isValid()){
    std::cout << "CRTTrack product " << fCRTTrackModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<sbn::crt::CRTTrack>> CRTTrackVec;
  art::fill_ptr_vector(CRTTrackVec, CRTTrackHandle);

  // Get CRTTrack -> CRTHit Assn
  art::FindManyP<sbn::crt::CRTHit> CRTTrackToCRTHits(CRTTrackHandle, e, fCRTTrackModuleLabel);

  // Fill CRTHit variables
  AnalyseCRTTracks(CRTTrackVec, CRTTrackToCRTHits);

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

  for(unsigned i = 0; i < nFEBData; ++i)
    {
      auto data = FEBDataVec[i];
      
      _feb_mac5[i]  = data->Mac5();
      _feb_flags[i] = data->Flags();
      _feb_ts0[i]   = data->Ts0();
      _feb_ts1[i]   = data->Ts1();
      _feb_unixs[i] = data->UnixS();
      _feb_coinc[i] = data->Coinc();

      for(unsigned j = 0; j < 32; ++j)
        _feb_adc[i][j] = data->ADC(j);
    }
}

void CRTSharpsAnalysis::AnalyseSPECTDCTimestamps(std::vector<art::Ptr<sbnd::timing::DAQTimestamp>> &DAQTimestampVec)
{
  unsigned nDAQTimestamps = DAQTimestampVec.size();

  _tdc_channel.resize(nDAQTimestamps);
  _tdc_timestamp.resize(nDAQTimestamps);
  _tdc_offset.resize(nDAQTimestamps);
  _tdc_name.resize(nDAQTimestamps);

  for(unsigned i = 0; i < nDAQTimestamps; ++i)
    {
      auto ts = DAQTimestampVec[i];

      _tdc_channel[i] = ts->Channel();
      _tdc_timestamp[i] = ts->Timestamp();
      _tdc_offset[i]    = ts->Offset();
      _tdc_name[i]      = ts->Name();
    }
}

void CRTSharpsAnalysis::AnalyseCRTHits(std::vector<art::Ptr<sbn::crt::CRTHit>> &CRTHitVec, art::FindManyP<sbnd::crt::FEBData> &CRTHitToFEBData)
{
  unsigned nCRTHits = CRTHitVec.size();

  _chit_x.resize(nCRTHits);
  _chit_y.resize(nCRTHits);
  _chit_z.resize(nCRTHits);
  _chit_ex.resize(nCRTHits);
  _chit_ey.resize(nCRTHits);
  _chit_ez.resize(nCRTHits);
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
  _chit_sipm_raw_adc.resize(nCRTHits);
  _chit_sipm_adc.resize(nCRTHits);
  _chit_sipm_corr_adc.resize(nCRTHits);
  _chit_sipm_channel_id.resize(nCRTHits);
  _chit_sipm_feb_mac5.resize(nCRTHits);

  for(unsigned i = 0; i < nCRTHits; ++i)
    {
      auto hit = CRTHitVec[i];
      
      _chit_x[i]       = hit->x_pos;
      _chit_y[i]       = hit->y_pos;
      _chit_z[i]       = hit->z_pos;
      _chit_ex[i]      = hit->x_err;
      _chit_ey[i]      = hit->y_err;
      _chit_ez[i]      = hit->z_err;
      _chit_t0[i]      = hit->ts0_ns;
      _chit_t1[i]      = hit->ts1_ns;
      _chit_t1_diff[i] = hit->ts0_ns_corr; // the variable name in the object is old and is just a placeholder for diff, don't worry!
      _chit_unix_s[i]  = hit->ts0_s;
      _chit_pes[i]     = hit->peshit;

      if(hit->tagger == "volTESTTaggerSouth_0")
        _chit_plane[i] = 0; // upstream
      else
        _chit_plane[i] = 1; // downstream
      
      _chit_sipm_raw_adc[i].resize(4);
      _chit_sipm_adc[i].resize(4);
      _chit_sipm_corr_adc[i].resize(4);
      const std::array<uint16_t, 4> raw_adcs  = hit->raw_adcs;
      const std::array<uint16_t, 4> adcs      = hit->adcs;
      const std::array<uint16_t, 4> corr_adcs = hit->corr_adcs;

      for(unsigned adc_i = 0; adc_i < 4; ++adc_i)
        {
          _chit_sipm_raw_adc[i][adc_i]      = raw_adcs[adc_i];
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

void CRTSharpsAnalysis::AnalyseCRTTracks(std::vector<art::Ptr<sbn::crt::CRTTrack>> &CRTTrackVec, art::FindManyP<sbn::crt::CRTHit> &CRTTrackToCRTHits)
{
  unsigned nCRTTracks = CRTTrackVec.size();

  _ct_pes.resize(nCRTTracks);
  _ct_time.resize(nCRTTracks);
  _ct_length.resize(nCRTTracks);
  _ct_tof.resize(nCRTTracks);
  _ct_hit1_x.resize(nCRTTracks);
  _ct_hit1_y.resize(nCRTTracks);
  _ct_hit1_z.resize(nCRTTracks);
  _ct_hit2_x.resize(nCRTTracks);
  _ct_hit2_y.resize(nCRTTracks);
  _ct_hit2_z.resize(nCRTTracks);
  _ct_hit1_ex.resize(nCRTTracks);
  _ct_hit1_ey.resize(nCRTTracks);
  _ct_hit1_ez.resize(nCRTTracks);
  _ct_hit2_ex.resize(nCRTTracks);
  _ct_hit2_ey.resize(nCRTTracks);
  _ct_hit2_ez.resize(nCRTTracks);
  _ct_hit1_t0.resize(nCRTTracks);
  _ct_hit1_t1.resize(nCRTTracks);
  _ct_hit2_t0.resize(nCRTTracks);
  _ct_hit2_t1.resize(nCRTTracks);
  _ct_hit1_nhits.resize(nCRTTracks);
  _ct_hit2_nhits.resize(nCRTTracks);
  _ct_hit1_sipm_raw_adc.resize(nCRTTracks);
  _ct_hit1_sipm_adc.resize(nCRTTracks);
  _ct_hit1_sipm_corr_adc.resize(nCRTTracks);
  _ct_hit2_sipm_raw_adc.resize(nCRTTracks);
  _ct_hit2_sipm_adc.resize(nCRTTracks);
  _ct_hit2_sipm_corr_adc.resize(nCRTTracks);

  for(unsigned i = 0; i < nCRTTracks; ++i)
    {
      auto track = CRTTrackVec[i];

      _ct_pes[i]     = track->peshit;
      _ct_time[i]    = track->ts1_ns;
      _ct_length[i]  = track->length;
      _ct_hit1_x[i]  = track->x1_pos;
      _ct_hit1_y[i]  = track->y1_pos;
      _ct_hit1_z[i]  = track->z1_pos;
      _ct_hit2_x[i]  = track->x2_pos;
      _ct_hit2_y[i]  = track->y2_pos;
      _ct_hit2_z[i]  = track->z2_pos;
      _ct_hit1_ex[i] = track->x1_err;
      _ct_hit1_ey[i] = track->y1_err;
      _ct_hit1_ez[i] = track->z1_err;
      _ct_hit2_ex[i] = track->x2_err;
      _ct_hit2_ey[i] = track->y2_err;
      _ct_hit2_ez[i] = track->z2_err;

      std::vector<art::Ptr<sbn::crt::CRTHit>> CRTHitVec = CRTTrackToCRTHits.at(track.key());
      _ct_hit1_nhits[i] = 0;
      _ct_hit2_nhits[i] = 0;

      for(art::Ptr<sbn::crt::CRTHit> hit : CRTHitVec)
        {
          if(std::signbit(hit->z_pos) == std::signbit(track->z1_pos))
            ++_ct_hit1_nhits[i];
          if(std::signbit(hit->z_pos) == std::signbit(track->z2_pos))
            ++_ct_hit2_nhits[i];
        }

      _ct_hit1_sipm_raw_adc[i].resize(4 * _ct_hit1_nhits[i]);
      _ct_hit1_sipm_adc[i].resize(4 * _ct_hit1_nhits[i]);
      _ct_hit1_sipm_corr_adc[i].resize(4 * _ct_hit1_nhits[i]);
      _ct_hit2_sipm_raw_adc[i].resize(4 * _ct_hit2_nhits[i]);
      _ct_hit2_sipm_adc[i].resize(4 * _ct_hit2_nhits[i]);
      _ct_hit2_sipm_corr_adc[i].resize(4 * _ct_hit2_nhits[i]);

      _ct_hit1_t0[i] = 0;
      _ct_hit1_t1[i] = 0;
      _ct_hit2_t0[i] = 0;
      _ct_hit2_t1[i] = 0;

      unsigned used_hits_1 = 0, used_hits_2 = 0;

      for(unsigned i_hit = 0; i_hit < CRTHitVec.size(); ++i_hit)
        {
          art::Ptr<sbn::crt::CRTHit> hit = CRTHitVec[i_hit];

          if(std::signbit(hit->z_pos) == std::signbit(track->z1_pos))
            {
              _ct_hit1_t0[i] += hit->ts0_ns;
              _ct_hit1_t1[i] += hit->ts1_ns;
            
              const std::array<uint16_t,4> raw_adcs  = hit->raw_adcs;
              const std::array<uint16_t,4> adcs      = hit->adcs;
              const std::array<uint16_t,4> corr_adcs = hit->corr_adcs;

              for(unsigned adc_i = 0; adc_i < 4; ++adc_i)
                {
                  _ct_hit1_sipm_raw_adc[i][4 * used_hits_1 + adc_i]  = raw_adcs[adc_i];
                  _ct_hit1_sipm_adc[i][4 * used_hits_1 + adc_i]      = adcs[adc_i];
                  _ct_hit1_sipm_corr_adc[i][4 * used_hits_1 + adc_i] = corr_adcs[adc_i];
                }
              ++used_hits_1;
            }
          else if(std::signbit(hit->z_pos) == std::signbit(track->z2_pos))
            {
              _ct_hit2_t0[i] += hit->ts0_ns;
              _ct_hit2_t1[i] += hit->ts1_ns;

              const std::array<uint16_t,4> raw_adcs  = hit->raw_adcs;
              const std::array<uint16_t,4> adcs      = hit->adcs;
              const std::array<uint16_t,4> corr_adcs = hit->corr_adcs;

              for(unsigned adc_i = 0; adc_i < 4; ++adc_i)
                {
                  _ct_hit2_sipm_raw_adc[i][4 * used_hits_2 + adc_i]  = raw_adcs[adc_i];
                  _ct_hit2_sipm_adc[i][4 * used_hits_2 + adc_i]      = adcs[adc_i];
                  _ct_hit2_sipm_corr_adc[i][4 * used_hits_2 + adc_i] = corr_adcs[adc_i];
                }
              ++used_hits_2;
            }
        }

      _ct_hit1_t0[i] /= _ct_hit1_nhits[i];
      _ct_hit1_t1[i] /= _ct_hit1_nhits[i];
      _ct_hit2_t0[i] /= _ct_hit2_nhits[i];
      _ct_hit2_t1[i] /= _ct_hit2_nhits[i];
      
      _ct_tof[i] = _ct_hit2_t1[i] - _ct_hit1_t1[i];
    }
}

DEFINE_ART_MODULE(CRTSharpsAnalysis)

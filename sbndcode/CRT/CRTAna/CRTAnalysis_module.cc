////////////////////////////////////////////////////////////////////////
// Class:       CRTAnalysis
// Plugin Type: analyzer
// File:        CRTAnalysis_module.cc
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
#include "sbnobj/SBND/CRT/CRTStripHit.hh"
#include "sbnobj/SBND/CRT/CRTCluster.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

class CRTAnalysis;


class CRTAnalysis : public art::EDAnalyzer {
public:
  explicit CRTAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTAnalysis(CRTAnalysis const&) = delete;
  CRTAnalysis(CRTAnalysis&&) = delete;
  CRTAnalysis& operator=(CRTAnalysis const&) = delete;
  CRTAnalysis& operator=(CRTAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  void AnalyseMCParticles(std::vector<art::Ptr<simb::MCParticle>> &MCParticleVec);

  void AnalyseSimDeposits(std::vector<art::Ptr<sim::AuxDetSimChannel>> &SimDepositVec);

  void AnalyseFEBDatas(std::vector<art::Ptr<sbnd::crt::FEBData>> &FEBDataVec);

  void AnalyseCRTStripHits(std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &CRTStripHitVec);

  void AnalyseCRTClusters(std::vector<art::Ptr<sbnd::crt::CRTCluster>> &CRTClusterVec);  

private:

  sbnd::CRTGeoAlg fCRTGeoAlg;

  std::string fMCParticleModuleLabel, fSimDepositModuleLabel, fFEBDataModuleLabel, fCRTStripHitModuleLabel,
    fCRTClusterModuleLabel;
  bool fDebug;

  TTree* fTree;

  // Tree variables

  int _run;
  int _subrun;
  int _event;

  std::vector<int16_t>              _mc_trackid;
  std::vector<int16_t>              _mc_pdg;
  std::vector<int16_t>              _mc_status;
  std::vector<uint16_t>             _mc_ndaughters;
  std::vector<std::vector<int16_t>> _mc_daughters;
  std::vector<double>               _mc_vx;
  std::vector<double>               _mc_vy;
  std::vector<double>               _mc_vz;
  std::vector<double>               _mc_vt;
  std::vector<double>               _mc_endx;
  std::vector<double>               _mc_endy;
  std::vector<double>               _mc_endz;
  std::vector<double>               _mc_endt;
  std::vector<double>               _mc_vpx;
  std::vector<double>               _mc_vpy;
  std::vector<double>               _mc_vpz;
  std::vector<double>               _mc_ve;
  std::vector<double>               _mc_endpx;
  std::vector<double>               _mc_endpy;
  std::vector<double>               _mc_endpz;
  std::vector<double>               _mc_ende;

  std::vector<int16_t> _ide_trackid;
  std::vector<float>   _ide_e;
  std::vector<float>   _ide_entryx;
  std::vector<float>   _ide_entryy;
  std::vector<float>   _ide_entryz;
  std::vector<float>   _ide_entryt;
  std::vector<float>   _ide_exitx;
  std::vector<float>   _ide_exity;
  std::vector<float>   _ide_exitz;
  std::vector<float>   _ide_exitt;

  std::vector<uint16_t>              _feb_mac5;
  std::vector<uint16_t>              _feb_flags;
  std::vector<uint32_t>              _feb_ts0;
  std::vector<uint32_t>              _feb_ts1;
  std::vector<uint32_t>              _feb_unixs;
  std::vector<std::vector<uint16_t>> _feb_adc;
  std::vector<uint32_t>              _feb_coinc;

  std::vector<uint32_t> _sh_channel;
  std::vector<uint32_t> _sh_ts0;
  std::vector<uint32_t> _sh_ts1;
  std::vector<uint32_t> _sh_unixs;
  std::vector<double>   _sh_pos;
  std::vector<double>   _sh_err;
  std::vector<uint16_t> _sh_adc1;
  std::vector<uint16_t> _sh_adc2;
  std::vector<bool>     _sh_saturated;

  std::vector<uint32_t>    _cl_ts0;
  std::vector<uint32_t>    _cl_ts1;
  std::vector<uint32_t>    _cl_unixs;
  std::vector<uint16_t>    _cl_nhits;
  std::vector<std::string> _cl_tagger;
  std::vector<bool>        _cl_threed;
};


CRTAnalysis::CRTAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fCRTGeoAlg(p.get<fhicl::ParameterSet>("CRTGeoAlg", fhicl::ParameterSet()))
  {
    fMCParticleModuleLabel  = p.get<std::string>("MCParticleModuleLabel", "largeant");
    fSimDepositModuleLabel  = p.get<std::string>("SimDepositModuleLabel", "genericcrt");
    fFEBDataModuleLabel     = p.get<std::string>("FEBDataModuleLabel", "crtsim");
    fCRTStripHitModuleLabel = p.get<std::string>("CRTStripHitModuleLabel", "crtstrips");
    fCRTClusterModuleLabel  = p.get<std::string>("CRTClusterModuleLabel", "crtclustering");
    fDebug                  = p.get<bool>("Debug", false);

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

    fTree->Branch("mc_trackid", "std::vector<int16_t>", &_mc_trackid);
    fTree->Branch("mc_pdg", "std::vector<int16_t>", &_mc_pdg);
    fTree->Branch("mc_status", "std::vector<int16_t>", &_mc_status);
    fTree->Branch("mc_ndaughters", "std::vector<uint16_t>", &_mc_ndaughters);
    fTree->Branch("mc_daughters", "std::vector<std::vector<int16_t>>", &_mc_daughters);
    fTree->Branch("mc_vx", "std::vector<double>", &_mc_vx);
    fTree->Branch("mc_vy", "std::vector<double>", &_mc_vy);
    fTree->Branch("mc_vz", "std::vector<double>", &_mc_vz);
    fTree->Branch("mc_vt", "std::vector<double>", &_mc_vt);
    fTree->Branch("mc_endx", "std::vector<double>", &_mc_endx);
    fTree->Branch("mc_endy", "std::vector<double>", &_mc_endy);
    fTree->Branch("mc_endz", "std::vector<double>", &_mc_endz);
    fTree->Branch("mc_endt", "std::vector<double>", &_mc_endt);
    fTree->Branch("mc_vpx", "std::vector<double>", &_mc_vpx);
    fTree->Branch("mc_vpy", "std::vector<double>", &_mc_vpy);
    fTree->Branch("mc_vpz", "std::vector<double>", &_mc_vpz);
    fTree->Branch("mc_ve", "std::vector<double>", &_mc_ve);
    fTree->Branch("mc_endpx", "std::vector<double>", &_mc_endpx);
    fTree->Branch("mc_endpy", "std::vector<double>", &_mc_endpy);
    fTree->Branch("mc_endpz", "std::vector<double>", &_mc_endpz);
    fTree->Branch("mc_ende", "std::vector<double>", &_mc_ende);

    fTree->Branch("ide_trackid", "std::vector<int16_t>", &_ide_trackid);
    fTree->Branch("ide_e", "std::vector<float>", &_ide_e);
    fTree->Branch("ide_entryx", "std::vector<float>", &_ide_entryx);
    fTree->Branch("ide_entryy", "std::vector<float>", &_ide_entryy);
    fTree->Branch("ide_entryz", "std::vector<float>", &_ide_entryz);
    fTree->Branch("ide_entryt", "std::vector<float>", &_ide_entryt);
    fTree->Branch("ide_exitx", "std::vector<float>", &_ide_exitx);
    fTree->Branch("ide_exity", "std::vector<float>", &_ide_exity);
    fTree->Branch("ide_exitz", "std::vector<float>", &_ide_exitz);
    fTree->Branch("ide_exitt", "std::vector<float>", &_ide_exitt);

    fTree->Branch("sh_channel", "std::vector<uint32_t>", &_sh_channel);
    fTree->Branch("sh_ts0", "std::vector<uint32_t>", &_sh_ts0);
    fTree->Branch("sh_ts1", "std::vector<uint32_t>", &_sh_ts1);
    fTree->Branch("sh_unixs", "std::vector<uint32_t>", &_sh_unixs);
    fTree->Branch("sh_pos", "std::vector<double>", &_sh_pos);
    fTree->Branch("sh_err", "std::vector<double>", &_sh_err);
    fTree->Branch("sh_adc1", "std::vector<uint16_t>", &_sh_adc1);
    fTree->Branch("sh_adc2", "std::vector<uint16_t>", &_sh_adc2);
    fTree->Branch("sh_saturated", "std::vector<bool>", &_sh_saturated);

    fTree->Branch("cl_ts0", "std::vector<uint32_t>", &_cl_ts0);
    fTree->Branch("cl_ts1", "std::vector<uint32_t>", &_cl_ts1);
    fTree->Branch("cl_unixs", "std::vector<uint32_t>", &_cl_unixs);
    fTree->Branch("cl_nhits", "std::vector<uint16_t>", &_cl_nhits);
    fTree->Branch("cl_tagger", "std::vector<std::string>", &_cl_tagger);
    fTree->Branch("cl_threed", "std::vector<bool>", &_cl_threed);

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

void CRTAnalysis::analyze(art::Event const& e)
{
  _run = e.id().run();
  _subrun = e.id().subRun();
  _event =  e.id().event();

  if(fDebug) std::cout << "This is event " << _run << "-" << _subrun << "-" << _event << std::endl;

  // Get MCParticles
  art::Handle<std::vector<simb::MCParticle>> MCParticleHandle;
  e.getByLabel(fMCParticleModuleLabel, MCParticleHandle);
  if(!MCParticleHandle.isValid()){
    std::cout << "MCParticle product " << fMCParticleModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<simb::MCParticle>> MCParticleVec;
  art::fill_ptr_vector(MCParticleVec, MCParticleHandle);

  // Fill MCParticle variables
  AnalyseMCParticles(MCParticleVec);

  // Get SimDeposits
  art::Handle<std::vector<sim::AuxDetSimChannel>> SimDepositHandle;
  e.getByLabel(fSimDepositModuleLabel, SimDepositHandle);
  if(!SimDepositHandle.isValid()){
    std::cout << "SimDeposit product " << fSimDepositModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<sim::AuxDetSimChannel>> SimDepositVec;
  art::fill_ptr_vector(SimDepositVec, SimDepositHandle);

  // Fill SimDeposit variables
  AnalyseSimDeposits(SimDepositVec);

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

  // Get CRTStripHits
  art::Handle<std::vector<sbnd::crt::CRTStripHit>> CRTStripHitHandle;
  e.getByLabel(fCRTStripHitModuleLabel, CRTStripHitHandle);
  if(!CRTStripHitHandle.isValid()){
    std::cout << "CRTStripHit product " << fCRTStripHitModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<sbnd::crt::CRTStripHit>> CRTStripHitVec;
  art::fill_ptr_vector(CRTStripHitVec, CRTStripHitHandle);

  // Fill CRTStripHit variables
  AnalyseCRTStripHits(CRTStripHitVec);

  // Get CRTClusters
  art::Handle<std::vector<sbnd::crt::CRTCluster>> CRTClusterHandle;
  e.getByLabel(fCRTClusterModuleLabel, CRTClusterHandle);
  if(!CRTClusterHandle.isValid()){
    std::cout << "CRTCluster product " << fCRTClusterModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<sbnd::crt::CRTCluster>> CRTClusterVec;
  art::fill_ptr_vector(CRTClusterVec, CRTClusterHandle);

  // Fill CRTCluster variables
  AnalyseCRTClusters(CRTClusterVec);

  // Fill the Tree
  fTree->Fill();
}

void CRTAnalysis::AnalyseMCParticles(std::vector<art::Ptr<simb::MCParticle>> &MCParticleVec)
{
  unsigned nMCParticles = MCParticleVec.size();

  _mc_trackid.resize(nMCParticles);
  _mc_pdg.resize(nMCParticles);
  _mc_status.resize(nMCParticles);
  _mc_ndaughters.resize(nMCParticles);
  _mc_daughters.resize(nMCParticles);
  _mc_vx.resize(nMCParticles);
  _mc_vy.resize(nMCParticles);
  _mc_vz.resize(nMCParticles);
  _mc_vt.resize(nMCParticles);
  _mc_endx.resize(nMCParticles);
  _mc_endy.resize(nMCParticles);
  _mc_endz.resize(nMCParticles);
  _mc_endt.resize(nMCParticles);
  _mc_vpx.resize(nMCParticles);
  _mc_vpy.resize(nMCParticles);
  _mc_vpz.resize(nMCParticles);
  _mc_ve.resize(nMCParticles);
  _mc_endpx.resize(nMCParticles);
  _mc_endpy.resize(nMCParticles);
  _mc_endpz.resize(nMCParticles);
  _mc_ende.resize(nMCParticles);
  
  for(unsigned i = 0; i < nMCParticles; ++i)
    {
      auto mcp = MCParticleVec[i];

      _mc_trackid[i]    = mcp->TrackId();
      _mc_pdg[i]        = mcp->PdgCode();
      _mc_status[i]     = mcp->StatusCode();
      _mc_ndaughters[i] = mcp->NumberDaughters();
      _mc_vx[i]         = mcp->Vx();
      _mc_vy[i]         = mcp->Vy();
      _mc_vz[i]         = mcp->Vz();
      _mc_vt[i]         = mcp->T();
      _mc_endx[i]       = mcp->EndX();
      _mc_endy[i]       = mcp->EndY();
      _mc_endz[i]       = mcp->EndZ();
      _mc_endt[i]       = mcp->EndT();
      _mc_vpx[i]        = mcp->Px();
      _mc_vpy[i]        = mcp->Py();
      _mc_vpz[i]        = mcp->Pz();
      _mc_ve[i]         = mcp->E();
      _mc_endpx[i]      = mcp->EndPx();
      _mc_endpy[i]      = mcp->EndPy();
      _mc_endpz[i]      = mcp->EndPz();
      _mc_ende[i]       = mcp->EndE();

      _mc_daughters[i].resize(_mc_ndaughters[i]);

      for(unsigned ii = 0; ii < _mc_ndaughters[i]; ++ii)
	_mc_daughters[i][ii] = mcp->Daughter(ii);
    }
}

void CRTAnalysis::AnalyseSimDeposits(std::vector<art::Ptr<sim::AuxDetSimChannel>> &SimDepositVec)
{
  unsigned nAuxDetSimChannels = SimDepositVec.size();
  unsigned nIDEs              = 0;

  for(unsigned i = 0; i < nAuxDetSimChannels; ++i)
    {
      auto auxDetSimChannel = SimDepositVec[i];
      nIDEs += auxDetSimChannel->AuxDetIDEs().size();
    }

  _ide_trackid.resize(nIDEs);
  _ide_e.resize(nIDEs);
  _ide_entryx.resize(nIDEs);
  _ide_entryy.resize(nIDEs);
  _ide_entryz.resize(nIDEs);
  _ide_entryt.resize(nIDEs);
  _ide_exitx.resize(nIDEs);
  _ide_exity.resize(nIDEs);
  _ide_exitz.resize(nIDEs);
  _ide_exitt.resize(nIDEs);

  unsigned ide_counter = 0;

  for(unsigned i = 0; i < nAuxDetSimChannels; ++i)
    {
      auto auxDetSimChannel = SimDepositVec[i];
      auto ideVec           = auxDetSimChannel->AuxDetIDEs();

      for(unsigned ii = 0; ii < ideVec.size(); ++ii)
	{
	  auto ide = ideVec[ii];

	  _ide_trackid[ide_counter] = ide.trackID;
	  _ide_e[ide_counter]       = ide.energyDeposited;
	  _ide_entryx[ide_counter]  = ide.entryX;
	  _ide_entryy[ide_counter]  = ide.entryY;
	  _ide_entryz[ide_counter]  = ide.entryZ;
	  _ide_entryt[ide_counter]  = ide.entryT;
	  _ide_exitx[ide_counter]   = ide.exitX;
	  _ide_exity[ide_counter]   = ide.exitY;
	  _ide_exitz[ide_counter]   = ide.exitZ;
	  _ide_exitt[ide_counter]   = ide.exitT;
	  
	  ++ide_counter;
	}
    }
}

void CRTAnalysis::AnalyseFEBDatas(std::vector<art::Ptr<sbnd::crt::FEBData>> &FEBDataVec)
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

void CRTAnalysis::AnalyseCRTStripHits(std::vector<art::Ptr<sbnd::crt::CRTStripHit>> &CRTStripHitVec)
{
  unsigned nStripHits = CRTStripHitVec.size();
  
  _sh_channel.resize(nStripHits);
  _sh_ts0.resize(nStripHits);
  _sh_ts1.resize(nStripHits);
  _sh_unixs.resize(nStripHits);
  _sh_pos.resize(nStripHits);
  _sh_err.resize(nStripHits);
  _sh_adc1.resize(nStripHits);
  _sh_adc2.resize(nStripHits);
  _sh_saturated.resize(nStripHits);

  for(unsigned i = 0; i < nStripHits; ++i)
    {
      auto hit = CRTStripHitVec[i];

      _sh_channel[i]   = hit->Channel();
      _sh_ts0[i]       = hit->Ts0();
      _sh_ts1[i]       = hit->Ts1();
      _sh_unixs[i]     = hit->UnixS();
      _sh_pos[i]       = hit->Pos();
      _sh_err[i]       = hit->Error();
      _sh_adc1[i]      = hit->ADC1();
      _sh_adc2[i]      = hit->ADC2();
      _sh_saturated[i] = hit->Saturated();
    }
}

void CRTAnalysis::AnalyseCRTClusters(std::vector<art::Ptr<sbnd::crt::CRTCluster>> &CRTClusterVec)
{
  unsigned nClusters = CRTClusterVec.size();

  _cl_ts0.resize(nClusters);
  _cl_ts1.resize(nClusters);
  _cl_unixs.resize(nClusters);
  _cl_nhits.resize(nClusters);
  _cl_tagger.resize(nClusters);
  _cl_threed.resize(nClusters);

  for(unsigned i = 0; i < nClusters; ++i)
    {
      auto cluster = CRTClusterVec[i];
      
      _cl_ts0[i]    = cluster->Ts0();
      _cl_ts1[i]    = cluster->Ts1();
      _cl_unixs[i]  = cluster->UnixS();
      _cl_nhits[i]  = cluster->NHits();
      _cl_tagger[i] = cluster->Tagger();
      _cl_threed[i] = cluster->ThreeD();
    }
}

DEFINE_ART_MODULE(CRTAnalysis)

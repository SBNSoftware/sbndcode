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
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"
#include "sbndcode/CRT/CRTBackTracker/CRTBackTrackerAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTCommonUtils.h"

namespace sbnd::crt {
  class CRTAnalysis;
}

class sbnd::crt::CRTAnalysis : public art::EDAnalyzer {
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

  void AnalyseFEBDatas(std::vector<art::Ptr<FEBData>> &FEBDataVec);

  void AnalyseCRTStripHits(const art::Event &e, const std::vector<art::Ptr<CRTStripHit>> &CRTStripHitVec);

  void AnalyseCRTClusters(const art::Event &e, const std::vector<art::Ptr<CRTCluster>> &CRTClusterVec, const art::FindManyP<CRTSpacePoint> &clustersToSpacePoints);

  void AnalyseTrueDeposits(const std::map<std::pair<int, CRTTagger>, bool> &recoStatusMap);

private:

  CRTGeoAlg fCRTGeoAlg;
  CRTBackTrackerAlg fCRTBackTrackerAlg;

  std::string fMCParticleModuleLabel, fSimDepositModuleLabel, fFEBDataModuleLabel, fCRTStripHitModuleLabel,
    fCRTClusterModuleLabel, fCRTSpacePointModuleLabel;
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
  std::vector<bool>     _sh_saturated1;
  std::vector<bool>     _sh_saturated2;
  std::vector<int>      _sh_truth_trackid;
  std::vector<double>   _sh_truth_completeness;
  std::vector<double>   _sh_truth_purity;
  std::vector<double>   _sh_truth_pos;
  std::vector<double>   _sh_truth_energy;
  std::vector<double>   _sh_truth_time;

  std::vector<uint32_t> _cl_ts0;
  std::vector<uint32_t> _cl_ts1;
  std::vector<uint32_t> _cl_unixs;
  std::vector<uint16_t> _cl_nhits;
  std::vector<int16_t>  _cl_tagger;
  std::vector<uint8_t>  _cl_composition;
  std::vector<int>      _cl_truth_trackid;
  std::vector<double>   _cl_truth_completeness;
  std::vector<double>   _cl_truth_purity;
  std::vector<double>   _cl_truth_hit_completeness;
  std::vector<double>   _cl_truth_hit_purity;
  std::vector<int>      _cl_truth_pdg;
  std::vector<double>   _cl_truth_x;
  std::vector<double>   _cl_truth_y;
  std::vector<double>   _cl_truth_z;
  std::vector<double>   _cl_truth_energy;
  std::vector<double>   _cl_truth_time;
  std::vector<bool>     _cl_has_sp;
  std::vector<double>   _cl_sp_x;
  std::vector<double>   _cl_sp_ex;
  std::vector<double>   _cl_sp_y;
  std::vector<double>   _cl_sp_ey;
  std::vector<double>   _cl_sp_z;
  std::vector<double>   _cl_sp_ez;
  std::vector<double>   _cl_sp_pe;
  std::vector<double>   _cl_sp_time;
  std::vector<double>   _cl_sp_etime;
  std::vector<bool>     _cl_sp_complete;

  std::vector<int>     _td_trackid;
  std::vector<int>     _td_pdg;
  std::vector<int16_t> _td_tagger;
  std::vector<double>  _td_x;
  std::vector<double>  _td_y;
  std::vector<double>  _td_z;
  std::vector<double>  _td_energy;
  std::vector<double>  _td_time;
  std::vector<bool>    _td_reco_status;
};

sbnd::crt::CRTAnalysis::CRTAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fCRTGeoAlg(p.get<fhicl::ParameterSet>("CRTGeoAlg", fhicl::ParameterSet()))
  , fCRTBackTrackerAlg(p.get<fhicl::ParameterSet>("CRTBackTrackerAlg", fhicl::ParameterSet()))
  {
    fMCParticleModuleLabel    = p.get<std::string>("MCParticleModuleLabel", "largeant");
    fSimDepositModuleLabel    = p.get<std::string>("SimDepositModuleLabel", "genericcrt");
    fFEBDataModuleLabel       = p.get<std::string>("FEBDataModuleLabel", "crtsim");
    fCRTStripHitModuleLabel   = p.get<std::string>("CRTStripHitModuleLabel", "crtstrips");
    fCRTClusterModuleLabel    = p.get<std::string>("CRTClusterModuleLabel", "crtclustering");
    fCRTSpacePointModuleLabel = p.get<std::string>("CRTSpacePointModuleLabel", "crtspacepoints");
    fDebug                    = p.get<bool>("Debug", false);

    art::ServiceHandle<art::TFileService> fs;

    fTree = fs->make<TTree>("tree","");
    fTree->Branch("run", &_run);
    fTree->Branch("subrun", &_subrun);
    fTree->Branch("event", &_event);

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

    fTree->Branch("feb_mac5", "std::vector<uint16_t>", &_feb_mac5);
    fTree->Branch("feb_flags", "std::vector<uint16_t>", &_feb_flags);
    fTree->Branch("feb_ts0", "std::vector<uint32_t>", &_feb_ts0);
    fTree->Branch("feb_ts1", "std::vector<uint32_t>", &_feb_ts1);
    fTree->Branch("feb_unixs", "std::vector<uint32_t>", &_feb_unixs);
    fTree->Branch("feb_adc", "std::vector<std::vector<uint16_t>>", &_feb_adc);
    fTree->Branch("feb_coinc", "std::vector<uint32_t>", &_feb_coinc);

    fTree->Branch("sh_channel", "std::vector<uint32_t>", &_sh_channel);
    fTree->Branch("sh_ts0", "std::vector<uint32_t>", &_sh_ts0);
    fTree->Branch("sh_ts1", "std::vector<uint32_t>", &_sh_ts1);
    fTree->Branch("sh_unixs", "std::vector<uint32_t>", &_sh_unixs);
    fTree->Branch("sh_pos", "std::vector<double>", &_sh_pos);
    fTree->Branch("sh_err", "std::vector<double>", &_sh_err);
    fTree->Branch("sh_adc1", "std::vector<uint16_t>", &_sh_adc1);
    fTree->Branch("sh_adc2", "std::vector<uint16_t>", &_sh_adc2);
    fTree->Branch("sh_saturated1", "std::vector<bool>", &_sh_saturated1);
    fTree->Branch("sh_saturated2", "std::vector<bool>", &_sh_saturated2);
    fTree->Branch("sh_truth_trackid", "std::vector<int>", &_sh_truth_trackid);
    fTree->Branch("sh_truth_completeness", "std::vector<double>", &_sh_truth_completeness);
    fTree->Branch("sh_truth_purity", "std::vector<double>", &_sh_truth_purity);
    fTree->Branch("sh_truth_pos", "std::vector<double>", &_sh_truth_pos);
    fTree->Branch("sh_truth_energy", "std::vector<double>", &_sh_truth_energy);
    fTree->Branch("sh_truth_time", "std::vector<double>", &_sh_truth_time);

    fTree->Branch("cl_ts0", "std::vector<uint32_t>", &_cl_ts0);
    fTree->Branch("cl_ts1", "std::vector<uint32_t>", &_cl_ts1);
    fTree->Branch("cl_unixs", "std::vector<uint32_t>", &_cl_unixs);
    fTree->Branch("cl_nhits", "std::vector<uint16_t>", &_cl_nhits);
    fTree->Branch("cl_tagger", "std::vector<int16_t>", &_cl_tagger);
    fTree->Branch("cl_composition", "std::vector<uint8_t>", &_cl_composition);
    fTree->Branch("cl_truth_trackid", "std::vector<int>", &_cl_truth_trackid);
    fTree->Branch("cl_truth_completeness", "std::vector<double>", &_cl_truth_completeness);
    fTree->Branch("cl_truth_purity", "std::vector<double>", &_cl_truth_purity);
    fTree->Branch("cl_truth_hit_completeness", "std::vector<double>", &_cl_truth_hit_completeness);
    fTree->Branch("cl_truth_hit_purity", "std::vector<double>", &_cl_truth_hit_purity);
    fTree->Branch("cl_truth_pdg", "std::vector<int>", &_cl_truth_pdg);
    fTree->Branch("cl_truth_x", "std::vector<double>", &_cl_truth_x);
    fTree->Branch("cl_truth_y", "std::vector<double>", &_cl_truth_y);
    fTree->Branch("cl_truth_z", "std::vector<double>", &_cl_truth_z);
    fTree->Branch("cl_truth_energy", "std::vector<double>", &_cl_truth_energy);
    fTree->Branch("cl_truth_time", "std::vector<double>", &_cl_truth_time);
    fTree->Branch("cl_has_sp", "std::vector<bool>", &_cl_has_sp);
    fTree->Branch("cl_sp_x", "std::vector<double>", &_cl_sp_x);
    fTree->Branch("cl_sp_ex", "std::vector<double>", &_cl_sp_ex);
    fTree->Branch("cl_sp_y", "std::vector<double>", &_cl_sp_y);
    fTree->Branch("cl_sp_ey", "std::vector<double>", &_cl_sp_ey);
    fTree->Branch("cl_sp_z", "std::vector<double>", &_cl_sp_z);
    fTree->Branch("cl_sp_ez", "std::vector<double>", &_cl_sp_ez);
    fTree->Branch("cl_sp_pe", "std::vector<double>", &_cl_sp_pe);
    fTree->Branch("cl_sp_time", "std::vector<double>", &_cl_sp_time);
    fTree->Branch("cl_sp_etime", "std::vector<double>", &_cl_sp_etime);
    fTree->Branch("cl_sp_complete", "std::vector<bool>", &_cl_sp_complete);

    fTree->Branch("td_trackid", "std::vector<int>", &_td_trackid);
    fTree->Branch("td_pdg", "std::vector<int>", &_td_pdg);
    fTree->Branch("td_tagger", "std::vector<int16_t>", &_td_tagger);
    fTree->Branch("td_x", "std::vector<double>", &_td_x);
    fTree->Branch("td_y", "std::vector<double>", &_td_y);
    fTree->Branch("td_z", "std::vector<double>", &_td_z);
    fTree->Branch("td_energy", "std::vector<double>", &_td_energy);
    fTree->Branch("td_time", "std::vector<double>", &_td_time);
    fTree->Branch("td_reco_status", "std::vector<bool>", &_td_reco_status);

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
            std::cout << "Module:  " << module.name << " (" << module.taggerName << ")" << '\n';
            if(module.minos)
              std::cout << "MINOS module" << std::endl;
            std::cout << "X - Min: " << module.minX << " Max: " << module.maxX << " Diff: " << module.maxX - module.minX << '\n' 
                      << "Y - Min: " << module.minY << " Max: " << module.maxY << " Diff: " << module.maxY - module.minY << '\n' 
                      << "Z - Min: " << module.minZ << " Max: " << module.maxZ << " Diff: " << module.maxZ - module.minZ << '\n' 
                      << "Orientation: " << module.orientation << '\n' << std::endl;
          }

        std::cout << std::endl;

        for(auto const &[name, sipm] : fCRTGeoAlg.GetSiPMs())
          {
            std::cout << "SiPM:  " << sipm.channel << " (" << sipm.channel/32 << " - " << sipm.channel%32 << ")" << '\n'
                      << "x: " << sipm.x << " y: " << sipm.y << " z: " << sipm.z << std::endl;
          }
      }
  }

void sbnd::crt::CRTAnalysis::analyze(art::Event const& e)
{
  fCRTBackTrackerAlg.SetupMaps(e);
  fCRTBackTrackerAlg.RunRecoStatusChecks(e);

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
  art::Handle<std::vector<FEBData>> FEBDataHandle;
  e.getByLabel(fFEBDataModuleLabel, FEBDataHandle);
  if(!FEBDataHandle.isValid()){
    std::cout << "FEBData product " << fFEBDataModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<FEBData>> FEBDataVec;
  art::fill_ptr_vector(FEBDataVec, FEBDataHandle);

  // Fill FEBData variables
  AnalyseFEBDatas(FEBDataVec);

  // Get CRTStripHits
  art::Handle<std::vector<CRTStripHit>> CRTStripHitHandle;
  e.getByLabel(fCRTStripHitModuleLabel, CRTStripHitHandle);
  if(!CRTStripHitHandle.isValid()){
    std::cout << "CRTStripHit product " << fCRTStripHitModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<CRTStripHit>> CRTStripHitVec;
  art::fill_ptr_vector(CRTStripHitVec, CRTStripHitHandle);

  // Fill CRTStripHit variables
  AnalyseCRTStripHits(e, CRTStripHitVec);

  // Get CRTClusters
  art::Handle<std::vector<CRTCluster>> CRTClusterHandle;
  e.getByLabel(fCRTClusterModuleLabel, CRTClusterHandle);
  if(!CRTClusterHandle.isValid()){
    std::cout << "CRTCluster product " << fCRTClusterModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<CRTCluster>> CRTClusterVec;
  art::fill_ptr_vector(CRTClusterVec, CRTClusterHandle);

  // Get CRTSpacePoint to CRTCluster Assns
  art::FindManyP<CRTSpacePoint> clustersToSpacePoints(CRTClusterHandle, e, fCRTSpacePointModuleLabel);

  // Fill CRTCluster variables
  AnalyseCRTClusters(e, CRTClusterVec, clustersToSpacePoints);

  // Get Map of TrueDeposits from BackTracker
  std::map<std::pair<int, CRTTagger>, bool> recoStatusMap = fCRTBackTrackerAlg.GetRecoStatusMap();
  
  // Fill TrueDeposit variables
  AnalyseTrueDeposits(recoStatusMap);

  // Fill the Tree
  fTree->Fill();
}

void sbnd::crt::CRTAnalysis::AnalyseMCParticles(std::vector<art::Ptr<simb::MCParticle>> &MCParticleVec)
{
  const unsigned nMCParticles = MCParticleVec.size();

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
      const auto mcp = MCParticleVec[i];

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

void sbnd::crt::CRTAnalysis::AnalyseSimDeposits(std::vector<art::Ptr<sim::AuxDetSimChannel>> &SimDepositVec)
{
  const unsigned nAuxDetSimChannels = SimDepositVec.size();
  unsigned nIDEs = 0;

  for(unsigned i = 0; i < nAuxDetSimChannels; ++i)
    {
      const auto auxDetSimChannel = SimDepositVec[i];
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
      const auto auxDetSimChannel = SimDepositVec[i];
      const auto ideVec           = auxDetSimChannel->AuxDetIDEs();

      for(unsigned ii = 0; ii < ideVec.size(); ++ii)
        {
          const auto ide = ideVec[ii];

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

void sbnd::crt::CRTAnalysis::AnalyseFEBDatas(std::vector<art::Ptr<FEBData>> &FEBDataVec)
{
  const unsigned nFEBData = FEBDataVec.size();

  _feb_mac5.resize(nFEBData);
  _feb_flags.resize(nFEBData);
  _feb_ts0.resize(nFEBData);
  _feb_ts1.resize(nFEBData);
  _feb_unixs.resize(nFEBData);
  _feb_adc.resize(nFEBData, std::vector<uint16_t>(32));
  _feb_coinc.resize(nFEBData);

  for(unsigned i = 0; i < nFEBData; ++i)
    {
      const auto data = FEBDataVec[i];
      
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

void sbnd::crt::CRTAnalysis::AnalyseCRTStripHits(const art::Event &e, const std::vector<art::Ptr<CRTStripHit>> &CRTStripHitVec)
{
  const unsigned nStripHits = CRTStripHitVec.size();
  
  _sh_channel.resize(nStripHits);
  _sh_ts0.resize(nStripHits);
  _sh_ts1.resize(nStripHits);
  _sh_unixs.resize(nStripHits);
  _sh_pos.resize(nStripHits);
  _sh_err.resize(nStripHits);
  _sh_adc1.resize(nStripHits);
  _sh_adc2.resize(nStripHits);
  _sh_saturated1.resize(nStripHits);
  _sh_saturated2.resize(nStripHits);
  _sh_truth_trackid.resize(nStripHits);
  _sh_truth_completeness.resize(nStripHits);
  _sh_truth_purity.resize(nStripHits);
  _sh_truth_pos.resize(nStripHits);
  _sh_truth_energy.resize(nStripHits);
  _sh_truth_time.resize(nStripHits);

  for(unsigned i = 0; i < nStripHits; ++i)
    {
      const auto hit = CRTStripHitVec[i];

      _sh_channel[i]    = hit->Channel();
      _sh_ts0[i]        = hit->Ts0();
      _sh_ts1[i]        = hit->Ts1();
      _sh_unixs[i]      = hit->UnixS();
      _sh_pos[i]        = hit->Pos();
      _sh_err[i]        = hit->Error();
      _sh_adc1[i]       = hit->ADC1();
      _sh_adc2[i]       = hit->ADC2();
      _sh_saturated1[i] = hit->Saturated1();
      _sh_saturated2[i] = hit->Saturated2();

      const CRTBackTrackerAlg::TruthMatchMetrics truthMatch = fCRTBackTrackerAlg.TruthMatching(e, hit);
      const std::vector<double> localpos = fCRTGeoAlg.StripWorldToLocalPos(hit->Channel(), truthMatch.deposit.x, truthMatch.deposit.y, truthMatch.deposit.z);
      const double width = fCRTGeoAlg.GetStrip(hit->Channel()).width;

      _sh_truth_trackid[i]      = truthMatch.trackid;
      _sh_truth_completeness[i] = truthMatch.completeness;
      _sh_truth_purity[i]       = truthMatch.purity;
      _sh_truth_pos[i]          = localpos[1] + width / 2.;
      _sh_truth_energy[i]       = truthMatch.deposit.energy;
      _sh_truth_time[i]         = truthMatch.deposit.time;
    }
}

void sbnd::crt::CRTAnalysis::AnalyseCRTClusters(const art::Event &e, const std::vector<art::Ptr<CRTCluster>> &CRTClusterVec, const art::FindManyP<CRTSpacePoint> &clustersToSpacePoints)
{
  const unsigned nClusters = CRTClusterVec.size();

  _cl_ts0.resize(nClusters);
  _cl_ts1.resize(nClusters);
  _cl_unixs.resize(nClusters);
  _cl_nhits.resize(nClusters);
  _cl_tagger.resize(nClusters);
  _cl_composition.resize(nClusters);
  _cl_truth_trackid.resize(nClusters);
  _cl_truth_completeness.resize(nClusters);
  _cl_truth_purity.resize(nClusters);
  _cl_truth_hit_completeness.resize(nClusters);
  _cl_truth_hit_purity.resize(nClusters);
  _cl_truth_pdg.resize(nClusters);
  _cl_truth_x.resize(nClusters);
  _cl_truth_y.resize(nClusters);
  _cl_truth_z.resize(nClusters);
  _cl_truth_energy.resize(nClusters);
  _cl_truth_time.resize(nClusters);
  _cl_has_sp.resize(nClusters);
  _cl_sp_x.resize(nClusters);
  _cl_sp_ex.resize(nClusters);
  _cl_sp_y.resize(nClusters);
  _cl_sp_ey.resize(nClusters);
  _cl_sp_z.resize(nClusters);
  _cl_sp_ez.resize(nClusters);
  _cl_sp_pe.resize(nClusters);
  _cl_sp_time.resize(nClusters);
  _cl_sp_etime.resize(nClusters);
  _cl_sp_complete.resize(nClusters);

  for(unsigned i = 0; i < nClusters; ++i)
    {
      const auto cluster = CRTClusterVec[i];

      _cl_ts0[i]         = cluster->Ts0();
      _cl_ts1[i]         = cluster->Ts1();
      _cl_unixs[i]       = cluster->UnixS();
      _cl_nhits[i]       = cluster->NHits();
      _cl_tagger[i]      = cluster->Tagger();
      _cl_composition[i] = cluster->Composition();

      const CRTBackTrackerAlg::TruthMatchMetrics truthMatch = fCRTBackTrackerAlg.TruthMatching(e, cluster);
      _cl_truth_trackid[i]          = truthMatch.trackid;
      _cl_truth_completeness[i]     = truthMatch.completeness;
      _cl_truth_purity[i]           = truthMatch.purity;
      _cl_truth_hit_completeness[i] = truthMatch.hitcompleteness;
      _cl_truth_hit_purity[i]       = truthMatch.hitpurity;
      _cl_truth_pdg[i]              = truthMatch.deposit.pdg;
      _cl_truth_x[i]                = truthMatch.deposit.x;
      _cl_truth_y[i]                = truthMatch.deposit.y;
      _cl_truth_z[i]                = truthMatch.deposit.z;
      _cl_truth_energy[i]           = truthMatch.deposit.energy;
      _cl_truth_time[i]             = truthMatch.deposit.time;

      const auto spacepoints = clustersToSpacePoints.at(cluster.key());
      if(spacepoints.size() == 1)
        { 
          const auto spacepoint = spacepoints[0];
          
          _cl_has_sp[i]      = true;
          _cl_sp_x[i]        = spacepoint->X();
          _cl_sp_ex[i]       = spacepoint->XErr();
          _cl_sp_y[i]        = spacepoint->Y();
          _cl_sp_ey[i]       = spacepoint->YErr();
          _cl_sp_z[i]        = spacepoint->Z();
          _cl_sp_ez[i]       = spacepoint->ZErr();
          _cl_sp_pe[i]       = spacepoint->PE();
          _cl_sp_time[i]     = spacepoint->Time();
          _cl_sp_etime[i]    = spacepoint->TimeErr();
          _cl_sp_complete[i] = spacepoint->Complete();
        }
      else
        {
          _cl_has_sp[i]      = false;
          _cl_sp_x[i]        = -999999.;
          _cl_sp_ex[i]       = -999999.;
          _cl_sp_y[i]        = -999999.;
          _cl_sp_ey[i]       = -999999.;
          _cl_sp_z[i]        = -999999.;
          _cl_sp_ez[i]       = -999999.;
          _cl_sp_pe[i]       = -999999.;
          _cl_sp_time[i]     = -999999.;
          _cl_sp_etime[i]    = -999999.;
          _cl_sp_complete[i] = false;
        }
    }
}

void sbnd::crt::CRTAnalysis::AnalyseTrueDeposits(const std::map<std::pair<int, CRTTagger>, bool> &recoStatusMap)
{
  const unsigned nTrueDeposits = recoStatusMap.size();

  _td_trackid.resize(nTrueDeposits);
  _td_pdg.resize(nTrueDeposits);
  _td_tagger.resize(nTrueDeposits);
  _td_x.resize(nTrueDeposits);
  _td_y.resize(nTrueDeposits);
  _td_z.resize(nTrueDeposits);
  _td_energy.resize(nTrueDeposits);
  _td_time.resize(nTrueDeposits);
  _td_reco_status.resize(nTrueDeposits);

  unsigned entry = 0;
  for(auto const& [category, status] : recoStatusMap)
    {
      const CRTBackTrackerAlg::TrueDeposit deposit = fCRTBackTrackerAlg.GetTrueDeposit(category);

      _td_trackid[entry]     = deposit.trackid;
      _td_pdg[entry]         = deposit.pdg;
      _td_tagger[entry]      = deposit.tagger;
      _td_x[entry]           = deposit.x;
      _td_y[entry]           = deposit.y;
      _td_z[entry]           = deposit.z;
      _td_energy[entry]      = deposit.energy;
      _td_time[entry]        = deposit.time;
      _td_reco_status[entry] = status;

      ++entry;
    }
}

DEFINE_ART_MODULE(sbnd::crt::CRTAnalysis)

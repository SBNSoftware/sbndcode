////////////////////////////////////////////////////////////////////////
// Class:       SBNDFlashAna
// Plugin Type: analyzer (art v3_04_00)
// File:        SBNDFlashAna_module.cc
//
// Generated at Sun Mar  1 15:08:24 2020 by Marco Del Tutto using cetskelgen
// from cetlib version v3_09_00.
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

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"

#include <TTree.h>

class SBNDFlashAna;


class SBNDFlashAna : public art::EDAnalyzer {
public:
  explicit SBNDFlashAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SBNDFlashAna(SBNDFlashAna const&) = delete;
  SBNDFlashAna(SBNDFlashAna&&) = delete;
  SBNDFlashAna& operator=(SBNDFlashAna const&) = delete;
  SBNDFlashAna& operator=(SBNDFlashAna&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  std::vector<std::string> _flash_label_v;
  std::vector<std::string> _ophit_label_v;
  std::vector<std::string> _mcflash_label_v;
  std::vector<int> _tpc_v;

  opdet::sbndPDMapAlg _pds_map;

  TTree *_flash_tree;
  int _run, _subrun, _event;
  int _tpc;
  double _flash_time;
  double _flash_total_pe;
  std::vector<double> _flash_pe_v;
  double _flash_y;
  double _flash_yerr ;
  double _flash_z;
  double _flash_zerr;
  std::vector<double> _flash_ophit_time;
  std::vector<double> _flash_ophit_amp;
  std::vector<double> _flash_ophit_area;
  std::vector<double> _flash_ophit_width;
  std::vector<double> _flash_ophit_pe;
  std::vector<int> _flash_ophit_ch;
  std::vector<std::string> _flash_ophit_chname;

  TTree *_mcflash_tree;
};


SBNDFlashAna::SBNDFlashAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  _flash_label_v   = p.get<std::vector<std::string>>("FlashProducer",   _flash_label_v);
  _mcflash_label_v = p.get<std::vector<std::string>>("MCFlashProducer", _mcflash_label_v);
  _ophit_label_v   = p.get<std::vector<std::string>>("OpHitProducer", _ophit_label_v);
  _tpc_v           = p.get<std::vector<int>>("TPC", _tpc_v);

  art::ServiceHandle<art::TFileService> fs;
  _flash_tree = fs->make<TTree>("FlashTree","");
  _flash_tree->Branch("run", &_run, "run/I");
  _flash_tree->Branch("subrun", &_subrun, "subrun/I");
  _flash_tree->Branch("event", &_event, "event/I");
  _flash_tree->Branch("tpc", &_tpc, "tpc/I");
  _flash_tree->Branch("flash_time", &_flash_time, "flash_time/D");
  _flash_tree->Branch("flash_total_pe", &_flash_total_pe, "flash_total_pe/D");
  _flash_tree->Branch("flash_y", &_flash_y, "flash_y/D");
  _flash_tree->Branch("flash_yerr", &_flash_yerr, "flash_yerr/D");
  _flash_tree->Branch("flash_z", &_flash_z, "flash_z/D");
  _flash_tree->Branch("flash_zerr", &_flash_zerr, "flash_zerr/D");
  _flash_tree->Branch("flash_pe_v", "std::vector<double>", &_flash_pe_v);
  _flash_tree->Branch("flash_ophit_time", "std::vector<double>", &_flash_ophit_time);
  _flash_tree->Branch("flash_ophit_amp", "std::vector<double>", &_flash_ophit_amp);
  _flash_tree->Branch("flash_ophit_area", "std::vector<double>", &_flash_ophit_area);
  _flash_tree->Branch("flash_ophit_width", "std::vector<double>", &_flash_ophit_width);
  _flash_tree->Branch("flash_ophit_pe", "std::vector<double>", &_flash_ophit_pe);
  _flash_tree->Branch("flash_ophit_ch", "std::vector<int>", &_flash_ophit_ch);
  _flash_tree->Branch("flash_ophit_chname", "std::vector<std::string>", &_flash_ophit_chname);

  _mcflash_tree = fs->make<TTree>("MCFlashTree","");
  _mcflash_tree->Branch("run", &_run, "run/I");
  _mcflash_tree->Branch("subrun", &_subrun, "subrun/I");
  _mcflash_tree->Branch("event", &_event, "event/I");
  _mcflash_tree->Branch("tpc", &_tpc, "tpc/I");
  _mcflash_tree->Branch("flash_time", &_flash_time, "flash_time/D");
  _mcflash_tree->Branch("flash_total_pe", &_flash_total_pe, "flash_total_pe/D");
  _mcflash_tree->Branch("flash_y", &_flash_y, "flash_y/D");
  _mcflash_tree->Branch("flash_yerr", &_flash_yerr, "flash_yerr/D");
  _mcflash_tree->Branch("flash_z", &_flash_z, "flash_z/D");
  _mcflash_tree->Branch("flash_zerr", &_flash_zerr, "flash_zerr/D");
  _mcflash_tree->Branch("flash_pe_v", "std::vector<double>", &_flash_pe_v);
}

void SBNDFlashAna::analyze(art::Event const& e)
{
  // if (_flash_label_v.size() != 2) {
  //   std::cout << "FlashProducer size has to be 2!" << std::endl;
  //   throw std::exception();
  // }
  // if (_mcflash_label_v.size() != 2) {
  //   std::cout << "MCFlashProducer size has to be 2!" << std::endl;
  //   throw std::exception();
  // }

  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();
  
  for (size_t l = 0; l < _flash_label_v.size(); l++) {
    art::Handle<std::vector<recob::OpFlash>> flash_h;
    e.getByLabel(_flash_label_v[l], flash_h);
    if(!flash_h.isValid()){
      std::cout << "Invalid producer for recob::OpFlash: " << _flash_label_v[l] << ". Ignoring." << std::endl;
    }

    art::FindManyP<recob::OpHit> flashToOpHitAssns(flash_h, e, _ophit_label_v[l]);

    for (size_t i = 0; i < flash_h->size(); i++) {
      auto const& f = (*flash_h)[i];
      std::cout << "Flash " << i << ", time " << f.AbsTime() << std::endl;
      _tpc = _tpc_v[l];
      _flash_time = f.AbsTime();
      _flash_total_pe = f.TotalPE();
      _flash_pe_v = f.PEs();
      _flash_y = f.YCenter();
      _flash_yerr = f.YWidth();
      _flash_z = f.ZCenter();
      _flash_zerr = f.ZWidth();

      _flash_ophit_time.clear();
      _flash_ophit_amp.clear();
      _flash_ophit_area.clear();
      _flash_ophit_width.clear();
      _flash_ophit_pe.clear();
      _flash_ophit_ch.clear();

      std::vector<art::Ptr<recob::OpHit>> ophit_v = flashToOpHitAssns.at(i);
      for (auto ophit : ophit_v) {
        _flash_ophit_time.push_back(ophit->PeakTimeAbs());
        _flash_ophit_amp.push_back(ophit->Amplitude());
        _flash_ophit_area.push_back(ophit->Area());
        _flash_ophit_width.push_back(ophit->Width());
        _flash_ophit_pe.push_back(ophit->PE());
        _flash_ophit_ch.push_back(ophit->OpChannel());
        _flash_ophit_chname.push_back(_pds_map.pdType(ophit->OpChannel()));
        
      }

      _flash_tree->Fill();
    }
  }

  for (size_t l = 0; l < _mcflash_label_v.size(); l++) {
    art::Handle<std::vector<recob::OpFlash>> flash_h;
    e.getByLabel(_mcflash_label_v[l], flash_h);
    if(!flash_h.isValid()){
      std::cout << "Invalid producer for MC recob::OpFlash: " << _mcflash_label_v[l] << ". Ignoring." << std::endl;
    }
    for (size_t i = 0; i < flash_h->size(); i++) {
      auto const& f = (*flash_h)[i];
      _tpc = _tpc_v[l];
      _flash_time = f.Time();
      _flash_total_pe = f.TotalPE();
      _flash_pe_v = f.PEs();
      _flash_y = f.YCenter();
      _flash_yerr = f.YWidth();
      _flash_z = f.ZCenter();
      _flash_zerr = f.ZWidth();
      _mcflash_tree->Fill();
    }
  }

}




DEFINE_ART_MODULE(SBNDFlashAna)






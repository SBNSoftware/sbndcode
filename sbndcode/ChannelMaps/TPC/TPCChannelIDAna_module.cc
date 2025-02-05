////////////////////////////////////////////////////////////////////////
// Class:       TPCChannelIDAna
// Plugin Type: analyzer
// File:        PDSMeanRMSLZ_tool.cc
// Author:      tjyang@fnal.gov
// 
// Analyze channel ID data 
// 
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

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larcore/Geometry/Geometry.h"

#include "sbndcode/ChannelMaps/TPC/TPCChannelMapService.h"

#include "TTree.h"

using namespace std;

namespace sbnd {
  class TPCChannelIDAna;
}


class sbnd::TPCChannelIDAna : public art::EDAnalyzer {
public:
  explicit TPCChannelIDAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TPCChannelIDAna(TPCChannelIDAna const&) = delete;
  TPCChannelIDAna(TPCChannelIDAna&&) = delete;
  TPCChannelIDAna& operator=(TPCChannelIDAna const&) = delete;
  TPCChannelIDAna& operator=(TPCChannelIDAna&&) = delete;

  // Required functions.
  void beginJob() override;
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.
  art::InputTag fRawDigitModuleLabel;
  TTree *chs;
  int ch;     //offline channel
  int chid;   //measurement from channel ID data
  int tpc;
  int plane;
  int wire;
  int wibcrate;
  int wib;
  int fembonwib;
  int fembch;
  int fembonwib_fromchid;
  int fembch_fromchid;
  int femcrate;
  int fem;
  int femch;
};


sbnd::TPCChannelIDAna::TPCChannelIDAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fRawDigitModuleLabel{p.get<art::InputTag>("RawDigitModuleLabel")}
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void sbnd::TPCChannelIDAna::beginJob(){
  art::ServiceHandle<art::TFileService> tfs;
  chs = tfs->make<TTree>("chs","Channel information");
  chs->Branch("ch",        &ch,        "ch/I");
  chs->Branch("chid",      &chid,      "chid/I");
  chs->Branch("tpc",       &tpc,       "tpc/I");
  chs->Branch("plane",     &plane,     "plane/I");
  chs->Branch("wire",      &wire,      "wire/I");
  chs->Branch("wibcrate",  &wibcrate,  "wibcrate/I");
  chs->Branch("wib",       &wib,       "wib/I");
  chs->Branch("fembonwib", &fembonwib, "fembonwib/I");
  chs->Branch("fembch",    &fembch,    "fembch/I");
  chs->Branch("fembonwib_fromchid", &fembonwib_fromchid, "fembonwib_fromchid/I");
  chs->Branch("fembch_fromchid", &fembch_fromchid, "fembch_fromchid/I");
  chs->Branch("femcrate",  &femcrate,  "femcrate/I");
  chs->Branch("fem",       &fem,       "fem/I");
  chs->Branch("femch",     &femch,     "femch/I");
}

void sbnd::TPCChannelIDAna::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<SBND::TPCChannelMapService> channelMap;

  // Get raw digits
  auto const& rawdigts = e.getProduct<std::vector<raw::RawDigit>>(fRawDigitModuleLabel);
  for (const auto & rd : rawdigts){
    std::vector<short> rawadc;      //UNCOMPRESSED ADC VALUES.
    rawadc.resize(rd.Samples());
    raw::Uncompress(rd.ADCs(), rawadc, rd.GetPedestal(), rd.Compression());
    //cout<<rd.Channel()<<" "<<rawadc[0]<<endl;
    ch = rd.Channel();
    chid = rawadc[0];
    auto const & chids = geo->ChannelToWire(ch);
    tpc = chids[0].TPC;
    plane = chids[0].Plane;
    wire = chids[0].Wire;
    auto const & chaninfo = channelMap->GetChanInfoFromOfflChan(ch);
    wibcrate = chaninfo.WIBCrate;
    wib = chaninfo.WIB;
    fembonwib = chaninfo.FEMBOnWIB;
    fembch = chaninfo.FEMBCh;
    femcrate = chaninfo.FEMCrate;
    fem = chaninfo.FEM;
    femch = chaninfo.FEMCh;    
    fembonwib_fromchid = chid>>8;
    fembch_fromchid = chid&0xFF;
    chs->Fill();
  }
}

DEFINE_ART_MODULE(sbnd::TPCChannelIDAna)

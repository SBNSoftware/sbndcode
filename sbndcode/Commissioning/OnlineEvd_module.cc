////////////////////////////////////////////////////////////////////////
// Class:       OnlineEvd
// Plugin Type: analyzer
// File:        OnlineEvd_module.cc
// Author:      tjyang@fnal.gov
//
// Analyzer module to make event display for DQM
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

#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TTimeStamp.h"

using namespace std;

namespace sbnd {
  class OnlineEvd;
}


class sbnd::OnlineEvd : public art::EDAnalyzer {
public:
  explicit OnlineEvd(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  OnlineEvd(OnlineEvd const&) = delete;
  OnlineEvd(OnlineEvd&&) = delete;
  OnlineEvd& operator=(OnlineEvd const&) = delete;
  OnlineEvd& operator=(OnlineEvd&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Declare member data here.
  art::InputTag fRawDigitModuleLabel;
  TH2D *hrawadc[6];
  int padx;
  int pady;
  int palette;
};


sbnd::OnlineEvd::OnlineEvd(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fRawDigitModuleLabel{p.get<art::InputTag>("RawDigitModuleLabel")}
  , padx{p.get<int>("padx", 1000)}
  , pady{p.get<int>("pady", 618)}
  , palette{p.get<int>("palette", 87)}
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void sbnd::OnlineEvd::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  int run = e.run();
  int event = e.id().event();
  art::Timestamp ts = e.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  
  //art::ServiceHandle<art::TFileService> tfs;

  int nbins[3] = {1984, 1984, 1664};
  for (int i = 0; i<6; ++i){
    int tpc = i/3;
    int plane = i%3;
    hrawadc[i] = new TH2D(Form("hrawadc%d",i),Form("Run %d Event %d TPC %d Plane %d;Wire;Tick;ADC", run, event, tpc, plane), nbins[plane], 0, nbins[plane]-1, 3400, 0, 3400);
    hrawadc[i]->SetMaximum(60);
    hrawadc[i]->SetMinimum(-20);
    hrawadc[i]->GetXaxis()->CenterTitle(true);
    hrawadc[i]->GetYaxis()->CenterTitle(true);
    hrawadc[i]->GetZaxis()->CenterTitle(true);
  }
  art::ServiceHandle<geo::Geometry> geo;

  // Get raw digits
  auto const& rawdigts = e.getProduct<std::vector<raw::RawDigit>>(fRawDigitModuleLabel);
  for (const auto & rd : rawdigts){
    std::vector<short> rawadc;      //UNCOMPRESSED ADC VALUES.
    rawadc.resize(rd.Samples());
    raw::Uncompress(rd.ADCs(), rawadc, rd.GetPedestal(), rd.Compression());
    int ch = rd.Channel();
    auto const & chids = geo->ChannelToWire(ch);
    int tpc = chids[0].TPC;
    int plane = chids[0].Plane;
    int wire = chids[0].Wire;
    //cout<<tpc<<" "<<plane<<" "<<wire<<endl;
    for (unsigned short i = 0; i<rawadc.size(); ++i){
      hrawadc[plane + 3*tpc]->Fill(wire, i, rawadc[i] - rd.GetPedestal());
    }
  }

  TLatex t;
  t.SetNDC();
  t.SetTextFont(132);
  t.SetTextSize(0.03);

  TCanvas *can = new TCanvas("can","can",padx,pady);
  for (int i = 0; i<6; ++i){
    int tpc = i/3;
    int plane = i%3;
    hrawadc[i]->Draw("colz");
    t.DrawLatex(0.01,0.01,Form("%d/%d/%d %d:%d:%d UTC",
                             tts.GetDate()/10000,
                             tts.GetDate()%10000/100,
                             tts.GetDate()%10000%100,
                             tts.GetTime()/10000,
                             tts.GetTime()%10000/100,
                             tts.GetTime()%10000%100));
    can->Print(Form("sbnd_run%d_event%d_tpc%d_plane%d.png",run,event,tpc,plane));
    delete hrawadc[i];
  }
  
  delete can;
}

void sbnd::OnlineEvd::beginJob()
{
  // Implementation of optional member function here.
  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(0.085);
  gStyle->SetPadTopMargin   (0.08);
  gStyle->SetPadLeftMargin  (0.08);
  gStyle->SetPadRightMargin (0.11);
  gStyle->SetLabelFont  ( 62   ,"XYZ");
  gStyle->SetTitleFont  ( 62   ,"XYZ");
  gStyle->SetTitleOffset( 1.15  , "x");
  gStyle->SetTitleOffset( 1.15  , "y");
  gStyle->SetTitleOffset( .7  , "z");
  gStyle->SetNumberContours(256);
  //gStyle->SetPalette(kRainBow);
  //gStyle->SetPalette(kLightTemperature);
  //gStyle->SetPalette(kGreenPink);
  gStyle->SetPalette(palette);
}

DEFINE_ART_MODULE(sbnd::OnlineEvd)

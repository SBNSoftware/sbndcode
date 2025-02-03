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

#include "sbndcode/ChannelMaps/TPC/TPCChannelMapService.h"

#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TTimeStamp.h"
#include "TLine.h"

#include "hep_hpc/hdf5/File.hpp"
#include "hep_hpc/hdf5/Ntuple.hpp"

using wire_nt_t = hep_hpc::hdf5::Ntuple<float>;

using namespace std;

double calcrms(const vector<short> & adcs){

  if (adcs.empty()) return -100;
  double x = 0, x2 = 0;
  for (auto const & adc: adcs){
    x += adc;
    x2 += adc*adc;
  }
  return sqrt(x2/adcs.size() - pow(x/adcs.size(),2));
}

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
  virtual ~OnlineEvd() noexcept;
private:

  // Declare member data here.
  art::InputTag fRawDigitModuleLabel;
  TH2D *hrawadc[6];
  int padx;
  int pady;
  int palette;
  TH2D *wibplot;
  TH2D *femplot;
  bool saveh5;
  hep_hpc::hdf5::File hdffile;

};


sbnd::OnlineEvd::OnlineEvd(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fRawDigitModuleLabel{p.get<art::InputTag>("RawDigitModuleLabel")}
  , padx{p.get<int>("padx", 1000*1.1)}
  , pady{p.get<int>("pady", 618*1.1)}
  , palette{p.get<int>("palette", 87)}
  , saveh5{p.get<bool>("saveh5", true)}
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}
sbnd::OnlineEvd::~OnlineEvd() noexcept
{
}

void sbnd::OnlineEvd::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  int run = e.run();
  int event = e.id().event();
  art::Timestamp ts = e.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  
  //art::ServiceHandle<art::TFileService> tfs;
  if (saveh5){
    hdffile = hep_hpc::hdf5::File(Form("r%d_e%d.h5",run,event), H5F_ACC_TRUNC);
  }
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

  wibplot = new TH2D("wibplot", Form("Run %d Event %d;FEMB;WIB;RMS", run, event), 8*128, 0, 8*128, 12, 0, 12);
  wibplot->SetMaximum(10);
  wibplot->SetMinimum(-1);
  wibplot->GetXaxis()->CenterTitle(true);
  wibplot->GetYaxis()->CenterTitle(true);
  wibplot->GetZaxis()->CenterTitle(true);
  for (int i = 1; i<=wibplot->GetNbinsX(); ++i){
    for (int j = 1; j<=wibplot->GetNbinsY(); ++j){
      wibplot->SetBinContent(i, j, -100);
    }
  }

  femplot = new TH2D("femplot", Form("Run %d Event %d;FEM;TPC Crate;RMS", run, event), 16*64, 0, 16*64, 11, 0, 11);
  femplot->SetMaximum(10);
  femplot->SetMinimum(-1);
  femplot->GetXaxis()->CenterTitle(true);
  femplot->GetYaxis()->CenterTitle(true);
  femplot->GetZaxis()->CenterTitle(true);
  for (int i = 1; i<=femplot->GetNbinsX(); ++i){
    for (int j = 1; j<=femplot->GetNbinsY(); ++j){
      femplot->SetBinContent(i, j, -100);
    }
  }

  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<SBND::TPCChannelMapService> channelMap;

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
    // Fill event display
    for (unsigned short i = 0; i<rawadc.size(); ++i){
      hrawadc[plane + 3*tpc]->Fill(wire, i, rawadc[i] - rd.GetPedestal());
    }
    // Calculate rms
    double rms = calcrms(rawadc);
    if (rms >10) rms = 10;
    auto const & chaninfo = channelMap->GetChanInfoFromOfflChan(ch);
    // Fill wib plot
    int binx = (chaninfo.WIBCrate-1)%2*128*4 +
      chaninfo.FEMBOnWIB * 128 +
      chaninfo.FEMBCh + 1;
    int biny = (1-(chaninfo.WIBCrate-1)/2)*6 +
      5 - (chaninfo.WIB-1) + 1;
    wibplot->SetBinContent(binx, biny, rms);
    // Fill fem plot
    binx = (chaninfo.FEM-1)*64 + chaninfo.FEMCh + 1;
    biny = 11-chaninfo.FEMCrate + 1;
    femplot->SetBinContent(binx, biny, rms);
//    if (ch == 7163 || ch == 8372){
//      cout<<ch<<" "<<binx<<" "<<biny<<" "<<rms<<" "<<chaninfo.FEMCrate<<" "<<chaninfo.FEM<<" "<<chaninfo.FEMCh<<endl;
//    }
  }

  if (saveh5){
    wire_nt_t rawwf(hdffile, "rawwf", {"adc"});
    for (int i = 0; i<6; ++i){
      for (int j = 1; j<=hrawadc[i]->GetNbinsX(); ++j){
        for (int k = 1; k<=hrawadc[i]->GetNbinsY(); ++k){
          rawwf.insert(hrawadc[i]->GetBinContent(j,k));
        }
      }
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

  wibplot->Draw("colz");
  wibplot->GetXaxis()->SetNdivisions(408,false);
  wibplot->GetYaxis()->SetNdivisions(112,false);
  wibplot->GetXaxis()->SetLabelSize(0);
  wibplot->GetYaxis()->SetLabelSize(0);
  TLine l1(0,6,128*8,6);
  l1.SetLineWidth(2);
  l1.Draw();
  TLine l2(128*4,0,128*4,12);
  l2.SetLineWidth(2);
  l2.Draw();

  TLine *l3[5];
  TLine *l4[5];
  for (int i = 0; i<5; ++i){
    l3[i] = new TLine(0, i+1, 128*8, i+1);
    l3[i]->SetLineWidth(1);
    l3[i]->Draw();
    l4[i] = new TLine(0, i+7, 128*8, i+7);
    l4[i]->SetLineWidth(1);
    l4[i]->Draw();
  }
  TLine *l5[8];
  for (int i = 0; i<8; ++i){
    if (i%4==3) continue;
    l5[i] = new TLine((i+1)*128,0,(i+1)*128,12);
    l5[i]->SetLineColor(0);
    l5[i]->SetLineStyle(2);
    l5[i]->SetLineWidth(1);
    l5[i]->Draw();
  }

  TLine *l6[32];
  for (int i = 0; i<32; ++i){
    if (i%4==3) continue;
    l6[i] = new TLine((i+1)*32,0,(i+1)*32,12);
    l6[i]->SetLineColor(0);
    l6[i]->SetLineStyle(3);
    l6[i]->SetLineWidth(1);
    l6[i]->Draw();
  }

  TLatex t2;
  //t2.SetTextFont(132);
  t2.SetTextSize(0.03);
  for (int i = 0; i<8; ++i){
    t2.DrawLatex(128*i+60,-0.4,Form("%d",i%4));
  }
  for (int i = 0; i<12; ++i){
    t2.DrawLatex(-30,i+0.4,Form("%d",6-i%6));
  }
  t2.DrawLatex(10,11.5,"SW1");
  t2.DrawLatex(10+128*4,11.5,"NW2");
  t2.DrawLatex(10,5.5,"SE3");
  t2.DrawLatex(10+128*4,5.5,"NE4");
  t.DrawLatex(0.01,0.01,Form("%d/%d/%d %d:%d:%d UTC",
                             tts.GetDate()/10000,
                             tts.GetDate()%10000/100,
                             tts.GetDate()%10000%100,
                             tts.GetTime()/10000,
                             tts.GetTime()%10000/100,
                             tts.GetTime()%10000%100));
  can->Print("wibrms.png");

  femplot->Draw("colz");
  femplot->GetXaxis()->SetNdivisions(116,false);
  femplot->GetYaxis()->SetNdivisions(111,false);
  femplot->GetXaxis()->SetLabelSize(0);
  femplot->GetYaxis()->SetLabelSize(0);
  TLine *l7[11];
  for (int i = 0; i<11; ++i){
    l7[i] = new TLine(0, i+1, 16*64, i+1);
    l7[i]->SetLineWidth(1);
    l7[i]->Draw();
  }
  TLine *l8[15];
  for (int i = 0; i<15; ++i){
    l8[i] = new TLine((i+1)*64, 0, (i+1)*64, 11);
    l8[i]->SetLineWidth(1);
    l8[i]->SetLineColor(0);
    l8[i]->SetLineStyle(3);
    l8[i]->Draw();
  }
  for (int i = 0; i<16; ++i){
    t2.DrawLatex(64*i+28,-0.4,Form("%d",i+1));
  }
  for (int i = 0; i<11; ++i){
    t2.DrawLatex(-30,i+0.4,Form("%d",11-i));
  }
  t.DrawLatex(0.01,0.01,Form("%d/%d/%d %d:%d:%d UTC",
                             tts.GetDate()/10000,
                             tts.GetDate()%10000/100,
                             tts.GetDate()%10000%100,
                             tts.GetTime()/10000,
                             tts.GetTime()%10000/100,
                             tts.GetTime()%10000%100));
  can->Print("femrms.png");

  delete wibplot;
  delete femplot;
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

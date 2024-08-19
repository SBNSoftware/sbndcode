// Plugin Type: analyzer (Unknown Unknown)
// File:        FlagChoppy_module.cc
//
// Generated at Thu Jul 18 16:37:03 2024 by Mun Jung Jung using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <numeric>
#include <getopt.h>
#include <chrono>
#include <float.h>


//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TFFTReal.h"
#include "TTimeStamp.h"
#include "TLatex.h"
#include "TLine.h"


// art includes
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h" 
#include "art/Framework/Principal/SubRun.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"  
#include "lardataobj/RecoBase/Hit.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardataobj/RawData/raw.h"
#include "larcore/Geometry/Geometry.h"
#include "canvas/Persistency/Provenance/Timestamp.h"

#include "art_root_io/TFileService.h"

#include "sbndcode/ChannelMaps/TPC/TPCChannelMapService.h"

//#include "Analysis.hh"
//#include "ChannelData.hh"
//#include "sbndqm/Decode/TPC/HeaderData.hh"
//#include "FFT.hh"
//#include "Noise.hh"
//#include "PeakFinder.hh"

class FlagChoppy;


class FlagChoppy : public art::EDAnalyzer {
public:
  explicit FlagChoppy(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FlagChoppy(FlagChoppy const&) = delete;
  FlagChoppy(FlagChoppy&&) = delete;
  FlagChoppy& operator=(FlagChoppy const&) = delete;
  FlagChoppy& operator=(FlagChoppy&&) = delete;

  // Required functions.
  void analyze(art::Event const& evt) override;

  void beginJob() override;
  void endJob() override;

  std::vector<art::Ptr<raw::RDTimeStamp>> raw_timestamps_handle;

private:

  // Declare member data here.
  TTree *fEventTree;
  int run;
  int subrun;
  int event;
  long delta_tmaxmin;
  std::vector<long> ts_vec;
  std::vector<long> ts_hist_vec;
  std::vector<int> evd_vec_bin[6];
  std::vector<double> evd_vec_q[6];
  int npeak;

  TH2D *hrawadc[6];
  int padx;
  int pady;
  int palette;
  double qmax;
  double qmin;
};


FlagChoppy::FlagChoppy(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void FlagChoppy::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs; // TTree's are created in the memory managed by ROOT (you don't delete them)


  fEventTree = tfs->make<TTree>("EventTree", "event by event info");
  fEventTree->Branch("run", &run);
  fEventTree->Branch("subrun", &subrun);
  fEventTree->Branch("event", &event);
  fEventTree->Branch("delta_tmaxmin", &delta_tmaxmin);
  fEventTree->Branch("ts_vec", "std::vector<long>", &ts_vec);

  for (int i = 0; i<6; ++i){
    fEventTree->Branch(Form("evd_vec_bins_%d",i), "std::vector<int>", &evd_vec_bin[i]);
    fEventTree->Branch(Form("evd_vec_q_%d", i), "std::vector<double>", &evd_vec_q[i]);
  }
}

void FlagChoppy::analyze(art::Event const& evt)
{
  // Implementation of required member function here.
  raw_timestamps_handle.clear();
  run = evt.id().run();
  subrun = evt.id().subRun();
  event = evt.id().event();
  delta_tmaxmin = 0;

  ts_vec.clear();
  ts_hist_vec.clear();
  for (int i = 0; i<6; ++i){
    evd_vec_bin[i].clear();
    evd_vec_q[i].clear();
  }

  npeak = 0;


  // get the timestamps
  art::Handle<std::vector<raw::RDTimeStamp>> timestamp_handle;
  evt.getByLabel("daq", timestamp_handle);

  // exit if the data isn't present
  if (!timestamp_handle.isValid()) {
    std::cerr << "Error: missing timestamps with producer (" << "daq" << ")" << std::endl;
    return;
  }
  art::fill_ptr_vector(raw_timestamps_handle, timestamp_handle);

  //get the association
//  art::FindManyP<raw::RawDigit> rawdigFMts(timestamp_handle, evt, "raw::RDTimeStamp");
  //art::Assns<raw::RawDigit,raw::RDTimeStamp>

  // analyze timestamp distributions
  // collect timestamps
  unsigned timeindex = 0;
  for (auto const& timestamps: raw_timestamps_handle) {
    long this_ts = timestamps->GetTimeStamp();
    ts_vec.push_back(this_ts);
    timeindex++;
  }
  // make binned hist
  long ts_min = ts_vec[0];
  long ts_max = ts_vec[0];
  for (long ts : ts_vec) {
    if (ts < ts_min) {
      ts_min = ts;
    }
    if (ts > ts_max) {
      ts_max = ts;
    }
  }

  delta_tmaxmin = std::abs(ts_max - ts_min);

  float nbin = 100.;
  float bin_width = (ts_max - ts_min)/nbin;
  for (int i = 0; i < 100; ++i) {
    float bin_low = ts_min + i*bin_width;
    float bin_high = ts_min + (i+1)*bin_width;
    int bin_count = 0;
    for (long ts : ts_vec) {
      if ((bin_low < ts) & (ts <= bin_high)) {
        bin_count += 1;
      }
    }
    ts_hist_vec.push_back(bin_count);
  }
  // find peak bins
  for (int bc : ts_hist_vec) {
    //TODO: make threshold fcl parameter
    if (bc > 1000) {
      npeak += 1;
    }
  }
  std::cout << "THIS EVENT HAS " << npeak << " PEAKS" << std::endl;


  // 2D hist event display
  art::Timestamp ts = evt.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  auto ttm = tts.GetSec();
  auto tlocal = std::localtime(&ttm);
  std::time_t ttt = std::mktime(tlocal);
  std::ostringstream oss;
  oss << std::put_time(std::localtime(&ttt), "%c %Z");
  //art::ServiceHandle<art::TFileService> tfs;


  int nbins[3] = {1984, 1984, 1664};
  for (int i = 0; i<6; ++i){
    int tpc = i/3;
    int plane = i%3;
    hrawadc[i] = new TH2D(Form("hrawadc%d",i),Form("Run %d Event %d TPC %d Plane %d;Wire;Tick;ADC", run, event, tpc, plane), nbins[plane], 0, nbins[plane]-1, 3000, 0, 3000);
    if (plane == 2){
         qmax = 150;
         qmin = -50;
    }
    else if (plane == 0) { 
        qmax = 60;
        qmin = -20;
    }
    else if (plane == 1) { 
        qmax = 60;
        qmin = -20;
    }
    hrawadc[i]->SetMaximum(qmax);
    hrawadc[i]->SetMinimum(qmin);
    hrawadc[i]->GetXaxis()->CenterTitle(true);
    hrawadc[i]->GetYaxis()->CenterTitle(true);
    hrawadc[i]->GetZaxis()->CenterTitle(true);
    for (int j = 1; j<=hrawadc[i]->GetNbinsX(); ++j){
      for (int k = 1; k<=hrawadc[i]->GetNbinsY(); ++k){
        hrawadc[i]->SetBinContent(j,k,-10000);
      }
    }
 //   TLatex t;
 //   t.SetNDC();
 //   t.SetTextFont(132);
 //   t.SetTextSize(0.03);
  
//    TCanvas *can = new TCanvas("can","can",padx,pady);
//    for (int i = 0; i<6; ++i){
//      int tpc = i/3;
//      int plane = i%3;
//      hrawadc[i]->Draw("colz");
  //    t.DrawLatex(0.01,0.01,Form("%04d/%02d/%02d %02d:%02d:%02d %s",
  //			       tlocal->tm_year+1900,
  //			       tlocal->tm_mon+1,
  //			       tlocal->tm_mday,
  //			       tlocal->tm_hour,
  //			       tlocal->tm_min,
  //			       tlocal->tm_sec,
  //			       tlocal->tm_isdst?"CDT":"CST"));
//      t.DrawLatex(0.01,0.01,Form("%s",oss.str().c_str()));
//      can->Print(Form("sbnd_tpc%d_plane%d.png",tpc,plane));
//      //can->Print(Form("./sbnd_tpc%d_plane%d.pdf",tpc,plane));
//      delete hrawadc[i];
//    }
  }

  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<SBND::TPCChannelMapService> channelMap;


  // Get raw digits
  std::vector<art::Ptr<raw::RawDigit>> raw_digits_handle;
  art::Handle<std::vector<raw::RawDigit>> digit_handle;
  evt.getByLabel("daq", digit_handle);
  // exit if the data isn't present
  if (!digit_handle.isValid()) {
    std::cerr << "Error: missing digits with producer (" << "daq" << ")" << std::endl;
    return;
  }
  art::fill_ptr_vector(raw_digits_handle, digit_handle);

  for (const auto & rd : raw_digits_handle){
    auto adc_vec = rd->ADCs();
    int ch = rd->Channel();
    auto const & chids = geo->ChannelToWire(ch);
    int tpc = chids[0].TPC;
    int plane = chids[0].Plane;
    int wire = chids[0].Wire;
    //cout<<tpc<<" "<<plane<<" "<<wire<<endl;
    for (unsigned short i = 0; i<adc_vec.size(); ++i){
      //std::cout << "[tpcAnalysis::OnlineEvd::analyze] adc_vec[" << i << "] : " << adc_vec[i] << ", rd->GetPedestal() : " << rd->GetPedestal() << std::endl;
      int bin = hrawadc[plane + 3*tpc]->GetBin(wire+0.1,i+0.1);
      double q = adc_vec[i] - rd->GetPedestal();
      if (plane == 2){
           qmax = 150;
           qmin = -50;
      }
      else if (plane == 0) { 
          qmax = 60;
          qmin = -20;
      }
      else if (plane == 1) { 
          qmax = 60;
          qmin = -20;
      }
      if (q>qmax) q = qmax;
      if (q<qmin) q = qmin;
      hrawadc[plane + 3*tpc]->SetBinContent(bin, q);
      evd_vec_bin[plane + 3*tpc].push_back(bin);
      evd_vec_q[plane + 3*tpc].push_back(q);
    }
  }
  if (delta_tmaxmin > 100000){
    std::cout << "Event " << event << " has a delta_tmaxmin of " << delta_tmaxmin << std::endl;
  }

  fEventTree->Fill();

}

void FlagChoppy::endJob()
{
  mf::LogVerbatim("TPCChoppyAna") << "TPCChoppyAna finished job";
}

DEFINE_ART_MODULE(FlagChoppy)

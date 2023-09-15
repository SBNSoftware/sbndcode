////////////////////////////////////////////////////////////////////////
// Class:       HitWaveformDisplay
// Plugin Type: analyzer
// File:        HitWaveformDisplay_module.cc
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

#include "TCanvas.h"
#include "TH1D.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLegend.h"

#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"

constexpr int def_int = std::numeric_limits<int>::min();

namespace sbnd {
  class HitWaveformDisplay;
}

class sbnd::HitWaveformDisplay : public art::EDAnalyzer {
public:
  explicit HitWaveformDisplay(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HitWaveformDisplay(HitWaveformDisplay const&) = delete;
  HitWaveformDisplay(HitWaveformDisplay&&) = delete;
  HitWaveformDisplay& operator=(HitWaveformDisplay const&) = delete;
  HitWaveformDisplay& operator=(HitWaveformDisplay&&) = delete;

  // Required functions.
  void analyze(const art::Event &e) override;
  void SetStyle();

private:
  art::ServiceHandle<cheat::BackTrackerService>       backTracker;
  art::InputTag fHitModuleLabel;
};

sbnd::HitWaveformDisplay::HitWaveformDisplay(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  {
    fHitModuleLabel = p.get<art::InputTag>("HitModuleLabel", "gaushit");
  }

void sbnd::HitWaveformDisplay::analyze(const art::Event &e)
{
  SetStyle();

  art::Handle<std::vector<recob::Hit>> hitHandle;
  e.getByLabel(fHitModuleLabel, hitHandle);
  if(!hitHandle.isValid()){
    std::cout << "Hit product " << fHitModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }

  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  std::vector<art::Ptr<recob::Hit>> hitVec;
  art::fill_ptr_vector(hitVec, hitHandle);

  for(auto const& hit : hitVec)
    {
      const int trackID = TruthMatchUtils::TrueParticleID(clockData,hit,true);
      if(trackID == def_int)
        {
          const art::Ptr<sim::SimChannel> sc = backTracker->FindSimChannel(hit->Channel());

          const double hitStart = hit->PeakTimeMinusRMS(1.);
          const double hitEnd   = hit->PeakTimePlusRMS(1.);

          const unsigned short hitStartTDC = std::max(0, (int) clockData.TPCTick2TDC(hitStart));
          const unsigned short hitEndTDC   = std::max(0, (int) clockData.TPCTick2TDC(hitEnd));

          auto const &map = sc->TDCIDEMap();
          unsigned short mintdc = std::numeric_limits<unsigned short>::max(), maxtdc = 0;

          for(auto const& [tdc, ides] : map)
            {
              mintdc = std::min(mintdc, tdc);
              maxtdc = std::max(maxtdc, tdc);
            }

          mintdc = std::min(mintdc, hitStartTDC);
          maxtdc = std::max(maxtdc, hitEndTDC);

          std::vector<std::pair<unsigned short, unsigned short>> extraHits;

          for(auto const &otherHit : hitVec)
            {
              if(otherHit == hit || otherHit->Channel() != hit->Channel())
                continue;

              const double otherHitStart = otherHit->PeakTimeMinusRMS(1.);
              const double otherHitEnd   = otherHit->PeakTimePlusRMS(1.);

              const unsigned short otherHitStartTDC = std::max(0, (int) clockData.TPCTick2TDC(otherHitStart));
              const unsigned short otherHitEndTDC   = std::max(0, (int) clockData.TPCTick2TDC(otherHitEnd));

              mintdc = std::min(mintdc, otherHitStartTDC);
              maxtdc = std::max(maxtdc, otherHitEndTDC);

              extraHits.emplace_back(otherHitStartTDC, otherHitEndTDC);
            }

          const unsigned short nBins = maxtdc - mintdc + 21;
          const float xLow   = mintdc - 10.5;
          const float xHigh  = maxtdc + 10.5;

          TH1D *hist = new TH1D("hist", Form("Channel %d;Tick (TDC);N Electrons", hit->Channel()), nBins, xLow, xHigh);

          for(unsigned short tdc = mintdc; tdc < maxtdc+1; ++tdc)
            {
              const unsigned short bin = tdc - mintdc + 11;
              hist->SetBinContent(bin, sc->Charge(tdc));
            }

          TCanvas *canvas = new TCanvas("canvas", "canvas");
          canvas->cd();

          const double max = hist->GetMaximum();
          hist->SetMaximum(1.2 * max);
          hist->Draw("hist");

          TLine *startLine = new TLine(hitStartTDC, 0, hitStartTDC, max);
          startLine->SetLineColor(kRed);
          startLine->SetLineWidth(4);
          TLine *endLine = new TLine(hitEndTDC, 0, hitEndTDC, max);
          endLine->SetLineColor(kRed);
          endLine->SetLineWidth(4);

          startLine->Draw();
          endLine->Draw();

          TLegend *legend = new TLegend(0.3, 0.85, 0.8, 0.9);
          legend->SetNColumns(3);
          legend->AddEntry(hist, "Sim Deposits", "l");
          legend->AddEntry(startLine, "#pm1#sigma ghost hit", "l");

          int hitN = 0;
          for(auto const& [startTDC, endTDC] : extraHits)
            {
              TLine *otherStartLine = new TLine(startTDC, 0, startTDC, max);
              otherStartLine->SetLineColor(kOrange);
              otherStartLine->SetLineWidth(4);
              TLine *otherEndLine = new TLine(endTDC, 0, endTDC, max);
              otherEndLine->SetLineColor(kOrange);
              otherEndLine->SetLineWidth(4);

              otherStartLine->Draw();
              otherEndLine->Draw();

              if(hitN == 0)
                legend->AddEntry(otherStartLine, "#pm1#sigma good hit", "l");

              ++hitN;
            }

          legend->Draw();

          canvas->SaveAs(Form("channel%d.png",hit->Channel()));
          canvas->SaveAs(Form("channel%d.pdf",hit->Channel()));

          delete canvas;
          delete hist;
        }
    }
}

void sbnd::HitWaveformDisplay::SetStyle()
{
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameLineWidth(3);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetLegendFont(42);
  gStyle->SetLegendTextSize(0.04);

  gStyle->SetPaperSize(30,50);
  gStyle->SetCanvasDefH(1000);
  gStyle->SetCanvasDefW(1700);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetPadLeftMargin(0.18);

  gStyle->SetTextFont(62);
  gStyle->SetTextSize(0.09);

  gStyle->SetLabelFont(62,"xyz");
  gStyle->SetLabelSize(0.07,"xyz");
  gStyle->SetTitleSize(0.07,"xyz");
  gStyle->SetTitleFont(62,"xyz");

  gStyle->SetTitleOffset(1.1,"x");
  gStyle->SetTitleOffset(1.1,"y");
  gStyle->SetTitleOffset(1,"z");

  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.7);
  gStyle->SetHistLineWidth(6);
  gStyle->SetLineStyleString(2,"[12 12]");

  gStyle->SetLegendBorderSize(0);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

DEFINE_ART_MODULE(sbnd::HitWaveformDisplay)

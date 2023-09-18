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
#include "TGaxis.h"
#include "TSystem.h"

#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
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
  art::InputTag fHitModuleLabel, fWireModuleLabel;

  std::string fSaveDir;

  bool fROIOnly;
};

sbnd::HitWaveformDisplay::HitWaveformDisplay(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  {
    fHitModuleLabel  = p.get<art::InputTag>("HitModuleLabel", "gaushit");
    fWireModuleLabel = p.get<art::InputTag>("WireModuleLabel", "simtpc2d:gauss");
    fSaveDir         = p.get<std::string>("SaveDir", ".");
    fROIOnly         = p.get<bool>("ROIOnly", true);
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
  std::vector<art::Ptr<recob::Hit>> hitVec;
  art::fill_ptr_vector(hitVec, hitHandle);

  art::Handle<std::vector<recob::Wire>> wireHandle;
  e.getByLabel(fWireModuleLabel, wireHandle);
  if(!wireHandle.isValid()){
    std::cout << "Wire product " << fWireModuleLabel << " not found..." << std::endl;
    throw std::exception();
  }
  std::vector<art::Ptr<recob::Wire>> wireVec;
  art::fill_ptr_vector(wireVec, wireHandle);

  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  const int run = e.id().run(), subrun = e.id().subRun(), ev = e.id().event();
  const TString saveLoc  = fSaveDir + Form("/run%dsubrun%devent%d/", run, subrun, ev);

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

          const double hitStartWide = hit->PeakTimeMinusRMS(10.);
          const double hitEndWide   = hit->PeakTimePlusRMS(10.);

          const unsigned short hitStartWideTDC = std::max(0, (int) clockData.TPCTick2TDC(hitStartWide));
          const unsigned short hitEndWideTDC   = std::max(0, (int) clockData.TPCTick2TDC(hitEndWide));

          auto const &map = sc->TDCIDEMap();
          unsigned short mintdc = hitStartWideTDC, maxtdc = hitEndWideTDC;

          if(!fROIOnly)
            {
              for(auto const& [tdc, ides] : map)
                {
                  mintdc = std::min(mintdc, tdc);
                  maxtdc = std::max(maxtdc, tdc);
                }

              mintdc = std::min(mintdc, hitStartTDC);
              maxtdc = std::max(maxtdc, hitEndTDC);
            }

          std::vector<std::pair<unsigned short, unsigned short>> extraGoodHits;
          std::vector<std::pair<unsigned short, unsigned short>> extraBadHits;

          for(auto const &wire : wireVec)
            {
              if(wire->Channel() != hit->Channel())
                continue;

              const recob::Wire::RegionsOfInterest_t& signalROI = wire->SignalROI();

              for(auto const &roi : signalROI.get_ranges())
                {
                  const unsigned short roiStart = roi.begin_index();
                  const unsigned short roiEnd   = roi.begin_index() + roi.size();

                  const unsigned short roiStartTDC = std::max(0, (int) clockData.TPCTick2TDC(roiStart));
                  const unsigned short roiEndTDC   = std::max(0, (int) clockData.TPCTick2TDC(roiEnd));

                  if(!fROIOnly || (roiStartTDC > mintdc && roiStartTDC < maxtdc) ||
                     (roiEndTDC > mintdc && roiEndTDC < maxtdc) || (roiStartTDC < mintdc && roiEndTDC > maxtdc))
                    {
                      mintdc = std::min(mintdc, roiStartTDC);
                      maxtdc = std::max(maxtdc, roiEndTDC);
                    }
                }
            }

          for(auto const &otherHit : hitVec)
            {
              if(otherHit == hit || otherHit->Channel() != hit->Channel())
                continue;

              const double otherHitStart = otherHit->PeakTimeMinusRMS(1.);
              const double otherHitEnd   = otherHit->PeakTimePlusRMS(1.);

              const unsigned short otherHitStartTDC = std::max(0, (int) clockData.TPCTick2TDC(otherHitStart));
              const unsigned short otherHitEndTDC   = std::max(0, (int) clockData.TPCTick2TDC(otherHitEnd));

              if(!fROIOnly)
                {
                  mintdc = std::min(mintdc, otherHitStartTDC);
                  maxtdc = std::max(maxtdc, otherHitEndTDC);
                }
              else
                {
                  if(otherHitStartTDC > maxtdc || otherHitEndTDC < mintdc)
                    continue;
                }

              const int otherTrackID = TruthMatchUtils::TrueParticleID(clockData,otherHit,true);
              if(otherTrackID == def_int)
                extraBadHits.emplace_back(otherHitStartTDC, otherHitEndTDC);
              else
                extraGoodHits.emplace_back(otherHitStartTDC, otherHitEndTDC);
            }

          const unsigned short nBins = maxtdc - mintdc + 21;
          const float xLow   = mintdc - 10.5;
          const float xHigh  = maxtdc + 10.5;

          TH1D *simHist  = new TH1D("simHist", Form("Channel %d;Tick (TDC);True energy deposition (MeV)", hit->Channel()), nBins, xLow, xHigh);

          for(unsigned short tdc = mintdc - 11; tdc < maxtdc+11; ++tdc)
            {
              const int bin = simHist->FindBin(tdc);
              simHist->SetBinContent(bin, sc->Energy(tdc));
            }

          if(simHist->Integral() == 0)
            continue;

          TCanvas *canvas = new TCanvas("canvas", "canvas");
          canvas->cd();

          const double max = simHist->GetMaximum();
          simHist->SetMaximum(1.5 * max);
          simHist->Draw("hist");
          simHist->GetXaxis()->SetNdivisions(505);
          simHist->GetYaxis()->SetNdivisions(507);

          TLine *startLine = new TLine(hitStartTDC, 0, hitStartTDC, 1.2 * max);
          startLine->SetLineColor(kRed);
          startLine->SetLineWidth(4);
          TLine *endLine = new TLine(hitEndTDC, 0, hitEndTDC, 1.2 * max);
          endLine->SetLineColor(kRed);
          endLine->SetLineWidth(4);

          startLine->Draw();
          endLine->Draw();

          TLegend *legend = new TLegend(0.3, 0.78, 0.8, 0.9);
          legend->SetNColumns(2);
          legend->AddEntry(simHist, "Sim Deposits", "l");
          legend->AddEntry(startLine, "#pm1#sigma ghost hit", "l");

          int hitN = 0;
          for(auto const& [startTDC, endTDC] : extraGoodHits)
            {
              TLine *otherStartLine = new TLine(startTDC, 0, startTDC, 1.2 * max);
              otherStartLine->SetLineColor(kOrange);
              otherStartLine->SetLineWidth(4);
              TLine *otherEndLine = new TLine(endTDC, 0, endTDC, 1.2 * max);
              otherEndLine->SetLineColor(kOrange);
              otherEndLine->SetLineWidth(4);

              otherStartLine->Draw();
              otherEndLine->Draw();

              if(hitN == 0)
                legend->AddEntry(otherStartLine, "#pm1#sigma good hit", "l");

              ++hitN;
            }

          for(auto const& [startTDC, endTDC] : extraBadHits)
            {
              TLine *otherStartLine = new TLine(startTDC, 0, startTDC, 1.2 * max);
              otherStartLine->SetLineColor(kRed);
              otherStartLine->SetLineWidth(4);
              TLine *otherEndLine = new TLine(endTDC, 0, endTDC, 1.2 * max);
              otherEndLine->SetLineColor(kRed);
              otherEndLine->SetLineWidth(4);

              otherStartLine->Draw();
              otherEndLine->Draw();
            }

          simHist->Draw("histsame");
          TH1D *wireHist = new TH1D("wireHist", Form("Channel %d;Tick (TDC);N Electrons", hit->Channel()), nBins, xLow, xHigh);

          for(auto const &wire : wireVec)
            {
              if(wire->Channel() != hit->Channel())
                continue;

              const recob::Wire::RegionsOfInterest_t& signalROI = wire->SignalROI();

              for(auto const &roi : signalROI.get_ranges())
                {
                  unsigned short roiStart = roi.begin_index();

                  for(auto const& value : roi)
                    {
                      unsigned short tdc = clockData.TPCTick2TDC(roiStart);
                      const int bin = simHist->FindBin(tdc);
                      wireHist->SetBinContent(bin, 50 * value);
                      ++roiStart;
                    }
                }

              const double wireMax = wireHist->GetMaximum();
              wireHist->Scale(max / wireMax);
              wireHist->Draw("same hist");
              wireHist->SetLineColor(kMagenta);

              TGaxis *wireAxis = new TGaxis(xHigh, 0, xHigh, 1.5*max,
                                            0, 1.5*wireMax, 507, "+L");
              wireAxis->SetLineWidth(1);
              wireAxis->SetLabelSize(0.05);
              wireAxis->SetTitleSize(0.05);
              wireAxis->SetTitleOffset(1);
              wireAxis->SetTitle("e^{-}");
              wireAxis->Draw();

              legend->AddEntry(wireHist, "Deconv. Waveform", "l");
            }

          legend->Draw();

          const TString fileName = fROIOnly ? Form("channel%d_hit%lu", hit->Channel(), hit.key()) : Form("channel%d", hit->Channel());

          gSystem->Exec("mkdir -p " + saveLoc);

          canvas->SaveAs(saveLoc + fileName + ".png");
          canvas->SaveAs(saveLoc + fileName + ".pdf");

          delete canvas;
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
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.12);

  gStyle->SetTextFont(62);
  gStyle->SetTextSize(0.09);

  gStyle->SetLabelFont(62,"xyz");
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleFont(62,"xyz");

  gStyle->SetTitleOffset(1.07,"x");
  gStyle->SetTitleOffset(1.12,"y");
  gStyle->SetTitleOffset(1,"z");

  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.7);
  gStyle->SetHistLineWidth(6);
  gStyle->SetLineStyleString(2,"[12 12]");

  gStyle->SetLegendBorderSize(0);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(0);
}

DEFINE_ART_MODULE(sbnd::HitWaveformDisplay)

#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

using namespace ana;

#include "HenryVars.h"
#include "HenryCuts.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"

void SetupStyle()
{
  gStyle->SetLineWidth(4);
  gStyle->SetHistLineWidth(3);
  gStyle->SetMarkerStyle(8);
  gStyle->SetMarkerSize(.5);
  gStyle->SetTitleOffset(.8, "y");

  gROOT->ForceStyle();
}                                                                                                                                                                                             
void Make_CAF_Plots()
{
  SetupStyle();

  // This is the input dataset we're going to use
  const std::string inputName = "/scratch/LAR25/caf/reco2_tutorial.flat.caf.root";

  // The directory we want to save plots in
  const TString saveDir = "./plots";

  // This object takes the dataset and fills all the 'spectra' you request from it
  SpectrumLoader loader(inputName);

  // Binning schemes we want to use for our plots
  Binning trackLengthBins = Binning::Simple(70, 0, 350);
  Binning dEdxBins        = Binning::Simple(180, 0, 30);
  Binning resRangeBins    = Binning::Simple(200, 0, 50);
  Binning opt0TimeBins    = Binning::Simple(50, 1.6, 1.61);

  // The 'spectra' we want to create from our dataset
  Spectrum *sChildTrackLength = new Spectrum("sChildTrackLength", trackLengthBins, loader, kChildTrackLengths,
                                             kNoSpillCut, kIsNeutrinoCandidateSlice);
  Spectrum *sChildTrackdEdxResRange = new Spectrum("sChildTrackdEdxResRange", loader, resRangeBins, kChildTrackResRange,
                                                   dEdxBins, kChildTrackdEdx, kNoSpillCut, kIsNeutrinoCandidateSlice);

  Spectrum *sChildTrackLengthLongestTrack = new Spectrum("sChildTrackLengthLongestTrack", trackLengthBins, loader,
                                                         kChildTrackLengthLongestTrack, kNoSpillCut, kIsNeutrinoCandidateSlice);
  Spectrum *sChildTrackLengthOtherTracks = new Spectrum("sChildTrackLengthOtherTracks", trackLengthBins, loader,
                                                        kChildTrackLengthOtherTracks, kNoSpillCut, kIsNeutrinoCandidateSlice);
  Spectrum *sChildTrackdEdxResRangeLongestTrack = new Spectrum("sChildTrackdEdxResRangeLongestTrack", loader, resRangeBins,
                                                               kChildTrackResRangeLongestTrack, dEdxBins, kChildTrackdEdxLongestTrack,
                                                               kNoSpillCut, kIsNeutrinoCandidateSlice);
  Spectrum *sChildTrackdEdxResRangeOtherTracks = new Spectrum("sChildTrackdEdxResRangeOtherTracks", loader, resRangeBins,
                                                              kChildTrackResRangeOtherTracks, dEdxBins, kChildTrackdEdxOtherTracks,
                                                              kNoSpillCut, kIsNeutrinoCandidateSlice);

  Spectrum *sOpT0Time = new Spectrum("sOpT0Time", opt0TimeBins, loader, kOpT0Time, kNoSpillCut, kIsNeutrinoCandidateSlice);

  // Tell the loader we've declared all our spectra, time to fill them!
  loader.Go();

  // Child track length plot
  TCanvas *cChildTrackLength = new TCanvas("cChildTrackLength", "cChildTrackLength");
  cChildTrackLength->cd();

  TH1D *hChildTrackLength = sChildTrackLength->ToTH1(sChildTrackLength->Livetime(), kLivetime);
  hChildTrackLength->SetTitle(";Track Length (cm);Primary Children");
  hChildTrackLength->Draw("histe");

  cChildTrackLength->SaveAs(saveDir + "/child_track_length.pdf");
  cChildTrackLength->SaveAs(saveDir + "/child_track_length.png");

  // Child dEdx vs. residual range plot
  TCanvas *cChildTrackdEdxResRange = new TCanvas("cChildTrackdEdxResRange", "cChildTrackdEdxResRange");
  cChildTrackdEdxResRange->cd();

  TH2D *hChildTrackdEdxResRange = (TH2D*) sChildTrackdEdxResRange->ToTH2(sChildTrackdEdxResRange->Livetime(), kLivetime);
  hChildTrackdEdxResRange->SetTitle(";Residual Range (cm);dE/dx (MeV/cm);Primary Children");
  hChildTrackdEdxResRange->Draw("colz");

  cChildTrackdEdxResRange->SaveAs(saveDir + "/child_track_dedx_resrange.pdf");
  cChildTrackdEdxResRange->SaveAs(saveDir + "/child_track_dedx_resrange.png");

  // Separated track length plot
  TCanvas *cChildTrackLengthSeparated = new TCanvas("cChildTrackLengthSeparated", "cChildTrackLengthSeparated");
  cChildTrackLengthSeparated->cd();

  TH1D *hChildTrackLengthLongestTrack = sChildTrackLengthLongestTrack->ToTH1(sChildTrackLengthLongestTrack->Livetime(), kLivetime);
  hChildTrackLengthLongestTrack->SetTitle(";Track Length (cm);Primary Children");
  hChildTrackLengthLongestTrack->SetLineColor(kMagenta+2);
  hChildTrackLengthLongestTrack->Draw("histe");

  TH1D *hChildTrackLengthOtherTracks = sChildTrackLengthOtherTracks->ToTH1(sChildTrackLengthOtherTracks->Livetime(), kLivetime);
  hChildTrackLengthOtherTracks->SetLineColor(kOrange+2);
  hChildTrackLengthOtherTracks->Draw("histesame");

  TLegend *lChildTrackLengthSeparated = new TLegend(.4, .65, .6, .8);
  lChildTrackLengthSeparated->AddEntry(hChildTrackLengthLongestTrack, "Longest Track", "le");
  lChildTrackLengthSeparated->AddEntry(hChildTrackLengthOtherTracks, "Other Tracks", "le");
  lChildTrackLengthSeparated->Draw();

  cChildTrackLengthSeparated->SaveAs(saveDir + "/child_track_length_separated.pdf");
  cChildTrackLengthSeparated->SaveAs(saveDir + "/child_track_length_separated.png");

  // Separated dEdx vs. residual range plot
  TCanvas *cChildTrackdEdxResRangeSeparated = new TCanvas("cChildTrackdEdxResRangeSeparated", "cChildTrackdEdxResRangeSeparated");
  cChildTrackdEdxResRangeSeparated->cd();

  TH2D *hChildTrackdEdxResRangeLongestTrack = (TH2D*) sChildTrackdEdxResRangeLongestTrack->ToTH2(sChildTrackdEdxResRangeLongestTrack->Livetime(), kLivetime);
  hChildTrackdEdxResRangeLongestTrack->SetTitle(";Residual Range (cm);dE/dx (MeV/cm);Primary Children");
  hChildTrackdEdxResRangeLongestTrack->SetMarkerColor(kMagenta+2);
  hChildTrackdEdxResRangeLongestTrack->Draw();

  TH2D *hChildTrackdEdxResRangeOtherTracks = (TH2D*) sChildTrackdEdxResRangeOtherTracks->ToTH2(sChildTrackdEdxResRangeOtherTracks->Livetime(), kLivetime);
  hChildTrackdEdxResRangeOtherTracks->SetTitle(";Residual Range (cm);dE/dx (MeV/cm);Primary Children");
  hChildTrackdEdxResRangeOtherTracks->SetMarkerColor(kOrange+2);
  hChildTrackdEdxResRangeOtherTracks->Draw("same");

  TLegend *lChildTrackdEdxResRangeSeparated = new TLegend(.6, .65, .8, .8);
  lChildTrackdEdxResRangeSeparated->AddEntry(hChildTrackdEdxResRangeLongestTrack, "Longest Track", "p");
  lChildTrackdEdxResRangeSeparated->AddEntry(hChildTrackdEdxResRangeOtherTracks, "Other Tracks", "p");
  lChildTrackdEdxResRangeSeparated->Draw();

  cChildTrackdEdxResRangeSeparated->SaveAs(saveDir + "/child_track_dedx_resrange_separated.pdf");
  cChildTrackdEdxResRangeSeparated->SaveAs(saveDir + "/child_track_dedx_resrange_separated.png");

  // Child track length plot
  TCanvas *cOpT0Time = new TCanvas("cOpT0Time", "cOpT0Time");
  cOpT0Time->cd();

  TH1D *hOpT0Time = sOpT0Time->ToTH1(sOpT0Time->Livetime(), kLivetime);
  hOpT0Time->SetTitle(";OpT0 Matched T0 (#mu s);Slices");
  hOpT0Time->SetLineColor(kGreen+2);
  hOpT0Time->Draw("histe");

  cOpT0Time->SaveAs(saveDir + "/opt0_time.pdf");
  cOpT0Time->SaveAs(saveDir + "/opt0_time.png");
}

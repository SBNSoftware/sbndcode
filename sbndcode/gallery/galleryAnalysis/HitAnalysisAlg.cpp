
#include "HitAnalysisAlg.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <cmath>
#include <algorithm>

namespace HitAnalysis
{
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
HitAnalysisAlg::HitAnalysisAlg(fhicl::ParameterSet const & pset) :
    fRootDirectory(0)

{
    
    reconfigure(pset);
    
    // Report.
    mf::LogInfo("HitAnalysisAlg") << "HitAnalysisAlg configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
HitAnalysisAlg::~HitAnalysisAlg()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void HitAnalysisAlg::reconfigure(fhicl::ParameterSet const & pset)
{
    fLocalDirName = pset.get<std::string>("LocalDirName", std::string("wow"));
}

//----------------------------------------------------------------------------
/// Begin job method.
void HitAnalysisAlg::setup(const geo::GeometryCore&           geometry,
                           const detinfo::DetectorProperties& detectorProperties,
                           TDirectory*                        rootDirectory)
{
    // Get geometry and detector properties
    fGeometry           = &geometry;
    fDetectorProperties = &detectorProperties;
    fRootDirectory      = rootDirectory->mkdir(fLocalDirName.c_str());

    // Make a directory for these histograms
//    art::TFileDirectory dir = tfs->mkdir(fLocalDirName.c_str());

    fHitsByWire[0]            = std::make_unique<TH1D>("HitsByWire0", ";Wire #", fGeometry->Nwires(0), 0., fGeometry->Nwires(0));
    fHitsByWire[1]            = std::make_unique<TH1D>("HitsByWire1", ";Wire #", fGeometry->Nwires(1), 0., fGeometry->Nwires(1));
    fHitsByWire[2]            = std::make_unique<TH1D>("HitsByWire2", ";Wire #", fGeometry->Nwires(2), 0., fGeometry->Nwires(2));
    
    fDriftTimes[0]            = std::make_unique<TH1D>("DriftTime0",  ";time(ticks)", 3200, 0., 9600.);
    fDriftTimes[1]            = std::make_unique<TH1D>("DriftTime1",  ";time(ticks)", 3200, 0., 9600.);
    fDriftTimes[2]            = std::make_unique<TH1D>("DriftTime2",  ";time(ticks)", 3200, 0., 9600.);
    
    fHitsByTime[0]            = std::make_unique<TH1D>("HitsByTime0", ";Tick",   1600, 0., 6400.);
    fHitsByTime[1]            = std::make_unique<TH1D>("HitsByTime1", ";Tick",   1600, 0., 6400.);
    fHitsByTime[2]            = std::make_unique<TH1D>("HitsByTime2", ";Tick",   1600, 0., 6400.);
    
    fPulseHeight[0]           = std::make_unique<TH1D>("PulseHeight0",  "PH (ADC)",  300,  0.,  150.);
    fPulseHeight[1]           = std::make_unique<TH1D>("PulseHeight1",  "PH (ADC)",  300,  0.,  150.);
    fPulseHeight[2]           = std::make_unique<TH1D>("PulseHeight2",  "PH (ADC)",  300,  0.,  150.);
    fPulseHeightSingle[0]     = std::make_unique<TH1D>("PulseHeightS0", "PH (ADC)",  300,  0.,  150.);
    fPulseHeightSingle[1]     = std::make_unique<TH1D>("PulseHeightS1", "PH (ADC)",  300,  0.,  150.);
    fPulseHeightSingle[2]     = std::make_unique<TH1D>("PulseHeightS2", "PH (ADC)",  300,  0.,  150.);
    fPulseHeightMulti[0]      = std::make_unique<TH1D>("PulseHeightM0", "PH (ADC)",  300,  0.,  150.);
    fPulseHeightMulti[1]      = std::make_unique<TH1D>("PulseHeightM1", "PH (ADC)",  300,  0.,  150.);
    fPulseHeightMulti[2]      = std::make_unique<TH1D>("PulseHeightM2", "PH (ADC)",  300,  0.,  150.);
    fChi2DOF[0]               = std::make_unique<TH1D>("Chi2DOF0",      "Chi2DOF",   502, -1.,  250.);
    fChi2DOF[1]               = std::make_unique<TH1D>("Chi2DOF1",      "Chi2DOF",   502, -1.,  250.);
    fChi2DOF[2]               = std::make_unique<TH1D>("Chi2DOF2",      "Chi2DOF",   502, -1.,  250.);
    fNumDegFree[0]            = std::make_unique<TH1D>("NumDegFree0",   "NDF",       100,  0.,  100.);
    fNumDegFree[1]            = std::make_unique<TH1D>("NumDegFree1",   "NDF",       100,  0.,  100.);
    fNumDegFree[2]            = std::make_unique<TH1D>("NumDegFree2",   "NDF",       100,  0.,  100.);
    fChi2DOFSingle[0]         = std::make_unique<TH1D>("Chi2DOFS0",     "Chi2DOF",   502, -1.,  250.);
    fChi2DOFSingle[1]         = std::make_unique<TH1D>("Chi2DOFS1",     "Chi2DOF",   502, -1.,  250.);
    fChi2DOFSingle[2]         = std::make_unique<TH1D>("Chi2DOFS2",     "Chi2DOF",   502, -1.,  250.);
    fHitMult[0]               = std::make_unique<TH1D>("HitMult0",      "# hits",     15,  0.,   15.);
    fHitMult[1]               = std::make_unique<TH1D>("HitMult1",      "# hits",     15,  0.,   15.);
    fHitMult[2]               = std::make_unique<TH1D>("HitMult2",      "# hits",     15,  0.,   15.);
    fHitCharge[0]             = std::make_unique<TH1D>("HitCharge0",    "Charge",   1000,  0., 2000.);
    fHitCharge[1]             = std::make_unique<TH1D>("HitCharge1",    "Charge",   1000,  0., 2000.);
    fHitCharge[2]             = std::make_unique<TH1D>("HitCharge2",    "Charge",   1000,  0., 2000.);
    fFitWidth[0]              = std::make_unique<TH1D>("FitWidth0",     "Width",     100,  0.,   10.);
    fFitWidth[1]              = std::make_unique<TH1D>("FitWidth1",     "Width",     100,  0.,   10.);
    fFitWidth[2]              = std::make_unique<TH1D>("FitWidth2",     "Width",     100,  0.,   10.);
    fHitSumADC[0]             = std::make_unique<TH1D>("SumADC0",       "Sum ADC",  1000,  0., 2000.);
    fHitSumADC[1]             = std::make_unique<TH1D>("SumADC1",       "Sum ADC",  1000,  0., 2000.);
    fHitSumADC[2]             = std::make_unique<TH1D>("SumADC2",       "Sum ADC",  1000,  0., 2000.);
    
    fNDFVsChi2[0]             = std::make_unique<TH2D>("NDFVsChi20",    ";NDF;Chi2",  50,  0.,   50., 101, -1., 100.);
    fNDFVsChi2[1]             = std::make_unique<TH2D>("NDFVsChi21",    ";NDF;Chi2",  50,  0.,   50., 101, -1., 100.);
    fNDFVsChi2[2]             = std::make_unique<TH2D>("NDFVsChi22",    ";NDF;Chi2",  50,  0.,   50., 101, -1., 100.);
    
    fPulseHVsWidth[0]         = std::make_unique<TH2D>("PHVsWidth0",    ";PH;Width", 100,  0.,  100., 100,  0., 20.);
    fPulseHVsWidth[1]         = std::make_unique<TH2D>("PHVsWidth1",    ";PH;Width", 100,  0.,  100., 100,  0., 20.);
    fPulseHVsWidth[2]         = std::make_unique<TH2D>("PHVsWidth2",    ";PH;Width", 100,  0.,  100., 100,  0., 20.);
    
    fPulseHVsCharge[0]        = std::make_unique<TH2D>("PHVsChrg0",     ";PH;Q",     100,  0.,  100., 100,  0., 2000.);
    fPulseHVsCharge[1]        = std::make_unique<TH2D>("PHVsChrg1",     ";PH;Q",     100,  0.,  100., 100,  0., 2000.);
    fPulseHVsCharge[2]        = std::make_unique<TH2D>("PHVsChrg2",     ";PH;Q",     100,  0.,  100., 100,  0., 2000.);
    
    fPulseHVsHitNo[0]         = std::make_unique<TProfile>("PHVsNo0",   ";Hit #;PH", 1000, 0., 1000., 0., 100.);
    fPulseHVsHitNo[1]         = std::make_unique<TProfile>("PHVsNo1",   ";Hit #;PH", 1000, 0., 1000., 0., 100.);
    fPulseHVsHitNo[2]         = std::make_unique<TProfile>("PHVsNo2",   ";Hit #;PH", 1000, 0., 1000., 0., 100.);
    
    fChargeVsHitNo[0]         = std::make_unique<TProfile>("QVsNo0",    ";Hit No;Q", 1000, 0., 1000., 0., 2000.);
    fChargeVsHitNo[1]         = std::make_unique<TProfile>("QVsNo1",    ";Hit No;Q", 1000, 0., 1000., 0., 2000.);
    fChargeVsHitNo[2]         = std::make_unique<TProfile>("QVsNo2",    ";Hit No;Q", 1000, 0., 1000., 0., 2000.);
    
    fChargeVsHitNoS[0]        = std::make_unique<TProfile>("QVsNoS0",   ";Hit No;Q", 1000, 0., 1000., 0., 2000.);
    fChargeVsHitNoS[1]        = std::make_unique<TProfile>("QVsNoS1",   ";Hit No;Q", 1000, 0., 1000., 0., 2000.);
    fChargeVsHitNoS[2]        = std::make_unique<TProfile>("QVsNoS2",   ";Hit No;Q", 1000, 0., 1000., 0., 2000.);
    
    fBadWPulseHeight          = std::make_unique<TH1D>("BWPulseHeight", "PH (ADC)",  300,  0.,  150.);
    fBadWPulseHVsWidth        = std::make_unique<TH2D>("BWPHVsWidth",   ";PH;Width", 100,  0.,  100., 100,  0., 10.);
    fBadWHitsByWire           = std::make_unique<TH1D>("BWHitsByWire",  ";Wire #", fGeometry->Nwires(2), 0., fGeometry->Nwires(2));
    
    fSPHvsIdx[0]              = std::make_unique<TH2D>("SPHVsIdx0",     ";PH;Idx", 30,  0.,  30., 100,  0., 100.);
    fSPHvsIdx[1]              = std::make_unique<TH2D>("SPHVsIdx1",     ";PH;Idx", 30,  0.,  30., 100,  0., 100.);
    fSPHvsIdx[2]              = std::make_unique<TH2D>("SPHVsIdx2",     ";PH;Idx", 30,  0.,  30., 100,  0., 100.);
    
    fSWidVsIdx[0]             = std::make_unique<TH2D>("SWidsIdx0",     ";Width;Idx", 30,  0.,  30., 100,  0., 10.);
    fSWidVsIdx[1]             = std::make_unique<TH2D>("SWidsIdx1",     ";Width;Idx", 30,  0.,  30., 100,  0., 10.);
    fSWidVsIdx[2]             = std::make_unique<TH2D>("SWidsIdx2",     ";Width;Idx", 30,  0.,  30., 100,  0., 10.);
    
    f1PPHvsWid[0]             = std::make_unique<TH2D>("1PPHVsWid0",    ";PH;Width", 100,  0.,  100., 100,  0., 20.);
    f1PPHvsWid[1]             = std::make_unique<TH2D>("1PPHVsWid1",    ";PH;Width", 100,  0.,  100., 100,  0., 20.);
    f1PPHvsWid[2]             = std::make_unique<TH2D>("1PPHVsWid2",    ";PH;Width", 100,  0.,  100., 100,  0., 20.);
    
    fSPPHvsWid[0]             = std::make_unique<TH2D>("SPPHVsWid0",    ";PH;Width", 100,  0.,  100., 100,  0., 20.);
    fSPPHvsWid[1]             = std::make_unique<TH2D>("SPPHVsWid1",    ";PH;Width", 100,  0.,  100., 100,  0., 20.);
    fSPPHvsWid[2]             = std::make_unique<TH2D>("SPPHVsWid2",    ";PH;Width", 100,  0.,  100., 100,  0., 20.);
    
    fSOPHvsWid[0]             = std::make_unique<TH2D>("SOPHVsWid0",    ";PH;Width", 100,  0.,  100., 100,  0., 20.);
    fSOPHvsWid[1]             = std::make_unique<TH2D>("SOPHVsWid1",    ";PH;Width", 100,  0.,  100., 100,  0., 20.);
    fSOPHvsWid[2]             = std::make_unique<TH2D>("SOPHVsWid2",    ";PH;Width", 100,  0.,  100., 100,  0., 20.);
    
    fPHRatVsIdx[0]            = std::make_unique<TH2D>("PHRatVsIdx0",   ";PHRat;Idx", 30,  0.,  30., 51,  0., 1.02);
    fPHRatVsIdx[1]            = std::make_unique<TH2D>("PHRatVsIdx1",   ";PHRat;Idx", 30,  0.,  30., 51,  0., 1.02);
    fPHRatVsIdx[2]            = std::make_unique<TH2D>("PHRatVsIdx2",   ";PHRat;Idx", 30,  0.,  30., 51,  0., 1.02);

    // Now attach them to the output file
    for(int idx = 0; idx < 3; idx++)
    {
        fHitsByWire[idx]->SetDirectory(fRootDirectory);
        fDriftTimes[idx]->SetDirectory(fRootDirectory);
        fHitsByTime[idx]->SetDirectory(fRootDirectory);
        fPulseHeight[idx]->SetDirectory(fRootDirectory);
        fPulseHeightSingle[idx]->SetDirectory(fRootDirectory);
        fPulseHeightMulti[idx]->SetDirectory(fRootDirectory);
        fChi2DOF[idx]->SetDirectory(fRootDirectory);
        fNumDegFree[idx]->SetDirectory(fRootDirectory);
        fChi2DOFSingle[idx]->SetDirectory(fRootDirectory);
        fHitMult[idx]->SetDirectory(fRootDirectory);
        fHitCharge[idx]->SetDirectory(fRootDirectory);
        fFitWidth[idx]->SetDirectory(fRootDirectory);
        fHitSumADC[idx]->SetDirectory(fRootDirectory);
        fNDFVsChi2[idx]->SetDirectory(fRootDirectory);
        fPulseHVsWidth[idx]->SetDirectory(fRootDirectory);
        fPulseHVsCharge[idx]->SetDirectory(fRootDirectory);
        fPulseHVsHitNo[idx]->SetDirectory(fRootDirectory);
        fChargeVsHitNo[idx]->SetDirectory(fRootDirectory);
        fChargeVsHitNoS[idx]->SetDirectory(fRootDirectory);
        
        fSPHvsIdx[idx]->SetDirectory(fRootDirectory);
        fSWidVsIdx[idx]->SetDirectory(fRootDirectory);
        f1PPHvsWid[idx]->SetDirectory(fRootDirectory);
        fSPPHvsWid[idx]->SetDirectory(fRootDirectory);
        fSOPHvsWid[idx]->SetDirectory(fRootDirectory);
        fPHRatVsIdx[idx]->SetDirectory(fRootDirectory);
    }
    
    // get the overall now
    fBadWPulseHeight->SetDirectory(fRootDirectory);
    fBadWPulseHVsWidth->SetDirectory(fRootDirectory);
    fBadWHitsByWire->SetDirectory(fRootDirectory);

    return;
}
    
void HitAnalysisAlg::fillHistograms(const TrackPlaneHitMap& trackPlaneHitMap) const
{
    int    longTrackID(0);
    size_t longTrackLen(0);
    
    for(const auto& trackHitVecMapItr : trackPlaneHitMap)
    {
        size_t numHits(0);
        
        for(const auto& planeHitPair : trackHitVecMapItr.second)
        {
            fillHistograms(planeHitPair.second);
            
            if (planeHitPair.first == 2 && planeHitPair.second.size() > numHits) numHits = planeHitPair.second.size();
        }
        
        if (numHits > longTrackLen)
        {
            longTrackID  = trackHitVecMapItr.first;
            longTrackLen = numHits;
        }
    }
    
    if (longTrackLen > 0)
    {
        for(const auto& planeHitPair : trackPlaneHitMap.find(longTrackID)->second)
        {
            int hitNo(0);
        
            for(const auto& hit : planeHitPair.second)
            {
                if (hit.Multiplicity() < 2) fChargeVsHitNoS[planeHitPair.first]->Fill(float(hitNo)+0.5, std::min(float(1999.),hit.Integral()), 1.);
                fPulseHVsHitNo[planeHitPair.first]->Fill(float(hitNo)+0.5, std::min(float(99.9),hit.PeakAmplitude()), 1.);
                fChargeVsHitNo[planeHitPair.first]->Fill(float(hitNo)+0.5, std::min(float(1999.),hit.Integral()), 1.);
                hitNo++;
            }
        }
    }
    
    return;
}
    
void HitAnalysisAlg::fillHistograms(const HitVec& hitVec) const
{
    // Keep track of number of hits per plane
    size_t nHitsPerPlane[] = {0,0,0};
    size_t negCount(0);
    
    std::vector<const recob::Hit*> hitSnippetVec;
    
    // Loop the hits and make some plots
    for(const auto& hit : hitVec)
    {
        // Extract interesting hit parameters
        const geo::WireID& wireID   = hit.WireID();
        float              chi2DOF  = std::min(hit.GoodnessOfFit(),float(249.8));
        int                numDOF   = hit.DegreesOfFreedom();
        int                hitMult  = hit.Multiplicity();
        float              peakTime = hit.PeakTime();
        float              charge   = hit.Integral();
        float              sumADC   = hit.SummedADC();
        float              hitPH    = std::min(hit.PeakAmplitude(),float(249.8));
        float              hitSigma = hit.RMS();
        
        size_t             plane    = wireID.Plane;
        size_t             wire     = wireID.Wire;
        
//        if (plane == 2 && (wire == 18 || wire == 527 || wire == 528)) continue;
        
//        if (charge < 0. || sumADC < 0. || hitPH < 0.)
        if (hitPH < 0.)
        {
            negCount++;
            std::cout << "Hit plane: " << plane << ", wire: " << wire << ", T: " << peakTime << ", PH: " << hitPH << ", charge: " << charge << ", sumADC: " << sumADC << std::endl;
        }
        
        nHitsPerPlane[plane]++;
        
        fHitsByWire[plane]->Fill(wire,1.);
        fHitsByTime[plane]->Fill(peakTime, 1.);
        fPulseHeight[plane]->Fill(hitPH, 1.);
        fChi2DOF[plane]->Fill(chi2DOF, 1.);
        fNumDegFree[plane]->Fill(numDOF, 1.);
        fHitMult[plane]->Fill(hitMult, 1.);
        fHitCharge[plane]->Fill(charge, 1.);
        fFitWidth[plane]->Fill(std::min(float(19.99),hitSigma), 1.);
        fHitSumADC[plane]->Fill(sumADC, 1.);
        fNDFVsChi2[plane]->Fill(numDOF, chi2DOF, 1.);
        fDriftTimes[plane]->Fill(peakTime, 1.);
        
        if (hitMult == 1)
        {
            fPulseHeightSingle[plane]->Fill(hitPH, 1.);
            fChi2DOFSingle[plane]->Fill(chi2DOF, 1.);
            fPulseHVsWidth[plane]->Fill(std::min(float(99.9),hitPH), std::min(float(19.99),hitSigma), 1.);
            fPulseHVsCharge[plane]->Fill(std::min(float(99.9),hitPH), std::min(float(1999.),charge), 1.);
            
            if (plane == 2 && hitPH < 5 && hitSigma < 2.2)
            {
                std::cout << "++> plane: " << plane << ", wire: " << wire << ", peakTime: " << peakTime << ", ph: " << hitPH << ", w: " << hitSigma << std::endl;
                
                fBadWPulseHeight->Fill(hitPH,1.);
                fBadWPulseHVsWidth->Fill(std::min(float(99.9),hitPH), std::min(float(19.99),hitSigma), 1.);
                fBadWHitsByWire->Fill(wire,1.);
            }
        }
        else
            fPulseHeightMulti[plane]->Fill(hitPH, 1.);
        
        // Look at hits on snippets
        if (!hitSnippetVec.empty() && hitSnippetVec.back()->LocalIndex() >= hit.LocalIndex())
        {
            // Only worried about multi hit snippets
            if (hitSnippetVec.size() > 1)
            {
                // Sort in order of largest to smallest pulse height
                std::sort(hitSnippetVec.begin(),hitSnippetVec.end(),[](const auto* left, const auto* right){return left->PeakAmplitude() > right->PeakAmplitude();});
                
                float maxPulseHeight = hitSnippetVec.front()->PeakAmplitude();
                
                for(size_t idx = 0; idx < hitSnippetVec.size(); idx++)
                {
                    float pulseHeight      = hitSnippetVec.at(idx)->PeakAmplitude();
                    float pulseWid         = hitSnippetVec.at(idx)->RMS();
                    float pulseHeightRatio = pulseHeight / maxPulseHeight;
                    
                    size_t plane = hitSnippetVec.at(idx)->WireID().Plane;
                    
                    fSPHvsIdx[plane]->Fill(idx, std::min(float(99.9),pulseHeight), 1.);
                    fSWidVsIdx[plane]->Fill(idx, std::min(float(19.99),pulseWid), 1.);
                    fPHRatVsIdx[plane]->Fill(idx, pulseHeightRatio, 1.);
                    
                    if (idx == 0) fSPPHvsWid[plane]->Fill(std::min(float(99.9),pulseHeight), std::min(float(19.99),pulseWid), 1.);
                    else          fSOPHvsWid[plane]->Fill(std::min(float(99.9),pulseHeight), std::min(float(19.99),pulseWid), 1.);
                }
            }
            else
            {
                float  pulseHeight = hitSnippetVec.front()->PeakAmplitude();
                float  pulseWid    = hitSnippetVec.front()->RMS();
                size_t plane        = hitSnippetVec.front()->WireID().Plane;
                
                f1PPHvsWid[plane]->Fill(std::min(float(99.9),pulseHeight), std::min(float(19.99),pulseWid), 1.);
            }
            
            hitSnippetVec.clear();
        }
        
        hitSnippetVec.push_back(&hit);
    }
    
    return;
}
    
// Useful for normalizing histograms
void HitAnalysisAlg::endJob(int numEvents)
{
    // Normalize wire profiles to be hits/event
    double normFactor(1./numEvents);
    
    for(size_t idx = 0; idx < 3; idx++) fHitsByWire[idx]->Scale(normFactor);
    
    // Now write out everything...
    fRootDirectory->cd();
    
    // Start by looping over the planes
    for(int idx = 0; idx < 3; idx++)
    {
        fHitsByWire[idx]->Write();
        fDriftTimes[idx]->Write();
        fHitsByTime[idx]->Write();
        fPulseHeight[idx]->Write();
        fPulseHeightSingle[idx]->Write();
        fPulseHeightMulti[idx]->Write();
        fChi2DOF[idx]->Write();
        fNumDegFree[idx]->Write();
        fChi2DOFSingle[idx]->Write();
        fHitMult[idx]->Write();
        fHitCharge[idx]->Write();
        fFitWidth[idx]->Write();
        fHitSumADC[idx]->Write();
        fNDFVsChi2[idx]->Write();
        fPulseHVsWidth[idx]->Write();
        fPulseHVsCharge[idx]->Write();
        fPulseHVsHitNo[idx]->Write();
        fChargeVsHitNo[idx]->Write();
        fChargeVsHitNoS[idx]->Write();
        
        fSPHvsIdx[idx]->Write();
        fSWidVsIdx[idx]->Write();
        f1PPHvsWid[idx]->Write();
        fSPPHvsWid[idx]->Write();
        fSOPHvsWid[idx]->Write();
        fPHRatVsIdx[idx]->Write();
    }
    
    // get the overall now
    fBadWPulseHeight->Write();
    fBadWPulseHVsWidth->Write();
    fBadWHitsByWire->Write();

    return;
}

}

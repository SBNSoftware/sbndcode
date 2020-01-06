#! /usr/bin/env python
import sys, os
# Prevent root from printing garbage on initialization.
if os.environ.has_key('TERM'):
    del os.environ['TERM']

# Hide command line arguments from ROOT module.
myargv = sys.argv
sys.argv = myargv[0:1]
#sys.argv.append( '-b' )

from ROOT import TFile, TCanvas, TH1F, TH2F, TProfile, TLegend, TF1
from ROOT import gDirectory, gROOT, gPad, gApplication
#ROOT.gErrorIgnoreLevel = ROOT.kError
sys.argv = myargv

#from validation_utilities import *

# Get the object with a given name from a list of objects in the directory
def GetObject(name,list):
    for i in list:
        if i.GetName()==name:
            return i.ReadObj()
    return 0

# Print help.

def makePlots(inputHitFileName):
    can      = {}
    canIdx   = 0

    fitTrackPH = {}
    allHitPH   = {}
    
    ###########################################################################################################
    # Open all the input root files.
    hitTFile = TFile(inputHitFileName)
    list1    = hitTFile.GetListOfKeys()

    print '==> In plotHits, using file name {0:s}'.format(inputHitFileName)
    
    # Go to directory for the TrackHitAna output
    for i in list1:
        if i.GetClassName() == 'TDirectoryFile':
            if i.GetName() != "mcassociations":
                continue
            hitTFile.cd(i.GetName())
            histogramList = gDirectory.GetListOfKeys()
            can[canIdx] = TCanvas("MC vs Reco", "MC vs Reco", 800, 600)
            can[canIdx].Divide(2,2)
            nHitsPerPrimary  = GetObject("NHitsPerPrimary",  histogramList)
            primaryRecoNHits = GetObject("PrimaryRecoNHits", histogramList)
            deltaNHits       = GetObject("DeltaNHits",       histogramList)
            primaryLength    = GetObject("PrimaryLength",    histogramList)
            primaryRecLength = GetObject("PrimaryRecLength", histogramList)
            deltaTrackLen    = GetObject("DeltaTrackLen",    histogramList)
            can[canIdx].cd(1)
            nHitsPerPrimary.SetLineColor(4)
            nHitsPerPrimary.SetStats(0)
            nHitsPerPrimary.Draw()
            primaryRecoNHits.SetLineColor(2)
            primaryRecoNHits.SetStats(0)
            primaryRecoNHits.Draw("SAMES")
            can[canIdx].cd(2)
            deltaNHits.SetLineColor(4)
            deltaNHits.SetStats(0)
            deltaNHits.Draw()
            can[canIdx].cd(3)
            primaryLength.SetLineColor(4)
            primaryLength.SetStats(0)
            primaryLength.Draw()
            primaryRecLength.SetLineColor(2)
            primaryRecLength.SetStats(0)
            primaryRecLength.Draw("SAMES")
            can[canIdx].cd(4)
            deltaTrackLen.SetLineColor(4)
            deltaTrackLen.SetStats(0)
            deltaTrackLen.Draw()
            can[canIdx].Update()
            
            canIdx += 1
    
            can[canIdx] = TCanvas("Completeness", "Completeness", 800, 600)
            can[canIdx].Divide(2,2)
            primCompleteness = GetObject("PrimaryCompleteness",  histogramList)
            primCompVsNHits  = GetObject("PrimaryCompVsLogHits", histogramList)
            primCompVsLen    = GetObject("PrimaryCompVsLen",     histogramList)
            primCompVsMom    = GetObject("PrimaryCompVsMom",     histogramList)
            can[canIdx].cd(1)
            primCompleteness.SetLineColor(4)
            primCompleteness.SetStats(0)
            primCompleteness.Draw()
            can[canIdx].cd(2)
            primCompVsNHits.SetLineColor(4)
            primCompVsNHits.SetStats(0)
            primCompVsNHits.Draw()
            can[canIdx].cd(3)
            primCompVsLen.SetLineColor(4)
            primCompVsLen.SetStats(0)
            primCompVsLen.Draw()
            can[canIdx].cd(4)
            primCompVsMom.SetLineColor(4)
            primCompVsMom.SetStats(0)
            primCompVsMom.Draw()
            can[canIdx].Update()
            
            canIdx += 1
                
            can[canIdx] = TCanvas("Efficiency", "Efficiency", 800, 600)
            can[canIdx].Divide(2,2)
            primEfficiency = GetObject("PrimaryEfficiency",   histogramList)
            primEffVsNHits = GetObject("PrimaryEffVsLogHits", histogramList)
            primEffVsLen   = GetObject("PrimaryEffVsLen",     histogramList)
            primEffVsMom   = GetObject("PrimaryEffVsMom",     histogramList)
            can[canIdx].cd(1)
            primEfficiency.SetLineColor(4)
            primEfficiency.SetStats(0)
            primEfficiency.Draw()
            can[canIdx].cd(2)
            primEffVsNHits.SetLineColor(4)
            primEffVsNHits.SetStats(0)
            primEffVsNHits.Draw()
            can[canIdx].cd(3)
            primEffVsLen.SetLineColor(4)
            primEffVsLen.SetStats(0)
            primEffVsLen.Draw()
            can[canIdx].cd(4)
            primEffVsMom.SetLineColor(4)
            primEffVsMom.SetStats(0)
            primEffVsMom.Draw()
            can[canIdx].Update()

            canIdx += 1

            can[canIdx] = TCanvas("Purity", "Purity", 800, 600)
            can[canIdx].Divide(2,2)
            primPurity     = GetObject("PrimaryPurity",          histogramList)
            primPurVsNHits = GetObject("PrimaryPurityVsLogHits", histogramList)
            primPurVsLen   = GetObject("PrimaryPurVsLen",        histogramList)
            primPurVsMom   = GetObject("PrimaryPurVsMom",        histogramList)
            can[canIdx].cd(1)
            primPurity.SetLineColor(4)
            primPurity.SetStats(0)
            primPurity.Draw()
            can[canIdx].cd(2)
            primPurVsNHits.SetLineColor(4)
            primPurVsNHits.SetStats(0)
            primPurVsNHits.Draw()
            can[canIdx].cd(3)
            primPurVsLen.SetLineColor(4)
            primPurVsLen.SetStats(0)
            primPurVsLen.Draw()
            can[canIdx].cd(4)
            primPurVsMom.SetLineColor(4)
            primPurVsMom.SetStats(0)
            primPurVsMom.Draw()
            can[canIdx].Update()

            canIdx += 1


    gApplication.Run()


def main(argv):
    infile='TrackAnalysis.root'
    outdir= ''
    release=''
    calorimetry = 0
    hit = 0
    tracking = 0
    momresol = 0
    flash = 0
    pid = 0

    args = argv[1:]

    if len(args) > 0 :
        if args[0] == '--input' and len(args) > 1 :
            infile = args[1]

    makePlots(infile)

    currdir = os.getcwd()
    if outdir != currdir:
    	os.chdir(currdir) 	    
	    	  
if __name__ == '__main__':
    rc = main(sys.argv)
    sys.exit(rc)

#! /usr/bin/env python
######################################################################
#
# Name: generate_simple_weighted_template.py
#
# Purpose: Create templates
#
# Created: 14-Dec-2020  Iker Loïc de Icaza Astiz (icaza@fnal.gov)
#
# Usage:
#
# stat.py <options> [@filelist] [file1 file2 ...]
#
# Options:
#
# [-h|--help] - Print help message.
# --level n   - Branch level (default 1).  Use --level 1 to see top
#               branches only.  Use --level 2 to also see subbranches.
# --nfile n   - Number of files to analyze (default all).
# --all       - Print analysis of each file (default is only summary).
#
# Arguments:
#
# @filelist       - File list containing one input file per line.
# file1 file2 ... - Input files.
#
######################################################################

import sys, os, string
from array import array
# import project_utilities
import larbatch_posix

# Import ROOT module.
# Globally turn off root warnings.
# Don't let root see our command line options.

myargv = sys.argv
sys.argv = myargv[0:1]
if 'TERM' in os.environ:
    del os.environ['TERM']
import ROOT
from ROOT import TH1F, TH2F, TProfile, TFile
from ROOT import TStyle, TCanvas, TGraph, TGraphErrors
ROOT.gErrorIgnoreLevel = ROOT.kError
sys.argv = myargv

# Print help.
def help():

    filename = sys.argv[0]
    file = open(filename)

    doprint = False

    for line in file.readlines():
        if line[2:9] == 'stat.py':
            doprint = True
        elif line[0:6] == '######' and doprint:
            doprint = False
        if doprint:
            if len(line) > 2:
                print(line[2:], end=' ')
            else:
                print()


def generator(input_file, rootfile, gtrees, gbranches):
    dir = rootfile.Get(input_file+":/fmatch")
    nuslice_tree = dir.Get("nuslicetree")#, nuslice_tree)
    # nuslice_tree.Print()
    n_bins = 20
    xvals = list(range(5, 205, 10)) # TODO: un hardcode
    xerrs = [5] * len(xvals)
    dist_to_anode_bins =  100
    dist_to_anode_low = 0. # TODO: un hardcode
    dist_to_anode_up = 200. # TODO: un hardcode
    profile_bins = 40
    profile_option = 's' # errors are the standard deviation

    dy_spreads = [None] * n_bins
    dy_means = [None] * n_bins
    y_bins =  100
    y_low = -200. # TODO: un hardcode
    y_up = 200. # TODO: un hardcode
    dy_hist = TH2F("dy_hist", "#Delta y ",
               dist_to_anode_bins, dist_to_anode_low, dist_to_anode_up,
               y_bins, y_low, y_up)
    dy_hist.GetXaxis().SetTitle("distance from anode (cm)");
    dy_hist.GetYaxis().SetTitle("y_flash - y_TPC (cm)");
    dy_prof = TProfile("dy_prof","Profile of dy_spreads in #Delta y",
                       profile_bins, dist_to_anode_low, dist_to_anode_up,
                       y_low, y_up, profile_option);
    dy_prof.GetXaxis().SetTitle("distance from anode (cm)");
    dy_prof.GetYaxis().SetTitle("y_flash - y_TPC (cm)");
    dy_h1 = TH1F("dy_h1", "",
                 profile_bins, dist_to_anode_low, dist_to_anode_up);
    dy_h1.GetXaxis().SetTitle("distance from anode (cm)");
    dy_h1.GetYaxis().SetTitle("y_flash - y_TPC (cm)");

    dz_spreads = [None] * n_bins
    dz_means = [None] * n_bins
    z_bins =  100
    z_low = -200. # TODO: un hardcode
    z_up = 200. # TODO: un hardcode
    dz_hist = TH2F("dz_hist", "#Delta z ",
               dist_to_anode_bins, dist_to_anode_low, dist_to_anode_up,
               z_bins, z_low, z_up)
    dz_hist.GetXaxis().SetTitle("distance from anode (cm)");
    dz_hist.GetYaxis().SetTitle("z_flash - z_TPC (cm)");
    dz_prof = TProfile("dz_prof","Profile of dz_spreads in #Delta z",
                       profile_bins, dist_to_anode_low, dist_to_anode_up,
                       z_low, z_up, profile_option);
    dz_prof.GetXaxis().SetTitle("distance from anode (cm)");
    dz_prof.GetYaxis().SetTitle("z_flash - z_TPC (cm)");
    dz_h1 = TH1F("dz_h1", "",
                 profile_bins, dist_to_anode_low, dist_to_anode_up);
    dz_h1.GetXaxis().SetTitle("distance from anode (cm)");
    dz_h1.GetYaxis().SetTitle("z_flash - z_TPC (cm)");

    rr_spreads = [None] * n_bins
    rr_means = [None] * n_bins
    rr_bins =  100
    rr_low = 0. # TODO: un hardcode
    rr_up = 200. # TODO: un hardcode
    rr_hist = TH2F("rr_hist", "RMS Flash",
               dist_to_anode_bins, dist_to_anode_low, dist_to_anode_up,
               rr_bins, rr_low, rr_up)
    rr_hist.GetXaxis().SetTitle("distance from anode (cm)");
    rr_hist.GetYaxis().SetTitle("RMS flash (cm)");
    rr_prof = TProfile("rr_prof","Profile of dy_spreads in #Delta y", # TODO: title
                       profile_bins, dist_to_anode_low, dist_to_anode_up,
                       rr_low, rr_up, profile_option);
    rr_prof.GetXaxis().SetTitle("distance from anode (cm)");
    rr_prof.GetYaxis().SetTitle("y_flash - y_TPC (cm)"); # TODO: title
    rr_h1 = TH1F("rr_h1", "",
                 profile_bins, dist_to_anode_low, dist_to_anode_up);
    rr_h1.GetXaxis().SetTitle("distance from anode (cm)");
    rr_h1.GetYaxis().SetTitle("y_flash - y_TPC (cm)"); # TODO: title

    pe_spreads = [None] * n_bins
    pe_means = [None] * n_bins
    pe_bins =  100
    pe_low = 0. # TODO: un hardcode
    pe_up = 50. # TODO: un hardcode
    pe_hist = TH2F("pe_hist", "#Delta z ",
               dist_to_anode_bins, dist_to_anode_low, dist_to_anode_up,
               pe_bins, pe_low, pe_up)
    pe_hist.GetXaxis().SetTitle("distance from anode (cm)");
    pe_hist.GetYaxis().SetTitle("z_flash - z_TPC (cm)");
    pe_prof = TProfile("pe_prof","Profile of dz_spreads in #Delta z",
                       profile_bins, dist_to_anode_low, dist_to_anode_up,
                       pe_low, pe_up, profile_option);
    pe_prof.GetXaxis().SetTitle("distance from anode (cm)");
    pe_prof.GetYaxis().SetTitle("z_flash - z_TPC (cm)");
    pe_h1 = TH1F("pe_h1", "",
                 profile_bins, dist_to_anode_low, dist_to_anode_up);
    pe_h1.GetXaxis().SetTitle("distance from anode (cm)");
    pe_h1.GetYaxis().SetTitle("z_flash - z_TPC (cm)");

    score_hist_bins = 100
    score_hist_low = 0.
    score_hist_up = 50.
    match_score_scatter = TH2F("match_score_scatter", "Scatter plot of match scores",
                               dist_to_anode_bins, dist_to_anode_low, dist_to_anode_up,
                               score_hist_bins, score_hist_low, score_hist_up*(3./5.))
    match_score_scatter.GetXaxis().SetTitle("distance from anode (cm)")
    match_score_scatter.GetYaxis().SetTitle("match score (arbitrary)");
    match_score_hist = TH1F("match_score", "Match Score",
                            score_hist_bins, score_hist_low, score_hist_up)
    match_score_hist.GetXaxis().SetTitle("match score (arbitrary)")

    for e in nuslice_tree:
        slice = e.nuvtx_x
        uncoated_coated_ratio = 100.*e.flash_unpe/e.flashpe

        dy_hist.Fill(slice, e.flash_y - e.nuvtx_y)
        dy_prof.Fill(slice, e.flash_y - e.nuvtx_y)
        dz_hist.Fill(slice, e.flash_z - e.nuvtx_z)
        dz_prof.Fill(slice, e.flash_z - e.nuvtx_z)
        rr_hist.Fill(slice, e.flash_r)
        rr_prof.Fill(slice, e.flash_r)
        pe_hist.Fill(slice, uncoated_coated_ratio)
        pe_prof.Fill(slice, uncoated_coated_ratio)

    # fill histograms for match score calculation from profile histograms
    for ib in list(range(0, profile_bins)):
        dy_h1.SetBinContent(ib, dy_prof.GetBinContent(ib))
        dy_h1.SetBinError(ib, dy_prof.GetBinError(ib))
        dy_means[int(ib/2)] = dy_prof.GetBinContent(ib)
        dy_spreads[int(ib/2)] = dy_prof.GetBinError(ib)
        dz_h1.SetBinContent(ib, dz_prof.GetBinContent(ib))
        dz_h1.SetBinError(ib, dz_prof.GetBinError(ib))
        dz_means[int(ib/2)] = dz_prof.GetBinContent(ib)
        dz_spreads[int(ib/2)] = dz_prof.GetBinError(ib)
        rr_h1.SetBinContent(ib, rr_prof.GetBinContent(ib))
        rr_h1.SetBinError(ib, rr_prof.GetBinError(ib))
        rr_means[int(ib/2)] = rr_prof.GetBinContent(ib)
        rr_spreads[int(ib/2)] = rr_prof.GetBinError(ib)
        pe_h1.SetBinContent(ib, pe_prof.GetBinContent(ib))
        pe_h1.SetBinError(ib, pe_prof.GetBinError(ib))
        pe_means[int(ib/2)] = pe_prof.GetBinContent(ib)
        pe_spreads[int(ib/2)] = pe_prof.GetBinError(ib)

    for e in nuslice_tree:
        slice = e.nuvtx_x
        uncoated_coated_ratio = 100.*e.flash_unpe/e.flashpe
        # calculate match score
        isl = int(slice/10.)
        # if (isl>) isl=19;
        # if (isl<0) isl=0;

        score = 0.
        if dy_spreads[isl] <= 1e-8:
            print(f"isl {isl}. dy_spreads[isl]{dy_spreads[isl]} ")
            dy_spreads[isl] = 1.
        if dz_spreads[isl] <= 1e-8:
            print(f"isl {isl}. dz_spreads[isl]{dz_spreads[isl]} ")
            dz_spreads[isl] = 1.
        if rr_spreads[isl] <= 1e-8:
            print(f"isl {isl}. rr_spreads[isl]{rr_spreads[isl]} ")
            rr_spreads[isl] = 1.
        if pe_spreads[isl] <= 1e-8:
            print(f"isl {isl}. pe_spreads[isl]{pe_spreads[isl]} ")
            pe_spreads[isl] = 1.
        score += abs(abs(e.flash_y-e.nuvtx_y)- dy_means[isl])/dy_spreads[isl]
        score += abs(abs(e.flash_z-e.nuvtx_z)- dz_means[isl])/dz_spreads[isl]
        score += abs(e.flash_r-rr_means[isl])/rr_spreads[isl]
        # score += abs(uncoated_coated_ratio-pe_means[isl])/pe_spreads[isl]
        match_score_scatter.Fill(slice, score)
        match_score_hist.Fill(score)

    hfile = ROOT.gROOT.FindObject('fm_metrics_sbnd.root')
    if hfile:
        hfile.Close()
    hfile = TFile('fm_metrics_sbnd.root', 'RECREATE', 'Simple flash matching metrics for SBND')
    dy_hist.Write();
    dy_prof.Write();
    dy_h1.Write();
    dz_hist.Write();
    dz_prof.Write();
    dz_h1.Write();
    rr_hist.Write();
    rr_prof.Write();
    rr_h1.Write();
    pe_hist.Write();
    pe_prof.Write();
    pe_h1.Write();
    match_score_scatter.Write();
    match_score_hist.Write();
    hfile.Close();

    canv = TCanvas("canv");

    canv.Update();
    dy_hist.Draw();
    crosses = TGraphErrors(n_bins, array('f',xvals), array('f',dy_means), array('f',xerrs), array('f',dy_spreads))
    crosses.SetLineColor(9)
    crosses.SetLineWidth(3)
    crosses.Draw("Psame")
    canv.Print("dy.pdf")

    canv.Update();
    dz_hist.Draw();
    crosses = TGraphErrors(n_bins, array('f',xvals), array('f',dz_means), array('f',xerrs), array('f',dz_spreads))
    crosses.SetLineColor(9)
    crosses.SetLineWidth(3)
    crosses.Draw("Psame")
    canv.Print("dz.pdf")

    canv.Update();
    rr_hist.Draw();
    crosses = TGraphErrors(n_bins, array('f',xvals), array('f',rr_means), array('f',xerrs), array('f',rr_spreads))
    crosses.SetLineColor(9)
    crosses.SetLineWidth(3)
    crosses.Draw("Psame")
    canv.Print("rr.pdf")

    canv.Update();
    pe_hist.Draw();
    crosses = TGraphErrors(n_bins, array('f',xvals), array('f',pe_means), array('f',xerrs), array('f',pe_spreads))
    crosses.SetLineColor(9)
    crosses.SetLineWidth(3)
    crosses.Draw("Psame")
    canv.Print("pe.pdf")

    canv.Update()
    match_score_scatter.Draw()
    canv.Print("match_score_scatter.pdf")

    canv.Update()
    match_score_hist.Draw()
    canv.Print("match_score.pdf")


# Main program.
def main(argv):

    # Parse arguments.
    input_files = []
    # level = 1
    nfilemax = 1
    # all = 0

    args = argv[1:]
    while len(args) > 0:
        if args[0] == '-h' or args[0] == '--help':
            help()
            return 0
        elif args[0][0] == '-':
            # Unknown option.
            print('Unknown option %s' % args[0])
            return 1

        elif args[0][0] == '@':
            # Read in file list to input files.
            filelistname = args[0][1:]
            if larbatch_posix.exists(filelistname):
                for filename in larbatch_posix.readlines(filelistname):
                    input_files.append(string.strip(filename))
            else:
                print('File list %s does not exist.' % filelistname)
                return 1
            del args[0]
        else:

            # Add single file to input files.
            input_files.append(args[0])
            del args[0]

    # Loop over input files.

    gtrees = {}
    gbranches = {}
    nfile = 0

    for input_file in input_files:

        if nfilemax > 0 and nfile >= nfilemax:
            break
        nfile = nfile + 1

        if not larbatch_posix.exists(input_file):
            print('Input file %s does not exist.' % input_file)
            return 1

        print('\nOpening %s' % input_file)
        rootfile = TFile.Open(input_file)
        if not rootfile.IsOpen() or rootfile.IsZombie():
            print('Failed to open %s' % input_file)
            return 1

        generator(input_file, rootfile, gtrees, gbranches)

    print('\n%d files analyzed.' % nfile)


if __name__ == '__main__':
    sys.exit(main(sys.argv))

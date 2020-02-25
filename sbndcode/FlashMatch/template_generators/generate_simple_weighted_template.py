#! /usr/bin/env python
######################################################################
#
# Name: generate_simple_weighted_template.py
#
# Purpose: Create templates
#
# Created: February-2020  Iker Loïc de Icaza Astiz (icaza@fnal.gov)
#
# Usage:
#
# generate_simple_weighted_template.py file
#
# Options:
#
# [-h|--help] - Print help message.
#
# Arguments:
#
# file  ... - Input file.
#
######################################################################

import sys
import os
import string
import numpy as np
from time import sleep
from array import array

from ROOT import TStyle, TCanvas, TGraph, TGraphErrors
from ROOT import TH1D, TH2D, TProfile, TFile
import ROOT
# import project_utilities
import larbatch_posix
import fhicl


class dotDict(dict):
    def __getattr__(self,val):
        return self[val]


# Globally turn off root warnings.
# Don't let root see our command line options.
myargv = sys.argv
sys.argv = myargv[0:1]
if 'TERM' in os.environ:
    del os.environ['TERM']
ROOT.gErrorIgnoreLevel = ROOT.kError
sys.argv = myargv

# # Print help
# def help():

#     filename = sys.argv[0]
#     file = open(filename)

#     doprint = False

#     for line in file.readlines():
#         if line[2:9] == 'stat.py':
#             doprint = True
#         elif line[0:6] == '######' and doprint:
#             doprint = False
#         if doprint:
#             if len(line) > 2:
#                 print(line[2:], end=' ')
#             else:
#                 print()


def generator(input_file, rootfile, pset):
    dir = rootfile.Get(input_file+":/fmatch")
    nuslice_tree = dir.Get("nuslicetree")  # , nuslice_tree)
    # nuslice_tree.Print()

    n_bins = pset.n_bins
    drift_distance = pset.DriftDistance
    bin_width = drift_distance/n_bins
    half_bin_width = bin_width/2.

    xvals = np.arange(half_bin_width, drift_distance, bin_width)
    xerrs = np.array([half_bin_width] * len(xvals))
    dist_to_anode_bins = n_bins
    dist_to_anode_low = 0.
    dist_to_anode_up = drift_distance
    profile_bins = n_bins
    profile_option = 's'  # errors are the standard deviation

    dy_spreads = [None] * n_bins
    dy_means = [None] * n_bins
    dy_hist = TH2D("dy_hist", "#Delta y",
                   dist_to_anode_bins, dist_to_anode_low, dist_to_anode_up,
                   pset.dy_bins, pset.dy_low, pset.dy_up)
    dy_hist.GetXaxis().SetTitle("distance from anode (cm)")
    dy_hist.GetYaxis().SetTitle("y_flash - y_TPC (cm)")
    dy_prof = TProfile("dy_prof", "Profile of dy_spreads in #Delta y",
                       profile_bins, dist_to_anode_low, dist_to_anode_up,
                       pset.dy_low*2, pset.dy_up*2, profile_option)
    dy_prof.GetXaxis().SetTitle("distance from anode (cm)")
    dy_prof.GetYaxis().SetTitle("y_flash - y_TPC (cm)")
    dy_h1 = TH1D("dy_h1", "",
                 profile_bins, dist_to_anode_low, dist_to_anode_up)
    dy_h1.GetXaxis().SetTitle("distance from anode (cm)")
    dy_h1.GetYaxis().SetTitle("y_flash - y_TPC (cm)")

    dz_spreads = [None] * n_bins
    dz_means = [None] * n_bins
    dz_hist = TH2D("dz_hist", "#Delta z",
                   dist_to_anode_bins, dist_to_anode_low, dist_to_anode_up,
                   pset.dz_bins, pset.dz_low, pset.dz_up)
    dz_hist.GetXaxis().SetTitle("distance from anode (cm)")
    dz_hist.GetYaxis().SetTitle("z_flash - z_TPC (cm)")
    dz_prof = TProfile("dz_prof", "Profile of dz_spreads in #Delta z",
                       profile_bins, dist_to_anode_low, dist_to_anode_up,
                       pset.dz_low*2.5, pset.dz_up*2.5, profile_option)
    dz_prof.GetXaxis().SetTitle("distance from anode (cm)")
    dz_prof.GetYaxis().SetTitle("z_flash - z_TPC (cm)")
    dz_h1 = TH1D("dz_h1", "",
                 profile_bins, dist_to_anode_low, dist_to_anode_up)
    dz_h1.GetXaxis().SetTitle("distance from anode (cm)")
    dz_h1.GetYaxis().SetTitle("z_flash - z_TPC (cm)")

    rr_spreads = [None] * n_bins
    rr_means = [None] * n_bins
    rr_hist = TH2D("rr_hist", "PE Spread",
                   dist_to_anode_bins, dist_to_anode_low, dist_to_anode_up,
                   pset.rr_bins, pset.rr_low, pset.rr_up)
    rr_hist.GetXaxis().SetTitle("distance from anode (cm)")
    rr_hist.GetYaxis().SetTitle("RMS flash (cm)")
    rr_prof = TProfile("rr_prof", "Profile of PE Spread",
                       profile_bins, dist_to_anode_low, dist_to_anode_up,
                       pset.rr_low, pset.rr_up, profile_option)
    rr_prof.GetXaxis().SetTitle("distance from anode (cm)")
    rr_prof.GetYaxis().SetTitle("RMS flash (cm)")
    rr_h1 = TH1D("rr_h1", "",
                 profile_bins, dist_to_anode_low, dist_to_anode_up)
    rr_h1.GetXaxis().SetTitle("distance from anode (cm)")
    rr_h1.GetYaxis().SetTitle("RMS flash (cm)")

    pe_spreads = [None] * n_bins
    pe_means = [None] * n_bins
    pe_hist = TH2D("pe_hist", "Uncoated/Coated Ratio",
                   dist_to_anode_bins, dist_to_anode_low, dist_to_anode_up,
                   pset.pe_bins, pset.pe_low, pset.pe_up)
    pe_hist.GetXaxis().SetTitle("distance from anode (cm)")
    pe_hist.GetYaxis().SetTitle("ratio_{uncoated/coated}")
    pe_prof = TProfile("pe_prof", "Profile of Uncoated/Coated Ratio",
                       profile_bins, dist_to_anode_low, dist_to_anode_up,
                       pset.pe_low, pset.pe_up, profile_option)
    pe_prof.GetXaxis().SetTitle("distance from anode (cm)")
    pe_prof.GetYaxis().SetTitle("ratio_{uncoated/coated}")
    pe_h1 = TH1D("pe_h1", "",
                 profile_bins, dist_to_anode_low, dist_to_anode_up)
    pe_h1.GetXaxis().SetTitle("distance from anode (cm)")
    pe_h1.GetYaxis().SetTitle("ratio_{uncoated/coated}")

    match_score_scatter = TH2D("match_score_scatter", "Scatter plot of match scores",
                               dist_to_anode_bins, dist_to_anode_low, dist_to_anode_up,
                               pset.score_hist_bins, pset.score_hist_low, pset.score_hist_up*(3./5.))
    match_score_scatter.GetXaxis().SetTitle("distance from anode (cm)")
    match_score_scatter.GetYaxis().SetTitle("match score (arbitrary)")
    match_score_hist = TH1D("match_score", "Match Score",
                            pset.score_hist_bins, pset.score_hist_low, pset.score_hist_up)
    match_score_hist.GetXaxis().SetTitle("match score (arbitrary)")

    for e in nuslice_tree:
        slice = e.charge_x
        uncoated_coated_ratio = 100.*e.flash_unpe/e.flash_pe # percentage

        dy_hist.Fill(slice, e.flash_y - e.charge_y)
        dy_prof.Fill(slice, e.flash_y - e.charge_y)
        dz_hist.Fill(slice, e.flash_z - e.charge_z)
        dz_prof.Fill(slice, e.flash_z - e.charge_z)
        rr_hist.Fill(slice, e.flash_r)
        rr_prof.Fill(slice, e.flash_r)
        pe_hist.Fill(slice, uncoated_coated_ratio)
        pe_prof.Fill(slice, uncoated_coated_ratio)

    # fill histograms for match score calculation from profile histograms
    for ib in list(range(0, profile_bins)):
        ibp = ib + 1
        dy_h1.SetBinContent(ibp, dy_prof.GetBinContent(ibp))
        dy_h1.SetBinError(ibp, dy_prof.GetBinError(ibp))
        dy_means[int(ib)] = dy_prof.GetBinContent(ibp)
        dy_spreads[int(ib)] = dy_prof.GetBinError(ibp)
        dz_h1.SetBinContent(ibp, dz_prof.GetBinContent(ibp))
        dz_h1.SetBinError(ibp, dz_prof.GetBinError(ibp))
        dz_means[int(ib)] = dz_prof.GetBinContent(ibp)
        dz_spreads[int(ib)] = dz_prof.GetBinError(ibp)
        rr_h1.SetBinContent(ibp, rr_prof.GetBinContent(ibp))
        rr_h1.SetBinError(ibp, rr_prof.GetBinError(ibp))
        rr_means[int(ib)] = rr_prof.GetBinContent(ibp)
        rr_spreads[int(ib)] = rr_prof.GetBinError(ibp)
        pe_h1.SetBinContent(ibp, pe_prof.GetBinContent(ibp))
        pe_h1.SetBinError(ibp, pe_prof.GetBinError(ibp))
        pe_means[int(ib)] = pe_prof.GetBinContent(ibp)
        pe_spreads[int(ib)] = pe_prof.GetBinError(ibp)

    for e in nuslice_tree:
        slice = e.charge_x
        uncoated_coated_ratio = 100.*e.flash_unpe/e.flash_pe # percentage
        # calculate match score
        isl = int(slice/bin_width)
        score = 0.
        if dy_spreads[isl] <= 1.e-8:
            print("Warning zero spread.\n",
                  f"slice: {slice}. isl: {isl}. dy_spreads[isl]: {dy_spreads[isl]} ")
            dy_spreads[isl] = dy_spreads[isl+1]
        if dz_spreads[isl] <= 1.e-8:
            print("Warning zero spread.\n",
                  f"slice: {slice}. isl: {isl}. dz_spreads[isl]: {dz_spreads[isl]} ")
            dz_spreads[isl] = dz_spreads[isl+1]
        if rr_spreads[isl] <= 1.e-8:
            print("Warning zero spread.\n",
                  f"slice: {slice}. isl: {isl}. rr_spreads[isl]: {rr_spreads[isl]} ")
            rr_spreads[isl] = rr_spreads[isl+1]
        if pe_spreads[isl] <= 1.e-8:
            print("Warning zero spread.\n",
                  f"slice: {slice}. isl: {isl}. pe_spreads[isl]: {pe_spreads[isl]} ")
            pe_spreads[isl] = pe_spreads[isl+1]
        score += abs(abs(e.flash_y-e.charge_y) - dy_means[isl])/dy_spreads[isl]
        score += abs(abs(e.flash_z-e.charge_z) - dz_means[isl])/dz_spreads[isl]
        score += abs(e.flash_r-rr_means[isl])/rr_spreads[isl]
        # score += abs(uncoated_coated_ratio-pe_means[isl])/pe_spreads[isl]
        match_score_scatter.Fill(slice, score)
        match_score_hist.Fill(score)
    hfile = ROOT.gROOT.FindObject('fm_metrics_sbnd.root')
    if hfile:
        hfile.Close()
    hfile = TFile('fm_metrics_sbnd.root', 'RECREATE',
                  'Simple flash matching metrics for SBND')
    dy_hist.Write()
    dy_prof.Write()
    dy_h1.Write()
    dz_hist.Write()
    dz_prof.Write()
    dz_h1.Write()
    rr_hist.Write()
    rr_prof.Write()
    rr_h1.Write()
    pe_hist.Write()
    pe_prof.Write()
    pe_h1.Write()
    match_score_scatter.Write()
    match_score_hist.Write()
    hfile.Close()

    canv = TCanvas("canv")

    dy_hist.Draw()
    crosses = TGraphErrors(n_bins,
                           array('f', xvals), array('f', dy_means),
                           array('f', xerrs), array('f', dy_spreads))
    crosses.SetLineColor(9)
    crosses.SetLineWidth(3)
    crosses.Draw("Psame")
    canv.Print("dy.pdf")
    canv.Update()

    dz_hist.Draw()
    crosses = TGraphErrors(n_bins,
                           array('f', xvals), array('f', dz_means),
                           array('f', xerrs), array('f', dz_spreads))
    crosses.SetLineColor(9)
    crosses.SetLineWidth(3)
    crosses.Draw("Psame")
    canv.Print("dz.pdf")
    canv.Update()

    rr_hist.Draw()
    crosses = TGraphErrors(n_bins,
                           array('f', xvals), array('f', rr_means),
                           array('f', xerrs), array('f', rr_spreads))
    crosses.SetLineColor(9)
    crosses.SetLineWidth(3)
    crosses.Draw("Psame")
    canv.Print("rr.pdf")
    canv.Update()

    pe_hist.Draw()
    crosses = TGraphErrors(n_bins,
                           array('f', xvals), array('f', pe_means),
                           array('f', xerrs), array('f', pe_spreads))
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
    canv.Update()
    sleep(20)

# Main program.
def main(argv):

    # Parse arguments.
    input_files = []
    nfilemax = 1

    args = argv[1:]
    while len(args) > 0:
        if args[0] == '-h' or args[0] == '--help':
            # help()
            print("To run do:\n"/
                  "generate_simple_weighted_template.py file.root\n"/
                  "where file.root has a fmatch/nuslicetree")
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

        fcl_params = fhicl.make_pset('flashmatch_sbnd.fcl')
        pset = dotDict(fcl_params['sbnd_simple_flashmatch'])
        generator(input_file, rootfile, pset)

    print('\n%d files analyzed.' % nfile)


if __name__ == '__main__':
    sys.exit(main(sys.argv))

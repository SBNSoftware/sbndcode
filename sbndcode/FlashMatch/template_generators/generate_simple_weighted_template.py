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
    # level = 2
    # trees = {}
    # events = None
    # keys = rootfile.GetListOfKeys()
    # for key in keys:
    #     objname = key.GetName()
    #     print(objname)
    #     if objname not in trees:
    #         obj = rootfile.Get(objname)
    #         print("obj:{}".format(obj))
    #         if obj and obj.InheritsFrom('TTree'):
    #             trees[objname] = obj
    #             print("trees:{}".format(trees))

    #             if objname == 'Events':
    #                 events = obj

    # doprint = True
    # # Print summary of trees.
    # if doprint:
    #     print('\nTrees:\n')
    # for key in sorted(trees.keys()):
    #     tree = trees[key]
    #     nentry = tree.GetEntriesFast()
    #     if doprint:
    #         print('%s has %d entries.' % (key, nentry))

    #     # Remember information about trees.

    #     if key in gtrees:
    #         gtrees[key] = gtrees[key] + nentry
    #     else:
    #         gtrees[key] = nentry

    # # Print summary of branches in Events tree.

    # if doprint:
    #     print('\nBranches of Events tree:\n')

    # # If level is zero, we are done (don't analyze branches).

    # # if level == 0:
    # #     return

    # if events:

    #     if doprint:
    #         print('   Total bytes  Zipped bytes   Comp.  Branch name')
    #         print('   -----------  ------------   -----  -----------')

    #     branches = events.GetListOfBranches()
    #     ntotall = 0
    #     nzipall = 0

    #     # Loop over branche of Events tree.

    #     for branch in sorted(branches):
    #         branch_class = branch.GetClass().GetName()

    #         # Only look at data products (class art::Wrapper<T>).

    #         if branch_class[0: 13] == 'art::Wrapper<':

    #             # Loop over subbranches.

    #             subbranches = branch.GetListOfBranches()
    #             for subbranch in sorted(subbranches):
    #                 name = subbranch.GetName()

    #                 # Only look at '.obj' subbranch (wrapped object).

    #                 if name[-4:] == '.obj':
    #                     ntot = subbranch.GetTotBytes("*")
    #                     nzip = subbranch.GetZipBytes("*")
    #                     ntotall = ntotall + ntot
    #                     nzipall = nzipall + nzip
    #                     if doprint:
    #                         if nzip != 0:
    #                             comp = float(ntot) / float(nzip)
    #                         else:
    #                             comp = 0.
    #                         print('%14d%14d%8.2f  %s' % (ntot, nzip, comp, name))

    #                     # Remember information about branches.
                        
    #                     if name in gbranches:
    #                         gbranches[name][0] = gbranches[name][0] + ntot
    #                         gbranches[name][1] = gbranches[name][1] + nzip
    #                     else:
    #                         gbranches[name] = [ntot, nzip]

    #                     # Loop over subsubbranches (attributes of wrapped object).
                        
    #                     if level > 1:
    #                         subsubbranches = subbranch.GetListOfBranches()
    #                         for subsubbranch in sorted(subsubbranches):
    #                             name = subsubbranch.GetName()
    #                             ntot = subsubbranch.GetTotBytes("*")
    #                             nzip = subsubbranch.GetZipBytes("*")
    #                             if doprint:
    #                                 if nzip != 0:
    #                                     comp = float(ntot) / float(nzip)
    #                                 else:
    #                                     comp = 0.
    #                                 print('%14d%14d%8.2f  %s' % (ntot, nzip, comp,
    #                                                              subsubbranch.GetName()))

    #                             # Remember information about branches.
                        
    #                             if name in gbranches:
    #                                 gbranches[name][0] = gbranches[name][0] + ntot
    #                                 gbranches[name][1] = gbranches[name][1] + nzip
    #                             else:
    #                                 gbranches[name] = [ntot, nzip]

    # # Do summary of all branches.


    
    dir = rootfile.Get(input_file+":/fmatch")
    # print("dir: "+str(type(dir)))
    # print(dir)
    # TTree nuslice_tree  = None
    nuslice_tree = dir.Get("nuslicetree")#, nuslice_tree)
    # nuslice_tree.Print()
    n_bins = 20

    y_spreads = [None] * n_bins
    y_means = [None] * n_bins

    xvals = list(range(5, 205, 10)) # TODO: un hardcode
    xerrs = [5] * len(xvals)

    dist_to_anode_bins =  100
    dist_to_anode_low = 0. # TODO: un hardcode
    dist_to_anode_up = 200. # TODO: un hardcode

    y_bins =  100
    y_low = -200. # TODO: un hardcode
    y_up = 200. # TODO: un hardcode

    dy_hist = TH2F("dy_hist", "#Delta y ",
               dist_to_anode_bins, dist_to_anode_low, dist_to_anode_up,
               y_bins, y_low, y_up)
    dy_hist.GetXaxis().SetTitle("distance from anode (cm)");
    dy_hist.GetYaxis().SetTitle("y_flash - y_TPC (cm)");

    # TODO: the rest of TH2Fs
    score_hist_bins = 100
    score_hist_low = 0.
    score_hist_up = 50.
    match_score_scatter = TH2F("match_score_scatter", "Scatter plot of match scores",
                               dist_to_anode_bins, dist_to_anode_low, dist_to_anode_up,
                               score_hist_bins, score_hist_low, score_hist_up)
    match_score_scatter.GetXaxis().SetTitle("distance from anode (cm)")
    match_score_scatter.GetYaxis().SetTitle("match score (arbitrary)");
    match_score_hist = TH1F("match_score", "Match Score",
                            score_hist_bins, score_hist_low, score_hist_up)
    match_score_hist.GetXaxis().SetTitle("match score (arbitrary)")

    profile_bins = 40
    profile_option = 's' # errors are the standard deviation
    dy_prof = TProfile("dy_prof","Profile of y_spreads in #Delta y",
                       profile_bins, dist_to_anode_low, dist_to_anode_up,
                       y_low, y_up, profile_option);
    dy_prof.GetXaxis().SetTitle("distance from anode (cm)");
    dy_prof.GetYaxis().SetTitle("y_flash - y_TPC (cm)");

    # TODO: the rest of TProfiles

    dy_h1 = TH1F("dy_h1", "",
                 profile_bins, dist_to_anode_low, dist_to_anode_up);
    dy_h1.GetXaxis().SetTitle("distance from anode (cm)");
    dy_h1.GetYaxis().SetTitle("y_flash - y_TPC (cm)");

    # exit(0)
    # counter = 1
    for e in nuslice_tree:
        # print("Going in nuslice_tree for", counter)
        # counter += 1
        slice = e.nuvtx_x

        dy_hist.Fill(slice, e.flash_y - e.nuvtx_y)
        dy_prof.Fill(slice, e.flash_y - e.nuvtx_y)

        # fill histograms for match score calculation from profile histograms
        for ib in list(range(0, profile_bins)):
            dy_h1.SetBinContent(ib, dy_prof.GetBinContent(ib))
            dy_h1.SetBinError(ib, dy_prof.GetBinError(ib))
            y_means[int(ib/2)] = 0.
            y_spreads[int(ib/2)] = dy_prof.GetBinError(ib)
            #TODO: fill rest of d*_h1 hists

    for e in nuslice_tree:
        slice = e.nuvtx_x
        # calculate match score
        isl = int(slice/10.)
        # if (isl>) isl=19;
        # if (isl<0) isl=0;

        score = 0.
        if y_spreads[isl] <= 1e-8:
            print(f"isl {isl}. y_spreads[isl]{y_spreads[isl]} ")
            y_spreads[isl] = 1.
        score += abs(abs(e.flash_y-e.nuvtx_y)- y_means[isl])/y_spreads[isl]
        # score += abs(abs(e.flash_z-e.nuvtx_z)-dzmean[isl])/dzsp[isl]
        # score += abs(flash_r-rrmean[isl])/rrsp[isl]
        # score += abs(myratio-pemean[isl])/pesp[isl]
        match_score_hist.Fill(score)

    hfile = ROOT.gROOT.FindObject('fm_metrics_sbnd.root')
    if hfile:
        hfile.Close()
    hfile = TFile('fm_metrics_sbnd.root', 'RECREATE', 'Simple flash matching metrics for SBND')
    match_score_scatter.Write();
    match_score_hist.Write();
    dy_hist.Write();
    dy_prof.Write();
    dy_h1.Write();

    hfile.Close();

    canv = TCanvas("canv");
    canv.Update();
    dy_hist.Draw();
    crosses = TGraphErrors(n_bins, array('f',xvals), array('f',y_means), array('f',xerrs), array('f',y_spreads))
    crosses.SetLineColor(9)
    crosses.SetLineWidth(3)
    crosses.Draw("Psame")
    canv.Print("dy.pdf")
    canv.Update()

    match_score_scatter.Draw()
    canv.Print("match_score_scatter.pdf")
    canv.Update()

    match_score_hist.Draw()
    canv.Print("match_score.pdf")
    canv.Update()



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

        # elif args[0] == '--level' and len(args) > 1:
        #     # Analyze level.
        #     level = int(args[1])
        #     del args[0:2]

        # elif args[0] == '--nfile' and len(args) > 1:
        #     # Number of files.
        #     nfilemax = int(args[1])
        #     del args[0:2]

        # elif args[0] == '--all':
            # All files flag.
            # all = 1
            # del args[0]

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

        # Analyze this file.
        # analyze(root, level, gtrees, gbranches, all)
        generator(input_file, rootfile, gtrees, gbranches)

    print('\n%d files analyzed.' % nfile)


if __name__ == '__main__':
    sys.exit(main(sys.argv))

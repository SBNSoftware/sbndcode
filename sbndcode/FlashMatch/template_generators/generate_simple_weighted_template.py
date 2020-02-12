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
# import project_utilities
# import larbatch_posix

# Import ROOT module.
# Globally turn off root warnings.
# Don't let root see our command line options.

myargv = sys.argv
sys.argv = myargv[0:1]
if 'TERM' in os.environ: 
    del os.environ['TERM']
import ROOT
ROOT.gErrorIgnoreLevel = ROOT.kError
sys.argv = myargv

# Print help.

def help():

    filename = sys.argv[0]
    file = open(filename)

    doprint=0

    for line in file.readlines():
        if line[2:9] == 'stat.py':
            doprint = 1
        elif line[0:6] == '######' and doprint:
            doprint = 0
        if doprint:
            if len(line) > 2:
                print(line[2:], end=' ')
            else:
                print()


# Main program.

def main(argv):

    # Parse arguments.

    input_files = []
    level = 1
    nfilemax = 0
    all = 0

    args = argv[1:]
    while len(args) > 0:
        if args[0] == '-h' or args[0] == '--help':

            # Help.

            help()
            return 0

        elif args[0] == '--level' and len(args) > 1:

            # Analyze level.

            level = int(args[1])
            del args[0:2]

        elif args[0] == '--nfile' and len(args) > 1:

            # Number of files.

            nfilemax = int(args[1])
            del args[0:2]

        elif args[0] == '--all':

            # All files flag.

            all = 1
            del args[0]

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

    # # Loop over input files.

    # gtrees = {}
    # gbranches = {}
    # nfile = 0

    # for input_file in input_files:

    #     if nfilemax > 0 and nfile >= nfilemax:
    #         break
    #     nfile = nfile + 1

    #     if not larbatch_posix.exists(input_file):
    #         print('Input file %s does not exist.' % input_file)
    #         return 1

    #     print('\nOpening %s' % input_file)
    #     root = ROOT.TFile.Open(input_file)
    #     if not root.IsOpen() or root.IsZombie():
    #         print('Failed to open %s' % input_file)
    #         return 1

    #     # Analyze this file.

    #     analyze(root, level, gtrees, gbranches, all)

    # print '\n%d files analyzed.' % nfile

# Invoke main program.

if __name__ == '__main__':
    sys.exit(main(sys.argv))

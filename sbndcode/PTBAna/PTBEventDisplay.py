################################################
# Python script to display LLT and HLT trigger information for a given event
# in a stacked histogram format.
#
# Author: fnicolas@gnal.fov
#
# Date: 09/11/2025
#
# Input: PTBAna TTree
# Output: Stacked histograms of LLT and HLT triggers with timing information
#
# Run as follows:
# python PTBEventDisplay.py -s <path_to_TFile>
# 
# Command optional arguments:
#  -t <path_to_TTree>, default="ptbana/tree"
#  -e <event_number>, default=-1 (all events)
#  -m <trigger_mode>, defines the reference time HLT ID, e,g.
#       1: BNB zero bias
#       2: BNB light + beam
#       3: Off beam zero bias
#       ...
#       default=2
#  -bs <bin_size>, default=2000, in number of 20 ns ticks
#  -tw <time_window>, default=1e5, in number of 20 ns ticks
################################################



################################################
# Import necessary libraries
################################################
import matplotlib.pyplot as plt
import ROOT
import argparse
import numpy as np

rcParams = { 'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'axes.titlesize': 14,
    'legend.fontsize': 10,
}
plt.rcParams.update(rcParams)


################################################
# Define LLT types with names and colors
################################################

LLTTypes = {
    0: {"name": "Post ET Flash", "color": "gray"},
    1: {"name": "Post ET Inhibit", "color": "slategray"},
    2: {"name": "Random trigger 1", "color": "magenta"},
    3: {"name": "Random trigger 2", "color": "magenta"},
    4: {"name": "CRT reset Beam", "color": "blue"},
    5: {"name": "CRT reset EXT", "color": "blue"},
    6: {"name": "CRTs all, i.e. reset veto", "color": "lightseagreen"},
    7: {"name": "CRT top high AND", "color": "deepskyblue"},
    8: {"name": "CRT east AND", "color": "darkcyan"},
    9: {"name": "CRT south AND", "color": "royalblue"},
    10: {"name": "CRT bottom AND", "color": "steelblue"},
    11: {"name": "CRT top low AND", "color": "skyblue"},
    12: {"name": "CRT north AND", "color": "navy"},
    13: {"name": "CRT west AND", "color": "cyan"},
    14: {"name": "Any MTCA", "color": "orange"},
    15: {"name": "CAEN 10", "color": "red"},
    16: {"name": "", "color": "red"},
    17: {"name": "Any CAEN", "color": "red"},
    18: {"name": "MTC/A LLT18", "color": "lightcoral"},
    19: {"name": "", "color": "green"},
    20: {"name": "MTC/A LLT20 - flash", "color": "firebrick"},
    21: {"name": "MTC/A LLT21", "color": "indianred"},
    22: {"name": "", "color": "green"},
    23: {"name": "", "color": "green"},
    24: {"name": "EXT Early Warning", "color": "forestgreen"},
    25: {"name": "EXT Pre-Arrival", "color": "darkgreen"},
    26: {"name": "EXT Beam Accept", "color": "green"},
    27: {"name": "EXT Post-Arrival", "color": "seagreen"},
    28: {"name": "Early Warning", "color": "forestgreen"},
    29: {"name": "Pre-Arrival", "color": "darkgreen"},
    30: {"name": "Beam Accept", "color": "green"},
    31: {"name": "Post-Arrival", "color": "seagreen"}
}

################################################
# Define HLT types with names and colors
################################################

HLTTypes = {
    0: {"name": "Random Trigger", "color": "plum"},
    1: {"name": "Zero Bias", "color": "indigo"},
    2: {"name": "Light + Beam", "color": "blue"},
    3: {"name": "Off Beam 0 Bias", "color": "green"},
    4: {"name": "Light + Off Beam", "color": "yellow"},
    5: {"name": "Crossing Muon Placeholder", "color": "magenta"},
    6: {"name": "Low light + beam", "color": "cyan"},
    7: {"name": "CRT north & bottom", "color": "royalblue"},
    8: {"name": "CRT east & north", "color": "royalblue"},
    9: {"name": "CRT west & east", "color": "royalblue"},
    10: {"name": "CRT west & north", "color": "royalblue"},
    11: {"name": "CRT south & west", "color": "royalblue"},
    12: {"name": "CRT south & east", "color": "royalblue"},
    13: {"name": "CRT south & north", "color": "royalblue"},
    14: {"name": "CRT SN only", "color": "royalblue"},
    15: {"name": "CRT WE only", "color": "royalblue"},
    16: {"name": "Light + Post-Beam", "color": "orange"},
    17: {"name": "Light + Post-Off-Beam", "color": "blue"},
    18: {"name": "", "color": "green"},
    19: {"name": "Sunset", "color": "yellow"},
    20: {"name": "CRT Beam Reset", "color": "magenta"},
    21: {"name": "CRT Off Beam Reset", "color": "cyan"},
    22: {"name": "", "color": "violet"},
    23: {"name": "", "color": "azure"},
    24: {"name": "", "color": "gray"},
    25: {"name": "", "color": "pink"},
    26: {"name": "Beam acceptance flash", "color": "orange"},
    27: {"name": "Off-beam acceptance flash", "color": "peru"},
    28: {"name": "Light + Prearrival", "color": "darkgoldenrod"},
    29: {"name": "Light + Post ET", "color": "gold"},
    30: {"name": "Light + Off Beam Prearrival", "color": "chocolate"},
    31: {"name": "Random Trigger 2", "color": "orange"}
}


################################################
# Parse command line arguments
################################################
parser = argparse.ArgumentParser()
parser.add_argument("-s", help="Path to the TFile", required=True)
parser.add_argument("-t", help="Path to the TTree", default="ptbana/tree")
parser.add_argument("-e", help="Event number", default=-1)
parser.add_argument("-m", help="TriggerMode", default=2, type=int)
parser.add_argument("-bs", help="BinSize", default=2000, type=int)
parser.add_argument("-tw", help="TimeWindow", default=1e5, type=int)
parserargs = parser.parse_args()


################################################
# Function to get reference time based on trigger mode
################################################
def getReferenceTime(hlt_id, hlt_timestamp, mode):
    time = 0
    for k in range(len(hlt_id)):
        if(mode == hlt_id[k]): time = hlt_timestamp[k]
        
    return time


################################################
# Read the TFile and TTree
################################################
file = ROOT.TFile.Open(parserargs.s)
tree = file.Get(parserargs.t)
tree.Print()

################################################
# Loop over events in the PTBAna TTree
################################################
for i in range(tree.GetEntries()):

    ## Get entry
    tree.GetEntry(i)

    ## Skip if not the desired event
    if(parserargs.e != -1 and tree.event != int(parserargs.e)):
        continue

    ## Event label
    eventLabel = "Run_" + str(tree.run) + "_Subrun_" + str(tree.subrun) + "_Event_" + str(tree.event)
    
    ## Load LLT variables
    llt_id = np.array(list(tree.ptb_llt_trunmask))
    llt_timestamp = np.array(list(tree.ptb_llt_unmask_timestamp))

    ## Load HLT variables
    hlt_id = np.array(list(tree.ptb_hlt_trunmask))
    hlt_timestamp = np.array(list(tree.ptb_hlt_unmask_timestamp))

    ## Get reference time
    referenceTime = getReferenceTime(hlt_id, hlt_timestamp, parserargs.m)
    if referenceTime == 0:
        print("No reference time HLT found for event : ", eventLabel, ", skipping...")
        continue
    print("Event: ", eventLabel, ", with reference time: ", referenceTime)

    ## Plot: two pads with common x-axis
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True, num=eventLabel)
    plt.subplots_adjust(hspace=.025, left=0.1, right=0.95, top=0.95, bottom=0.1)

    ## Binning
    bins = np.arange(-1.*parserargs.tw, +1.*parserargs.tw, parserargs.bs)
    
    ## Sort and stack histograms for LLT
    lltSorted = {}
    for j in range(len(llt_id)):
        # skip CRT LLTs
        #if(llt_id[j] >= 6 and llt_id[j] < 15): continue
        if(llt_id[j] in lltSorted):
            lltSorted[llt_id[j]].append(llt_timestamp[j]-referenceTime)
        else:
            lltSorted[llt_id[j]] = [llt_timestamp[j]-referenceTime]
    ax1.hist( [ lltSorted[key] for key in lltSorted ], bins=bins, color=[LLTTypes[key]["color"] for key in lltSorted], label=[LLTTypes[key]["name"] for key in lltSorted], stacked=True)
        

    ## Sort and stack histograms for HLT
    hltSorted = {}
    for j in range(len(hlt_id)):
        if(hlt_id[j] in hltSorted):
            hltSorted[hlt_id[j]].append(hlt_timestamp[j]-referenceTime)
        else:
            hltSorted[hlt_id[j]] = [hlt_timestamp[j]-referenceTime]
    ax2.hist( [ hltSorted[key] for key in hltSorted ], bins=bins, color=[HLTTypes[key]["color"] for key in hltSorted], label=[HLTTypes[key]["name"] for key in hltSorted], stacked=True)


    ## Vertical line at reference time
    ax1.axvline(x=0, color='k', linestyle='--')
    ax2.axvline(x=0, color='k', linestyle='--')
    ax2.set_xlabel("TimeTick [20 ns]")
    ax1.set_ylabel("# LLTs")
    ax2.set_ylabel("# HLTs")

    # Shrink current axis by 20%
    box1 = ax1.get_position()
    ax1.set_position([box1.x0, box1.y0, box1.width, box1.height * 0.9])
    ax1.legend(loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.2))
    box2 = ax2.get_position()
    ax2.set_position([box2.x0, box2.y0, box2.width, box2.height * 0.85])
    ax2.legend(loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.15))

    fig.savefig("EventDisplay_" + eventLabel + ".pdf")
    plt.show()
    plt.close()

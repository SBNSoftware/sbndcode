from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import argparse, glob, os
import ROOT


VarsList = ["Eta_C", "Delta_C", "FitScore_C", "Alpha_C", "IntMode", "IntType"]
#VarsList = ["Delta_C", "FitScore_C", "Alpha_C"]

outputFileName = "hyperon_tmva.root"

# input dataframes
UseReco=False; UseReco=True
BGName="QE"; #BGLabel="RES"
BGLabel="BG_"+BGName
inputfile_dict = {}
ParentPath="/Users/franciscojaviernicolas/Work/HyperonsAna/SelectionCuts/DFAna"
if(UseReco==True):
    ParentPath+="Reco"
inputfile_dict["BG"] = ParentPath+"/HyperonSelectionVars_BG"+BGName+".csv"
inputfile_dict["S"] = ParentPath+"/HyperonSelectionVars_Signal.csv"
labelsV = ["S", "BG"]

optionsRDF = ROOT.RDF.RSnapshotOptions();
optionsRDF.fMode = "UPDATE";
optionsRDF.fOverwriteIfExists = 1;
# input dataframes
for l in labelsV:
    data = pd.read_csv(inputfile_dict[l], index_col=0)
    data = data[ data["Eta_C"].gt(0) ]
    data = {key: data[key].values for key in VarsList}


    if(l=="S"):
        rdfS = ROOT.RDF.MakeNumpyDataFrame(data)
        rdfS.Snapshot("TreeS", outputFileName,"", optionsRDF)
        rdfS.Display().Print()
    elif(l=="BG"):
        rdfBG = ROOT.RDF.MakeNumpyDataFrame(data)
        rdfBG.Snapshot("TreeBG", outputFileName,"", optionsRDF)
        rdfBG.Display().Print()

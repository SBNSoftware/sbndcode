#!/usr/bin/env python
# PandoraMVA.py

from sklearn import preprocessing
from datetime import datetime

import numpy as np
import sys
import time
import pickle

def LoadData(trainingFileName, delimiter=','):
    # Use the first example to get the number of columns
    with open(trainingFileName) as file:
        ncols = len(file.readline().split(delimiter))
        
    # First column is a datestamp, so skip it
    trainingSet = np.genfromtxt(trainingFileName, delimiter=delimiter, usecols=range(1,ncols), 
                                dtype=None)
                                
    nExamples = trainingSet.size
    nFeatures = ncols - 2 # last column is the response
    
    return np.array(trainingSet), nFeatures, nExamples
    
#--------------------------------------------------------------------------------------------------

def SplitTrainingSet(trainingSet, nFeatures):
    X=[] # features sets
    Y=[] # responses

    for example in trainingSet:
        Y.append(int(example[nFeatures])) # type of Y should be bool or int
        features = []
        for i in range(0, nFeatures):
            features.append(float(example[i])) # features in this SVM must be Python float
            
        X.append(features)

    return np.array(X).astype(np.float64), np.array(Y).astype(np.int)
    
#--------------------------------------------------------------------------------------------------

def Randomize(X, Y, setSameSeed=False):
    if setSameSeed:
        np.random.seed(0)

    order = np.random.permutation(Y.size)
    return X[order], Y[order]
    
#--------------------------------------------------------------------------------------------------

def Sample(X, Y, testFraction=0.1):
    trainSize = int((1.0 - testFraction) * Y.size)
    
    X_train = X[:trainSize]
    Y_train = Y[:trainSize]
    X_test  = X[trainSize:]
    Y_test  = Y[trainSize:]
    
    return X_train, Y_train, X_test, Y_test
    
#--------------------------------------------------------------------------------------------------

def ValidateModel(model, X_test, Y_test):               
    return model.score(X_test, Y_test)

#--------------------------------------------------------------------------------------------------

def OverwriteStdout(text):
    sys.stdout.write('\x1b[2K\r' + text)
    sys.stdout.flush()

#--------------------------------------------------------------------------------------------------

def OpenXmlTag(modelFile, tag, indentation):
    modelFile.write((' ' * indentation) + '<' + tag + '>\n')
    return indentation + 4
    
#--------------------------------------------------------------------------------------------------

def CloseXmlTag(modelFile, tag, indentation):
    indentation = max(indentation - 4, 0)
    modelFile.write((' ' * indentation) + '</' + tag + '>\n')
    return indentation

#--------------------------------------------------------------------------------------------------

def WriteXmlFeatureVector(modelFile, featureVector, tag, indentation):
    modelFile.write((' ' * indentation) + '<' + tag + '>')

    firstTime=True
    for feature in featureVector:
        if firstTime:
            modelFile.write(str(feature))
            firstTime=False
        else:
            modelFile.write(' ' + str(feature))
            
    modelFile.write('</' + tag + '>\n')
    
#--------------------------------------------------------------------------------------------------

def WriteXmlFeature(modelFile, feature, tag, indentation):
    modelFile.write((' ' * indentation) + '<' + tag + '>')
    modelFile.write(str(feature))     
    modelFile.write('</' + tag + '>\n')


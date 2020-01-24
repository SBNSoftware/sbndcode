#!/usr/bin/env python
# PandoraSVM.py

from sklearn import svm
from sklearn import preprocessing
from datetime import datetime

import numpy as np
import sys
import time
import pickle

from PandoraMVA import *

def StandardizeFeatures(X):
    muValues    = np.mean(X, axis=0)
    sigmaValues = np.std(X, axis=0)
    return np.divide((X - muValues), sigmaValues), muValues, sigmaValues
    
#--------------------------------------------------------------------------------------------------

def TrainModel(X_train, Y_train, kernelString, kernelDegree=2, gammaValue=0.05, coef0Value=1.0, 
               cValue=1.0, tol=0.001, cache_size=1000, enableProbability=False, shrinking=True):
    # Load the SVC object
    svmModel = svm.SVC(C=cValue, cache_size=cache_size, class_weight=None, coef0=coef0Value,
                       decision_function_shape=None, degree=kernelDegree, gamma=gammaValue, 
                       kernel=kernelString, max_iter=-1, probability=enableProbability, 
                       random_state=None, shrinking=shrinking, tol=tol, verbose=False)
    
    # Train the model   
    startTime = time.time() 
    svmModel.fit(X_train, Y_train)
    
    endTime = time.time()
    nSupportVectors = svmModel.support_vectors_.shape[0]

    return svmModel, endTime - startTime, nSupportVectors
    
#--------------------------------------------------------------------------------------------------

def QuickTest(X_train, Y_train, X_test, Y_test, kernelString, kernelDegree=2, gammaValue=0.05, 
              coef0Value=1.0, cValue=1.0):
    # Train and validate the model
    svmModel, trainingTime, nSupportVectors = TrainModel(X_train, Y_train, kernelString, 
                                                         kernelDegree, gammaValue, coef0Value, 
                                                         cValue)
                                                         
    modelScore = ValidateModel(svmModel, X_test, Y_test)
    
    # Write validation output to screen
    stdoutString = '[' + kernelString
    if kernelString == 'poly':
        stdoutString += ',deg=' + str(kernelDegree)
        
    elif kernelString == 'rbf':
        stdoutString += ',gamma=' + str(gammaValue)
        
    if kernelString != 'rbf' and kernelString != 'linear':
        stdoutString += ',coef0=' + str(coef0Value)

    stdoutString += (',C=' + str(cValue) + '] : %.1f%% (%d seconds, %d SVs)\n' % 
                     (modelScore * 100, trainingTime, nSupportVectors))
                     
    OverwriteStdout(stdoutString)
    
#--------------------------------------------------------------------------------------------------

def WriteXmlFile(filePath, svmName, datetimeString, yAlpha, bias, kernel, mu, scale, sigma, 
                 supportVectors, standardize=True, enableProbability=False, probAParam=0.0, 
                 probBParam=0.0):
    with open(filePath, "a") as modelFile:
        standStr = str(standardize).lower()
        probStr = str(enableProbability).lower()
        indentation = OpenXmlTag(modelFile,    'SupportVectorMachine', 0)
        WriteXmlFeature(modelFile, svmName,        'Name', indentation)
        WriteXmlFeature(modelFile, datetimeString, 'Timestamp', indentation)
        
        indentation = OpenXmlTag(modelFile,        'Machine', indentation)
        WriteXmlFeature(modelFile, kernel,             'KernelType', indentation)
        WriteXmlFeature(modelFile, bias,               'Bias', indentation)
        WriteXmlFeature(modelFile, scale,              'ScaleFactor', indentation)
        WriteXmlFeature(modelFile, standStr,           'Standardize', indentation)
        WriteXmlFeature(modelFile, probStr,            'EnableProbability', indentation)
        WriteXmlFeature(modelFile, probAParam,         'ProbAParameter', indentation)
        WriteXmlFeature(modelFile, probBParam,         'ProbBParameter', indentation)
        indentation = CloseXmlTag(modelFile,       'Machine', indentation)
        
        indentation = OpenXmlTag(modelFile,        'Features', indentation)
        WriteXmlFeatureVector(modelFile, mu,           'MuValues', indentation)
        WriteXmlFeatureVector(modelFile, sigma,        'SigmaValues', indentation)
        indentation = CloseXmlTag(modelFile,       'Features', indentation)
        
        for supVec, yAlphaValue in zip(supportVectors, yAlpha):
            indentation = OpenXmlTag(modelFile,    'SupportVector', indentation)
            WriteXmlFeature(modelFile, yAlphaValue,    'AlphaY', indentation)
            WriteXmlFeatureVector(modelFile, supVec,   'Values', indentation)
            indentation = CloseXmlTag(modelFile,   'SupportVector', indentation)
        
        CloseXmlTag(modelFile,                 'SupportVectorMachine', indentation)

#--------------------------------------------------------------------------------------------------

def GetKernelInt(kernelType, kernelDegree=2):   
    if kernelType == 'linear':
        return 1
        
    if kernelType == 'poly' and kernelDegree == 2:
        return 2
        
    if kernelType == 'poly' and kernelDegree == 3:
        return 3
        
    if kernelType == 'rbf':
        return 4

    raise ValueError('Unknown kernel type for Pandora kernel enum: ' + kernelType)

#--------------------------------------------------------------------------------------------------

def SerializeToPkl(fileName, svmModel, mu, sigma):
    with open(fileName, 'w') as f:
        pickle.dump(svmModel, f)
        pickle.dump(mu, f)
        pickle.dump(sigma, f)

#--------------------------------------------------------------------------------------------------

def LoadFromPkl(fileName, svmModel, mu, sigma):
    with open(fileName, 'r') as f:
        svmModel = pickle.load(f)
        mu       = pickle.load(f)
        sigma    = pickle.load(f)

        return svmModel, mu, sigma


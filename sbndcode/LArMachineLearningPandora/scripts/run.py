#!/usr/bin/env python
# example.py

from PandoraSVM import *
import sys

if __name__=="__main__":

    # Settings ------------------------------------------------------------------------------------
    
    trainingFile      = '../to_use/cut_7600_300238.txt'
    svmName           = 'SVMVertexSelection'
    kernelType        = 'rbf'       # poly, rbf, linear
    kernelDegree      = 2           # only for poly kernel
    cValue            = float(sys.argv[1])
    
    gammaValue        = float(sys.argv[2])
    coef0Value        = 1.0         # default 1.0
    enableProbability = False
    
    serializeToPkl    = True
    serializeToXml    = True
    loadFromPkl       = False
    xmlFileName       = 'max_14650_300075_rbf_' + str(cValue) + '_' + str(gammaValue) + '.xml'
    pklFileName       = 'max_14650_300075_rbf_' + str(cValue) + '_' + str(gammaValue) + '.pkl'
    
    tol       = 0.001
    shrinking = False 
    
    #----------------------------------------------------------------------------------------------
    
    if loadFromPkl:
        OverwriteStdout('Loading model from file ' + pklFileName + '\n')
        svmModel, mu, sigma = LoadFromPkl(pklFileName)
    
    else:
        # Load the data
        OverwriteStdout('Loading training set data from file ' + trainingFile + '\n')
        trainSet, nFeatures, nExamples = LoadData(trainingFile, ',')

        # Standardize the data and hold onto the means and stddevs for later
        OverwriteStdout(('Preprocessing ' + str(nExamples) + ' training examples of ' + 
                         str(nFeatures) + ' features'))
        
        X_org, Y_org     = SplitTrainingSet(trainSet, nFeatures)
        X_org, mu, sigma = StandardizeFeatures(X_org)
        
        # Train the SVM
        X, Y = Randomize(X_org, Y_org)
        #X_train, Y_train, X_test, Y_test = Sample(X, Y, 0.1)
        
        OverwriteStdout('Training SVM...')
        svmModel, trainingTime, nSupportVectors = TrainModel(X, Y, kernelType, gammaValue=gammaValue, 
                                                             cValue=cValue, kernelDegree=kernelDegree, 
                                                             coef0Value=coef0Value, shrinking=shrinking,
                                                             tol=tol, enableProbability=enableProbability)
        
        OverwriteStdout(('Trained SVM with ' + str(nFeatures) + ' features and ' + str(nExamples) + 
                         ' examples (%d seconds, %d SVs)\n' % (trainingTime, nSupportVectors)))
                  
        # Validate the model 
        #modelScore = ValidateModel(svmModel, X_test, Y_test)
        #OverwriteStdout('Model score: %.2f%%\n' % (modelScore * 100))
        
        # Save the model                  
        modelParams    = svmModel.get_params(deep=True)
        supportVectors = svmModel.support_vectors_
        yAlpha         = svmModel.dual_coef_[0]
        bias           = svmModel.intercept_[0]

        if enableProbability:
            probAParam     = svmModel.probA_[0]
            probBParam     = svmModel.probB_[0]

        else:
            probAParam = 0.0
            probBParam = 0.0
        
        if serializeToXml:
            OverwriteStdout('Writing model to xml file ' + xmlFileName + '\n')
            datetimeString = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
            WriteXmlFile(xmlFileName, svmName, datetimeString, yAlpha, bias, GetKernelInt(kernelType), mu, 
                         gammaValue, sigma, supportVectors, standardize=True, enableProbability=enableProbability, 
                         probAParam=probAParam, probBParam=probBParam) 
                     
        if serializeToPkl:
            OverwriteStdout('Writing model to pkl file ' + pklFileName + '\n')
            SerializeToPkl(pklFileName, svmModel, mu, sigma)
            
    # Do other stuff with your trained/loaded model
    # ...
    

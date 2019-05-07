#!/usr/bin/env python
# example.py

from PandoraSVM import *

if __name__=="__main__":

    # Settings ------------------------------------------------------------------------------------
    
    trainingFile    = '/sbnd/data/users/dbarker/showers/PandoraTuning/Pandora_SVM/VertexSelection100k.txt'
    svmName         = 'VertexSVM'
    kernelType      = 'poly' # poly, rbf, linear
    kernelDegree    = 3      # only for poly kernel
    cValue          = 10    # default 1.0
    gammaValue      = 0.046415888336127725  # default 1/n_features
    coef0Value      = 1.0    # default 1.0
    
    serializeToPkl  = True
    serializeToXml  = True
    loadFromPkl     = False
    xmlFileName     = 'SVM_SBND.xml'
    pklFileName     = 'SVM_SBND.pkl'
    
    #----------------------------------------------------------------------------------------------
    
    if loadFromPkl:
        print 1
        OverwriteStdout('Loading model from file ' + pklFileName + '\n')
        svmModel, mu, sigma = LoadFromPkl(pklFileName)
    
    else:
        # Load the data
        OverwriteStdout('Loading training set data from file ' + trainingFile + '\n')
        trainSet, nFeatures, nExamples = LoadData(trainingFile, ',')

        # Standardize the data and hold onto the means and stddevs for later
        OverwriteStdout(('Preprocessing ' + str(nExamples) + ' training examples of ' + 
                         str(nFeatures) + ' features'))
        
        X, Y         = SplitTrainingSet(trainSet, nFeatures)
        X, mu, sigma = StandardizeFeatures(X)
        X, Y         = Randomize(X, Y)
        X_train, Y_train, X_test, Y_test = Sample(X, Y, 0.1)
        
        # Train the SVM
        OverwriteStdout('Training SVM...')
        svmModel, trainingTime, nSupportVectors = TrainModel(X_train, Y_train, kernelType, gammaValue=gammaValue, 
                                                             cValue=coef0Value, kernelDegree=kernelDegree, 
                                                             coef0Value=coef0Value)
        
        OverwriteStdout(('Trained SVM with ' + str(nFeatures) + ' features and ' + str(nExamples) + 
                         ' examples (%d seconds, %d SVs)\n' % (trainingTime, nSupportVectors)))
                  
        # Validate the model 
        modelScore = ValidateModel(svmModel, X_test, Y_test)
        OverwriteStdout('Model score: %.2f%%\n' % (modelScore * 100))
        
        # Save the model                  
        modelParams    = svmModel.get_params(deep=True)
        supportVectors = svmModel.support_vectors_
        yAlpha         = svmModel.dual_coef_[0]
        bias           = svmModel.intercept_[0]
        
        if serializeToXml:
            OverwriteStdout('Writing model to xml file ' + xmlFileName + '\n')
            datetimeString = datetime.now().strftime('%Y-%m-%d_%H:%M:%S')
            WriteXmlFile(xmlFileName, svmName, datetimeString, yAlpha, bias, GetKernelInt(kernelType), mu, 
                         gammaValue, sigma, supportVectors, standardize=True) 
                     
        if serializeToPkl:
            OverwriteStdout('Writing model to pkl file ' + pklFileName + '\n')
            SerializeToPkl(pklFileName, svmModel, mu, sigma)
            
    # Do other stuff with your trained/loaded model
    # ...
    

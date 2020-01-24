#!/usr/bin/env python
# example.py

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import load_iris
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import GridSearchCV

from PandoraSVM import *

# Utility function to move the midpoint of a colormap to be around
# the values of interest.

class MidpointNormalize(Normalize):

    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

if __name__=="__main__":

    trainingFile    = '/sbnd/data/users/dbarker/showers/PandoraTuning/Pandora_SVM/VertexSelection.txt'

    # Load the data
    OverwriteStdout('Loading training set data from file ' + trainingFile + '\n')
    trainSet, nFeatures, nExamples = LoadData(trainingFile, ',')

    # Standardize the data and hold onto the means and stddevs for later
    OverwriteStdout(('Preprocessing ' + str(nExamples) + ' training examples of ' + 
                     str(nFeatures) + ' features'))
    
    X_org, Y_org     = SplitTrainingSet(trainSet, nFeatures)
    X_org, mu, sigma = StandardizeFeatures(X_org)
    
    C_range = np.logspace(-3, 9, 13)
    gamma_range = np.logspace(-10, 3, 13)
    #C_range = np.linspace(700000, 1400000, num=15)
    #gamma_range = np.linspace(0.4642, 0.4642, num=1)
    
    param_grid = dict(gamma=gamma_range, C=C_range)
    cv = StratifiedShuffleSplit(n_splits=5, test_size=0.2, random_state=42)
    grid = GridSearchCV(SVC(), param_grid=param_grid, cv=cv, n_jobs=8)
    grid.fit(X_org, Y_org)
    
    print("The best parameters are %s with a score of %0.2f"
      % (grid.best_params_, grid.best_score_))
    
    scores = grid.cv_results_['mean_test_score'].reshape(len(C_range),
                                                     len(gamma_range))
    
    plt.figure(figsize=(8, 6))
    plt.subplots_adjust(left=.2, right=0.95, bottom=0.15, top=0.95)
    plt.imshow(scores, interpolation='nearest', cmap=plt.cm.hot,
               norm=MidpointNormalize(vmin=0.2, midpoint=0.88))
    plt.xlabel('gamma')
    plt.ylabel('C')
    plt.colorbar()
    plt.xticks(np.arange(len(gamma_range)), gamma_range, rotation=45)
    plt.yticks(np.arange(len(C_range)), C_range)
    plt.title('Validation accuracy')
    plt.show()

### LArSvmVertexSelection algorithm

1. Run the reconstruction over the training set samples with this algorithm in training set mode; e.g.
```xml
<algorithm type = "LArSvmVertexSelection">
	<InputCaloHitListNames>CaloHitListU CaloHitListV CaloHitListW</InputCaloHitListNames>
	<InputClusterListNames>ClustersU ClustersV ClustersW</InputClusterListNames>
	<OutputVertexListName>NeutrinoVertices3D</OutputVertexListName>
	<ReplaceCurrentVertexList>true</ReplaceCurrentVertexList>
	<TrainingSetMode>true</TrainingSetMode>
	<TrainingOutputFileRegion>VertexSelection_Region_MCC7</TrainingOutputFileRegion>
	<TrainingOutputFileVertex>VertexSelection_Vertex_MCC7</TrainingOutputFileVertex>
	<MCParticleListName>MCParticleList3D</MCParticleListName>
	<MCVertexXCorrection>0.495694</MCVertexXCorrection>
	<CaloHitListName>CaloHitList2D</CaloHitListName>
	<FeatureTools>
    		<tool type = "LArEnergyKickFeature"/>
    		<tool type = "LArLocalAsymmetryFeature"/>
    		<tool type = "LArGlobalAsymmetryFeature"/>
    		<tool type = "LArShowerAsymmetryFeature"/>
    		<tool type = "LArRPhiFeature"/>
	</FeatureTools>
</algorithm>
```
2. The above produces sets of training data files of the format `VertexSelection_Region_MCC7_[interaction_type].txt` and `VertexSelection_Vertex_MCC7_[interaction_type].txt`. Randomly permute the vertex SVM data files and sample from different interaction types as desired to produce a vertex SVM training set&mdash;and similarly for the region SVM data.

3. The vertex and region SVMs are treated independently. Analogously to the below, use the sklearn-based Python script `scripts/rbf_gridsearch_test.py` to find graphically the optimal values of `C` and `gamma`.

4. Use the `C` and `gamma` values to train the SVM and produce XML files using the script at `scripts/run.py`. Ensure that the vertex SVM is named `VertexSelectionVertex` and the region SVM `VertexSelectionRegion` (or whatever names the algorithm is pointed to when not in training set mode).

5. Concatenate the two XML files and point the algorithm towards the combined file when running via the `SvmFileName` parameter; e.g. 
```xml
<algorithm type = "LArSvmVertexSelection">
	<InputCaloHitListNames>CaloHitListU CaloHitListV CaloHitListW</InputCaloHitListNames>
	<InputClusterListNames>ClustersU ClustersV ClustersW</InputClusterListNames>
	<OutputVertexListName>NeutrinoVertices3D</OutputVertexListName>
	<ReplaceCurrentVertexList>true</ReplaceCurrentVertexList>
	<SvmFileName>PandoraSvm_VertexSelection_MicroBooNE_mcc7.xml</SvmFileName>
	<RegionSvmName>VertexSelectionRegion</RegionSvmName>
	<VertexSvmName>VertexSelectionVertex</VertexSvmName>
	<FeatureTools>
	    <tool type = "LArEnergyKickFeature"/>
	    <tool type = "LArLocalAsymmetryFeature"/>
	    <tool type = "LArGlobalAsymmetryFeature"/>
	    <tool type = "LArShowerAsymmetryFeature"/>
	    <tool type = "LArRPhiFeature"/>
	</FeatureTools>
    </algorithm>
```

********************************************************

Instructions for creating Pfo Characterisation SVM model, 
implemented through the algorithm LArContent/larpandoracontent/LArTrackShowerId/SvmPfoCharacterisationAlgorithm.cc


0 - Separate your sample of events in two, one for training and one for testing (not necessarily 50% each, but representative of the spectrum of events for your problem)

1 - Run over your training subsample with the option TrainingSetMode = true. To do this, edit the corresponding PandoraSettings file (e.g. LArReco/scripts/uboone/PandoraSettings_MicroBooNE_Neutrino.xml) and provide the output training file name like: 
     <algorithm type = "LArSVMClusterCharacterisation">
     ...
     <TrainingSetMode>true</TrainingSetMode>
     <TrainingOutputFileName>OUTPUT_NAME</TrainingOutputFileName>
This will create a txt file (note that the .txt will be appended to the OUTPUT_NAME provided) with the features calculated for the events in your input file. You can see an example in the file SVM_training_data_pfocharacterisation_example.txt, containing lines like: 
     05/30/17_16:10:06,98.7135,0.00741252,0.00355123,0.00949093,0.00658028,0.0324565,0.0307609,0.000491683,1
which are a list of the features starting with a timestamp and finishing with the true value to train for latter classification. In this example, the true value is 1 for a track and 0 for a shower, and the features are the variables computed by the tools:
      <FeatureTools>                                                                                                            
      	<tool type = "LArLinearFitFeatureTool"/>
	<tool type = "LArShowerFitFeatureTool"/>
	<tool type = "LArVertexDistanceFeatureTool"/>
      </FeatureTools>	
which are added in this order: 1) straight line length, 2) mean of difference with straight line, 3) sigma (standard deviation) of difference with straight line, 4) dTdL width, 5) max gap lenght, 6) RMS of linear fit, 7) shower fit width, 8) vertex distance. Expect the first feature (straight line lenght) the other ones are normalized divided by the straight line lenght if the option RatioVariables = true. 
The implementation of the tools and available variables can be found in LArContent/larpandoracontent/LArTrackShowerId/TrackShowerIdFeatureTool.cc

*** Note: the next two steps are specific for the rbf (radial basis function) kernel option. For other options, check http://scikit-learn.org/stable/

2 - Use the python script rbf_gridsearch_test.py to search for the optimal values of C and gamma for your training data. Edit rbf_gridsearch_test.py and give the text file calculated in step 1 as trainingFile. This script will do a grid search which is time and memory intense, consider sampling your training data accordingly (for example, the SVM_training_data_pfocharacterisation_example.txt contains randomly selected 1000 training examples from the entire training data for this step). This script will report at the end that "The best parameters are C: and gamma: " with a given score. The score is a measurement of the classification, for example in the track-shower characterisation it would be: ntracks*tracks_eff + nshowers*showers_eff. The python script produces also a plot like Example_rbf_output.png, with indicative values of the score in the searched grid, which is to be checked to ensure that it is smooth and the selected grid was enough to find a reliable best score (otherwise, if the selected point is at an edge of the grid, consider extending the grid extremes and running again this step).

3 - With the values of C and gamma obtained in the previous step, run example.py. Edit it and change C and gamma, and give the appropriate trainingFile name.
This step is less time and memory consuming, so the input data can be scaled (for example using 100k training examples). This will give another score, which is to be checked against the one in the previous step. If it is very different, it could mean that the sampled training examples used in step 2 were not representative enough of the entire training data, and you might consider running again from step 2 with a larger training sampled input. The output of step 3 will be a .xml file (as well as a .pkl file) with the model, i.e. the SVs, to be used for solving the problem afterwards in your testing data (separated in setp 0), i.e. the input file to be given to the algorithm using it, like in this case: 
     <algorithm type = "LArSvmPfoCharacterisation">
     		...
                <SvmFileName>/r05/dune/MachineLearningData/PandoraSvm_PfoCharacterisation_MicroBooNE_mcc7.xml</SvmFileName>
                <SvmName>FinalPfoCharacterisation</SvmName>
		...		
The SvmName is the one added in the output .xml file, and can be changed in example.py 
      			 

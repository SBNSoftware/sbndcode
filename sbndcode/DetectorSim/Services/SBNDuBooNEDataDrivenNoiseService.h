// SBNDuBooNEDataDrivenNoiseService.h
// Andrew Scarff (University of Sheffield)
// July 2019

// Based upon SPhaseChannelNoiseService.h created by Jingbo Wang (UC Davis) for ProtoDUNE.

// Implementation of a general TPC channel noise model with:
// (1) white noise
// (2) Inherent Gaussian noise in frequency
// (3) MicroBooNE noise in frequency
// (4) Coherent noise (exponential + Gaussian) in frequency 
//     (Note a: phase at each frequency bin is randamized at the moment. Will be updated soon
//      Note b: Currently, consecutive offline channels (configurable) are grouped together and 
//              the same coherent noise waveform is assigned to channels within the same group. )
//
// The default parameters are obtained from the ProtoDUNE-SP data (run 4096)
// fcl file: sbndcode/DetectorSim/Services/noiseservices_sbnd.fcl
//

#ifndef SBNDuBooNEDataDrivenNoiseService_H
#define SBNDuBooNEDataDrivenNoiseService_H

#include "sbndcode/DetectorSim/Services/ChannelNoiseService.h"

#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/LArFFT.h"
#include "larcore/Geometry/Geometry.h"
#include "nurandom/RandomUtils/NuRandomService.h"

#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "TH1F.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TMath.h"

#include <vector>
#include <iostream>
#include <sstream>

class TH1;
namespace CLHEP {
class HepRandomEngine;
}

class SBNDuBooNEDataDrivenNoiseService : public ChannelNoiseService {

public:

  // Ctor.
  SBNDuBooNEDataDrivenNoiseService(fhicl::ParameterSet const& pset);

  // Ctor.
  SBNDuBooNEDataDrivenNoiseService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Dtor.
  ~SBNDuBooNEDataDrivenNoiseService();

  // Add noise to a signal array.
  int addNoise(detinfo::DetectorClocksData const& clockData, Channel chan, AdcSignalVector& sigs) const override;

  void generateNoise(detinfo::DetectorClocksData const& clockData) override;
 
  // Print the configuration.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const override;

private:
 
  // Fill the noise vectors.
  //void generateNoise();
  
  // Fill a noise vector.
  // Input vector contents are lost.
  // The size of the vector is obtained from the FFT service.
  void generateMicroBooNoise(float wirelength, float ENOB, 
                     AdcSignalVector& noise, TH1* aNoiseHist) const;
  void generateGaussianNoise(detinfo::DetectorClocksData const& clockData,
                             AdcSignalVector& noise, std::vector<float> gausNorm,
	                    std::vector<float> gausMean, std::vector<float> gausSigma,
	                    TH1* aNoiseHist) const;
  void generateCoherentNoise(detinfo::DetectorClocksData const& clockData,
                             AdcSignalVector& noise, std::vector<float> gausNorm,
	                    std::vector<float> gausMean, std::vector<float> gausSigma,
	                    float cohExpNorm, float cohExpWidth, float cohExpOffset, 
	                    TH1* aNoiseHist) const;
  
  // Make coherent groups
  void makeCoherentGroupsByOfflineChannel(unsigned int nchpergroup);
  std::vector<unsigned int> fChannelGroupMap;   ///< assign each channel a group number
  std::vector<int> fGroupCoherentNoiseMap; ///< assign each group a noise 
  unsigned int getGroupNumberFromOfflineChannel(unsigned int offlinechan) const;
  unsigned int getCohNoiseChanFromGroup(unsigned int cohgroup) const;
  
  // General parameters
  unsigned int fNoiseArrayPoints;  ///< number of points in randomly generated noise array
  int          fRandomSeed;        ///< Seed for random number service. If absent or zero, use SeedSvc.
  int          fLogLevel;          ///< Log message level: 0=quiet, 1=init only, 2+=every event
  
  // Inherent White noise parameters
  bool         fEnableWhiteNoise;
  float        fWhiteNoiseZ;       ///< Level (per freq bin) for white noise for Z.
  float        fWhiteNoiseU;       ///< Level (per freq bin) for white noise for U.
  float        fWhiteNoiseV;       ///< Level (per freq bin) for white noise for V.  
  
  // Inherent Gaussian noise
  bool         fEnableGaussianNoise;
  std::vector<float>	fGausNormU;		///<  noise scale factor for the gaussian component in coherent noise
  std::vector<float>	fGausMeanU;		///<  mean of the gaussian component in coherent noise
  std::vector<float>	fGausSigmaU;	///<  sigma of the gaussian component in coherent noise
  std::vector<float>	fGausNormV;		///<  noise scale factor for the gaussian component in coherent noise
  std::vector<float>	fGausMeanV;		///<  mean of the gaussian component in coherent noise
  std::vector<float>	fGausSigmaV;	///<  sigma of the gaussian component in coherent noise
  std::vector<float>	fGausNormZ;		///<  noise scale factor for the gaussian component in coherent noise
  std::vector<float>	fGausMeanZ;		///<  mean of the gaussian component in coherent noise
  std::vector<float>	fGausSigmaZ;	///<  sigma of the gaussian component in coherent noise
  
  // Inherent MicroBoone noise parameters
  bool         fEnableMicroBooNoise;  ///< enable MicroBooNE noise model
  float        fENOB;                 ///< Effective number of bits
  bool         fIncludeJumpers;       ///< Include jumper term
  float        fJumperCapacitance;    ///< Capacitance of jumper cables in pF. Defaults to 0 if not included.
  float        fUFirstJumper;         ///< Wire number of first wire on U layer to include a jumper cable. Defaults to 0 if not included.
  float        fULastJumper;          ///< Wire number of last wire on U layer to include a jumper cable. Defaults to 0 if not included.
  float        fVFirstJumper;         ///< Wire number of first wire on V layer to include a jumper cable. Defaults to 0 if not included.
  float        fVLastJumper;          ///< Wire number of last wire on V layer to include a jumper cable. Defaults to 0 if not included.
  std::vector<float>  fNoiseFunctionParameters;  ///< Parameters in the MicroBooNE noise model
  
  // Coherent Noise parameters
  bool         fEnableCoherentNoise;
  std::vector<unsigned int> fNChannelsPerCoherentGroup;
  unsigned int fExpNoiseArrayPoints;  ///< number of points in randomly generated noise array
  unsigned int fCohNoiseArrayPoints;  ///< number of points in randomly generated noise array
  float        fCohExpNorm;           ///< noise scale factor for the exponential component component in coherent noise
  float        fCohExpWidth;          ///< width of the exponential component in coherent noise
  float        fCohExpOffset;         ///< Amplitude offset of the exponential background component in coherent noise
  std::vector<float>	fCohGausNorm;		///<  noise scale factor for the gaussian component in coherent noise
  std::vector<float>	fCohGausMean;		///<  mean of the gaussian component in coherent noise
  std::vector<float>	fCohGausSigma;	///<  sigma of the gaussian component in coherent noise 
  
  
  // Inherent Gausian noise arrays
  AdcSignalVectorVector fGausNoiseZ;
  AdcSignalVectorVector fGausNoiseU;
  AdcSignalVectorVector fGausNoiseV;
  
  // Inherent MicroBoo noise arrays
  AdcSignalVectorVector fMicroBooNoiseZ;  
  AdcSignalVectorVector fMicroBooNoiseU; 
  AdcSignalVectorVector fMicroBooNoiseV;
  
  // Coherent Noise array.
  AdcSignalVectorVector fCohNoiseZ;  ///< noise on each channel for each time for all planes  
  AdcSignalVectorVector fCohNoiseU;  ///< noise on each channel for each time for all planes
  AdcSignalVectorVector fCohNoiseV;  ///< noise on each channel for each time for all planes


  // Histograms.
  
  TH1* fGausNoiseHistZ;     ///< distribution of noise counts for Z
  TH1* fGausNoiseHistU;     ///< distribution of noise counts for U
  TH1* fGausNoiseHistV;     ///< distribution of noise counts for V
  TH1* fGausNoiseChanHist;  ///< distribution of accessed noise samples

  TH1* fMicroBooNoiseHistZ;     ///< distribution of noise counts for Z
  TH1* fMicroBooNoiseHistU;     ///< distribution of noise counts for U
  TH1* fMicroBooNoiseHistV;     ///< distribution of noise counts for V
  TH1* fMicroBooNoiseChanHist;  ///< distribution of accessed noise samples
  
  TH1* fCohNoiseHist;      ///< distribution of noise counts
  TH1* fCohNoiseChanHist;  ///< distribution of accessed noise samples

  TF1* _wld_f;
  double wldparams[2];

  TF1* _poisson;


  // Randomisation.
  bool haveSeed;
  CLHEP::HepRandomEngine* m_pran;
  CLHEP::HepRandomEngine* ConstructRandomEngine(const bool haveSeed);
  double GetRandomTF1(TF1* func) const;
  TRandom3* fTRandom3;


};

DECLARE_ART_SERVICE_INTERFACE_IMPL(SBNDuBooNEDataDrivenNoiseService, ChannelNoiseService, LEGACY)

#endif

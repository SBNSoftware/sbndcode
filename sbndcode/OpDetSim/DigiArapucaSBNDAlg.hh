////////////////////////////////////////////////////////////////////////
// File:       DigiArapucaSBNDAlg.hh
//
// This algorithm is used for the electronic response of arapucas
// Created by L. Paulucci, F. Marinho and I.L. de Icaza
// Based on OpDetDigitizerDUNE_module.cc and SimPMTSBND_module.cc
////////////////////////////////////////////////////////////////////////

#ifndef SBND_OPDETSIM_DIGIARAPUCASBNDALG_HH
#define SBND_OPDETSIM_DIGIARAPUCASBNDALG_HH

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/DelegatedParameter.h"
#include "fhiclcpp/types/OptionalDelegatedParameter.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandExponential.h"
#include "art/Utilities/make_tool.h"

#include <algorithm>
#include <memory>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include <sstream>
#include <fstream>

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

#include "sbndcode/OpDetSim/HDWvf/HDOpticalWaveforms.hh"

#include "TFile.h"

namespace opdet {

  class DigiArapucaSBNDAlg {

  public:

    struct ConfigurationParameters_t {
      double PeakTime;       // maximum time of each component pulse in ns
      double PulseLength;    // in ns
      double MeanAmplitude;  // Mean amplitude for pulses in mV
      double RiseTime;       // pulse rising time constant (exponential)
      double FallTime;       // pulse decaying time constant (exponential)
      double ADC;            //voltage to ADC convertion scale
      double Baseline;       //waveform baseline
      double BaselineRMS;    //Pedestal RMS in ADC counts
      double DarkNoiseRate;  //in Hz
      double CrossTalk;      //probability for producing a signal of 2 PE in response to 1 photon
      double SaturationHigh;     //Saturation in number of p.e.
      double SaturationLow;     //Saturation in number of p.e.
      double XArapucaVUVEffVis;    //XArapucaVUV efficiency to visible light
      double XArapucaVUVEffVUV;    //XArapucaVUV efficiency to VUV light
      double XArapucaVISEff;    //XArapucaVIS efficiency (optical window + cavity)
      double DecayTXArapucaVIS;// Decay time of EJ280 in ns
      std::string ArapucaDataFile; //File containing timing structure for arapucas
      bool ArapucaSinglePEmodel; //Model for single pe response, false for ideal, true for test bench meas
      bool MakeAmpFluctuations;  //Add amplitude fluctuations to the simulation
      double AmpFluctuation; //Model for single pe response, false for ideal, true for test bench meas

      detinfo::LArProperties const* larProp = nullptr; ///< LarProperties service provider.
      double frequency; ///< Optical-clock frequency
      double frequency_Daphne; ///< Optical-clock frequency for daphne readouts	

      CLHEP::HepRandomEngine* engine = nullptr;
      fhicl::ParameterSet HDOpticalWaveformParams;
    };// ConfigurationParameters_t

    //Default constructor
    DigiArapucaSBNDAlg(ConfigurationParameters_t const& config);
    //Default destructor
    ~DigiArapucaSBNDAlg();

    double Baseline()
      {
        return fParams.Baseline;
      }

    void ConstructWaveform(int ch,
                           sim::SimPhotons const& simphotons,
                           std::vector<short unsigned int>& waveform,
                           std::string pdtype,
                           bool is_daphne,
                           double start_time,
                           unsigned n_samples);
    void ConstructWaveformVUVXA(int ch,
                                    std::vector<short unsigned int>& waveform,
                                    std::unordered_map<int, sim::SimPhotons>& DirectPhotonsMap,
                                    std::unordered_map<int, sim::SimPhotons>& ReflectedPhotonsMap,
                                    double start_time,
                                    unsigned n_samples);
    void ConstructWaveformLite(int ch,
                               sim::SimPhotonsLite const& litesimphotons,
                               std::vector<short unsigned int>& waveform,
                               std::string pdtype,
                               bool is_daphne,
                               double start_time,
                               unsigned n_samples);
    void ConstructWaveformLiteVUVXA(int ch,
                                    std::vector<short unsigned int>& waveform,
                                    std::unordered_map<int, sim::SimPhotonsLite>& DirectPhotonsMap,
                                    std::unordered_map<int, sim::SimPhotonsLite>& ReflectedPhotonsMap,
                                    double start_time,
                                    unsigned n_samples);

  private:

    // Declare member data here.
    const ConfigurationParameters_t fParams;

    const double fSampling;        //wave sampling frequency (GHz)
    const double fSampling_Daphne;        //wave sampling frequency (GHz)
    const double fXArapucaVUVEffVis; //VUV XArapuca efficiency to visible light
    const double fXArapucaVUVEffVUV; //VUV XArapuca efficiency to VUV light
    const double fXArapucaVISEff;
    const double fADCSaturationHigh;
    const double fADCSaturationLow;

    CLHEP::HepRandomEngine* fEngine; //!< Reference to art-managed random-number engine
    CLHEP::RandFlat fFlatGen;
    CLHEP::RandPoissonQ fPoissonQGen;
    CLHEP::RandGaussQ fGaussQGen;
    CLHEP::RandExponential fExponentialGen;
    std::unique_ptr<CLHEP::RandGeneral> fTimeXArapucaVUV;// histogram for getting the photon time distribution inside the XArapuca VUV box (considering the optical window)
    std::unique_ptr<CLHEP::RandGeneral> fTimeTPB; // histogram for getting the TPB emission time for visible (x)arapucas

    std::vector<double> fWaveformSP; //single photon pulse vector
    std::vector<double> fWaveformSP_Daphne; //single photon pulse vector
    std::vector<std::vector<double>> fWaveformSP_Daphne_HD; //single photon pulse vector
    
    std::unordered_map< raw::Channel_t, std::vector<double> > fFullWaveforms;

    //HDWaveforms
    std::unique_ptr<opdet::HDOpticalWaveform> fPMTHDOpticalWaveformsPtr;


    void CreatePDWaveform(sim::SimPhotons const& SimPhotons,
                          double t_min,
                          std::vector<double>& wave,
                          std::string pdtype,
                          bool is_daphne);
    void CreatePDWaveformLite(std::map<int, int> const& photonMap,
                              double t_min,
                              std::vector<double>& wave,
                              std::string pdtype,
                              bool is_daphne);
    void SinglePDWaveformCreatorLite(double effT,
                                     std::unique_ptr<CLHEP::RandGeneral>& timeHisto,
                                     std::vector<double>& wave,
                                     std::map<int, int> const& photonMap,
                                     double const& t_min,
                                     bool is_daphne);
    void SinglePDWaveformCreatorLite(double effT,
                                     std::vector<double>& wave,
                                     std::map<int, int> const& photonMap,
                                     double const& t_min,
                                     bool is_daphne);
    void AddSPE(size_t time_bin, std::vector<double>& wave, const std::vector<double>& fWaveformSP, int nphotons); // add single pulse to auxiliary waveform
    void Pulse1PE(std::vector<double>& wave,const double sampling);
    // void produceSER_HD(std::vector<double> *SER_HD, std::vector<double>& SER);
    void AddLineNoise(std::vector<double>& wave);
    void AddDarkNoise(std::vector<double>& wave , std::vector<double>& WaveformSP);
    double FindMinimumTime(sim::SimPhotons const& simphotons);
    double FindMinimumTimeLite(std::map< int, int > const& photonMap);
    void CreateSaturation(std::vector<double>& wave);//Including saturation effects
  };//class DigiArapucaSBNDAlg

  class DigiArapucaSBNDAlgMaker {

  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;

      fhicl::Atom<double> peakTime {
        Name("ArapucaPeakTime"),
        Comment("Single pe: Time of maximum amplitude in the SiPM pulse in ns")
      };

      fhicl::Atom<double> pulseLength {
        Name("ArapucaPulseLength"),
        Comment("Single pe: Time length of SiPM pulse")
      };

      fhicl::Atom<double> meanAmplitude {
        Name("ArapucaMeanAmplitude"),
        Comment("Single pe: mean amplitude of SiPM pulse in mV")
      };

      fhicl::Atom<double> riseTime {
        Name("ArapucaRiseTime"),
        Comment("Single pe: Pulse rise time constant (exponential), from 0.1 to 0.9 of maximum amplitude")
      };

      fhicl::Atom<double> fallTime {
        Name("ArapucaFallTime"),
        Comment("Single pe: Pulse decay time constant (exponential), from 0.1 to 0.9 of maximum amplitude")
      };

      fhicl::Atom<double> baseline {
        Name("ArapucaBaseline"),
        Comment("Waveform baseline in ADC")
      };

      fhicl::Atom<double> baselineRMS {
        Name("ArapucaBaselineRMS"),
        Comment("RMS of the electronics noise fluctuations in ADC counts")
      };

      fhicl::Atom<double> darkNoiseRate {
        Name("ArapucaDarkNoiseRate"),
        Comment("Dark noise rate in Hz")
      };

      fhicl::Atom<double> crossTalk {
        Name("CrossTalk"),
        Comment("Probability for producing a signal of 2 PE in response to 1 photon")
      };

      fhicl::Atom<double> saturationHigh {
        Name("ArapucaSaturationHigh"),
        Comment("Saturation in number of ADC counts")
      };

      fhicl::Atom<double> saturationLow {
        Name("ArapucaSaturationLow"),
        Comment("Saturation in number of ADC counts")
      };

      fhicl::Atom<double> xArapucaVUVEffVis {
        Name("XArapucaVUVEffVis"),
        Comment("XArapuca VUV efficiency (optical window + cavity)")
      };

      fhicl::Atom<double> xArapucaVUVEffVUV {
        Name("XArapucaVUVEffVUV"),
        Comment("XArapuca VUV efficiency (optical window + cavity)")
      };

      fhicl::Atom<double> xArapucaVISEff {
        Name("XArapucaVISEff"),
        Comment("XArapuca VIS efficiency (optical window + cavity)")
      };

      fhicl::Atom<double> decayTXArapucaVIS {
        Name("DecayTXArapucaVIS"),
        Comment("X-Arapuca VIS decay time of EJ280 in ns")
      };

      fhicl::Atom<std::string> arapucaDataFile {
        Name("ArapucaDataFile"),
        Comment("File containing timing distribution for ArapucaVUV (optical window + cavity), ArapucaVIS (optical window + cavity), XArapuca VUV (optical window)")
      };
      
      fhicl::Atom<bool> ArapucasinglePEmodel {
        Name("ArapucaSinglePEmodel"),
        Comment("Model used for single PE response of PMT. =0 is ideal, =1 is from X-TDBoard data (with overshoot)")
      };

      fhicl::Atom<double> DaphneFrequency {
        Name("DaphneFrequency"),
        Comment("Sampling Frequency of the XArapucas with Daphne readouts (SBND Light detection system has 2 readout frequencies). Apsaia readouts read the frec value from LArSoft.")
      };

      fhicl::Atom<bool> makeAmpFluctuations {
        Name("MakeAmpFluctuations"),
        Comment("Include amplitude fluctuations to the simulated wvfs. For each PE, roll a rand number(gaussian distributed) with std=AmpFluctuation.")
      };

      fhicl::Atom<double> ampFluctuation {
        Name("AmpFluctuation"),
        Comment("Std of the gaussian fit (from data) to the 1PE distribution. Relative units i.e.: AmpFluctuation=0.1-> Amp(1PE)= 18+-1.8 ADC counts")
      };

      fhicl::OptionalDelegatedParameter hdOpticalWaveformParams {
        Name("HDOpticalWaveformParamsXARAPUCA"),
        Comment("Parameters used for high definition waveform")
      };


    };    //struct Config

    DigiArapucaSBNDAlgMaker(Config const& config); //Constructor

    std::unique_ptr<DigiArapucaSBNDAlg> operator()(
      detinfo::LArProperties const& larProp,
      detinfo::DetectorClocksData const& clockData,
      CLHEP::HepRandomEngine* engine
      ) const;

  private:
    /// Part of the configuration learned from configuration files.
    DigiArapucaSBNDAlg::ConfigurationParameters_t fBaseConfig;
  }; //class DigiArapucaSBNDAlgMaker

} //namespace

#endif // SBND_OPDETSIM_DIGIARAPUCASBNDALG_HH

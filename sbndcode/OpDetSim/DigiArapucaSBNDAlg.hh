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
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandExponential.h"

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
      double Saturation;     //Saturation in number of p.e.
      double ArapucaVUVEff;   //ArapucaVUV efficiency (optical window + cavity)
      double ArapucaVISEff;   //ArapucaVIS efficiency (optical window + cavity)
      double XArapucaVUVEff;    //XArapucaVUV efficiency (optical window + cavity)
      double XArapucaVISEff;    //XArapucaVIS efficiency (optical window + cavity)
      double DecayTXArapucaVIS;// Decay time of EJ280 in ns
      std::string ArapucaDataFile; //File containing timing structure for arapucas

      detinfo::LArProperties const* larProp = nullptr; ///< LarProperties service provider.
      detinfo::DetectorClocks const* timeService = nullptr; ///< DetectorClocks service provider.
      CLHEP::HepRandomEngine* engine = nullptr;
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
                           double start_time,
                           unsigned n_samples);
    void ConstructWaveformLite(int ch,
                               sim::SimPhotonsLite const& litesimphotons,
                               std::vector<short unsigned int>& waveform,
                               std::string pdtype,
                               double start_time,
                               unsigned n_samples);

  private:

    // Declare member data here.
    ConfigurationParameters_t fParams;

    double fSampling;        //wave sampling frequency (GHz)
    int pulsesize;
    double fArapucaVUVEff;
    double fArapucaVISEff;
    double fXArapucaVUVEff;
    double fXArapucaVISEff;

    double saturation;

    CLHEP::HepRandomEngine* fEngine; //!< Reference to art-managed random-number engine

    CLHEP::RandGeneral* fTimeArapucaVUV; //histogram for getting the photon time distribution inside the Arapuca VUV box (considering the optical window)
    CLHEP::RandGeneral* fTimeArapucaVIS; //histogram for getting the photon time distribution inside the Arapuca VIS box (considering the optical window)
    CLHEP::RandGeneral* fTimeXArapucaVUV; //histogram for getting the photon time distribution inside the XArapuca VUV box (considering the optical window)

    std::vector<double> wsp; //single photon pulse vector
    std::unordered_map< raw::Channel_t, std::vector<double> > fFullWaveforms;

    void CreatePDWaveform(sim::SimPhotons const& SimPhotons,
                          double t_min,
                          std::vector<double>& wave,
                          std::string pdtype);
    void CreatePDWaveformLite(std::map<int, int> const& photonMap,
                              double t_min,
                              std::vector<double>& wave,
                              std::string pdtype);
    void SinglePDWaveformCreatorLite(double effT,
                                     CLHEP::RandGeneral** timeHisto,
                                     std::vector<double>& wave,
                                     std::map<int, int> const& photonMap,
                                     double const& t_min);
    void SinglePDWaveformCreatorLite(double effT,
                                     std::vector<double>& wave,
                                     std::map<int, int> const& photonMap,
                                     double const& t_min);
    void AddSPE(size_t time_bin, std::vector<double>& wave, int nphotons); // add single pulse to auxiliary waveform
    void Pulse1PE(std::vector<double>& wave);
    void AddLineNoise(std::vector<double>& wave);
    void AddDarkNoise(std::vector<double>& wave);
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

      fhicl::Atom<double> voltageToADC {
        Name("ArapucaVoltageToADC"),
        Comment("Voltage to ADC convertion factor")
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

      fhicl::Atom<double> saturation {
        Name("ArapucaSaturation"),
        Comment("Saturation in number of p.e.")
      };

      fhicl::Atom<double> arapucaVUVEff {
        Name("ArapucaVUVEff"),
        Comment("Arapuca VUV efficiency (optical window + cavity)")
      };

      fhicl::Atom<double> arapucaVISEff {
        Name("ArapucaVISEff"),
        Comment("Arapuca VIS efficiency (optical window + cavity)")
      };

      fhicl::Atom<double> xArapucaVUVEff {
        Name("XArapucaVUVEff"),
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
    };    //struct Config

    DigiArapucaSBNDAlgMaker(Config const& config); //Constructor

    std::unique_ptr<DigiArapucaSBNDAlg> operator()(
      detinfo::LArProperties const& larProp,
      detinfo::DetectorClocks const& detClocks,
      CLHEP::HepRandomEngine* engine
      ) const;

  private:
    /// Part of the configuration learned from configuration files.
    DigiArapucaSBNDAlg::ConfigurationParameters_t fBaseConfig;
  }; //class DigiArapucaSBNDAlgMaker

} //namespace

#endif // SBND_OPDETSIM_DIGIARAPUCASBNDALG_HH

////////////////////////////////////////////////////////////////////////
// Class:       opDetDigitizerSBND
// Module Type: producer
// File:        opDetDigitizerSBND_module.cc
//
// This module produces digitized waveforms of the optical detectors
// Created by L. Paulucci, F. Marinho, and I.L. de Icaza
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/TableFragment.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nurandom/RandomUtils/NuRandomService.h"
#include "CLHEP/Random/JamesRandom.h"

#include <memory>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include <sstream>
#include <fstream>
#include <thread>
#include <cstdlib>
#include <stdexcept>
#include <iostream>

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"

#include "TMath.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TF1.h"

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "sbndcode/OpDetSim/DigiArapucaSBNDAlg.hh"
#include "sbndcode/OpDetSim/DigiPMTSBNDAlg.hh"
#include "sbndcode/OpDetSim/opDetSBNDTriggerAlg.hh"
#include "sbndcode/OpDetSim/opDetDigitizerWorker.hh"

#include "TriggerEmulationService.h"

namespace opdet {

  /*
  * This module simulates the digitization of SBND photon detectors response.
  *
  * The module has an interface to the simulation algorithms for PMTs and Arapucas,
  * opdet::DigiPMTSBNDAlg e opdet::DigiArapucaSBNDAlg.
  *
  * Input
  * ======
  * The module utilizes as input a collection of `sim::SimPhotons` or `sim::SimPhotonsLite`, each
  * containing the photons propagated to a single optical detector channel.
  *
  * Output
  * =======
  * A collection of optical detector waveforms (`std::vector<raw::OpDetWaveform>`) is produced.
  *
  * Requirements
  * =============
  * This module currently requires LArSoft services:
  * * `DetectorClocksService` for timing conversions and settings
  * * `LArPropertiesService` for the scintillation yield(s)
  *
  */

  class opDetDigitizerSBND;

  class opDetDigitizerSBND : public art::EDProducer {
  public:
    struct Config {
      using Comment = fhicl::Comment;
      using Name = fhicl::Name;

      fhicl::Atom<art::InputTag> InputModuleName {
        Name("InputModule"),
        Comment("Simulated photons to be digitized")
      };
      fhicl::Atom<double> WaveformSize {
        Name("WaveformSize"),
        Comment("Value to initialize the waveform vector in ns. It is resized in the algorithms according to readout window of PDs")
      };
      fhicl::Atom<bool> UseSimPhotonsLite {
        Name("UseSimPhotonsLite"),
        Comment("Whether SimPhotonsLite or SimPhotons will be used")
      };

      fhicl::Atom<bool> ApplyTriggers {
        Name("ApplyTriggers"),
        Comment("Whether to apply trigger algorithm to waveforms"),
        true
      };

      fhicl::Atom<unsigned> NThreads {
        Name("NThreads"),
        Comment("Number of threads to split waveform process into. Defaults to 1.\
                     Set 0 to autodetect. Autodection will first check $SBNDCODE_OPDETSIM_NTHREADS for number of threads. \
                     If this is not set, then NThreads is set to the number of hardware cores on the host machine."),
        1
      };

      fhicl::Atom<int> ticksPerSlice{
        Name("ticksPerSlice"),
        Comment("Number of ticks for width of sliced waveform. Width of waveform in data is 5000 samples = 10us, since 1 tick of PMT data is 2ns of data.")
      };
      
      fhicl::Atom<float> PercentTicksBeforeCross{
        Name("PercentTicksBeforeCross"),
        Comment("Given a certain width of the waveform, how much of the waveform width is before the crossing point. To match data, this should be 0.2.")
      };

      fhicl::Atom<int> MonThreshold{
        Name("MonThreshold"),
        Comment("Threshold")
      };

      fhicl::Atom<int> MonPulseThresh{
        Name("MonPulseThresh"),
        Comment("Threshold for MonPulse (to determine interesting trigger)")
      };

      fhicl::Atom<bool> SaveNewPlots{
        Name("SaveNewPlots"),
        Comment("Save plots of triggered waveforms and MonPulse with new trigger logic.")
      };

      fhicl::Atom<bool> SaveOldPlots{
        Name("SaveOldPlots"),
        Comment("Save plots of triggered waveforms with old trigger logic.")
      };

      fhicl::TableFragment<opdet::DigiPMTSBNDAlgMaker::Config> pmtAlgoConfig;
      fhicl::TableFragment<opdet::DigiArapucaSBNDAlgMaker::Config> araAlgoConfig;
      fhicl::TableFragment<opdet::opDetSBNDTriggerAlg::Config> trigAlgoConfig;
    }; // struct Config

    using Parameters = art::EDProducer::Table<Config>;

    explicit opDetDigitizerSBND(Parameters const& config);
    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.
    // Add a destructor to deal with random number generator pointer
    ~opDetDigitizerSBND();

    // Plugins should not be copied or assigned.
    opDetDigitizerSBND(opDetDigitizerSBND const &) = delete;
    opDetDigitizerSBND(opDetDigitizerSBND &&) = delete;
    opDetDigitizerSBND & operator = (opDetDigitizerSBND const &) = delete;
    opDetDigitizerSBND & operator = (opDetDigitizerSBND &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;
    std::vector<raw::OpDetWaveform> sliceWaveforms(std::vector<raw::OpDetWaveform> fWaveforms,
                                                        int WaveIndex,
                                                        std::vector<int> *MonPulse,
                                                        int MonPulseThresh,
                                                        double tickPeriod,
                                                        int ticksPerSlice,
                                                        float PercentTicksBeforeCross,
                                                        int PMTPerBoard); 
    std::vector<std::vector<int>> sliceMonPulse(std::vector<int> *MonPulse,
                                                        int MonPulseThresh,
                                                        double tickPeriod,
                                                        int ticksPerSlice,
                                                        float PercentTicksBeforeCross); 
    void PlotWaveforms(const std::vector<raw::OpDetWaveform>& waveforms,
                                           const std::string& basename);
    opdet::sbndPDMapAlg map; //map for photon detector types
    unsigned int nChannels = map.size();
    std::vector<raw::OpDetWaveform> fWaveforms; // holder for un-triggered waveforms

  private:
    bool fApplyTriggers;
    art::InputTag fInputModuleName;
    std::unordered_map< raw::Channel_t, std::vector<double> > fFullWaveforms;

    bool fUseSimPhotonsLite;
    unsigned fPMTBaseline;
    unsigned fArapucaBaseline;
    unsigned fNThreads;
    // digitizer workers
    std::vector<opdet::opDetDigitizerWorker> fWorkers;
    std::vector<std::vector<raw::OpDetWaveform>> fTriggeredWaveforms;
    std::vector<std::thread> fWorkerThreads;
    std::vector<std::vector<raw::OpDetWaveform>> fSlicedWaveformsAll;
    //std::vector<std::vector<int>> fMonPulsesAll;

    // product containers
    std::vector<art::Handle<std::vector<sim::SimPhotonsLite>>> fPhotonLiteHandles;
    std::vector<art::Handle<std::vector<sim::SimPhotons>>> fPhotonHandles;

    // sync stuff
    opdet::opDetDigitizerWorker::Semaphore fSemStart;
    opdet::opDetDigitizerWorker::Semaphore fSemFinish;
    bool fFinished;

    // trigger algorithm
    opdet::opDetSBNDTriggerAlg fTriggerAlg;

    int ticksPerSlice;
    float PercentTicksBeforeCross; 
    int MonThreshold;
    int MonPulseThresh;
    bool SaveNewPlots;
    bool SaveOldPlots;
  };

  opDetDigitizerSBND::opDetDigitizerSBND(Parameters const& config)
    : EDProducer{config}
    , fApplyTriggers(config().ApplyTriggers())
    , fUseSimPhotonsLite(config().UseSimPhotonsLite())
    , fPMTBaseline(config().pmtAlgoConfig().pmtbaseline())
    , fArapucaBaseline(config().araAlgoConfig().baseline())
    , fTriggerAlg(config().trigAlgoConfig())
    , ticksPerSlice(config().ticksPerSlice())
    , PercentTicksBeforeCross(config().PercentTicksBeforeCross())
    , MonThreshold(config().MonThreshold())
    , MonPulseThresh(config().MonPulseThresh())
    , SaveNewPlots(config().SaveNewPlots())
    , SaveOldPlots(config().SaveOldPlots())
  {
    opDetDigitizerWorker::Config wConfig( config().pmtAlgoConfig(), config().araAlgoConfig());

    fNThreads = config().NThreads();
    if (fNThreads == 0) { // autodetect -- first check env var
      const char *env = std::getenv("SBNDCODE_OPDETSIM_NTHREADS");
      // try to parse into positive integer
      if (env != NULL) {
        try {
          int n_threads = std::stoi(env);
          if (n_threads <= 0) {
            throw std::invalid_argument("Expect positive integer");
          }
          fNThreads = n_threads;
        }
        catch (...) {
          mf::LogError("OpDetDigitizer") << "Unable to parse number of threads "
                                         << "in environment variable (SBNDCODE_OPDETSIM_NTHREADS): (" << env << ").\n"
                                         << "Setting Number opdet threads to 1." << std::endl;
          fNThreads = 1;
        }
      }
    }

    if (fNThreads == 0) { // autodetect -- now try to get number of cpu's
      fNThreads = std::thread::hardware_concurrency();
    }
    if (fNThreads == 0) { // autodetect failed
      fNThreads = 1;
    }
    mf::LogInfo("OpDetDigitizer") << "Digitizing on n threads: " << fNThreads << std::endl;

    wConfig.nThreads = fNThreads;

    wConfig.UseSimPhotonsLite = config().UseSimPhotonsLite();
    wConfig.InputModuleName = config().InputModuleName();

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clockData);
    wConfig.Sampling = (clockData.OpticalClock().Frequency()) / 1000.0; //in GHz
    wConfig.Sampling_Daphne =  config().araAlgoConfig().DaphneFrequency() / 1000.0; //in GHz
    wConfig.EnableWindow = fTriggerAlg.TriggerEnableWindow(clockData, detProp); // us
    wConfig.Nsamples = (wConfig.EnableWindow[1] - wConfig.EnableWindow[0]) * 1000. /*us -> ns*/ * wConfig.Sampling /* GHz */;
    wConfig.Nsamples_Daphne = (wConfig.EnableWindow[1] - wConfig.EnableWindow[0]) * 1000. /*us -> ns*/ * wConfig.Sampling_Daphne /* GHz */;

    fFinished = false;

    fWorkers.reserve(fNThreads);
    fTriggeredWaveforms.reserve(fNThreads);
    for (unsigned i = 0; i < fNThreads; i++) {
      // Set random number gen seed from the NuRandomService
      art::ServiceHandle<rndm::NuRandomService> seedSvc;
      CLHEP::HepJamesRandom *engine = new CLHEP::HepJamesRandom;
      seedSvc->registerEngine(rndm::NuRandomService::CLHEPengineSeeder(engine), "opDetDigitizerSBND" + std::to_string(i));

      fTriggeredWaveforms.emplace_back();

      // setup worker
      fWorkers.emplace_back(i, wConfig, engine, fTriggerAlg);
      fWorkers[i].SetPhotonLiteHandles(&fPhotonLiteHandles);
      fWorkers[i].SetPhotonHandles(&fPhotonHandles);
      fWorkers[i].SetWaveformHandle(&fWaveforms);
      fWorkers[i].SetTriggeredWaveformHandle(&fTriggeredWaveforms[i]);

      // start worker thread
      fWorkerThreads.emplace_back(opdet::opDetDigitizerWorkerThread,
                                  std::cref(fWorkers[i]),
                                  clockData,
                                  std::ref(fSemStart),
                                  std::ref(fSemFinish),
                                  fApplyTriggers,
                                  &fFinished);
    }

    // Call appropriate produces<>() functions here.
    produces< std::vector< raw::OpDetWaveform > >();
    produces<bool>("triggerEmulation");
    produces<int>("pairsOverThreshold");
    produces< std::vector< raw::OpDetWaveform > >("slicedWaveforms");
    //produces< std::vector< std::vector<int> > >("MonPulses");
  }

  opDetDigitizerSBND::~opDetDigitizerSBND()
  {
    // cleanup all of the workers
    fFinished = true;
    opdet::StartopDetDigitizerWorkers(fNThreads, fSemStart);

    // join the threads
    for (std::thread &thread : fWorkerThreads) thread.join();

  }

  void opDetDigitizerSBND::produce(art::Event & e)
  {
    std::unique_ptr< std::vector< raw::OpDetWaveform > > pulseVecPtr(std::make_unique< std::vector< raw::OpDetWaveform > > ());
    // Implementation of required member function here.
    mf::LogInfo("opDetDigitizer") << "Event: " << e.id().event() << std::endl;

    // setup the waveforms
    fWaveforms = std::vector<raw::OpDetWaveform> (nChannels);

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e, clockData);

    if (fUseSimPhotonsLite) {
      fPhotonLiteHandles.clear();
      //Get *ALL* SimPhotonsCollectionLite from Event
      fPhotonLiteHandles = e.getMany<std::vector<sim::SimPhotonsLite>>();
      if (fPhotonLiteHandles.size() == 0)
        mf::LogError("OpDetDigitizer") << "sim::SimPhotonsLite not found -> No Optical Detector Simulation!\n";
    }
    else {
      fPhotonHandles.clear();
      //Get *ALL* SimPhotonsCollection from Event
      fPhotonHandles = e.getMany<std::vector<sim::SimPhotons>>();
      if (fPhotonHandles.size() == 0)
        mf::LogError("OpDetDigitizer") << "sim::SimPhotons not found -> No Optical Detector Simulation!\n";
    }
    // Start the workers!
    // Run the digitizer over the full readout window
    opdet::StartopDetDigitizerWorkers(fNThreads, fSemStart);
    opdet::WaitopDetDigitizerWorkers(fNThreads, fSemFinish);

    if (fApplyTriggers) {

      // clear previous
      fSlicedWaveformsAll.clear();
      //fMonPulsesAll.clear();
      // find the trigger locations for the waveforms using the LArService
      if (!fWaveforms.empty()) {
          // Implement service
          art::ServiceHandle<art::TFileService> tfs;
          art::ServiceHandle<calib::TriggerEmulationService> fTriggerService;

          int PMTPerBoard = fTriggerService->getPMTPerBoard();
          int fTotalCAENBoards = fTriggerService->getTotalCAENBoards();

          int TotalFlash = fWaveforms.size() / (fTotalCAENBoards * PMTPerBoard);

          int numPairsOverThreshold = 0;

          // Loop through by flash -> compatible with ConstructMonPulse logic
          for (int FlashCounter = 0; FlashCounter < TotalFlash; ++FlashCounter) {
              int WaveIndex = FlashCounter*PMTPerBoard;
              int WaveformSize = fWaveforms[WaveIndex].size();

              std::vector<int>* MonPulse = new std::vector<int>(WaveformSize, 0);
    
              int pairThisFlash = 0;
              // Send 3ms waveforms to ConstructMonPulse
              fTriggerService->ConstructMonPulse(fWaveforms, MonThreshold, MonPulse, false, FlashCounter, &pairThisFlash);
              numPairsOverThreshold = numPairsOverThreshold + pairThisFlash;

              double tickPeriod = sampling_rate(clockData);

              std::vector<raw::OpDetWaveform> fSlicedWaveforms = sliceWaveforms(fWaveforms, WaveIndex, MonPulse, MonPulseThresh, tickPeriod, ticksPerSlice, PercentTicksBeforeCross, PMTPerBoard);
              std::vector<std::vector<int>> fSlicedMonPulse = sliceMonPulse(MonPulse, MonPulseThresh, tickPeriod, ticksPerSlice, PercentTicksBeforeCross);

              fSlicedWaveformsAll.push_back(std::move(fSlicedWaveforms));
              //fMonPulsesAll.push_back(*MonPulse);

              if (SaveNewPlots) {
                  // Save histograms
                  // Sliced waveforms
                  for (size_t j; j < fSlicedWaveformsAll.size(); ++j) { 
                      std::stringstream plotname2; 
                      plotname2 << "Sliced_waveforms_" << e.id().event() << "_Mon_" << MonThreshold << "_" << FlashCounter << "_slice" << j;
                      PlotWaveforms(fSlicedWaveformsAll[j], plotname2.str());
                  }
                  // Long MonPulse
                  std::stringstream histname;
                  histname << "Long_event_" << e.id().event() << "_Mon_" << MonThreshold << "_" << FlashCounter;
                  TH1D* MonHist = tfs->make<TH1D>(histname.str().c_str(), histname.str().c_str(),
                                                  MonPulse->size(), 0.0, MonPulse->size() - 1);
                  for (size_t i = 0; i < MonPulse->size(); i++) {
                      MonHist->SetBinContent(i + 1, (*MonPulse)[i]);
                  }
                  // Sliced MonPulse
                  for (size_t idx = 0; idx < fSlicedMonPulse.size(); ++idx) {
                      auto const& vec = fSlicedMonPulse[idx];
                      std::stringstream histname;
                      histname << "Sliced_event_" << e.id().event() << "_Mon_" << MonThreshold << "_" << FlashCounter << "_slice" << idx;

                      TH1D* MonHist = tfs->make<TH1D>(histname.str().c_str(), histname.str().c_str(),
                                                      vec.size(), 0.0, vec.size() - 1);
                      for (size_t i = 0; i < vec.size(); i++) {
                          MonHist->SetBinContent(i + 1, vec[i]);
                      }
                  }
              }   
          delete MonPulse;
          }

          // find the trigger locations for the waveforms - old version, keeping for validation
          for (const raw::OpDetWaveform &waveform : fWaveforms) {
            raw::Channel_t ch = waveform.ChannelNumber();
            // skip light channels which don't correspond to readout channels
            if (ch == std::numeric_limits<raw::Channel_t>::max() /* "NULL" value*/) {
              continue;
            }
            raw::ADC_Count_t baseline = (map.isPDType(ch, "pmt_uncoated") || map.isPDType(ch, "pmt_coated")) ?
                                        fPMTBaseline : fArapucaBaseline;
            fTriggerAlg.FindTriggerLocations(clockData, detProp, waveform, baseline);
          }

          // combine the triggers
          fTriggerAlg.MergeTriggerLocations();
          // Start the workers!
          // Apply the trigger locations
          opdet::StartopDetDigitizerWorkers(fNThreads, fSemStart);
          opdet::WaitopDetDigitizerWorkers(fNThreads, fSemFinish);

          // plot fTriggeredWaveforms
          if (SaveOldPlots) {
              for (size_t j; j < fTriggeredWaveforms.size(); ++j) { 
                  std::stringstream plotnameTW; 
                  plotnameTW << "Triggered_waveforms_" << e.id().event() << "_Mon_" << MonThreshold;
                  PlotWaveforms(fTriggeredWaveforms[j], plotnameTW.str());
              }
          }

          // put triggered waveforms in the event (old trigger logic)
          for (std::vector<raw::OpDetWaveform> &waveforms : fTriggeredWaveforms) {
            // move these waveforms into the pulseVecPtr
            pulseVecPtr->reserve(pulseVecPtr->size() + waveforms.size());
            std::move(waveforms.begin(), waveforms.end(), std::back_inserter(*pulseVecPtr));
          }
          // clean up the vector
          for (unsigned i = 0; i < fTriggeredWaveforms.size(); i++) {
            fTriggeredWaveforms[i] = std::vector<raw::OpDetWaveform>();
          }
          // put the waveforms in the event
          e.put(std::move(pulseVecPtr));
          // clear out the triggers
          fTriggerAlg.ClearTriggerLocations();


          // put boolean trigger result in the event
          bool passedTrigger = false;
          // passes trigger if any of the fSlicedWaveforms have size > 0 
          for (auto wav : fSlicedWaveformsAll) if (wav.size() > 0) passedTrigger = true;
          auto triggerFlag = std::make_unique<bool>(passedTrigger);
          e.put(std::move(triggerFlag), "triggerEmulation");


          // put trigger pair result in the event
          auto pairFlag = std::make_unique<int>(numPairsOverThreshold);
          e.put(std::move(pairFlag), "pairsOverThreshold");


          // put sliced waveforms in the event
          std::unique_ptr< std::vector< raw::OpDetWaveform > > SlicedWaveformsPtr(std::make_unique< std::vector< raw::OpDetWaveform > > ());
          for (std::vector<raw::OpDetWaveform> &waveforms : fSlicedWaveformsAll) {
            // move sliced waveforms into the SlicedWaveformsPtr
            SlicedWaveformsPtr->reserve(SlicedWaveformsPtr->size() + waveforms.size());
            std::move(waveforms.begin(), waveforms.end(), std::back_inserter(*SlicedWaveformsPtr));
          }
          // clean up the vector
          for (unsigned i = 0; i < fSlicedWaveformsAll.size(); i++) {
            fSlicedWaveformsAll[i] = std::vector<raw::OpDetWaveform>();
          }
          // put the waveforms in the event
          e.put(std::move(SlicedWaveformsPtr), "slicedWaveforms");


          // put MonPulses in the event
          /*auto MonPulsesPtr = std::make_unique<std::vector<int>>();
          for (auto& pulses : fMonPulsesAll) { 
              MonPulsesPtr->reserve(MonPulsesPtr->size() + pulses.size());
              std::move(pulses.begin(), pulses.end(), std::back_inserter(*MonPulsesPtr));
          }
          // clean up the vector
          for (auto& pulses : fMonPulsesAll) pulses.clear();
          // put the waveforms in the event
          e.put(std::move(MonPulsesPtr), "MonPulses");*/

      } else std::cout << "Empty waveforms found on event " << e.id().event() << "  " << fWaveforms.empty() << std::endl; 
    }
    else {
      std::cout<<"ApplyTriggers is false"<<std::endl;
      // put the full waveforms in the event
      for (const raw::OpDetWaveform &waveform : fWaveforms) {
        if (waveform.ChannelNumber() == std::numeric_limits<raw::Channel_t>::max() /* "NULL" value*/) {
          continue;
        }
        pulseVecPtr->push_back(waveform);
      }
      e.put(std::move(pulseVecPtr));
    }

    // clear out the full waveforms
    fWaveforms.clear();

  }//produce end


  // sliced MonPulse function: same logic as sliceWaveforms function
  std::vector<std::vector<int>> opDetDigitizerSBND::sliceMonPulse(std::vector<int> *MonPulse, 
                                                        int MonPulseThresh, 
                                                        double tickPeriod, 
                                                        int ticksPerSlice, 
                                                        float PercentTicksBeforeCross
  )
  { 
              // Slice up each waveform into 10us chunks based on if "interesting" or not
              // before and after crossing point (default is ~20% and ~80%)
              int ticksBeforeCross = static_cast<int>(std::round(PercentTicksBeforeCross*ticksPerSlice));
              int ticksAfterCross = ticksPerSlice-ticksBeforeCross;
              // find interesting area
              std::vector<std::pair<int,int>> interestIntervals;
              std::vector<int> crossingPoints;
              // clear
              interestIntervals.clear();
              crossingPoints.clear();

              bool interest = false;
              for (int i = 0; i < (int)MonPulse->size(); ++i) {
                  // find crossing point
                  if ((*MonPulse)[i] > MonPulseThresh && interest == false) {
                      crossingPoints.push_back(i);
                      interest = true;
                  } else if ((*MonPulse)[i] <= MonPulseThresh) interest = false;
              }    

              // create 10us slices around crossingPoints
              for (int j = 0; j < (int)crossingPoints.size(); ++j) {
    
                  // if near end of full waveform
                  if (crossingPoints[j]+ticksAfterCross > static_cast<int>(MonPulse->size())) {
                      interestIntervals.emplace_back(static_cast<int>(MonPulse->size())-ticksPerSlice, static_cast<int>(MonPulse->size()));
                  } 
                  // if near beginning of full waveform
                  else if (crossingPoints[j]-ticksBeforeCross < 0) {
                      interestIntervals.emplace_back(0, ticksPerSlice);
                  }  
                  else {
                      // check if overlaps with previous interval
                      if (!interestIntervals.empty()) {
                          if (crossingPoints[j]-ticksBeforeCross < interestIntervals.back().second) {
                              // if overlaps, extend interval
                              interestIntervals.back() = {interestIntervals.back().first, crossingPoints[j]+ticksAfterCross};
                          // if does not overlap, use typical interval length
                          } else interestIntervals.emplace_back(crossingPoints[j]-ticksBeforeCross, crossingPoints[j]+ticksAfterCross);
                      // if first, use typical interval length
                      } else interestIntervals.emplace_back(crossingPoints[j]-ticksBeforeCross, crossingPoints[j]+ticksAfterCross);
                  }
              }

              std::vector<std::vector<int>> fSlicedMonPulses;
              fSlicedMonPulses.clear();
              // loop through intervals
              for (auto [start, end] : interestIntervals) { 
                  std::vector<int> sliceMonPulse((*MonPulse).begin() + start, (*MonPulse).begin() + end);
                  if (!sliceMonPulse.empty()) {
                      fSlicedMonPulses.push_back(std::move(sliceMonPulse));
                  }
              }

      return fSlicedMonPulses;
  }


  // sliceWaveforms function
  std::vector<raw::OpDetWaveform> opDetDigitizerSBND::sliceWaveforms(std::vector<raw::OpDetWaveform> fWaveforms, 
                                                        int WaveIndex, 
                                                        std::vector<int> *MonPulse, 
                                                        int MonPulseThresh, 
                                                        double tickPeriod, 
                                                        int ticksPerSlice, 
                                                        float PercentTicksBeforeCross,
                                                        int PMTPerBoard
  )
  { 
              // Slice up each waveform into 10us chunks based on if "interesting" or not
              // how many ticks correspond to 10us
              // before and after crossing point (default is ~20% and ~80%)
              int ticksBeforeCross = static_cast<int>(std::round(PercentTicksBeforeCross*ticksPerSlice));
              int ticksAfterCross = ticksPerSlice-ticksBeforeCross;
              // find interesting area
              std::vector<std::pair<int,int>> interestIntervals;
              std::vector<int> crossingPoints;
              // clear initial variables
              crossingPoints.clear();
              interestIntervals.clear();

              bool interest = false;
              for (int i = 0; i < (int)MonPulse->size(); ++i) {
                  // find crossing point
                  if ((*MonPulse)[i] > MonPulseThresh && interest == false) {
                      crossingPoints.push_back(i);
                      interest = true;
                  } else if ((*MonPulse)[i] <= MonPulseThresh) interest = false;
              }    

              // create 10us slices around crossingPoints
              for (int j = 0; j < (int)crossingPoints.size(); ++j) {
    
                  // if near end of full waveform
                  if (crossingPoints[j]+ticksAfterCross > static_cast<int>(MonPulse->size())) {
                      interestIntervals.emplace_back(static_cast<int>(MonPulse->size())-ticksPerSlice, static_cast<int>(MonPulse->size()));
                  } 
                  // if near beginning of full waveform
                  else if (crossingPoints[j]-ticksBeforeCross < 0) {
                      interestIntervals.emplace_back(0, ticksPerSlice);
                  }  
                  else {
                      // check if overlaps with previous interval
                      if (!interestIntervals.empty()) {
                          if (crossingPoints[j]-ticksBeforeCross < interestIntervals.back().second) {
                              // if overlaps, extend interval
                              interestIntervals.back() = {interestIntervals.back().first, crossingPoints[j]+ticksAfterCross};
                          // if does not overlap, use typical interval length
                          } else interestIntervals.emplace_back(crossingPoints[j]-ticksBeforeCross, crossingPoints[j]+ticksAfterCross);
                      // if first, use typical interval length
                      } else interestIntervals.emplace_back(crossingPoints[j]-ticksBeforeCross, crossingPoints[j]+ticksAfterCross);
                  }
              }

              std::vector<raw::OpDetWaveform> fSlicedWaveforms;
              fSlicedWaveforms.clear();
              // loop through channels
              for (int chan = 0; chan < PMTPerBoard; ++chan) {
                  const raw::OpDetWaveform& wf = fWaveforms[WaveIndex + chan];

                  for (auto [start, end] : interestIntervals) { 
                      double sliceTime = wf.TimeStamp() + start * tickPeriod;
                      std::vector<uint16_t> sliceData(wf.begin() + start, wf.begin() + end);

                      if (!sliceData.empty()) {
                          raw::OpDetWaveform slice(sliceTime, wf.ChannelNumber(), sliceData);
                          fSlicedWaveforms.push_back(std::move(slice));
                      }
                  }
              }
              
      return fSlicedWaveforms;
  }

  void opDetDigitizerSBND::PlotWaveforms(const std::vector<raw::OpDetWaveform>& waveforms,
                                           const std::string& basename)
  {
      art::ServiceHandle<art::TFileService> tfs;

      for (size_t i = 0; i < waveforms.size(); ++i) {
          const auto& wf = waveforms[i];

          // unique name per event + waveform index + channel
          std::stringstream histName;
          histName << basename << "_wf" << i << "_chan" << wf.ChannelNumber();

          TH1F* h = tfs->make<TH1F>(histName.str().c_str(),
                                   Form("OpDet waveform (chan %d, wf %zu)", wf.ChannelNumber(), i),
                                   wf.size(), 0, wf.size());

          for (size_t tick = 0; tick < wf.size(); ++tick) {
              h->SetBinContent(tick + 1, wf[tick]);
          }
      }
  }

  DEFINE_ART_MODULE(opdet::opDetDigitizerSBND)

}//closing namespace

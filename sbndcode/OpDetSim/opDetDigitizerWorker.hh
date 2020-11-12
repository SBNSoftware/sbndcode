////////////////////////////////////////////////////////////////////////
// Class:       opDetDigitizerWorker
//
// This module handles the calls to the digitization functions
// Created by G. Putnam and I.L. de Icaza
////////////////////////////////////////////////////////////////////////

#ifndef SBND_OPDETSIM_OPDETDIGITIZERWORKER_HH
#define SBND_OPDETSIM_OPDETDIGITIZERWORKER_HH

#include <unordered_map>
#include <vector>
#include <mutex>
#include <condition_variable>
// #include <memory>

#include "sbndcode/OpDetSim/sbndPDMapAlg.hh"
#include "sbndcode/OpDetSim/DigiArapucaSBNDAlg.hh"
#include "sbndcode/OpDetSim/DigiPMTSBNDAlg.hh"
#include "sbndcode/OpDetSim/opDetSBNDTriggerAlg.hh"
namespace detinfo {
  class DetectorClocksData;
}

namespace opdet {

  class opDetDigitizerWorker {
  public:
    class Config {
    public:
      //arapuca and PMT digitization algorithms
      opdet::DigiPMTSBNDAlgMaker makePMTDigi;
      opdet::DigiArapucaSBNDAlgMaker makeArapucaDigi;

      opdet::sbndPDMapAlg pdsMap;  //map for photon detector types
      unsigned int nChannels = pdsMap.size();

      unsigned nThreads;

      art::InputTag InputModuleName;
      bool UseSimPhotonsLite; // SimPhotons have more information that SimPhotonsLite

      std::array<double, 2> EnableWindow;
      double Sampling;       //wave sampling frequency (GHz)
      unsigned int Nsamples; //Samples per waveform

      Config(const opdet::DigiPMTSBNDAlgMaker::Config &pmt_config, const opdet::DigiArapucaSBNDAlgMaker::Config &arapuca_config);
    };

    class Semaphore {
    public:
      Semaphore(unsigned count_ = 0): count(count_) {}
      void increment(unsigned n = 1);
      void decrement(unsigned n = 1);

    private:
      std::mutex mtx;
      std::condition_variable cv;
      unsigned count;
    };

    opDetDigitizerWorker(unsigned no, const Config &config, CLHEP::HepRandomEngine *Engine, const opDetSBNDTriggerAlg &trigger_alg);
    ~opDetDigitizerWorker();

    void SetPhotonLiteHandles(const std::vector<art::Handle<std::vector<sim::SimPhotonsLite>>> *PhotonLiteHandles)
    {
      fPhotonLiteHandles = PhotonLiteHandles;
    }
    void SetPhotonHandles(const std::vector<art::Handle<std::vector<sim::SimPhotons>>> *PhotonHandles)
    {
      fPhotonHandles = PhotonHandles;
    }
    void SetWaveformHandle(std::vector<raw::OpDetWaveform> *Waveforms)
    {
      fWaveforms = Waveforms;
    }
    void SetTriggeredWaveformHandle(std::vector<raw::OpDetWaveform> *Waveforms)
    {
      fTriggeredWaveforms = Waveforms;
    }

    void Start(detinfo::DetectorClocksData const& clockData) const;
    void ApplyTriggerLocations(detinfo::DetectorClocksData const& clockData) const;

  private:
    unsigned NChannelsToProcess(unsigned n) const;
    unsigned StartChannelToProcess(unsigned n) const;
    void CreateDirectPhotonMap(
      std::unordered_map<int, sim::SimPhotons>& directPhotonsOnPMTS,
      std::vector<art::Handle<std::vector<sim::SimPhotons>>> photon_handles) const;
    void CreateDirectPhotonMapLite(
      std::unordered_map<int, sim::SimPhotonsLite>& directPhotonsOnPMTS,
      std::vector<art::Handle<std::vector<sim::SimPhotonsLite>>> photon_handles) const;
    void MakeWaveforms(
      opdet::DigiPMTSBNDAlg *pmtDigitizer,
      opdet::DigiArapucaSBNDAlg *arapucaDigitizer) const;

    Config fConfig;
    unsigned fThreadNo;
    CLHEP::HepRandomEngine *fEngine;
    const opDetSBNDTriggerAlg &fTriggerAlg;

    const std::vector<art::Handle<std::vector<sim::SimPhotonsLite>>> *fPhotonLiteHandles;
    const std::vector<art::Handle<std::vector<sim::SimPhotons>>> *fPhotonHandles;
    std::vector<raw::OpDetWaveform> *fWaveforms;
    std::vector<raw::OpDetWaveform> *fTriggeredWaveforms;
  };

  void StartopDetDigitizerWorkers(unsigned n_workers, opDetDigitizerWorker::Semaphore &sem_start);
  void WaitopDetDigitizerWorkers(unsigned n_workers, opDetDigitizerWorker::Semaphore &sem_finish);
  void opDetDigitizerWorkerThread(const opDetDigitizerWorker &worker,
                                  detinfo::DetectorClocksData const& clockData,
                                  opDetDigitizerWorker::Semaphore &sem_start,
                                  opDetDigitizerWorker::Semaphore &sem_finish,
                                  bool ApplyTriggerLocations, bool *finished);

} // end namespace opdet

#endif // SBND_OPDETSIM_OPDETDIGITIZERWORKER_HH

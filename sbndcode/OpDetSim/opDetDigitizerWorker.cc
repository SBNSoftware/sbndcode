// TODO: plenty of refactoring potential in here! ~icaza

#include "larcore/CoreUtils/ServiceUtil.h"
#include "sbndcode/OpDetSim/opDetDigitizerWorker.hh"

opdet::opDetDigitizerWorker::Config::Config(const opdet::DigiPMTSBNDAlgMaker::Config &pmt_config,
                                            const opdet::DigiArapucaSBNDAlgMaker::Config &arapuca_config):
  makePMTDigi(pmt_config),
  makeArapucaDigi(arapuca_config) {}

opdet::opDetDigitizerWorker::opDetDigitizerWorker(unsigned no,
                                                  const Config &config,
                                                  CLHEP::HepRandomEngine *Engine,
                                                  const opDetSBNDTriggerAlg &trigger_alg):
  fConfig(config),
  fThreadNo(no),
  fEngine(Engine),
  fTriggerAlg(trigger_alg)
{}

void opdet::opDetDigitizerWorkerThread(const opdet::opDetDigitizerWorker &worker,
                                       detinfo::DetectorClocksData const& clockData,
                                       opdet::opDetDigitizerWorker::Semaphore &sem_start,
                                       opdet::opDetDigitizerWorker::Semaphore &sem_finish,
                                       bool ApplyTriggerLocations,
                                       bool *finished)
{

  bool do_apply_trigger_locations = false;
  while (1) {
    sem_start.decrement();

    if (*finished) break;

    if (do_apply_trigger_locations) {
      worker.ApplyTriggerLocations(clockData);
      do_apply_trigger_locations = false;
    }
    else {
      worker.Start(clockData);
      do_apply_trigger_locations = ApplyTriggerLocations;
    }

    sem_finish.increment();
  }

}

void opdet::StartopDetDigitizerWorkers(unsigned n_workers,
                                       opdet::opDetDigitizerWorker::Semaphore &sem_start)
{
  sem_start.increment(n_workers);
}

void opdet::WaitopDetDigitizerWorkers(unsigned n_workers,
                                      opdet::opDetDigitizerWorker::Semaphore &sem_finish)
{
  sem_finish.decrement(n_workers);
}

void opdet::opDetDigitizerWorker::Semaphore::increment(unsigned n)
{
  std::unique_lock<std::mutex> lock(mtx);
  count += n;
  cv.notify_all();
}

void opdet::opDetDigitizerWorker::Semaphore::decrement(unsigned n)
{
  std::unique_lock<std::mutex> lock(mtx);

  while (count < n) {
    cv.wait(lock);
  }
  count -= n;
}

unsigned opdet::opDetDigitizerWorker::NChannelsToProcess(unsigned n) const
{
  return (n + fConfig.nThreads - fThreadNo - 1) / fConfig.nThreads;
}

unsigned opdet::opDetDigitizerWorker::StartChannelToProcess(unsigned n) const
{
  unsigned n_per_job = n / fConfig.nThreads;
  unsigned leftover = std::min(fThreadNo, n % fConfig.nThreads);
  return n_per_job * fThreadNo + leftover;
}

void opdet::opDetDigitizerWorker::Start(detinfo::DetectorClocksData const& clockData) const
{
  auto arapucaDigitizer = fConfig.makeArapucaDigi(
                            *(lar::providerFrom<detinfo::LArPropertiesService>()),
                            clockData,
                            fEngine
                          );

  auto pmtDigitizer = fConfig.makePMTDigi(
                        *(lar::providerFrom<detinfo::LArPropertiesService>()),
                        clockData,
                        fEngine
                      );
  MakeWaveforms(pmtDigitizer.get(), arapucaDigitizer.get());

}

opdet::opDetDigitizerWorker::~opDetDigitizerWorker()
{
  delete fEngine;
}

void opdet::opDetDigitizerWorker::ApplyTriggerLocations(detinfo::DetectorClocksData const& clockData) const
{
  unsigned start = StartChannelToProcess(fConfig.nChannels);
  unsigned n = NChannelsToProcess(fConfig.nChannels);

  fTriggeredWaveforms->clear();

  // apply the triggers and save the output
  for (const raw::OpDetWaveform &waveform : *fWaveforms) {
    if (waveform.ChannelNumber() == std::numeric_limits<raw::Channel_t>::max() /* "NULL" value*/) {
      continue;
    }
    // only work on the prescribed channels
    if (waveform.ChannelNumber() < start || waveform.ChannelNumber() >= start + n) continue;

    std::vector<raw::OpDetWaveform> waveforms = fTriggerAlg.ApplyTriggerLocations(clockData, waveform);

    std::move(waveforms.begin(), waveforms.end(), std::back_inserter(*fTriggeredWaveforms));
  }
}

void opdet::opDetDigitizerWorker::MakeWaveforms(opdet::DigiPMTSBNDAlg *pmtDigitizer,
                                                opdet::DigiArapucaSBNDAlg *arapucaDigitizer) const
{
  if(fConfig.UseSimPhotonsLite) {
    const std::vector<art::Handle<std::vector<sim::SimPhotonsLite>>> &photon_handles = *fPhotonLiteHandles;

    // TODO: Instead of looping and evaluating if/else through all the photon_handles
    // we should get smaller containers with only the relevant handles for each case
    // ~icaza
    // std::vector<art::Handle<std::vector<sim::SimPhotonsLite>> const*> ptr_photon_handles(photon_handles.size());
    // std::transform(photon_handles.begin(), photon_handles.end(), ptr_photon_handles.begin(),
    //                [](auto& p) {return std::addressof(p);});

    // to temporarily store channel and combine PMT (direct and converted) time profiles
    std::unordered_map<int, sim::SimPhotonsLite> directPhotonsOnPMTS;
    CreateDirectPhotonMapLite(directPhotonsOnPMTS, photon_handles);
    const double startTime = fConfig.EnableWindow[0] * 1000 /*ns for digitizer*/;

    const unsigned start = StartChannelToProcess(fConfig.nChannels);
    const unsigned n = NChannelsToProcess(fConfig.nChannels);
    for (const art::Handle<std::vector<sim::SimPhotonsLite>> &opdetHandle : photon_handles) {
      // this now tells you if light collection is reflected
      const bool Reflected = (opdetHandle.provenance()->productInstanceName() == "Reflected");
      for (auto const& litesimphotons : (*opdetHandle)) {
        std::vector<short unsigned int> waveform;
        waveform.reserve(fConfig.Nsamples);
        const unsigned ch = litesimphotons.OpChannel;
        const std::string pdtype = fConfig.pdsMap.pdType(ch);
        // only work on the prescribed channels
        if (ch < start || ch >= start + n) continue;

        if((Reflected) &&
           ( (pdtype == "pmt_uncoated") || (pdtype == "pmt_coated")) ) { //All PMT channels
          pmtDigitizer->ConstructWaveformLite(ch,
                                              litesimphotons,
                                              waveform,
                                              pdtype,
                                              directPhotonsOnPMTS,
                                              startTime,
                                              fConfig.Nsamples);
          // including pre trigger window and transit time
          fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0],
                                                  (unsigned int)ch,
                                                  waveform);
        }
        // getting only xarapuca channels with appropriate type of light
        else if((pdtype == "xarapuca_vuv" && !Reflected) ||
                (pdtype == "xarapuca_vis" && Reflected) ) {
          arapucaDigitizer->ConstructWaveformLite(ch,
                                                  litesimphotons,
                                                  waveform,
                                                  pdtype,
                                                  startTime,
                                                  fConfig.Nsamples);
          // including pre trigger window and transit time
          fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0],
                                                  (unsigned int)ch,
                                                  waveform);
        }
        // getting only arapuca channels with appropriate type of light
        else if((pdtype == "arapuca_vuv" && !Reflected) ||
                (pdtype == "arapuca_vis" && Reflected) ) {
          arapucaDigitizer->ConstructWaveformLite(ch,
                                                  litesimphotons,
                                                  waveform,
                                                  pdtype,
                                                  startTime,
                                                  fConfig.Nsamples);
          // including pre trigger window and transit time
          fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0],
                                                  (unsigned int)ch,
                                                  waveform);
        }
      }
    }  //end loop on simphoton lite collections
  }
  else { // for SimPhotons
    // to temporarily store channel and direct light distribution
    std::unordered_map<int, sim::SimPhotons> directPhotonsOnPMTS;
    const std::vector<art::Handle<std::vector<sim::SimPhotons>>> &photon_handles = *fPhotonHandles;
    CreateDirectPhotonMap(directPhotonsOnPMTS, photon_handles);
    const double startTime = fConfig.EnableWindow[0] * 1000 /*ns for digitizer*/;

    const unsigned start = StartChannelToProcess(fConfig.nChannels);
    const unsigned n = NChannelsToProcess(fConfig.nChannels);
    for (const art::Handle<std::vector<sim::SimPhotons>> &opdetHandle : photon_handles) {
      const bool Reflected = (opdetHandle.provenance()->productInstanceName() == "Reflected");
      for (auto const& simphotons : (*opdetHandle)) {
        std::vector<short unsigned int> waveform;
        const unsigned ch = simphotons.OpChannel();
        const std::string pdtype = fConfig.pdsMap.pdType(ch);
        // only work on the prescribed channels
        if (ch < start || ch >= start + n) continue;
        // all PMTs
        if((Reflected) && (pdtype == "pmt_uncoated" || pdtype == "pmt_coated")) {
          pmtDigitizer->ConstructWaveform(ch,
                                          simphotons,
                                          waveform,
                                          pdtype,
                                          directPhotonsOnPMTS,
                                          startTime,
                                          fConfig.Nsamples);
          // including pre trigger window and transit time
          fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0],
                                                  (unsigned int)ch,
                                                  waveform);
        }
        // getting only arapuca channels with appropriate type of light
        if((pdtype == "arapuca_vuv" && !Reflected) ||
           (pdtype == "arapuca_vis" && Reflected)) {
          arapucaDigitizer->ConstructWaveform(ch,
                                              simphotons,
                                              waveform,
                                              pdtype,
                                              startTime,
                                              fConfig.Nsamples);
          // including pre trigger window and transit time
          fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0],
                                                  (unsigned int)ch,
                                                  waveform);
        }
        // getting only arapuca channels with appropriate type of light
        if((pdtype == "xarapuca_vuv" && !Reflected) ||
           (pdtype == "xarapuca_vis" && Reflected)) {
          arapucaDigitizer->ConstructWaveform(ch,
                                              simphotons,
                                              waveform,
                                              pdtype,
                                              startTime,
                                              fConfig.Nsamples);
          // including pre trigger window and transit time
          fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0],
                                                  (unsigned int)ch,
                                                  waveform);
        }
      }//optical channel loop
    }//type of light loop
  }//simphotons end
}


void opdet::opDetDigitizerWorker::CreateDirectPhotonMap(
  std::unordered_map<int, sim::SimPhotons>& directPhotonsOnPMTS,
  std::vector<art::Handle<std::vector<sim::SimPhotons>>> photon_handles) const
{
  // Loop over direct/reflected photons
  for (auto pmtHandle : photon_handles) {
    // Do some checking before we proceed
    if (!pmtHandle.isValid()) continue;
    if (pmtHandle.provenance()->moduleLabel() != fConfig.InputModuleName) continue;   //not the most efficient way of doing this, but preserves the logic of the module. Andrzej
    // this now tells you if light collection is reflected
    const bool Reflected = (pmtHandle.provenance()->productInstanceName() == "Reflected");
    for (auto const& simphotons : (*pmtHandle)) {
      const unsigned ch = simphotons.OpChannel();
      if(fConfig.pdsMap.isPDType(ch, "pmt_coated") && !Reflected)
        directPhotonsOnPMTS.insert(std::make_pair(ch, simphotons));
    }
  }
}


void opdet::opDetDigitizerWorker::CreateDirectPhotonMapLite(
  std::unordered_map<int, sim::SimPhotonsLite>& directPhotonsOnPMTS,
  std::vector<art::Handle<std::vector<sim::SimPhotonsLite>>> photon_handles) const
{
  // Loop over direct/reflected photons
  for (auto pmtHandle : photon_handles) {
    // Do some checking before we proceed
    if (!pmtHandle.isValid()) continue;
    if (pmtHandle.provenance()->moduleLabel() != fConfig.InputModuleName) continue;   //not the most efficient way of doing this, but preserves the logic of the module. Andrzej
    // this now tells you if light collection is reflected
    if (pmtHandle.provenance()->productInstanceName() == "Reflected") continue;
    for (auto const& litesimphotons : (*pmtHandle)) {
      const unsigned ch = litesimphotons.OpChannel;
      if(fConfig.pdsMap.isPDType(ch, "pmt_coated"))
        directPhotonsOnPMTS.insert(std::make_pair(ch, litesimphotons));
    }
  }
}

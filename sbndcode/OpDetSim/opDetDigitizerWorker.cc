// TODO: plenty of refactoring potential in here! ~icaza

#include "sbndcode/OpDetSim/opDetDigitizerWorker.h"

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
      worker.ApplyTriggerLocations();
      do_apply_trigger_locations = false;
    }
    else {
      worker.Start();
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

void opdet::opDetDigitizerWorker::Start() const
{
  auto arapucaDigitizer = fConfig.makeArapucaDigi(
                            *(lar::providerFrom<detinfo::LArPropertiesService>()),
                            *(lar::providerFrom<detinfo::DetectorClocksService>()),
                            fEngine
                          );

  auto pmtDigitizer = fConfig.makePMTDigi(
                        *(lar::providerFrom<detinfo::LArPropertiesService>()),
                        *(lar::providerFrom<detinfo::DetectorClocksService>()),
                        fEngine
                      );
  MakeWaveforms(pmtDigitizer.get(), arapucaDigitizer.get());

}

opdet::opDetDigitizerWorker::~opDetDigitizerWorker()
{
  delete fEngine;
}

void opdet::opDetDigitizerWorker::ApplyTriggerLocations() const
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

    std::vector<raw::OpDetWaveform> waveforms = fTriggerAlg.ApplyTriggerLocations(waveform);

    std::move(waveforms.begin(), waveforms.end(), std::back_inserter(*fTriggeredWaveforms));
  }
}

void opdet::opDetDigitizerWorker::MakeWaveforms(opdet::DigiPMTSBNDAlg *pmtDigitizer,
                                                opdet::DigiArapucaSBNDAlg *arapucaDigitizer) const
{
  unsigned ch, channel;
  std::string pdtype;
  if(fConfig.UseSimPhotonsLite) {
    const std::vector<art::Handle<std::vector<sim::SimPhotonsLite>>> &photon_handles = *fPhotonLiteHandles;
    // to temporarily store channel and combine PMT (direct and converted) time profiles
    std::unordered_map<int, sim::SimPhotonsLite> directPhotonsOnPMTS;
    CreateDirectPhotonMapLite(directPhotonsOnPMTS, photon_handles);

    unsigned start = StartChannelToProcess(fConfig.nChannels);
    unsigned n = NChannelsToProcess(fConfig.nChannels);
    for (const art::Handle<std::vector<sim::SimPhotonsLite>> &opdetHandle : photon_handles) {
      // this now tells you if light collection is reflected
      bool Reflected = (opdetHandle.provenance()->productInstanceName() == "Reflected");
      for (auto const& litesimphotons : (*opdetHandle)) {
        std::vector<short unsigned int> waveform;
        ch = litesimphotons.OpChannel;
        pdtype = fConfig.pdsMap.pdType(ch);
        // only work on the prescribed channels
        if (ch < start || ch >= start + n) continue;
        if((Reflected) &&
           ( (pdtype == "pmt_uncoated") || (pdtype == "pmt_coated")) ) { //All PMT channels
          pmtDigitizer->ConstructWaveformLite(ch,
                                              litesimphotons,
                                              waveform,
                                              fConfig.EnableWindow[0] * 1000 /*ns for digitizer*/,
                                              pdtype,
                                              directPhotonsOnPMTS,
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
                                                  fConfig.EnableWindow[0] * 1000 /*ns for digitizer*/,
                                                  pdtype,
                                                  fConfig.Nsamples);
          // including pre trigger window and transit time
          fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0],
                                                  (unsigned int)ch,
                                                  waveform);
        }
        // getting only xarapuca channels with appropriate type of light
        // (this separation is needed because xarapucas are set as
        // two different optical channels but are actually only one readout channel)
        else if((pdtype == "xarapuca_vuv" && !Reflected)) {
          sim::SimPhotonsLite auxLite;
          for (auto const& litesimphotons : (*opdetHandle)) {
            channel = litesimphotons.OpChannel;
            if(channel == ch) auxLite = (litesimphotons);
            if(channel == (ch + 2)) auxLite += (litesimphotons);
          }
          arapucaDigitizer->ConstructWaveformLite(ch,
                                                  auxLite,
                                                  waveform,
                                                  fConfig.EnableWindow[0] * 1000 /*ns for digitizer*/,
                                                  pdtype,
                                                  fConfig.Nsamples);
          // including pre trigger window and transit time
          fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0],
                                                  (unsigned int)ch,
                                                  waveform);
        }
        // getting only xarapuca channels with appropriate type of light
        // (this separation is needed because xarapucas are set as
        // two different optical channels but are actually only one readout channel)
        else if(pdtype == "xarapuca_vis" && Reflected) {
          sim::SimPhotonsLite auxLite;
          for (auto const& litesimphotons : (*opdetHandle)) {
            channel = litesimphotons.OpChannel;
            if(channel == ch) auxLite = (litesimphotons);
            if(channel == (ch + 2)) auxLite += (litesimphotons);
          }
          arapucaDigitizer->ConstructWaveformLite(ch,
                                                  auxLite,
                                                  waveform,
                                                  fConfig.EnableWindow[0] * 1000 /*ns for digitizer*/,
                                                  pdtype,
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

    unsigned start = StartChannelToProcess(fConfig.nChannels);
    unsigned n = NChannelsToProcess(fConfig.nChannels);
    for (const art::Handle<std::vector<sim::SimPhotons>> &opdetHandle : photon_handles) {
      bool Reflected = (opdetHandle.provenance()->productInstanceName() == "Reflected");
      for (auto const& simphotons : (*opdetHandle)) {
        std::vector<short unsigned int> waveform;
        ch = simphotons.OpChannel();
        pdtype = fConfig.pdsMap.pdType(ch);
        // only work on the prescribed channels
        if (ch < start || ch >= start + n) continue;
        // all PMTs
        if((Reflected) && (pdtype == "pmt_uncoated" || pdtype == "pmt_coated")) {
          pmtDigitizer->ConstructWaveform(ch,
                                          simphotons,
                                          waveform,
                                          pdtype,
                                          fConfig.EnableWindow[0] * 1000 /*ns for digitizer*/,
                                          directPhotonsOnPMTS,
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
                                              fConfig.EnableWindow[0] * 1000 /*ns for digitizer*/,
                                              fConfig.Nsamples);
          // including pre trigger window and transit time
          fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0],
                                                  (unsigned int)ch,
                                                  waveform);
        }
        // getting only xarapuca channels with appropriate type of light
        // (this separation is needed because xarapucas are set as
        // two different optical channels but are actually only one readout channel)
        if(pdtype == "xarapuca_vuv" && !Reflected) {
          sim::SimPhotons auxPhotons;
          for (auto const& simphotons : (*opdetHandle)) {
            channel = simphotons.OpChannel();
            if(channel == ch) auxPhotons = (simphotons);
            if(channel == (ch + 2)) auxPhotons += (simphotons);
          }
          arapucaDigitizer->ConstructWaveform(ch,
                                              auxPhotons,
                                              waveform,
                                              pdtype,
                                              fConfig.EnableWindow[0] * 1000 /*ns for digitizer*/,
                                              fConfig.Nsamples);
          // including pre trigger window and transit time
          fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0],
                                                  (unsigned int)ch,
                                                  waveform);
        }
        // getting only xarapuca channels with appropriate type of light
        // (this separation is needed because xarapucas are set as
        // two different optical channels but are actually only one readout channel)
        if(pdtype == "xarapuca_vis" && Reflected) {
          sim::SimPhotons auxPhotons;
          for (auto const& simphotons : (*opdetHandle)) {
            channel = simphotons.OpChannel();
            if(channel == ch) auxPhotons = (simphotons);
            if(channel == (ch + 2)) auxPhotons += (simphotons);
          }
          arapucaDigitizer->ConstructWaveform(ch,
                                              auxPhotons,
                                              waveform,
                                              pdtype,
                                              fConfig.EnableWindow[0] * 1000 /*ns for digitizer*/,
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
  int ch;
  // Loop over direct/reflected photons
  for (auto pmtHandle : photon_handles) {
    // Do some checking before we proceed
    if (!pmtHandle.isValid()) continue;
    if (pmtHandle.provenance()->moduleLabel() != fConfig.InputModuleName) continue;   //not the most efficient way of doing this, but preserves the logic of the module. Andrzej
    // this now tells you if light collection is reflected
    bool Reflected = (pmtHandle.provenance()->productInstanceName() == "Reflected");
    for (auto const& simphotons : (*pmtHandle)) {
      ch = simphotons.OpChannel();
      if(fConfig.pdsMap.isPDType(ch, "pmt_coated") && !Reflected)
        directPhotonsOnPMTS.insert(std::make_pair(ch, simphotons));
    }
  }
}


void opdet::opDetDigitizerWorker::CreateDirectPhotonMapLite(
  std::unordered_map<int, sim::SimPhotonsLite>& directPhotonsOnPMTS,
  std::vector<art::Handle<std::vector<sim::SimPhotonsLite>>> photon_handles) const
{
  int ch;
  // Loop over direct/reflected photons
  for (auto pmtHandle : photon_handles) {
    // Do some checking before we proceed
    if (!pmtHandle.isValid()) continue;
    if (pmtHandle.provenance()->moduleLabel() != fConfig.InputModuleName) continue;   //not the most efficient way of doing this, but preserves the logic of the module. Andrzej
    // this now tells you if light collection is reflected
    bool Reflected = (pmtHandle.provenance()->productInstanceName() == "Reflected");
    for (auto const& litesimphotons : (*pmtHandle)) {
      ch = litesimphotons.OpChannel;
      if(fConfig.pdsMap.isPDType(ch, "pmt_coated"))
        directPhotonsOnPMTS.insert(std::make_pair(ch, litesimphotons));
    }
  }
}

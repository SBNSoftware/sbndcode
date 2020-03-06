#include <map>
// TODO: plenty of refactoring potential in here! ~icaza
// TODO: map.isPDType and map.pdType are called a bunch of times, maybe it could be saved into a variable?

#include "sbndcode/OpDetSim/opDetDigitizerWorker.h"

opdet::opDetDigitizerWorker::Config::Config(const opdet::DigiPMTSBNDAlgMaker::Config &pmt_config, const opdet::DigiArapucaSBNDAlgMaker::Config &arapuca_config):
  makePMTDigi(pmt_config),
  makeArapucaDigi(arapuca_config) {}

opdet::opDetDigitizerWorker::opDetDigitizerWorker(unsigned no, const Config &config, CLHEP::HepRandomEngine *Engine, const opDetSBNDTriggerAlg &trigger_alg):
  fConfig(config),
  fThreadNo(no),
  fEngine(Engine),
  fTriggerAlg(trigger_alg)
{}

void opdet::opDetDigitizerWorkerThread(const opdet::opDetDigitizerWorker &worker,
                                       opdet::opDetDigitizerWorker::Semaphore &sem_start, opdet::opDetDigitizerWorker::Semaphore &sem_finish,
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

void opdet::StartopDetDigitizerWorkers(unsigned n_workers, opdet::opDetDigitizerWorker::Semaphore &sem_start)
{
  sem_start.increment(n_workers);
}

void opdet::WaitopDetDigitizerWorkers(unsigned n_workers, opdet::opDetDigitizerWorker::Semaphore &sem_finish)
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
  unsigned start = StartChannelToProcess(fConfig.map.size());
  unsigned n = NChannelsToProcess(fConfig.map.size());

  fTriggeredWaveforms->clear();

  // apply the triggers and save the output
  for (const raw::OpDetWaveform &waveform : *fWaveforms) {
    if (waveform.ChannelNumber() == std::numeric_limits<raw::Channel_t>::max() /* "NULL" value*/) {
      continue;
    }
    // only work on the perscribed channels
    if (waveform.ChannelNumber() < start || waveform.ChannelNumber() >= start + n) continue;

    std::vector<raw::OpDetWaveform> waveforms = fTriggerAlg.ApplyTriggerLocations(waveform);

    std::move(waveforms.begin(), waveforms.end(), std::back_inserter(*fTriggeredWaveforms));
  }
}

void opdet::opDetDigitizerWorker::MakeWaveforms(opdet::DigiPMTSBNDAlg *pmtDigitizer, opdet::DigiArapucaSBNDAlg *arapucaDigitizer) const
{
  unsigned ch, channel;
  if(fConfig.UseLitePhotons == 1) { //using SimPhotonsLite
    const std::vector<art::Handle<std::vector<sim::SimPhotonsLite>>> &photon_handles = *fPhotonLiteHandles;

    std::map<int, sim::SimPhotonsLite> auxmap;  // to temporarily store channel and combine PMT (direct and converted) time profiles
    CreateDirectPhotonMapLite(auxmap, photon_handles);

    unsigned start = StartChannelToProcess(fConfig.map.size());
    unsigned n = NChannelsToProcess(fConfig.map.size());
    for (const art::Handle<std::vector<sim::SimPhotonsLite>> &opdetHandle : photon_handles) {
      //this now tells you if light collection is reflected
      bool Reflected = (opdetHandle.provenance()->productInstanceName() == "Reflected");

      for (auto const& litesimphotons : (*opdetHandle)) {
        std::vector<short unsigned int> waveform;
        ch = litesimphotons.OpChannel;

        // only work on the perscribed channels
        if (ch < start || ch >= start + n) continue;

        if((Reflected) && (fConfig.map.isPDType(ch, "uncoatedpmt") || fConfig.map.isPDType(ch, "coatedpmt") )) { //All PMT channels
          pmtDigitizer->ConstructWaveformLite(ch, litesimphotons, waveform, fConfig.map.pdType(ch), auxmap, fConfig.EnableWindow[0] * 1000 /*ns for digitizer*/, fConfig.Nsamples);
          fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0], (unsigned int)ch, waveform);//including pre trigger window and transit time
        }
        else if((fConfig.map.isPDType(ch, "arapucaT1") && !Reflected) || (fConfig.map.isPDType(ch, "arapucaT2") && Reflected) ) { //getting only arapuca channels with appropriate type of light
          arapucaDigitizer->ConstructWaveformLite(ch, litesimphotons, waveform, fConfig.map.pdType(ch), fConfig.EnableWindow[0] * 1000 /*ns for digitizer*/, fConfig.Nsamples);
          fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0], (unsigned int)ch, waveform);//including pre trigger window and transit time
        }
        else if((fConfig.map.isPDType(ch, "xarapucaT1") && !Reflected)) { //getting only xarapuca channels with appropriate type of light (this separation is needed because xarapucas are set as two different optical channels but are actually only one readout channel)
          sim::SimPhotonsLite auxLite;
          for (auto const& litesimphotons : (*opdetHandle)) {
            channel = litesimphotons.OpChannel;
            if(channel == ch) auxLite = (litesimphotons);
            if(channel == (ch + 2)) auxLite += (litesimphotons);
          }
          arapucaDigitizer->ConstructWaveformLite(ch, auxLite, waveform, fConfig.map.pdType(ch), fConfig.EnableWindow[0] * 1000 /*ns for digitizer*/, fConfig.Nsamples);
          fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0], (unsigned int)ch, waveform);//including pre trigger window and transit time
        }
        else if((fConfig.map.isPDType(ch, "xarapucaT2") && Reflected)) { //getting only xarapuca channels with appropriate type of light (this separation is needed because xarapucas are set as two different optical channels but are actually only one readout channel)
          sim::SimPhotonsLite auxLite;
          for (auto const& litesimphotons : (*opdetHandle)) {
            channel = litesimphotons.OpChannel;
            if(channel == ch) auxLite = (litesimphotons);
            if(channel == (ch + 2)) auxLite += (litesimphotons);
          }
          arapucaDigitizer->ConstructWaveformLite(ch, auxLite, waveform, fConfig.map.pdType(ch), fConfig.EnableWindow[0] * 1000 /*ns for digitizer*/, fConfig.Nsamples);
          fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0], (unsigned int)ch, waveform);//including pre trigger window and transit time
        }
      }
    }  //end loop on simphoton lite collections
  }
  else { //for SimPhotons
    std::map<int, sim::SimPhotons> auxmap;  // to temporarily store channel and direct light distribution

    const std::vector<art::Handle<std::vector<sim::SimPhotons>>> &photon_handles = *fPhotonHandles;
    CreateDirectPhotonMap(auxmap, photon_handles);

    unsigned start = StartChannelToProcess(fConfig.map.size());
    unsigned n = NChannelsToProcess(fConfig.map.size());
    for (const art::Handle<std::vector<sim::SimPhotons>> &opdetHandle : photon_handles) {
      bool Reflected = (opdetHandle.provenance()->productInstanceName() == "Reflected");

      for (auto const& simphotons : (*opdetHandle)) {
        std::vector<short unsigned int> waveform;
        ch = simphotons.OpChannel();

        // only work on the perscribed channels
        if (ch < start || ch >= start + n) continue;

        if((Reflected) && (fConfig.map.isPDType(ch, "uncoatedpmt") || fConfig.map.isPDType(ch, "coatedpmt"))) { //all PMTs
          pmtDigitizer->ConstructWaveform(ch, simphotons, waveform, fConfig.map.pdType(ch), auxmap, fConfig.EnableWindow[0] * 1000 /*ns for digitizer*/, fConfig.Nsamples);
          fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0], (unsigned int)ch, waveform);//including pre trigger window and transit time
        }
        if((fConfig.map.isPDType(ch, "arapucaT1") && !Reflected) || (fConfig.map.isPDType(ch, "arapucaT2") && Reflected) ) { //getting only arapuca channels with appropriate type of light
          arapucaDigitizer->ConstructWaveform(ch, simphotons, waveform, fConfig.map.pdType(ch), fConfig.EnableWindow[0] * 1000 /*ns for digitizer*/, fConfig.Nsamples);
          fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0], (unsigned int)ch, waveform);//including pre trigger window and transit time
        }
        if((fConfig.map.isPDType(ch, "xarapucaT1") && !Reflected)) { //getting only xarapuca channels with appropriate type of light (this separation is needed because xarapucas are set as two different optical channels but are actually only one readout channel)
          sim::SimPhotons auxPhotons;
          for (auto const& simphotons : (*opdetHandle)) {
            channel = simphotons.OpChannel();
            if(channel == ch) auxPhotons = (simphotons);
            if(channel == (ch + 2)) auxPhotons += (simphotons);
          }
          arapucaDigitizer->ConstructWaveform(ch, auxPhotons, waveform, fConfig.map.pdType(ch), fConfig.EnableWindow[0] * 1000 /*ns for digitizer*/, fConfig.Nsamples);
          fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0], (unsigned int)ch, waveform);//including pre trigger window and transit time
        }
        if((fConfig.map.isPDType(ch, "xarapucaT2") && !Reflected)) { //getting only xarapuca channels with appropriate type of light (this separation is needed because xarapucas are set as two different optical channels but are actually only one readout channel)
          sim::SimPhotons auxPhotons;
          for (auto const& simphotons : (*opdetHandle)) {
            channel = simphotons.OpChannel();
            if(channel == ch) auxPhotons = (simphotons);
            if(channel == (ch + 2)) auxPhotons += (simphotons);
          }
          arapucaDigitizer->ConstructWaveform(ch, auxPhotons, waveform, fConfig.map.pdType(ch), fConfig.EnableWindow[0] * 1000 /*ns for digitizer*/, fConfig.Nsamples);
          fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0], (unsigned int)ch, waveform);//including pre trigger window and transit time
        }
      }//optical channel loop
    }//type of light loop
  }//simphotons end
}
void opdet::opDetDigitizerWorker::CreateDirectPhotonMapLite(std::map<int, sim::SimPhotonsLite>& auxmap, std::vector< art::Handle< std::vector< sim::SimPhotonsLite > > > photon_handles) const
{
  int ch;
  // Loop over direct/reflected photons
  for (auto pmtHandle : photon_handles) {
    // Do some checking before we proceed
    if (!pmtHandle.isValid()) continue;
    if (pmtHandle.provenance()->moduleLabel() != fConfig.InputModuleName) continue;   //not the most efficient way of doing this, but preserves the logic of the module. Andrzej
    //this now tells you if light collection is reflected
    bool Reflected = (pmtHandle.provenance()->productInstanceName() == "Reflected");

    for (auto const& litesimphotons : (*pmtHandle)) {
      ch = litesimphotons.OpChannel;
      if(fConfig.map.isPDType(ch, "coatedpmt") && !Reflected)
        auxmap.insert(std::make_pair(ch, litesimphotons));
    }
  }
}

void opdet::opDetDigitizerWorker::CreateDirectPhotonMap(std::map<int, sim::SimPhotons>& auxmap, std::vector< art::Handle< std::vector< sim::SimPhotons > > > photon_handles) const
{
  int ch;
  // Loop over direct/reflected photons
  for (auto pmtHandle : photon_handles) {
    // Do some checking before we proceed
    if (!pmtHandle.isValid()) continue;
    if (pmtHandle.provenance()->moduleLabel() != fConfig.InputModuleName) continue;   //not the most efficient way of doing this, but preserves the logic of the module. Andrzej
    //this now tells you if light collection is reflected
    bool Reflected = (pmtHandle.provenance()->productInstanceName() == "Reflected");
    for (auto const& simphotons : (*pmtHandle)) {
      ch = simphotons.OpChannel();
      if(fConfig.map.isPDType(ch, "coatedpmt") && !Reflected)
        auxmap.insert(std::make_pair(ch, simphotons));
    }
  }
}

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
  for (const raw::OpDetWaveform &waveform : *fWaveforms){
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


    //need to combine direct and reflected photons
    std::unordered_map<std::string, std::unordered_map<int, sim::SimPhotonsLite> > AllPhotonsMap;
    AllPhotonsMap["Direct"] = std::unordered_map<int, sim::SimPhotonsLite>();
    AllPhotonsMap["Reflected"] = std::unordered_map<int, sim::SimPhotonsLite>();
    for (const art::Handle<std::vector<sim::SimPhotonsLite>> &opdetHandle : photon_handles) {
      const bool Reflected = (opdetHandle.provenance()->productInstanceName() == "Reflected");
      for (auto const& litesimphotons : (*opdetHandle)){

        if(Reflected){
          auto it = AllPhotonsMap["Reflected"].find(litesimphotons.OpChannel);
          if(it==AllPhotonsMap["Reflected"].end())
            AllPhotonsMap["Reflected"][litesimphotons.OpChannel] = litesimphotons;
          else
            AllPhotonsMap["Reflected"][litesimphotons.OpChannel] += litesimphotons;
        }

        else{ //Direct
          auto it = AllPhotonsMap["Direct"].find(litesimphotons.OpChannel);
          if(it==AllPhotonsMap["Direct"].end())
            AllPhotonsMap["Direct"][litesimphotons.OpChannel] = litesimphotons;
          else
            AllPhotonsMap["Direct"][litesimphotons.OpChannel] += litesimphotons;
        }
      }
    }

    const double startTime = fConfig.EnableWindow[0] * 1000. /*ns for digitizer*/;

    std::unordered_map<int, sim::SimPhotonsLite> DirectPhotonsMap;
    std::unordered_map<int, sim::SimPhotonsLite> ReflectedPhotonsMap;
    std::unordered_set<short unsigned int> coatedpmts_todigitize;
    // XARAPUCAs
    std::unordered_map<int, sim::SimPhotonsLite> XADirectPhotonsMap;
    std::unordered_map<int, sim::SimPhotonsLite> XAReflectedPhotonsMap;
    std::unordered_set<short unsigned int>   vuvxarapucas_todigitize; //sensible to both VUV and visible light


    const unsigned start = StartChannelToProcess(fConfig.nChannels);
    const unsigned n = NChannelsToProcess(fConfig.nChannels);
    for (auto const& simphotons_cols : AllPhotonsMap){
      // this now tells you if light collection is reflected
      const bool Reflected = (simphotons_cols.first=="Reflected");
      for (auto const& simphotons_map : simphotons_cols.second) {

        auto const& litesimphotons = simphotons_map.second;

        std::vector<short unsigned int> waveform;
        waveform.reserve(fConfig.Nsamples);
        const unsigned ch = litesimphotons.OpChannel;
        const std::string pdtype = fConfig.pdsMap.pdType(ch);

        // only work on the prescribed channels
        if (ch < start || ch >= start + n) continue;


        if( pdtype == "pmt_coated" ){
          if(Reflected)
            ReflectedPhotonsMap.insert(std::make_pair(ch, litesimphotons));
          else
            DirectPhotonsMap.insert(std::make_pair(ch, litesimphotons));

          coatedpmts_todigitize.insert(ch);
        }
        else if( pdtype == "xarapuca_vuv" ){
          if(Reflected)
            XAReflectedPhotonsMap.insert(std::make_pair(ch, litesimphotons));
          else
            XADirectPhotonsMap.insert(std::make_pair(ch, litesimphotons));

          vuvxarapucas_todigitize.insert(ch);
        }
        else if( (Reflected) && (pdtype == "pmt_uncoated") ) { //Uncoated PMT channels
          pmtDigitizer->ConstructWaveformLiteUncoatedPMT(ch,
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
        // getting only xarapuca channels with appropriate type of light
        else if( (pdtype == "xarapuca_vis" && Reflected) ) {
          const bool is_daphne= fConfig.pdsMap.isElectronics(ch,"daphne");
          if(is_daphne){
            arapucaDigitizer->ConstructWaveformLite(ch,
                                                  litesimphotons,
                                                  waveform,
                                                  pdtype,
                                                  is_daphne,
                                                  startTime,
                                                  fConfig.Nsamples_Daphne);
            // including pre trigger window and transit time
            fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0],
                                                  (unsigned int)ch,
                                                  waveform);
            }
            else{
            arapucaDigitizer->ConstructWaveformLite(ch,
                                                  litesimphotons,
                                                  waveform,
                                                  pdtype,
                                                  is_daphne,
                                                  startTime,
                                                  fConfig.Nsamples);
            // including pre trigger window and transit time
            fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0],
                                                  (unsigned int)ch,
                                                  waveform);
          }
        }
      }
    }  //end loop on simphoton lite collections

    //Constructing Waveforms for hybrid OpChannels (coated pmts)
    for(auto ch : coatedpmts_todigitize){
      std::vector<short unsigned int> waveform;
      waveform.reserve(fConfig.Nsamples);
      pmtDigitizer->ConstructWaveformLiteCoatedPMT(ch, waveform, DirectPhotonsMap, ReflectedPhotonsMap, startTime, fConfig.Nsamples);
      fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0],
                                              (unsigned int)ch,
                                              waveform);
    }
    //VUV XAs, sensible to VUV and visible light
    for(auto ch : vuvxarapucas_todigitize){
      std::vector<short unsigned int> waveform;
      waveform.reserve(fConfig.Nsamples_Daphne);
      arapucaDigitizer->ConstructWaveformLiteVUVXA(ch, waveform, XADirectPhotonsMap, XAReflectedPhotonsMap, startTime, fConfig.Nsamples_Daphne);
      fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0],
                                              (unsigned int)ch,
                                              waveform);
    }
  }
  else { // for SimPhotons
    const std::vector<art::Handle<std::vector<sim::SimPhotons>>> &photon_handles = *fPhotonHandles;

    //need to combine direct and reflected photons
    std::unordered_map<std::string, std::unordered_map<int, sim::SimPhotons> > AllPhotonsMap;
    AllPhotonsMap["Direct"] = std::unordered_map<int, sim::SimPhotons>();
    AllPhotonsMap["Reflected"] = std::unordered_map<int, sim::SimPhotons>();
    for (const art::Handle<std::vector<sim::SimPhotons>> &opdetHandle : photon_handles) {
      const bool Reflected = (opdetHandle.provenance()->productInstanceName() == "Reflected");
      for (auto const& simphotons : (*opdetHandle)){

        if(Reflected){
          auto it = AllPhotonsMap["Reflected"].find(simphotons.OpChannel());
          if(it==AllPhotonsMap["Reflected"].end())
            AllPhotonsMap["Reflected"][simphotons.OpChannel()] = simphotons;
          else
            AllPhotonsMap["Reflected"][simphotons.OpChannel()] += simphotons;
        }

        else{ //Direct
          auto it = AllPhotonsMap["Direct"].find(simphotons.OpChannel());
          if(it==AllPhotonsMap["Direct"].end())
            AllPhotonsMap["Direct"][simphotons.OpChannel()] = simphotons;
          else
            AllPhotonsMap["Direct"][simphotons.OpChannel()] += simphotons;
        }
      }
    }


    // to temporarily store channel and direct light distribution
    std::unordered_map<int, sim::SimPhotons> DirectPhotonsMap;
    std::unordered_map<int, sim::SimPhotons> ReflectedPhotonsMap;
    std::unordered_set<short unsigned int> coatedpmts_todigitize;
    // XARAPUCAs
    std::unordered_map<int, sim::SimPhotons> XADirectPhotonsMap;
    std::unordered_map<int, sim::SimPhotons> XAReflectedPhotonsMap;
    std::unordered_set<short unsigned int>   vuvxarapucas_todigitize; //sensible to both VUV and visible light

    const double startTime = fConfig.EnableWindow[0] * 1000. /*ns for digitizer*/;

    const unsigned start = StartChannelToProcess(fConfig.nChannels);
    const unsigned n = NChannelsToProcess(fConfig.nChannels);
    for (auto const& simphotons_cols : AllPhotonsMap){
      // this now tells you if light collection is reflected
      const bool Reflected = (simphotons_cols.first=="Reflected");
      for (auto const& simphotons_map : simphotons_cols.second) {

        auto const& simphotons = simphotons_map.second;

        std::vector<short unsigned int> waveform;
        const unsigned ch = simphotons.OpChannel();
        const std::string pdtype = fConfig.pdsMap.pdType(ch);
        const bool is_daphne = fConfig.pdsMap.isElectronics(ch,"daphne");
        // only work on the prescribed channels
        if (ch < start || ch >= start + n) continue;
        //coated PMTs
        if( pdtype == "pmt_coated" ){
          if(Reflected)
            ReflectedPhotonsMap.insert(std::make_pair(ch, simphotons));
          else
            DirectPhotonsMap.insert(std::make_pair(ch, simphotons));

          coatedpmts_todigitize.insert(ch);
        }
        // uncoated PMTs
        else if(Reflected && pdtype == "pmt_uncoated") {
          pmtDigitizer->ConstructWaveformUncoatedPMT(ch,
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
        if( pdtype == "xarapuca_vuv" ){
          if(Reflected)
            XAReflectedPhotonsMap.insert(std::make_pair(ch, simphotons));
          else
            XADirectPhotonsMap.insert(std::make_pair(ch, simphotons));

          vuvxarapucas_todigitize.insert(ch);
        }
        // getting only xarapuca channels with appropriate type of light
        if((pdtype == "xarapuca_vis" && Reflected)) {
          if(is_daphne){
          arapucaDigitizer->ConstructWaveform(ch,
                                              simphotons,
                                              waveform,
                                              pdtype,
                                              is_daphne,
                                              startTime,
                                              fConfig.Nsamples_Daphne);
            // including pre trigger window and transit time
            fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0],
                                                    (unsigned int)ch,
                                                  waveform);
            }
            else{
              arapucaDigitizer->ConstructWaveform(ch,
                                              simphotons,
                                              waveform,
                                              pdtype,
                                              is_daphne,
                                              startTime,
                                              fConfig.Nsamples);
              // including pre trigger window and transit time
              fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0],
                                                  (unsigned int)ch,
                                                  waveform);
          }
        }
      }//optical channel loop
    }//type of light loop
    //Constructing Waveforms for hybrid OpChannels (coated pmts and VUV XAs)
    for(auto ch : coatedpmts_todigitize){
      std::vector<short unsigned int> waveform;
      waveform.reserve(fConfig.Nsamples);
      pmtDigitizer->ConstructWaveformCoatedPMT(ch, waveform, DirectPhotonsMap, ReflectedPhotonsMap, startTime, fConfig.Nsamples);
      fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0],
                                              (unsigned int)ch,
                                              waveform);
    }
    for(auto ch : vuvxarapucas_todigitize){
      std::vector<short unsigned int> waveform;
      waveform.reserve(fConfig.Nsamples_Daphne);
      arapucaDigitizer->ConstructWaveformVUVXA(ch, waveform, XADirectPhotonsMap, XAReflectedPhotonsMap, startTime, fConfig.Nsamples_Daphne);
      fWaveforms->at(ch) = raw::OpDetWaveform(fConfig.EnableWindow[0],
                                              (unsigned int)ch,
                                              waveform);
    }
  }//simphotons end
}

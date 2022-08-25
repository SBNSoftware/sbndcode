#ifndef SBND_CRTDETSIMALG_CC
#define SBND_CRTDETSIMALG_CC

#include "sbndcode/CRT/CRTSimulation/CRTDetSimAlg.h"

namespace sbnd {
namespace crt {

    CRTDetSimAlg::CRTDetSimAlg(const Parameters & params, CLHEP::HepRandomEngine& engine, double g4RefTime)
    : fParams(params())
    , fEngine(engine)
    , fG4RefTime(g4RefTime)
    {
        ConfigureWaveform();
        ConfigureTimeOffset();

        fTaggers.clear();
        fData.clear();
        fAuxData.clear();
    }


    void CRTDetSimAlg::ConfigureWaveform() {

        std::vector<double> wvf_x = fParams.WaveformX();
        std::vector<double> wvf_y = fParams.WaveformY();

        // Normalize the waveform
        float max_y = *std::max_element(std::begin(wvf_y), std::end(wvf_y));
        for (auto & y : wvf_y) y /= max_y;

        // Reverse the waveform, as we are only interested to use it to
        // estimate its effect on time delays.
        std::reverse(wvf_y.begin(),wvf_y.end());

        fInterpolator = std::make_unique<ROOT::Math::Interpolator>
            (wvf_y.size(), ROOT::Math::Interpolation::kLINEAR);

        fInterpolator->SetData(wvf_x, wvf_y);

    }

    void CRTDetSimAlg::ConfigureTimeOffset()
    {
        if (fParams.UseG4RefTimeOffset())
        {
            fTimeOffset = fG4RefTime;
            mf::LogInfo("CRTDetSimAlg") << "Configured to use G4 ref time as time offset: "
                                        << fTimeOffset << " ns." << std::endl;
        }
        else
        {
            fTimeOffset = fParams.GlobalT0Offset();
            mf::LogInfo("CRTDetSimAlg") << "Configured to use GlobalT0Offset as time offset: "
                                        << fTimeOffset << " ns." << std::endl;
        }
    }

    std::vector<std::pair<sbnd::crt::FEBData, std::vector<AuxDetIDE>>> CRTDetSimAlg::GetData()
    {
        return fData;
    }

    std::vector<std::vector<int>> CRTDetSimAlg::GetAuxData()
    {
        return fAuxData;
    }



    uint16_t CRTDetSimAlg::WaveformEmulation(const uint32_t & time_delay, const double & adc)
    {

        if (!fParams.DoWaveformEmulation())
        {
            return static_cast<uint16_t>(adc);
        }

        if (time_delay < 0)
        {
            throw art::Exception(art::errors::LogicError)
                << "Time delay cannot be negative for waveform emulation to happen." << std::endl;
        }

        if (time_delay > fParams.WaveformX().back())
        {
            // If the time delay is more than the waveform rise time, we
            // will never be able to see this signal. So return a 0 ADC value.
            return 0;
        }

        // Evaluate the waveform
        double wf = fInterpolator->Eval(time_delay);
        wf *= adc;

        if (fParams.DebugTrigger()) std::cout << "WaveformEmulation, time_delay " << time_delay
                                              << ", adc " << adc
                                              << ", waveform " << wf << std::endl;

        return static_cast<uint16_t>(wf);
    }


    void CRTDetSimAlg::AddADC(sbnd::crt::FEBData & feb_data,
                              const int & sipmID, const uint16_t & adc)
    {
        uint16_t original_adc = feb_data.ADC(sipmID);
        uint16_t new_adc = original_adc + adc;

        if (new_adc > fParams.AdcSaturation())
        {
            new_adc = fParams.AdcSaturation();
        }

        feb_data.SetADC(sipmID, new_adc);

        if (fParams.DebugTrigger()) std::cout << "Updating ADC value for FEB " << feb_data.Mac5()
                                              << ", sipmID " << sipmID
                                              << " with adc " << adc
                                              << ": was " << original_adc
                                              << ", now is " << feb_data.ADC(sipmID) << std::endl;

    }




    void CRTDetSimAlg::ProcessStrips(const std::vector<StripData> & strips)
    {
        std::map<uint16_t, sbnd::crt::FEBData> mac_to_febdata;
        std::map<uint16_t, std::vector<AuxDetIDE>> mac_to_ides;
        std::map<uint16_t, std::vector<int>> mac_to_sipmids;

        // TODO Add pedestal fluctuations
        std::array<uint16_t, 32> adc_pedestal = {static_cast<uint16_t>(fParams.QPed())};

        for (auto & strip : strips)
        {

            // First time we encounter this FEB...
            if (!mac_to_febdata.count(strip.mac5))
            {
                // Construct a new FEBData object with only pedestal values (will be filled later)
                mac_to_febdata[strip.mac5] = sbnd::crt::FEBData(strip.mac5,          // FEB ID
                                                                strip.sipm0.t0,      // Ts0
                                                                strip.sipm0.t1,      // Ts1
                                                                adc_pedestal,        // ADCs
                                                                strip.sipm0.sipmID); // Coinc

                mac_to_ides[strip.mac5] = std::vector<AuxDetIDE>();
                mac_to_sipmids[strip.mac5] = std::vector<int>();
            }
            // ... all the other times we encounter this FEB
            else
            {
                // We want to save the earliest t1 and t0 for each FEB.
                if (strip.sipm0.t1 < mac_to_febdata[strip.mac5].Ts1())
                {
                    mac_to_febdata[strip.mac5].SetTs1(strip.sipm0.t1);
                    mac_to_febdata[strip.mac5].SetTs0(strip.sipm0.t0);
                    mac_to_febdata[strip.mac5].SetCoinc(strip.sipm0.sipmID);
                }
            }
        }


        for (auto & strip : strips)
        {

            auto &feb_data = mac_to_febdata[strip.mac5];
            uint32_t trigger_time = feb_data.Ts1();

            uint16_t adc_sipm0 = WaveformEmulation(strip.sipm0.t1 - trigger_time, strip.sipm0.adc);
            uint16_t adc_sipm1 = WaveformEmulation(strip.sipm1.t1 - trigger_time, strip.sipm1.adc);

            AddADC(feb_data, strip.sipm0.sipmID, adc_sipm0);
            AddADC(feb_data, strip.sipm1.sipmID, adc_sipm1);

            auto &ides = mac_to_ides[strip.mac5];
            ides.push_back(strip.ide);

            auto &sipmids = mac_to_sipmids[strip.mac5];
            sipmids.push_back(std::min(strip.sipm0.sipmID, strip.sipm1.sipmID));
        }

        for (auto const& [mac, feb_data] : mac_to_febdata)
        {
            fData.push_back(std::make_pair(feb_data, mac_to_ides[mac]));
            fAuxData.push_back(mac_to_sipmids[mac]);
        }

        if (fParams.DebugTrigger()) std::cout << "Constructed " << mac_to_febdata.size()
                                              << " FEBData object(s)." << std::endl << std::endl;
    }



    void CRTDetSimAlg::CreateData()
    {

        /** A struct to temporarily store information on a CRT Tagger trigger.
         */
        struct Trigger {
            bool _is_bottom;
            double _dead_time;
            bool _planeX;
            bool _planeY;
            std::set<int> _mac5s;
            std::vector<StripData> _strips;
            uint32_t _trigger_time;
            bool _debug;

            Trigger(bool is_bottom, double dead_time, bool debug) {
                _planeX = _planeY = false;
                _is_bottom = is_bottom;
                _dead_time = dead_time;
                _debug = debug;
            }

            /** \brief Resets this trigger object */
            void reset(uint32_t trigger_time) {
                _planeX = _planeY = false;
                _strips.clear();
                _mac5s.clear();
                _trigger_time = trigger_time;

                if(_debug) std::cout << "TRIGGER TIME IS " << _trigger_time << std::endl;
            }

            /** \brief Add a strip belonging to a particular trigger */
            void add_strip(StripData strip) {
                _strips.push_back(strip);
                _mac5s.insert(strip.mac5);

                if (strip.sipm_coinc) {
                    if (strip.planeID == 0) _planeX = true;
                    if (strip.planeID == 1) _planeY = true;
                }

                if(_debug) std::cout << "\tAdded strip with mac " << strip.mac5
                                     << " on plane " << strip.planeID
                                     << ", with time " << strip.sipm0.t1 << std::endl;
            }

            /** \brief Tells is a tagger is triggering or not */
            bool tagger_triggered() {
                if (_is_bottom) {
                    return _planeX or _planeY;
                }
                return _planeX and _planeY;
            }

            /** \brief Returns true is the strip is in dead time */
            bool is_in_dead_time(int mac5, int time) {
                if (!_mac5s.count(mac5)) { return false; }
                return (time <= _dead_time);
            }

            void print_no_coinc(StripData strip) {
                if (_debug) std::cout << "\tStrip with mac " << strip.mac5
                                     << " on plane " << strip.planeID
                                     << ", with time " << strip.sipm0.t1
                                     << " -> didn't have SiPMs coincidence" << std::endl;
            }

            void print_dead_time(StripData strip) {
                if (_debug) std::cout << "\tStrip with mac " << strip.mac5
                                     << " on plane " << strip.planeID
                                     << ", with time " << strip.sipm0.t1
                                     << " -> happened during dead time." << std::endl;
            }
        };


        // Loop over all the CRT Taggers and simulate triggering, dead time, ...
        for (auto iter : fTaggers)
        {
            auto & name = iter.first;
            auto & tagger = iter.second;

           mf::LogInfo("CRTDetSimAlg") << "Simulating trigger for tagger " << name << std::endl;

            bool is_bottom = name.find("Bottom") != std::string::npos;
            Trigger trigger(is_bottom, fParams.DeadTime(), fParams.DebugTrigger());

            auto & strip_data_v = tagger.data;

            // Time order the data
            std::sort(strip_data_v.begin(), strip_data_v.end(),
                      [](const StripData& strip1,
                         const StripData& strip2) {
                         return strip1.sipm0.t1 < strip2.sipm0.t1;
                         });

            uint32_t trigger_ts1 = 0, current_time = 0;
            bool first_trigger = true;

            // Loop over all the strips in this tagger
            for (size_t i = 0; i < strip_data_v.size(); i++)
            {
                auto & strip_data = strip_data_v[i];

                current_time = strip_data.sipm0.t1;

                // Save the first trigger
                if (strip_data.sipm_coinc and first_trigger)
                {
                    first_trigger = false;
                    trigger_ts1 = current_time;
                    trigger.reset(trigger_ts1);
                }

                // Save strips belonging to this trigger
                if (current_time - trigger_ts1 < fParams.TaggerPlaneCoincidenceWindow())
                {
                    trigger.add_strip(strip_data);
                }
                // Create a new trigger if either the current tagger is not trigger, or,
                // if it is triggered, if we are past the dead time. Also, always require
                // both sipm coincidence.
                else if ((!trigger.tagger_triggered() or
                         (trigger.tagger_triggered() and
                          !trigger.is_in_dead_time(strip_data.mac5, current_time - trigger_ts1))) and
                         strip_data.sipm_coinc)
                {
                    if (trigger.tagger_triggered()) {
                        ProcessStrips(trigger._strips);
                    }

                    // Set the current, new, trigger
                    trigger_ts1 = current_time;

                    // Reset the trigger object
                    trigger.reset(trigger_ts1);

                    // Add this strip, which created the trigger
                    trigger.add_strip(strip_data);
                }
                else if (!strip_data.sipm_coinc)
                {
                    trigger.print_no_coinc(strip_data);
                }
                else
                {
                    trigger.print_dead_time(strip_data);
                }

            } // loop over strips

            if (trigger.tagger_triggered()) {
                ProcessStrips(trigger._strips);
            }

        } // loop over taggers

        mf::LogInfo("CRTDetSimAlg") << "There are " << fData.size() << " FEBData objects." << std::endl;

    }


    void CRTDetSimAlg::FillTaggers(const uint32_t adid, const uint32_t adsid,
                                   vector<sim::AuxDetIDE> ides) {

        art::ServiceHandle<geo::Geometry> geoService;

        const geo::AuxDetGeo& adGeo = geoService->AuxDet(adid);
        const geo::AuxDetSensitiveGeo& adsGeo = adGeo.SensitiveVolume(adsid);

        // Return the vector of IDEs
        std::sort(ides.begin(), ides.end(),
                  [](const sim::AuxDetIDE & a, const sim::AuxDetIDE & b) -> bool{
                    return ((a.entryT + a.exitT)/2) < ((b.entryT + b.exitT)/2);
                  });

        // std::set<std::string> volNames = { adsGeo.TotalVolume()->GetName() };
        std::set<std::string> volNames = { adGeo.TotalVolume()->GetName() };
        std::vector<std::vector<TGeoNode const*> > paths =
          geoService->FindAllVolumePaths(volNames);

        std::string path = "";
        for (size_t inode=0; inode<paths.at(0).size(); inode++) {
          path += paths.at(0).at(inode)->GetName();
          if (inode < paths.at(0).size() - 1) {
            path += "/";
          }
        }

        TGeoManager* manager = geoService->ROOTGeoManager();
        manager->cd(path.c_str());

        // We get the array of strips first, which is the AuxDet,
        // then from the AuxDet, we get the strip by picking the
        // daughter with the ID of the AuxDetSensitive, and finally
        // from the AuxDet, we go up and pick the module and tagger
        TGeoNode* nodeArray = manager->GetCurrentNode();
        TGeoNode* nodeStrip = nodeArray->GetDaughter(adsid);
        TGeoNode* nodeModule = manager->GetMother(1);
        TGeoNode* nodeTagger = manager->GetMother(2);

        std::string volumeName = nodeStrip->GetVolume()->GetName();

        mf::LogDebug("CRTDetSimAlg") << "Strip name: " << nodeStrip->GetName()
                                     << ", number = " << nodeStrip->GetNumber()
                                     << "\n Array name: " << nodeArray->GetName()
                                     << ", number = " << nodeArray->GetNumber()
                                     << "\n Module name: " << nodeModule->GetName()
                                     << ", number = " << nodeModule->GetNumber()
                                     << "\n Tagger name: " << nodeTagger->GetName()
                                     << ", number = " << nodeTagger->GetNumber()
                                     << "\n Strip volume name: " << volumeName << std::endl;

        // Retrive the ID of this CRT module
        uint16_t mac5 = static_cast<uint16_t>(nodeModule->GetNumber());

        // Module position in parent (tagger) frame
        double origin[3] = {0, 0, 0};
        double modulePosMother[3];
        nodeModule->LocalToMaster(origin, modulePosMother);

        // Determine plane ID (1 for z > 0, 0 for z < 0 in local coordinates)
        unsigned planeID = (modulePosMother[2] > 0);

        // Determine module orientation: which way is the top (readout end)?
        bool top = (planeID == 1) ? (modulePosMother[1] > 0) : (modulePosMother[0] < 0);

        // Simulate the CRT response for each hit
        mf::LogInfo("CRTDetSimAlg") << "We have " << ides.size() << " IDE for this SimChannel." << std::endl;
        for (size_t ide_i = 0; ide_i < ides.size(); ide_i++) {

            sim::AuxDetIDE ide = ides[ide_i];

            // Finally, what is the distance from the hit (centroid of the entry
            // and exit points) to the readout end?
            double x = (ide.entryX + ide.exitX) / 2;
            double y = (ide.entryY + ide.exitY) / 2;
            double z = (ide.entryZ + ide.exitZ) / 2;

            double tTrue = (ide.entryT + ide.exitT) / 2 + fTimeOffset; // ns
            double eDep = ide.energyDeposited;

            mf::LogInfo("CRTDetSimAlg") << "True IDE with time " << tTrue
                                      << ", energy " << eDep << std::endl;

            if (tTrue < 0) {
                throw art::Exception(art::errors::LogicError)
                    << "Time cannot be negative. Check the time offset used for the CRT simulation.\n"
                    << "True time: " << (ide.entryT + ide.exitT) / 2 << "\n"
                    << "TimeOffset: " << fTimeOffset << std::endl;
            }

            double world[3] = {x, y, z};
            double svHitPosLocal[3];
            adsGeo.WorldToLocal(world, svHitPosLocal);

            // Calculate distance to the readout
            double distToReadout;
            if (top) {
                distToReadout = abs( adsGeo.HalfWidth1() - svHitPosLocal[0]);
            }
            else {
                distToReadout = abs(-adsGeo.HalfWidth1() - svHitPosLocal[0]);
            }

            // Calculate distance to fibers
            double d0 = abs(-adsGeo.HalfHeight() - svHitPosLocal[1]);  // L
            double d1 = abs( adsGeo.HalfHeight() - svHitPosLocal[1]);  // R

            // Simulate time response
            // Waveform emulation is added later, because
            // it depends on trigger time
            long npe0, npe1;
            double q0, q1;
            ChargeResponse(eDep, d0, d1, distToReadout,
                           npe0, npe1, q0, q1);

            // Time relative to trigger, accounting for propagation delay and 'walk'
            // for the fixed-threshold discriminator
            uint32_t ts1_ch0 =
              getChannelTriggerTicks(/*trigClock,*/ tTrue, npe0, distToReadout);
            uint32_t ts1_ch1 =
              getChannelTriggerTicks(/*trigClock,*/ tTrue, npe1, distToReadout);

            if (fParams.EqualizeSiPMTimes()) {
                mf::LogWarning("CRTDetSimAlg") << "EqualizeSiPMTimes is on." << std::endl;
                ts1_ch1 = ts1_ch0;
            }

            // Time relative to PPS: Random for now! (FIXME)
            uint32_t ppsTicks =
              CLHEP::RandFlat::shootInt(&fEngine, /*trigClock.Frequency()*/ fParams.ClockSpeedCRT() * 1e6);

            // Adjacent channels on a strip are numbered sequentially.
            //
            // In the AuxDetChannelMapAlg methods, channels are identified by an
            // AuxDet name (retrievable given the hit AuxDet ID) which specifies a
            // module, and a channel number from 0 to 32.
            // uint32_t moduleID = adid;
            uint32_t moduleID = mac5;
            uint32_t stripID = adsid;
            uint32_t channel0ID = 32 * moduleID + 2 * stripID + 0;
            uint32_t channel1ID = 32 * moduleID + 2 * stripID + 1;
            uint32_t sipm0ID = stripID * 2 + 0;
            uint32_t sipm1ID = stripID * 2 + 1;

            if (volumeName.find("MINOS") != std::string::npos) {continue;} // Ignoring MINOS modules for now.

            // Apply ADC threshold and strip-level coincidence (both fibers fire)
            double threshold = static_cast<double>(fParams.QThreshold());
            bool sipm_coinc = false;

            if (q0 > threshold &&
                q1 > threshold &&
                util::absDiff(ts1_ch0, ts1_ch1) < fParams.StripCoincidenceWindow())
            {
                sipm_coinc = true;
            }

            // Time ordered
            if (ts1_ch0 > ts1_ch1)
            {
                std::swap(ts1_ch0, ts1_ch1);
                std::swap(channel0ID, channel1ID);
                std::swap(q0, q1);
                std::swap(sipm0ID, sipm1ID);
            }

            SiPMData sipm0 = SiPMData(sipm0ID,
                                      channel0ID,
                                      ppsTicks,
                                      ts1_ch0,
                                      q0);
            SiPMData sipm1 = SiPMData(sipm1ID,
                                      channel1ID,
                                      ppsTicks,
                                      ts1_ch1,
                                      q1);

            StripData strip_data = StripData(mac5,
                                             planeID,
                                             sipm0,
                                             sipm1,
                                             sipm_coinc,
                                             ide);

            // Retrive the Tagger object
            Tagger& tagger = fTaggers[nodeTagger->GetName()];
            tagger.data.push_back(strip_data);

            double poss[3];
            adsGeo.LocalToWorld(origin, poss);
            mf::LogInfo("CRTDetSimAlg")
                << "CRT HIT in adid/adsid " << adid << "/" << adsid << "\n"
                << "MAC5 " << mac5 << "\n"
                << "TRUE TIME  " << tTrue << "\n"
                << "TRACK ID  " << ide.trackID << "\n"
                << "CRT HIT POS " << x << " " << y << " " << z << "\n"
                << "CRT STRIP POS " << poss[0] << " " << poss[1] << " " << poss[2] << "\n"
                << "CRT MODULE POS " << modulePosMother[0] << " "
                                     << modulePosMother[1] << " "
                                     << modulePosMother[2] << " "
                                    << "\n"
                << "CRT PATH: " << path << "\n"
                << "CRT level 0 (strip): " << nodeStrip->GetName() << "\n"
                << "CRT level 1 (array): " << nodeArray->GetName() << "\n"
                << "CRT level 2 (module): " << nodeModule->GetName() << "\n"
                << "CRT level 3 (tagger): " << nodeTagger->GetName() << "\n"
                << "CRT PLANE ID: " << planeID << "\n"
                << "CRT distToReadout: " << distToReadout << " " << (top ? "top" : "bot") << "\n"
                << "CRT Q SiPM 0: " << q0 << ", SiPM 1: " << q1
                << "CRT Ts1 SiPM 0: " << ts1_ch0 << " SiPM 1: " << ts1_ch1 << "\n";
        }
    } //end FillTaggers


    void CRTDetSimAlg::ChargeResponse(double eDep, double d0, double d1, double distToReadout, // input
                                      long & npe0, long & npe1, double & q0, double &q1) // output
    {
        std::cout << "[ChargeResponse] eDep = " << eDep << std::endl;
        // The expected number of PE, using a quadratic model for the distance
        // dependence, and scaling linearly with deposited energy.
        double qr = fParams.UseEdep() ? 1.0 * eDep / fParams.Q0() : 1.0;
        std::cout << "[ChargeResponse] qr = " << qr << std::endl;

        double npeExpected =
            fParams.NpeScaleNorm() / pow(distToReadout - fParams.NpeScaleShift(), 2) * qr;
        std::cout << "[ChargeResponse] npeExpected = " << npeExpected << std::endl;

        // Put PE on channels weighted by transverse distance across the strip,
        // using an exponential model

        double abs0 = exp(-d0 / fParams.AbsLenEff());
        double abs1 = exp(-d1 / fParams.AbsLenEff());
        double npeExp0 = npeExpected * abs0 / (abs0 + abs1);
        double npeExp1 = npeExpected * abs1 / (abs0 + abs1);

        // Observed PE (Poisson-fluctuated)
        npe0 = CLHEP::RandPoisson::shoot(&fEngine, npeExp0);
        npe1 = CLHEP::RandPoisson::shoot(&fEngine, npeExp1);

        // SiPM and ADC response: Npe to ADC counts, pedestal is added later
        q0 = CLHEP::RandGauss::shoot(&fEngine, /*fQPed + */fParams.QSlope() * npe0,
                                     fParams.QRMS() * sqrt(npe0));
        q1 = CLHEP::RandGauss::shoot(&fEngine, /*fQPed + */fParams.QSlope() * npe1,
                                     fParams.QRMS() * sqrt(npe1));
        std::cout << "[ChargeResponse] npe0 = " << npe0 << " -> q0 = " << q0 << std::endl;
        std::cout << "[ChargeResponse] npe1 = " << npe1 << " -> q1 = " << q1 << std::endl;

        // Apply saturation
        double saturation = static_cast<double>(fParams.AdcSaturation());
        if (q0 > saturation) q0 = saturation;
        if (q1 > saturation) q1 = saturation;

        mf::LogInfo("CRTSetSimAlg")
            << "CRT CHARGE RESPONSE: eDep = " << eDep
            << ", npeExpected = " << npeExpected
            << ", npe0 = " << npe0 << " -> q0 = " << q0
            << ", npe1 = " << npe1 << " -> q1 = " << q1 << std::endl;
    }

    uint32_t CRTDetSimAlg::getChannelTriggerTicks(/*detinfo::ElecClock& clock,*/
                                             float t0, float npeMean, float r)
    {
        // Hit timing, with smearing and NPE dependence
        double tDelayMean =
          fParams.TDelayNorm() *
            exp(-0.5 * pow((npeMean - fParams.TDelayShift()) / fParams.TDelaySigma(), 2)) +
          fParams.TDelayOffset();

        double tDelayRMS =
          fParams.TDelayRMSGausNorm() *
            exp(-pow(npeMean - fParams.TDelayRMSGausShift(), 2) / fParams.TDelayRMSGausSigma()) +
          fParams.TDelayRMSExpNorm() *
            exp(-(npeMean - fParams.TDelayRMSExpShift()) / fParams.TDelayRMSExpScale());

        double tDelay = CLHEP::RandGauss::shoot(&fEngine, tDelayMean, tDelayRMS);

        // Time resolution of the interpolator
        tDelay += CLHEP::RandGauss::shoot(&fEngine, 0, fParams.TResInterpolator());

        // Propagation time
        double tProp = CLHEP::RandGauss::shoot(fParams.PropDelay(), fParams.PropDelayError()) * r;

        double t = t0 + tProp + tDelay;

        // Get clock ticks
        // FIXME no clock available for CRTs, have to do it by hand
        //clock.SetTime(t / 1e3);  // SetTime takes microseconds
        float time = (t / 1e3) * fParams.ClockSpeedCRT();

        if (time < 0) {
          mf::LogWarning("CRTSetSimAlg") << "Time is negative. Check the time offset used." << std::endl;
        }

        uint32_t time_int = static_cast<uint32_t>(time);

        mf::LogInfo("CRTSetSimAlg")
            << "CRT TIMING: t0 = " << t0 << " (true G4 time)"
            << ", tDelayMean = " << tDelayMean
            << ", tDelayRMS = " << tDelayRMS
            << ", tDelay = " << tDelay
            << ", tProp = " << tProp
            << ", t = " << t
            << ", time = " << time
            << ", time_int = " << time_int << " (time in uint32_t)" << std::endl;

        return time_int; // clock.Ticks();
    }


    void CRTDetSimAlg::ClearTaggers()
    {

        fTaggers.clear();
        fData.clear();
        fAuxData.clear();
    }

} // namespace crt
} // namespace sbnd

#endif
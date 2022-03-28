#ifndef SBND_CRTDETSIMALG_CC
#define SBND_CRTDETSIMALG_CC

#include "sbndcode/CRT/CRTSimulation/CRTDetSimAlg.h"

namespace sbnd {
namespace crt {

    CRTDetSimAlg::CRTDetSimAlg(const Parameters & params, CLHEP::HepRandomEngine& engine)
    : fParams(params())
    , fEngine(engine)
    {
        ConfigureWaveform();
        fTaggers.clear();
        fFEBDatas.clear();
    }

    void CRTDetSimAlg::ConfigureWaveform()
    {

        // Normalize the waveform
        std::vector<double> wvf_x = fParams.WaveformX();
        std::vector<double> wvf_y = fParams.WaveformY();

        float max_y = *std::max_element(std::begin(wvf_y), std::end(wvf_y));
        for (auto & y : wvf_y) y /= max_y;
        // Reverse
        std::reverse(wvf_y.begin(),wvf_y.end());

        // Construct the interpolator
        fInterpolator = std::make_unique<ROOT::Math::Interpolator>
            (wvf_x.size(), ROOT::Math::Interpolation::kLINEAR);

        fInterpolator->SetData(wvf_x, wvf_y);

    }

    std::vector<std::pair<sbnd::crt::FEBData, std::vector<AuxDetIDE>>> CRTDetSimAlg::GetData()
    {
        return fData;
    }



    uint16_t CRTDetSimAlg::WaveformEmulation(uint32_t time_delay, uint16_t adc)
    {

        if (time_delay < 0) {
            throw art::Exception(art::errors::LogicError)
                << "time_delay cannot be negative." << std::endl;;
        }

        if (time_delay > fParams.WaveformX().back()) {
            // If the time delay is more than the waveform rise time, we
            // will never be able to see this signal. So return a 0 ADC value.
            return 0;
        }

        // Evaluate the waveform
        double wf = fInterpolator->Eval(time_delay);

        std::cout << "[WaveformEmulation] time_delay " << time_delay
                  << ", waveform " << wf << std::endl;

        return adc * static_cast<uint16_t>(wf);
    }


    void CRTDetSimAlg::AddADC(sbnd::crt::FEBData & feb_data, int sipmID, uint16_t adc)
    {
        uint16_t original_adc = feb_data.ADC(sipmID);
        uint16_t new_adc = original_adc + adc;

        if (new_adc > fParams.AdcSaturation()) {
            new_adc = fParams.AdcSaturation();
        }

        mf::LogDebug("CRTDetSimAlg") << "Updating ADC value for FEB " << feb_data.Mac5()
                                     << "sipmID " << sipmID << ": was " << original_adc
                                     << ", now is " << new_adc << std::endl;

        feb_data.SetADC(sipmID, new_adc);
    }




    void CRTDetSimAlg::ProcessStrips(uint32_t trigger_time, uint32_t coinc, std::vector<StripData> strips)
    {
        std::map<uint16_t, sbnd::crt::FEBData> mac_to_febdata;
        std::map<uint16_t, std::vector<AuxDetIDE>> mac_to_ides;

        bool plane0_hit, plane1_hit = false;

        // TODO Add pedestal fluctuations
        std::array<uint16_t, 32> adc_pedestal = {static_cast<uint16_t>(fParams.QPed())};

        for (auto & strip : strips) {
            if (strip.planeID == 0) plane0_hit = true;
            if (strip.planeID == 1) plane1_hit = true;

            if (!mac_to_febdata.count(strip.mac5)) {
                // Construct a new FEBData object with only pedestal values (will be filled later)
                mac_to_febdata[strip.mac5] = sbnd::crt::FEBData(strip.mac5,   // FEB ID
                                                                trigger_time, // Ts0
                                                                trigger_time, // Ts1
                                                                adc_pedestal, // ADCs
                                                                coinc);       // Coinc

                mac_to_ides[strip.mac5] = std::vector<AuxDetIDE>();
            }
        }

        // Emulate the tagger-level trigger
        if (!(plane0_hit and plane1_hit)) {
            return;
        }


        for (auto & strip : strips)
        {

            std::cout << "sipm 0 time " << strip.sipm0.t0 << ", time diff " << strip.sipm0.t0 - trigger_time << std::endl;
            std::cout << "sipm 1 time " << strip.sipm1.t0 << ", time diff " << strip.sipm1.t0 - trigger_time << std::endl;

            uint16_t adc_sipm0 = WaveformEmulation(strip.sipm0.t0 - trigger_time, strip.sipm0.adc);
            uint16_t adc_sipm1 = WaveformEmulation(strip.sipm1.t0 - trigger_time, strip.sipm1.adc);

            auto &feb_data = mac_to_febdata[strip.mac5];

            AddADC(feb_data, strip.sipm0.sipmID, adc_sipm0);
            AddADC(feb_data, strip.sipm1.sipmID, adc_sipm1);

            auto &ides = mac_to_ides[strip.mac5];
            ides.push_back(strip.ide);

        }

        for (auto const& [mac, feb_data] : mac_to_febdata) {
            fFEBDatas.push_back(feb_data);
            fData.push_back(std::make_pair(feb_data, mac_to_ides[mac]));
        }

        std::cout << "Constructed " << mac_to_febdata.size()
                  << "FEBData objects for trigger time " << trigger_time << std::endl;
    }



    void CRTDetSimAlg::CreateData()
    {

        // Loop over all the CRT Taggers and simulate triggering, dead time, ...
        for (auto iter : fTaggers) {
            auto & name = iter.first;
            auto & tagger = iter.second;

            std::cout << "[CreateData] This is tagger " << name << std::endl;

            auto & strip_data_v = tagger.data;

            // Time order the data
            std::sort(strip_data_v.begin(), strip_data_v.end(),
                      [](const StripData& strip1,
                         const StripData& strip2) {
                         return strip1.sipm0.t0 < strip2.sipm0.t0;
                         });

            uint32_t trigger_time = 0, current_time = 0, coinc = 0;
            int trigger_index = -1;
            std::vector<StripData> strips;

            // Loop over all the strips in this tagger
            for (size_t i = 0; i < strip_data_v.size(); i++)
            {
                auto & strip_data = strip_data_v[i];

                current_time = strip_data.sipm0.t0;

                if (i == 0) {
                    trigger_time = current_time;
                    coinc = strip_data.sipm0.sipmID;
                    trigger_index ++;
                    std::cout << "[CreateData]   TRIGGER TIME IS " << trigger_time << std::endl;
                }

                std::cout << "[CreateData]   This is strip " << i << " mac " << strip_data.mac5 << " on plane " << strip_data.planeID << ", with time " << current_time << std::endl;

                if (current_time - trigger_time < fParams.TaggerPlaneCoincidenceWindow())
                {
                    std::cout << "[CreateData]   -> Is in" << std::endl;
                    strips.push_back(strip_data);
                }
                else if (current_time > trigger_time + fParams.DeadTime())
                {
                    std::cout << "[CreateData]   -> Created new trigger. " << std::endl;
                    ProcessStrips(trigger_time, coinc, strips);

                    trigger_time = current_time;
                    coinc = strip_data.sipm0.sipmID;
                    trigger_index ++;

                    // Clear the current strips as we are working on a new trigger
                    strips.clear();

                    // Add this strip, which created the trigger
                    strips.push_back(strip_data);
                    std::cout << "[CreateData]   TRIGGER TIME IS " << trigger_time << std::endl;
                }
                else
                {
                    // mf::LogDebug("CRTDetSimAlg") << "Strip happened during dead time." << std::endl;
                    std::cout << "[CreateData]   -> Strip happened during dead time. " << std::endl;
                }

            }

            ProcessStrips(trigger_time, coinc, strips);

        } // loop over taggers

        std::cout << "We have " << fFEBDatas.size() << " FEBData objects." << std::endl;

        for (auto & feb_data : fFEBDatas) {
            std::cout << "FEBData mac5: " << feb_data.Mac5() << std::endl;
        }


        // Logic: For normal taggers, require at least one hit in each perpendicular
        // plane. For the bottom tagger, any hit triggers read out.
        // for (auto trg : fTaggers) {
        //     bool trigger = false;

        //     auto name = trg.first;
        //     auto tagger = trg.second;

        //     // Loop over pairs of hits
        //     for (size_t i = 0; i < tagger.size(); i++) {
        //       auto planeID1 =  tagger.planeID[i];
        //       auto time1 = tagger.time[i];

        //       for (size_t j = 0; j < tagger.size(); j++) {
        //         auto planeID2 =  tagger.planeID[j];
        //         auto time2 = tagger.time[j];

        //         if (planeID1 == planeID2 && time1 == time2) {
        //           continue;
        //         }

        //         // Two hits on different planes with proximal t0 times
        //         if (planeID1 != planeID2 && util::absDiff(time1, time2) < fTaggerPlaneCoincidenceWindow) {
        //           trigger = true;
        //           break;
        //         }
        //       }
        //     }

        //     if (trigger || name.find("TaggerBot") != std::string::npos) {
        //         // Write out all hits on a tagger when there is any coincidence FIXME this reads out everything!
        //         for (size_t d_i = 0; d_i < tagger.data.size(); d_i++) {

        //             dataCol.push_back(std::make_pair(tagger.data[d_i], tagger.ides[d_i]));

        //         }
        //     }
        // }

        // std::vector<std::pair<sbnd::crt::FEBData, std::vector<AuxDetIDE>>> dataCol;
        // return dataCol;

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

        std::cout << "Strip name: " << nodeStrip->GetName() << ", number = " << nodeStrip->GetNumber() << std::endl;
        std::cout << "Array name: " << nodeArray->GetName() << ", number = " << nodeArray->GetNumber() << std::endl;
        std::cout << "Module name: " << nodeModule->GetName() << ", number = " << nodeModule->GetNumber() << std::endl;
        std::cout << "Tagger name: " << nodeTagger->GetName() << ", number = " << nodeTagger->GetNumber() << std::endl;

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
        std::cout << "We have " << ides.size() << " IDE for this SimChannel." << std::endl;
        for (size_t ide_i = 0; ide_i < ides.size(); ide_i++) {

          sim::AuxDetIDE ide = ides[ide_i];

          // Finally, what is the distance from the hit (centroid of the entry
          // and exit points) to the readout end?
          double x = (ide.entryX + ide.exitX) / 2;
          double y = (ide.entryY + ide.exitY) / 2;
          double z = (ide.entryZ + ide.exitZ) / 2;

          double tTrue = (ide.entryT + ide.exitT) / 2 + fParams.GlobalT0Offset(); // ns
          double eDep = ide.energyDeposited;

          double world[3] = {x, y, z};
          double svHitPosLocal[3];
          adsGeo.WorldToLocal(world, svHitPosLocal);

          double distToReadout;
          if (top) {
            distToReadout = abs( adsGeo.HalfWidth1() - svHitPosLocal[0]);
          }
          else {
            distToReadout = abs(-adsGeo.HalfWidth1() - svHitPosLocal[0]);
          }

          // The expected number of PE, using a quadratic model for the distance
          // dependence, and scaling linearly with deposited energy.
          double qr = fParams.UseEdep() ? 1.0 * eDep / fParams.Q0() : 1.0;

          double npeExpected =
            fParams.NpeScaleNorm() / pow(distToReadout - fParams.NpeScaleShift(), 2) * qr;

          // Put PE on channels weighted by transverse distance across the strip,
          // using an exponential model
          double d0 = abs(-adsGeo.HalfHeight() - svHitPosLocal[1]);  // L
          double d1 = abs( adsGeo.HalfHeight() - svHitPosLocal[1]);  // R
          double abs0 = exp(-d0 / fParams.AbsLenEff());
          double abs1 = exp(-d1 / fParams.AbsLenEff());
          double npeExp0 = npeExpected * abs0 / (abs0 + abs1);
          double npeExp1 = npeExpected * abs1 / (abs0 + abs1);

          // Observed PE (Poisson-fluctuated)
          long npe0 = CLHEP::RandPoisson::shoot(&fEngine, npeExp0);
          long npe1 = CLHEP::RandPoisson::shoot(&fEngine, npeExp1);

          // Time relative to trigger, accounting for propagation delay and 'walk'
          // for the fixed-threshold discriminator
          uint32_t t0 =
            getChannelTriggerTicks(&fEngine, /*trigClock,*/ tTrue, npe0, distToReadout);
          uint32_t t1 =
            getChannelTriggerTicks(&fEngine, /*trigClock,*/ tTrue, npe1, distToReadout);

          // Time relative to PPS: Random for now! (FIXME)
          uint32_t ppsTicks =
            CLHEP::RandFlat::shootInt(&fEngine, /*trigClock.Frequency()*/ fParams.ClockSpeedCRT() * 1e6);

          // SiPM and ADC response: Npe to ADC counts, pedestal is added later
          uint32_t q0 =
            CLHEP::RandGauss::shoot(&fEngine, /*fQPed + */fParams.QSlope() * npe0, fParams.QRMS() * sqrt(npe0));
          if(q0 > fParams.AdcSaturation()) q0 = fParams.AdcSaturation();
          uint32_t q1 =
            CLHEP::RandGauss::shoot(&fEngine, /*fQPed + */fParams.QSlope() * npe1, fParams.QRMS() * sqrt(npe1));
          if(q1 > fParams.AdcSaturation()) q1 = fParams.AdcSaturation();

          // Adjacent channels on a strip are numbered sequentially.
          //
          // In the AuxDetChannelMapAlg methods, channels are identified by an
          // AuxDet name (retrievable given the hit AuxDet ID) which specifies a
          // module, and a channel number from 0 to 32.
          uint32_t moduleID = adid;
          uint32_t stripID = adsid;
          uint32_t channel0ID = 32 * moduleID + 2 * stripID + 0;
          uint32_t channel1ID = 32 * moduleID + 2 * stripID + 1;
          uint32_t sipm0ID = stripID * 2 + 0;
          uint32_t sipm1ID = stripID * 2 + 1;

          if (moduleID >= 127) {continue;}        //Ignoring MINOS modules for now.

          // Apply ADC threshold and strip-level coincidence (both fibers fire)
          if (q0 > fParams.QThreshold() &&
              q1 > fParams.QThreshold() &&
              util::absDiff(t0, t1) < fParams.StripCoincidenceWindow()) {

            // Time ordered
            if (t0 > t1) {
                std::swap(t0, t1);
                std::swap(channel0ID, channel1ID);
                std::swap(q0, q1);
                std::swap(sipm0ID, sipm1ID);
            }

            SiPMData sipm0 = SiPMData(sipm0ID,
                                      channel0ID,
                                      t0,
                                      ppsTicks,
                                      q0);
            SiPMData sipm1 = SiPMData(sipm1ID,
                                      channel1ID,
                                      t1,
                                      ppsTicks,
                                      q1);

            StripData strip_data = StripData(mac5,
                                             planeID,
                                             sipm0,
                                             sipm1,
                                             ide);

            std::cout << "Constructed StripData for mac " << mac5 << " with true time " << tTrue << " ns. and for trackID " << ide.trackID << std::endl;
            std::cout << std::endl;

            // Retrive the Tagger object
            Tagger& tagger = fTaggers[nodeTagger->GetName()];
            tagger.data.push_back(strip_data);
          }

          double poss[3];
          adsGeo.LocalToWorld(origin, poss);
          mf::LogInfo("CRTDetSimAlg")
            << "CRT HIT in " << adid << "/" << adsid << "\n"
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
            << "CRT q0: " << q0 << ", q1: " << q1 << ", t0: " << t0 << ", t1: " << t1 << ", dt: " << util::absDiff(t0,t1) << "\n";
        }


    } //end FillTaggers



    uint32_t CRTDetSimAlg::getChannelTriggerTicks(CLHEP::HepRandomEngine* engine,
                                             /*detinfo::ElecClock& clock,*/
                                             float t0, float npeMean, float r) {
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

      double tDelay = CLHEP::RandGauss::shoot(engine, tDelayMean, tDelayRMS);

      // Time resolution of the interpolator
      tDelay += CLHEP::RandGauss::shoot(engine, 0, fParams.TResInterpolator());

      // Propagation time
      double tProp = CLHEP::RandGauss::shoot(fParams.PropDelay(), fParams.PropDelayError()) * r;

      double t = t0 + tProp + tDelay;

      // Get clock ticks
      // FIXME no clock available for CRTs, have to do it by hand
      //clock.SetTime(t / 1e3);  // SetTime takes microseconds
      int time = (t / 1e3) * fParams.ClockSpeedCRT();

      mf::LogInfo("CRT")
        << "CRT TIMING: t0=" << t0
        << ", tDelayMean=" << tDelayMean << ", tDelayRMS=" << tDelayRMS
        << ", tDelay=" << tDelay << ", tDelay(interp)="
        << tDelay << ", tProp=" << tProp << ", t=" << t << /*", ticks=" << clock.Ticks() <<*/ "\n";

      return time;//clock.Ticks();
    }


    void CRTDetSimAlg::ClearTaggers() {

        fTaggers.clear();
        fFEBDatas.clear();
    }

 }//namespace crt
}//namespace sbnd

#endif
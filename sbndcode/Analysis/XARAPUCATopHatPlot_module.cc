////////////////////////////////////////////////////////////////////////
// Class:       XARAPUCATopHatPlot
// Plugin Type: analyzer (Unknown Unknown)
// File:        XARAPUCATopHatPlot_module.cc
//
// Generated at Mon Jan 27 04:25:33 2025 by Alicia Vazquez Ramos using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RawData/OpDetWaveform.h"

#include "art_root_io/TFileService.h"

#include <TTree.h>


namespace thp {
  class XARAPUCATopHatPlot;
}


class thp::XARAPUCATopHatPlot : public art::EDAnalyzer {
public:
  explicit XARAPUCATopHatPlot(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  XARAPUCATopHatPlot(XARAPUCATopHatPlot const&) = delete;
  XARAPUCATopHatPlot(XARAPUCATopHatPlot&&) = delete;
  XARAPUCATopHatPlot& operator=(XARAPUCATopHatPlot const&) = delete;
  XARAPUCATopHatPlot& operator=(XARAPUCATopHatPlot&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:
  constexpr static int BOARD_SLOT = 7; // Board slot.
  constexpr static double MIN_WINDOW = -2.5; // us.
  constexpr static double MAX_WINDOW = -2.4; // us.
  constexpr static double NS_PER_SAMPLE = 16.0; // ns.
  constexpr static double NS_TO_US = 1E-3; // Conversion factor.
  constexpr static double BASELINE_PERCETAGE = 0.1; // (10 %) 

  int fdebug;

  TTree* fTree;
  std::vector<int> events;                            // Event number for each summed waveform.
  std::vector<double> flash_times;                    // Flash time for each summed waveform.
  std::vector<int> max_ADCs;                          // ADC max for each summed waveform.
  std::vector<double> init_time_stamps;               // Init time stamp for each summed waveform.
  std::vector<int> time_stamp_numbers;                // Time stamp ID number (for this event).
  std::vector<int> summed_wvfms_baselines;            // Baseline of each summed waveform.
  std::vector <std::vector <int> > summed_wvfms;      // Summed waveform constructed (N samples).
  std::vector <std::vector <int> > rm_baselines;      // Vector with the baselines removed before constructing the summed waveform (64 baselines).

  std::string fwaveforms_module_label;

};


thp::XARAPUCATopHatPlot::XARAPUCATopHatPlot(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  fwaveforms_module_label = p.get<std::string> ("waveforms_module_label", "xarapucadecoder");
  fdebug = p.get<int> ("debug", 0);
}

void thp::XARAPUCATopHatPlot::analyze(art::Event const& e)
{
  events.clear();
  flash_times.clear();
  max_ADCs.clear();
  init_time_stamps.clear();
  time_stamp_numbers.clear();
  summed_wvfms_baselines.clear();
  summed_wvfms.clear();
  rm_baselines.clear();

  int current_event = e.id().event();

  if (fdebug) std::cout << "\nANALYZING " << e.id() << std::endl;

  if (fdebug) std::cout << "Accessing to " << fwaveforms_module_label << " label." << std::endl;
  
  art::Handle <std::vector <raw::OpDetWaveform> > waveforms_handle;
  e.getByLabel(fwaveforms_module_label, waveforms_handle);

  // The art::Handle object is not valid.
  if (!waveforms_handle.isValid()) {
    if (fdebug) std::cout << "\nWaveforms handle not valid for " << fwaveforms_module_label << "." << std::endl;
    if (fdebug) std::cout << "\n" << waveforms_handle.whyFailed() << std::endl;
  
  // The art::Handle object is empty.
  } else if (waveforms_handle->empty()) {
    if (fdebug) std::cout << "\nWaveforms handle is empty." << std::endl;
  
  // The art::Handle object is valid and not empty
  } else {
    size_t num_opdetwvfm = waveforms_handle->size();
    if (fdebug) std::cout << "Valid ( " << num_opdetwvfm << " objects)." << std::endl;

    std::map <double, std::vector<int> > summed_wvfms_map;
    std::map <double, std::vector<int> > summed_wvfm_rm_baselines_map;
    
    // Navigates over all the waveforms of this event.
    for (size_t w = 0; w < num_opdetwvfm; w++) {
      const raw::OpDetWaveform wvfm = waveforms_handle->at(w);
      unsigned int channel_number = wvfm.ChannelNumber();
      unsigned int board_slot = channel_number / 100;
      double time_stamp = wvfm.TimeStamp();

      // Filters waveforms near to the zero.
      if (time_stamp > MIN_WINDOW && time_stamp < MAX_WINDOW && board_slot == BOARD_SLOT) { // us.
        if (fdebug) std::cout << "[WVFM_ID: " << w << "]\tChannel number: " << channel_number << "\tTimeStamp: " << time_stamp << " us\tBoard slot: " << board_slot << std::endl;
        
        if (summed_wvfms_map.find(time_stamp) == summed_wvfms_map.end()) {
          summed_wvfms_map[time_stamp] = std::vector<int>(wvfm.size(), 0);
          summed_wvfm_rm_baselines_map[time_stamp] = std::vector<int>();
        }

        // Compute baseline.
        double rm_baseline = 0;
        size_t num_baseline_samples = wvfm.size() * BASELINE_PERCETAGE;
        for (size_t s = 0; s < num_baseline_samples; s++) {
          rm_baseline += wvfm[s];
        }
        rm_baseline /= num_baseline_samples;
        if (fdebug) std::cout << "REMOVING BASELINE [" << num_baseline_samples << " samples (" << BASELINE_PERCETAGE << "%)]: " << rm_baseline << " ADCs." << std::endl; 
        summed_wvfm_rm_baselines_map[time_stamp].push_back(rm_baseline);

        // Sum waveforms element by element removing the baseline.
        for (size_t s = 0; s < wvfm.size(); s++) {
          summed_wvfms_map[time_stamp][s] += std::round(wvfm[s] - rm_baseline);
          //if (fdebug) std::cout << "\t\t\t" << summed_wvfms_map[time_stamp][s] << ":::\t\t" << std::round(wvfm[s] - rm_baseline) <<  "\t\t = " << wvfm[s] << "\t - \t" << rm_baseline << std::endl;
        }
      }
    } // End sum waveforms loop. 

    
    // ========== COMPUTES STATISTICS FOR EACH SUMMED WAVEFORM ==========
    int num_time_stamp = 0;
    for (const auto& [time_stamp, summed_wvfm] : summed_wvfms_map) {
      num_time_stamp++;

      // Computes the baseline for the summed waveform.
      double summed_wvfm_baseline = 0;
      size_t num_baseline_samples = summed_wvfm.size() * BASELINE_PERCETAGE;
      for (size_t s = 0; s < num_baseline_samples; s++) {
          summed_wvfm_baseline += summed_wvfm[s];
      }
      summed_wvfm_baseline /= num_baseline_samples;
      if (fdebug) std::cout << "SUMMED BASELINE [" << num_baseline_samples << " samples (" << BASELINE_PERCETAGE << "%)]: " << std::round(summed_wvfm_baseline) << " ADCs." << std::endl; 
      
      // Computes the maximum ADC in the summed waveform and its flash peak time.
      int max_ADC = -99999;
      int max_ADC_pos = -1;
      for (size_t s = 0; s < summed_wvfm.size(); s++) {
        if (summed_wvfm[s] > max_ADC) {
          max_ADC = summed_wvfm[s];
          max_ADC_pos = s;
        }
      }
      double flash_time = max_ADC_pos * NS_PER_SAMPLE * NS_TO_US + time_stamp; // us.
      
      events.push_back(current_event);
      flash_times.push_back(flash_time);
      max_ADCs.push_back(max_ADC);
      init_time_stamps.push_back(time_stamp);
      time_stamp_numbers.push_back(num_time_stamp);
      summed_wvfms_baselines.push_back(summed_wvfm_baseline);
      
      summed_wvfms.push_back(summed_wvfm);
      rm_baselines.push_back(summed_wvfm_rm_baselines_map[time_stamp]);

  
      if (fdebug) std::cout << "DATA: ["<< num_time_stamp << "] Init time stamp " << time_stamp << " us.\tMaximum ADC: " << max_ADC << ".\tPos: " << max_ADC_pos << "\tFlash time: " << flash_time << " us." << std::endl;
    } // End maximum search loop. 
    
  }
  
  fTree->Fill();
}

void thp::XARAPUCATopHatPlot::beginJob()
{
  // Implementation of optional member function here.
  if (fdebug) std::cout << "Time window of " << MAX_WINDOW - MIN_WINDOW << " us: [" << MIN_WINDOW << ", " << MAX_WINDOW << "] us." << std::endl;
  if (fdebug) std::cout << "Nanoseconds per sample: " << NS_PER_SAMPLE << " ns." << std::endl;

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("thp", "Top Hat Plot data");

  fTree->Branch("event", &events);
  fTree->Branch("flash_time", &flash_times);
  fTree->Branch("max_ADC", &max_ADCs);
  fTree->Branch("init_time_stamp", &init_time_stamps);
  fTree->Branch("num_time_stamp", &time_stamp_numbers);
  fTree->Branch("summed_wvfm_baseline", &summed_wvfms_baselines);

  fTree->Branch("summed_wvfm", &summed_wvfms);
  fTree->Branch("rm_baselines", &rm_baselines);

}

void thp::XARAPUCATopHatPlot::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(thp::XARAPUCATopHatPlot)

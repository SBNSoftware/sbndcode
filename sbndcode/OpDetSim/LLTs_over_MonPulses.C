//#include "/exp/sbnd/app/users/npallat/NC_disappearance/TriggerValidation/simplePlottingFunctions.h"

void makeComparisonHist(TH1D* h1, TH1D* h2,const std::string& label1, const std::string& label2, const std::string& save_filename, const std::string& plot_title, const std::string& x_label, const std::string& y_label) {
    TCanvas* c = new TCanvas("c_compare", "c_compare", 800, 600);
    gStyle->SetOptStat(0);

    h1->SetLineColor(kBlue);
    h2->SetLineColor(kRed);

    h1->SetTitle(plot_title.c_str());
    h1->GetXaxis()->SetTitle(x_label.c_str());
    h1->GetYaxis()->SetTitle(y_label.c_str());

    // Make both visible
    double max_val = std::max(h1->GetMaximum(), h2->GetMaximum());
    h1->SetMaximum(max_val * 1.2);
    
    h1->Draw("HIST E1");
    h2->Draw("HIST E1 SAME");
    
    // Legend
    TLegend* leg = new TLegend(0.7, 0.75, 0.88, 0.88);
    leg->AddEntry(h1, label1.c_str(), "l");
    leg->AddEntry(h2, label2.c_str(), "l");
    leg->Draw();
    c->SaveAs(save_filename.c_str());
    
    delete leg;
    delete c;
};

void makeAreaNormComparisonHist(TH1D* h1, TH1D* h2, const std::string& label1, const std::string& label2, const std::string& save_filename, const std::string& plot_title, const std::string& x_label, const std::string& y_label) {
    TCanvas* c = new TCanvas("c_compare", "c_compare", 800, 600);
    gStyle->SetOptStat(0);

    TH1D* h1n = (TH1D*)h1->Clone("h1n_norm");
    TH1D* h2n = (TH1D*)h2->Clone("h2n_norm");
    h1n->SetDirectory(0);
    h2n->SetDirectory(0);

    if (h1n->Integral() > 0) h1n->Scale(1.0 / h1n->Integral());
    if (h2n->Integral() > 0) h2n->Scale(1.0 / h2n->Integral());

    h1n->SetLineColor(kBlue);
    h2n->SetLineColor(kRed);

    h1n->SetTitle(plot_title.c_str());
    h1n->GetXaxis()->SetTitle(x_label.c_str());
    h1n->GetYaxis()->SetTitle(y_label.c_str());

    double max_val = std::max(h1n->GetMaximum(), h2n->GetMaximum());
    h1n->SetMaximum(max_val * 1.2);

    h1n->Draw("HIST E1");
    h2n->Draw("HIST E1 SAME");

    TLegend* leg = new TLegend(0.7, 0.75, 0.88, 0.88);
    leg->AddEntry(h1n, label1.c_str(), "l");
    leg->AddEntry(h2n, label2.c_str(), "l");
    leg->Draw();
    c->SaveAs(save_filename.c_str());

    delete leg;
    delete h1n;
    delete h2n;
    delete c;
};

int findPeak(const std::vector<int>& monpulse, int start=-1, int end=-1) {
  if (monpulse.empty()) return 0;

  const int s = static_cast<int>(monpulse.size());
  const int first = std::max(0, start);
  const int last = (end < 0) ? s - 1 : std::min(end, s - 1);

  if (first > last) return 0;

  return *std::max_element(monpulse.begin() + first,
                           monpulse.begin() + last + 1);
};


int findPeakTime(const std::vector<int>& monpulse, int start=-1, int end=-1) {
  if (monpulse.empty()) return 0;

  const int s = static_cast<int>(monpulse.size());
  const int first = std::max(0, start);
  const int last = (end < 0) ? s - 1 : std::min(end, s - 1);

  if (first > last) return 0;

  return static_cast<int>(
      std::distance(monpulse.begin(),
                    std::max_element(monpulse.begin() + first,
                                     monpulse.begin() + last + 1)));
};


auto makeHistSave = [](const std::vector<int>& data, const char* filename, int nbins, int min_bin, int max_bin, const char* xlabel="Missing Flash Trigger", const char* ylabel="Number of MonPulses") {

    static int hist_counter = 0;
    TH1D* h = new TH1D(Form("h_%d", hist_counter++), " ;Value;Counts", nbins, min_bin, max_bin);
    for (double x : data) h->Fill(x);

    TCanvas* c = new TCanvas("c", "canvas", 800, 600);
    h->GetXaxis()->SetTitle(xlabel);
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->SetTitle(ylabel);
    h->Draw();
    c->SaveAs(filename);

    delete c;
    return h;

};


void makeHist(const std::vector<int>& data, const char* filename, int nbins, int min_bin, int max_bin, const char* xlabel="Missing Flash Trigger", const char* ylabel="Number of MonPulses") {
    
    TH1D* h = new TH1D("h", " ;Value;Counts", nbins, min_bin, max_bin);
    for (double x : data) h->Fill(x);

    TCanvas* c = new TCanvas("c", "canvas", 800, 600);
    h->GetXaxis()->SetTitle(xlabel);
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->SetTitle(ylabel);
    h->Draw();
    c->SaveAs(filename);

    delete c;
    delete h;

};

auto makeFloatHist = [](const std::vector<float> vecToPlot, const char* title, const char* x_label, const char* y_label, int nbins, double xlow, double xhigh, const char* savename=nullptr) {

  // Histogram 
  TH1F *h = new TH1F("h", title, nbins, xlow, xhigh);
  h->SetDirectory(0);

  // Fill hist
  for (auto i = 0; i < vecToPlot.size(); ++i) h->Fill(vecToPlot[i]);

  if (savename) {
    TCanvas *c2 = new TCanvas("c2", title, 800, 600);
    h->Draw();
    h->GetXaxis()->SetTitle(x_label);
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->SetTitle(y_label);
    c2->SaveAs(savename);
    delete c2;
  }

  return h;
};


auto makeLLTsMonPulsePlot = [](const std::vector<int>& monpulse, std::vector<int> llts_18, std::vector<int> llts_20, std::vector<int> llts_21, int event, int num_pulse, const char* tag, int start=-1, int end=-1) {
  // Make histogram for this monpulse
  int s = int(monpulse.size());
  TH1I *h = new TH1I("h", "MonPulse", s, 0, s-1);
  h->SetDirectory(0);
  // loop through this monpulse
  for (int m = 0; m < s; ++m) {
    h->SetBinContent(m + 1, monpulse[m]);
  }
  // Draw histogram
  TCanvas *c0 = new TCanvas("c0", "MonPulse", 800, 600);
  h->Draw();
  int num_llt18_this_monpulse = 0;
  int num_llt20_this_monpulse = 0;
  int num_llt21_this_monpulse = 0;
  // Add vertical lines for LLTs
  for (int llt18 : llts_18) { 
    Double_t x_val = llt18;
    Double_t y_min = 0;
    Double_t y_max = h->GetMaximum(); 
    if (!(start == -9999 && end == -9999) && end >= start) {
      if (x_val >= start && x_val <= end) {
        TLine *line = new TLine(x_val, y_min, x_val, y_max); 
        line->SetLineColor(kRed);
        line->SetLineWidth(1);
        line->Draw();
        ++num_llt18_this_monpulse;
      }
    } else {
        TLine *line = new TLine(x_val, y_min, x_val, y_max); 
        line->SetLineColor(kRed);
        line->SetLineWidth(1);
        line->Draw();
        ++num_llt18_this_monpulse;
    }
  }
  for (int llt20 : llts_20) { 
    Double_t x_val = llt20;
    Double_t y_min = 0;
    Double_t y_max = h->GetMaximum(); 
    if (!(start == -9999 && end == -9999) && end >= start) {
      if (x_val >= start && x_val <= end) {
        TLine *line = new TLine(x_val, y_min, x_val, y_max); 
        line->SetLineColor(kGreen+1);
        line->SetLineWidth(1);
        line->Draw();
        ++num_llt20_this_monpulse;
      } 
    } else {
        TLine *line = new TLine(x_val, y_min, x_val, y_max); 
        line->SetLineColor(kGreen+1);
        line->SetLineWidth(1);
        line->Draw();
        ++num_llt20_this_monpulse;
    }
  }
  for (int llt21 : llts_21) { 
    Double_t x_val = llt21;
    Double_t y_min = 0;
    Double_t y_max = h->GetMaximum(); 
    if (!(start == -9999 && end == -9999) && end >= start) {
      if (x_val >= start && x_val <= end) {
        TLine *line = new TLine(x_val, y_min, x_val, y_max); 
        line->SetLineColor(kMagenta);
        line->SetLineWidth(1);
        line->Draw();
        ++num_llt21_this_monpulse;
      } 
    } else {
        TLine *line = new TLine(x_val, y_min, x_val, y_max); 
        line->SetLineColor(kMagenta);
        line->SetLineWidth(1);
        line->Draw();
        ++num_llt21_this_monpulse;
    }
  }
  char title[100];
  if (strcmp(tag, "PTB") == 0) sprintf(title, "PTB_monpulse_event_%d_numpulse_%d_llt18_%d_llt20_%d_llt21_%d.png", event, num_pulse, 
                                          num_llt18_this_monpulse, num_llt20_this_monpulse, num_llt21_this_monpulse);
  else if (strcmp(tag, "MonPulse") == 0) sprintf(title, "Sim_monpulse_event_%d_numpulse_%d_llt18_%d_llt20_%d_llt21_%d.png", event, num_pulse, 
                                          num_llt18_this_monpulse, num_llt20_this_monpulse, num_llt21_this_monpulse);
  else sprintf(title, "monpulse_event_%d_numpulse_%d_llt18_%d_llt20_%d_llt21_%d.png", event, num_pulse, num_llt18_this_monpulse, num_llt20_this_monpulse, num_llt21_this_monpulse);
  h->GetXaxis()->SetTitle("Ticks");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->SetTitle("Number of PMT Pairs Above Threshold");
  if (!(start == -9999 && end == -9999) && end >= start) h->GetXaxis()->SetRangeUser(start, end);
  else std::cout<<"Plotting whole MonPulse"<<std::endl;

  c0->SaveAs(title);

  delete c0;
  delete h;
};


// --------------------------------------------------------------------------
// MAIN FUNCTION
// --------------------------------------------------------------------------

void LLTs_over_MonPulses() {

  int missing_peak_3 = 0;
  int missing_peak_more = 0;
  int missing_peak_less = 0;
  int ptb_missing_too = 0;
  int ptb_missing_new = 0;

  int num_monpulses = 0;

  std::vector<int> passed_events = {};
  std::vector<int> passed_events_PTB = {};

  std::vector<int> first_llt_diff = {};
  std::vector<int> peak_llt_diff = {};
  std::vector<int> deltaT_A = {}; // A: DeltaT distribution between MonPulse waveform time stamp and time tick with MaxMonPulse
  std::vector<int> deltaT_B = {}; // B: DeltaT distribution between MonPulse waveform time stamp and first LLT20 in the waveform
  std::vector<int> deltaT_B_21 = {}; // B: DeltaT distribution between MonPulse waveform time stamp and first **LLT21** in the waveform

  std::vector<int> mon_vals_monpulse_21_direct = {};
  std::vector<int> mon_vals_monpulse_20_direct = {};
  std::vector<int> mon_vals_monpulse_18_direct = {};

  std::vector<int> mon_vals_monpulse_21 = {};
  std::vector<int> mon_vals_ptb_21 = {};
  std::vector<int> mon_vals_ptb_21_tolerance = {};
  std::vector<int> mon_vals_monpulse_20 = {};
  std::vector<int> mon_vals_ptb_20 = {};
  std::vector<int> mon_vals_ptb_20_tolerance = {};
  std::vector<int> mon_vals_monpulse_18 = {};
  std::vector<int> mon_vals_ptb_18 = {};
  std::vector<int> mon_vals_ptb_18_tolerance = {};

  // Open file
  TFile *f = TFile::Open("/exp/sbnd/data/users/npallat/TriggerRuns/Run19737/BeamRateRes19737.root", "READ");
  //TFile *f = TFile::Open("BeamRateResRun19663.root", "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Error opening file!" << std::endl;
    return;
  }

  // Get tree
  TTree *tree = dynamic_cast<TTree*>(f->Get("BeamCalib/TriggerMetrics"));
  if (!tree) {
    std::cerr << "TriggerMetrics tree not found!" << std::endl;
    f->ls();
    return;
  }

  // Decode / analysis constants
  constexpr int LLT20 = 20;
  constexpr int LLT21 = 21;
  constexpr int LLT18 = 18;
  constexpr int ACCEPTANCE_WINDOW = 100;
  constexpr int PTB_TOLERANCE = 60;
  constexpr double PTB_OFFSET = -222.1;

  // Declare branch variables
  std::vector<int> *LLT = nullptr;
  std::vector<int> *HLT = nullptr;
  std::vector<ULong64_t> *LLTTimes = nullptr;
  std::vector<ULong64_t> *HLTTimes = nullptr;
  std::vector<std::vector<int>> *SimLLT = nullptr;
  std::vector<std::vector<int>> *SimLLTTimes = nullptr;
  int numMonPulses = 0;
  int numMergedMonPulses = 0;
  int MissingFlashTriggers = 0;
  std::vector<double> *MonPulseTimes = nullptr;
  std::vector<std::vector<int>> *MonPulses = nullptr;
  int passedTrigger = 0;
 
  // Set branch addresses
  tree->SetBranchAddress("LLT", &LLT);
  tree->SetBranchAddress("HLT", &HLT);
  tree->SetBranchAddress("LLTTimes", &LLTTimes);
  tree->SetBranchAddress("HLTTimes", &HLTTimes);
  tree->SetBranchAddress("SimLLT", &SimLLT);
  tree->SetBranchAddress("SimLLTTimes", &SimLLTTimes);
  tree->SetBranchAddress("numMonPulses", &numMonPulses);
  tree->SetBranchAddress("numMergedMonPulses", &numMergedMonPulses);
  tree->SetBranchAddress("MonPulseTimes", &MonPulseTimes);
  tree->SetBranchAddress("MonPulses", &MonPulses);
  tree->SetBranchAddress("passedTrigger", &passedTrigger);

  std::vector<int> Sim_LLTs18_this_event;
  std::vector<int> Sim_LLTs20_this_event;
  std::vector<int> Sim_LLTs21_this_event;

  std::vector<int> LLT21_Times;
  std::vector<int> LLT20_Times;
  std::vector<int> LLT18_Times;
  std::vector<int> HLT26_Times;

  // new
  std::vector<int> PTB_LLTs;
  std::vector<int> PTB_LLTTimes;

  std::vector<int> All_Sim_LLTs;
  int All_LLTs_this_event = 0;

  std::vector<int> LLTs_per_event;
  std::vector<int> SimLLTs_per_event;

  std::vector<int> LLTs_18_per_event;
  std::vector<int> SimLLTs_18_per_event;
  int LLTs_18_this_event = 0;
  int SimLLTs_18_this_event = 0;

  std::vector<int> LLTs_20_per_event;
  std::vector<int> SimLLTs_20_per_event;
  int LLTs_20_this_event = 0;
  int SimLLTs_20_this_event = 0;

  std::vector<int> LLTs_21_per_event;
  std::vector<int> SimLLTs_21_per_event;
  int LLTs_21_this_event = 0;
  int SimLLTs_21_this_event = 0;

  std::vector<int> events_0;
  std::vector<int> Simevents_0;

  std::vector<float> percent_missing;
  std::vector<float> percent_missing_merged;

  double max_time = 0;
  double min_time = 0;
  std::vector<float> pull_window_lengths;


  // Loop over entries
  Long64_t nEntries = tree->GetEntries();
 
  std::vector<int> missing_flash_this_event_monpulses = {};
  std::vector<int> missing_flash_this_event_ptb = {};
  std::vector<int> missing_flash_this_monpulse = {};
  std::vector<int> missing_event_this_monpulse = {};
  std::vector<int> missing_18_this_monpulse = {};

  std::vector<int> missing_llt18_ptb = {};
  std::vector<int> missing_llt20_ptb = {};
  std::vector<int> missing_llt21_ptb = {};
  std::vector<int> missing_llts_ptb = {};

  std::vector<int> peak_height_missing = {};
  std::vector<int> peak_height_not_missing = {};

  std::vector<int> time_bn_LLT18 = {};
  std::vector<int> time_bn_LLT18_PTB = {};
  std::vector<int> time_bn_LLT20 = {};
  std::vector<int> time_bn_LLT20_PTB = {};
  std::vector<int> time_bn_LLT21 = {};
  std::vector<int> time_bn_LLT21_PTB = {};

  for (Long64_t i = 0; i < nEntries; ++i) {

    tree->GetEntry(i);

    if (passedTrigger > 1) { std::cout<<"over 1"<<std::endl; }

    if (passedTrigger > 0) { passed_events.push_back(1); }
    else passed_events.push_back(0);
 
    // clear and zero per entry variables 
    Sim_LLTs18_this_event.clear();
    Sim_LLTs20_this_event.clear();
    Sim_LLTs21_this_event.clear();
    LLT21_Times.clear();
    LLT20_Times.clear();
    LLT18_Times.clear();
    HLT26_Times.clear();
    int passed_PTB_curr_event = 0; 
    PTB_LLTs.clear();
    PTB_LLTTimes.clear();
    All_LLTs_this_event = 0;
    LLTs_18_this_event = 0;
    SimLLTs_18_this_event = 0;
    LLTs_20_this_event = 0;
    SimLLTs_20_this_event = 0;
    LLTs_21_this_event = 0;
    SimLLTs_21_this_event = 0;

    // Extract HLTs from the PTB
    std::vector<ULong64_t> HLT1_Times;
    for (size_t h = 0; h < HLT->size(); ++h) {
      if (HLT->at(h) == 1 || HLT->at(h) == 3) HLT1_Times.push_back(HLTTimes->at(h));
      if (HLT->at(h) == 2 || HLT->at(h) == 4) passed_PTB_curr_event = 1;
    }
    if (passed_PTB_curr_event == 1) passed_events_PTB.push_back(1);
    else passed_events_PTB.push_back(0);

    // Check
    if (HLT1_Times.size() != 1) {
      std::cout << "Error: HLT1_Times.size() is "
                << HLT1_Times.size() << std::endl;
      continue;
    }

    // HLT 26 check
    HLT26_Times.clear();
    for (size_t k = 0; k < HLT->size(); ++k) {
      if (HLT->at(k) == 26) HLT26_Times.push_back((((int)(HLTTimes->at(k) - HLT1_Times[0]))/2 - (int)(500.0*MonPulseTimes->at(0))));
    }


    double curr_time_PTB = 0;
    const double PTB_offset = PTB_OFFSET;
    double prev_LLT18_time_PTB = 0;
    double prev_LLT20_time_PTB = 0;
    double prev_LLT21_time_PTB = 0;

    // Extract LLTs from the PTB
    for (size_t j = 0; j < LLT->size(); ++j) {
      if (HLT1_Times[0]-1500000 < LLTTimes->at(j) && HLT1_Times[0]+1500000 > LLTTimes->at(j)) { // only save if within 3ms pull window
        // in ticks from event timestamp (can be - or +)
        curr_time_PTB = ((double)(int)(LLTTimes->at(j) - HLT1_Times[0])) / 2.0 - 500.0*MonPulseTimes->at(0) + PTB_offset; 

        // LLT 18, 20, 21
        if (LLT->at(j) == LLT20 || LLT->at(j) == LLT21 || LLT->at(j) == LLT18) {
          PTB_LLTs.push_back(LLT->at(j));
          PTB_LLTTimes.push_back((int)curr_time_PTB);
        }
        // LLT 20
        if (LLT->at(j) == LLT20) { 
          LLT20_Times.push_back((int)curr_time_PTB);
          ++LLTs_20_this_event;
          if (prev_LLT20_time_PTB != 0) time_bn_LLT20_PTB.push_back((int)(curr_time_PTB - prev_LLT20_time_PTB));
          prev_LLT20_time_PTB = curr_time_PTB;
        }
        // LLT 21
        if (LLT->at(j) == LLT21) { 
          LLT21_Times.push_back((int)curr_time_PTB);
          ++LLTs_21_this_event;
          if (prev_LLT21_time_PTB != 0) time_bn_LLT21_PTB.push_back((int)(curr_time_PTB - prev_LLT21_time_PTB));
          prev_LLT21_time_PTB = curr_time_PTB;
        }
        // LLT 18
        if (LLT->at(j) == LLT18) { 
          LLT18_Times.push_back((int)curr_time_PTB);
          ++LLTs_18_this_event;
          if (prev_LLT18_time_PTB != 0) time_bn_LLT18_PTB.push_back((int)(curr_time_PTB - prev_LLT18_time_PTB));
          prev_LLT18_time_PTB = curr_time_PTB;
        }
      }
    }


    // Check
    if (SimLLT->size() != MonPulses->size()) { 
      std::cout<<"Error: Missing LLTs for at least one entire MonPulse: SimLLT->size() = "<<SimLLT->size()<<" MonPulses->size(): "<<MonPulses->size()<<std::endl; 
      continue; 
    }

    // for padding MonPulses
    std::vector<int> MonPulses_this_event;
    std::vector<int> padding;
    int curr_padding;
    size_t total_size = 0;
    auto append = [&](const std::vector<int>& v, size_t zeros) {
        MonPulses_this_event.insert(MonPulses_this_event.end(), v.begin(), v.end());
        MonPulses_this_event.insert(MonPulses_this_event.end(), zeros, 0);
    };


    // Extract LLTs from MonPulses
    std::vector<int> Sim_LLTs18_this_pulse;
    std::vector<int> Sim_LLTs20_this_pulse;
    std::vector<int> Sim_LLTs21_this_pulse;

    std::vector<int> starts;
    std::vector<int> ends;
    int max_t = 0;
    int min_t = 0;
    double curr_time = 0;
    int prev_LLT18_time = 0;
    int prev_LLT20_time = 0;
    int prev_LLT21_time = 0;

    for (size_t k = 0; k < SimLLT->size(); ++k) {
      Sim_LLTs18_this_pulse.clear();
      Sim_LLTs20_this_pulse.clear();
      Sim_LLTs21_this_pulse.clear();
      max_t = 0;
      min_t = -1;
      // padding
      curr_padding = 0;
      if (k < SimLLT->size()-1) {
        //curr_padding = 500*(MonPulseTimes->at(k+1)) - (500*(MonPulseTimes->at(k)) + MonPulses->at(k).size());
        curr_padding = 500*(MonPulseTimes->at(k+1)) - (500*(MonPulseTimes->at(k)) + MonPulses->at(k).size());
      }
      if (curr_padding < 0) throw std::runtime_error("Negative padding computed");
      padding.push_back(curr_padding);
      size_t prev_total_size = total_size;
      total_size = total_size + (MonPulses->at(k)).size() + curr_padding;

      for (size_t l = 0; l < (SimLLT->at(k)).size(); ++l) {

        if ((SimLLT->at(k)).at(l) == LLT20) { 
          curr_time = (SimLLTTimes->at(k)).at(l) + prev_total_size;

          Sim_LLTs20_this_pulse.push_back((SimLLTTimes->at(k)).at(l));
          Sim_LLTs20_this_event.push_back( curr_time ); // in ticks!
          ++SimLLTs_20_this_event;

          All_Sim_LLTs.push_back(LLT20);
          ++All_LLTs_this_event;

          if (curr_time > max_t) max_t = curr_time;
          if (min_t == -1) min_t = curr_time;
          if (prev_LLT20_time != 0) time_bn_LLT20.push_back(curr_time-prev_LLT20_time);
          prev_LLT20_time = curr_time;
        } 

        if ((SimLLT->at(k)).at(l) == LLT21) { 
          curr_time = (SimLLTTimes->at(k)).at(l) + prev_total_size;

          Sim_LLTs21_this_pulse.push_back((SimLLTTimes->at(k)).at(l));
          Sim_LLTs21_this_event.push_back( curr_time ); // in ticks!
          ++SimLLTs_21_this_event;

          All_Sim_LLTs.push_back(LLT21);
          ++All_LLTs_this_event;

          if (curr_time > max_t) max_t = curr_time;
          if (min_t == -1) min_t = curr_time;
          if (prev_LLT21_time != 0) time_bn_LLT21.push_back(curr_time-prev_LLT21_time);
          prev_LLT21_time = curr_time;
        } 

        if ((SimLLT->at(k)).at(l) == LLT18) { 
          curr_time = (SimLLTTimes->at(k)).at(l) + prev_total_size;

          Sim_LLTs18_this_pulse.push_back((SimLLTTimes->at(k)).at(l));
          Sim_LLTs18_this_event.push_back( curr_time ); // in ticks!
          ++SimLLTs_18_this_event;

          All_Sim_LLTs.push_back(LLT18);
          ++All_LLTs_this_event;

          if (curr_time > max_t) max_t = curr_time;
          if (min_t == -1) min_t = curr_time;
          if (prev_LLT18_time != 0) time_bn_LLT18.push_back(curr_time-prev_LLT18_time);
          prev_LLT18_time = curr_time;
        } 

      }
      //makeLLTsMonPulsePlot(MonPulses->at(k), Sim_LLTs18_this_pulse, Sim_LLTs20_this_pulse, Sim_LLTs21_this_pulse, i, k, "MonPulse");
      starts.push_back(min_t-1000);
      ends.push_back(max_t+1000);

      for (int llt_time : Sim_LLTs21_this_pulse) mon_vals_monpulse_21_direct.push_back((MonPulses->at(k))[llt_time]);
      for (int llt_time : Sim_LLTs20_this_pulse) mon_vals_monpulse_20_direct.push_back((MonPulses->at(k))[llt_time]);
      for (int llt_time : Sim_LLTs18_this_pulse) mon_vals_monpulse_18_direct.push_back((MonPulses->at(k))[llt_time]);

      // DELTA T A
      int curr_peak = findPeakTime(MonPulses->at(k), 0, (MonPulses->at(k)).size());
      deltaT_A.push_back(curr_peak); // - MonPulseTimes->at(k));

      // DELTA T B 
      if (Sim_LLTs20_this_pulse.size() > 0) { deltaT_B.push_back(Sim_LLTs20_this_pulse.at(0)); missing_event_this_monpulse.push_back(0); }
      else missing_event_this_monpulse.push_back(1);
      //else deltaT_B.push_back(0);
      if (Sim_LLTs21_this_pulse.size() > 0) { deltaT_B_21.push_back(Sim_LLTs21_this_pulse.at(0)); missing_flash_this_monpulse.push_back(0); }
      else missing_flash_this_monpulse.push_back(1);
      //else deltaT_B_21.push_back(0);
      if (Sim_LLTs18_this_pulse.size() > 0) { missing_18_this_monpulse.push_back(0); }
      else missing_18_this_monpulse.push_back(1);
    }
    // Missing flash LLTs
    if (SimLLTs_21_this_event == 0) missing_flash_this_event_monpulses.push_back(1);
    else missing_flash_this_event_monpulses.push_back(0);
    if (LLTs_21_this_event == 0) missing_flash_this_event_ptb.push_back(1);
    else missing_flash_this_event_ptb.push_back(0);

    // From MonPulses
    MonPulses_this_event.reserve(total_size);
    for (size_t m = 0; m < SimLLT->size(); ++m) append(MonPulses->at(m), padding[m]);
    //if (i < 10) makeLLTsMonPulsePlot(MonPulses_this_event, Sim_LLTs18_this_event, Sim_LLTs20_this_event, Sim_LLTs21_this_event, i, 1000, "MonPulse");
    // From PTB
    //if (i < 10) makeLLTsMonPulsePlot(MonPulses_this_event, LLT18_Times, LLT20_Times, LLT21_Times, i, 2000, "PTB");

    int start = -1;
    int end = -1;

    int ptb_21 = 0;
    int monpulse_21 = 0;

    int peak = 0;

    int max_events_save = nEntries;
    bool saveMissingOnly = true;
    int num_llt18_this_monpulse_sim = 0;
    int num_llt20_this_monpulse_sim = 0;
    int num_llt21_this_monpulse_sim = 0;
    int num_llt18_this_monpulse_ptb = 0;
    int num_llt20_this_monpulse_ptb = 0;
    int num_llt21_this_monpulse_ptb = 0;
    int max_val_monpulse = 0;
    if (i <= max_events_save) { 
      // FULL plot for first x events
      //makeLLTsMonPulsePlot(MonPulses_this_event, Sim_LLTs18_this_event, Sim_LLTs20_this_event, Sim_LLTs21_this_event, i, 1000, "MonPulse");
      // From PTB
      //makeLLTsMonPulsePlot(MonPulses_this_event, LLT18_Times, LLT20_Times, LLT21_Times, i, 2000, "PTB");
      // Sliced plots for first x events
      for (size_t m = 0; m < MonPulseTimes->size(); ++m) {
        start = starts.at(m);
        end = ends.at(m);

        num_llt18_this_monpulse_sim = 0;
        num_llt20_this_monpulse_sim = 0;
        num_llt21_this_monpulse_sim = 0;
        num_llt18_this_monpulse_ptb = 0;
        num_llt20_this_monpulse_ptb = 0;
        num_llt21_this_monpulse_ptb = 0;
        max_val_monpulse = 0;

        for (int llt : Sim_LLTs18_this_event) if (llt >= start && llt <= end) ++num_llt18_this_monpulse_sim;
        for (int llt : Sim_LLTs20_this_event) if (llt >= start && llt <= end) ++num_llt20_this_monpulse_sim;
        for (int llt : Sim_LLTs21_this_event) if (llt >= start && llt <= end) ++num_llt21_this_monpulse_sim;
        for (int llt : LLT18_Times) if (llt >= start && llt <= end) ++num_llt18_this_monpulse_ptb;
        for (int llt : LLT20_Times) if (llt >= start && llt <= end) ++num_llt20_this_monpulse_ptb;
        for (int llt : LLT21_Times) if (llt >= start && llt <= end)  ++num_llt21_this_monpulse_ptb; 

        for (int m_i = start; m_i <= end; ++m_i) if (max_val_monpulse < MonPulses_this_event[m_i]) max_val_monpulse = MonPulses_this_event[m_i];

        if ((num_llt18_this_monpulse_sim + num_llt20_this_monpulse_sim + num_llt21_this_monpulse_sim == 0) 
             || (num_llt18_this_monpulse_ptb + num_llt20_this_monpulse_ptb + num_llt21_this_monpulse_ptb == 0)) {
          //makeLLTsMonPulsePlot(MonPulses_this_event, Sim_LLTs18_this_event, Sim_LLTs20_this_event, Sim_LLTs21_this_event, i, 1000+m, "MonPulse", start, end);
          //makeLLTsMonPulsePlot(MonPulses_this_event, LLT18_Times, LLT20_Times, LLT21_Times, i, 2000+m, "PTB", start, end);
          peak_height_missing.push_back(max_val_monpulse);
        } else {
          peak_height_not_missing.push_back(max_val_monpulse);
        }

        if (num_llt18_this_monpulse_ptb + num_llt20_this_monpulse_ptb + num_llt21_this_monpulse_ptb == 0) missing_llts_ptb.push_back(1);
        else missing_llts_ptb.push_back(0);
        if (num_llt18_this_monpulse_ptb == 0) missing_llt18_ptb.push_back(1);
        else missing_llt18_ptb.push_back(0);
        if (num_llt20_this_monpulse_ptb == 0) missing_llt20_ptb.push_back(1);
        else missing_llt20_ptb.push_back(0);
        if (num_llt21_this_monpulse_ptb == 0) missing_llt21_ptb.push_back(1);
        else missing_llt21_ptb.push_back(0);
        
      }
      num_monpulses = num_monpulses + MonPulses->size();
    } // i == 5


    // first llt of each monpulse
    std::vector<int> first_llt_diff_this_event = {};
    for (size_t m = 0; m < MonPulseTimes->size(); ++m) {

      int first_llt_monpulse = -1;
      int first_llt_ptb = -1;
      bool first =  true;
      bool first_ptb =  true;

      start = starts.at(m);
      end = ends.at(m);
      if (start >= 0 && end > start) {
        for (size_t n = 0; n < Sim_LLTs21_this_event.size(); ++n) {
          if (Sim_LLTs21_this_event[n] >= start && Sim_LLTs21_this_event[n] <= end && first) {
            first_llt_monpulse = Sim_LLTs21_this_event[n];
            first = false;
          }
        }
        for (size_t p = 0; p < LLT21_Times.size(); ++p) {
          if (LLT21_Times[p] >= start && LLT21_Times[p] <= end && first_ptb) {
            first_llt_ptb = LLT21_Times[p];
            first_ptb = false;
          }
        }
        if (first_llt_monpulse >= 0 && first_llt_ptb >= 0) { 
          first_llt_diff.push_back(first_llt_monpulse - first_llt_ptb); 
          first_llt_diff_this_event.push_back(first_llt_monpulse - first_llt_ptb); 
        }
        else {
          first_llt_diff.push_back(-499);
          first_llt_diff_this_event.push_back(-499);
        }
      }
    }

    // peak llt of each monpulse
    int peak_time = -1;
    for (size_t m = 0; m < MonPulseTimes->size(); ++m) {

      int peak_time_llt_monpulse = -1;
      int peak_time_llt_ptb = -1;

      start = starts.at(m);
      end = ends.at(m);
      if (start >= 0 && end > start) {
        peak_time = findPeakTime(MonPulses_this_event, start, end);

        const int acceptance_window = ACCEPTANCE_WINDOW;

        /*for (size_t n = 0; n < Sim_LLTs21_this_event.size(); ++n) {
          if (Sim_LLTs21_this_event[n] >= peak_time-acceptance_window && Sim_LLTs21_this_event[n] <= peak_time+acceptance_window) peak_time_llt_monpulse = Sim_LLTs21_this_event[n];
        }
        for (size_t p = 0; p < LLT21_Times.size(); ++p) {
          if (LLT21_Times[p] >= peak_time-acceptance_window && LLT21_Times[p] <= peak_time+acceptance_window) peak_time_llt_ptb = LLT21_Times[p];
        }*/

        int best_dist_monpulse = INT_MAX;
        for (size_t n = 0; n < Sim_LLTs21_this_event.size(); ++n) {
          int d = std::abs(Sim_LLTs21_this_event[n] - peak_time);
          if (d <= acceptance_window && d < best_dist_monpulse) {
            peak_time_llt_monpulse = Sim_LLTs21_this_event[n];
            best_dist_monpulse = d;
          }
        }

        int best_dist_ptb = INT_MAX;
        for (size_t p = 0; p < LLT21_Times.size(); ++p) {
          int d = std::abs(LLT21_Times[p] - peak_time);
          if (d <= acceptance_window && d < best_dist_ptb) {
            peak_time_llt_ptb = LLT21_Times[p];
            best_dist_ptb = d;
          }
        }

        if (peak_time_llt_monpulse >= 0 && peak_time_llt_ptb >= 0) peak_llt_diff.push_back(peak_time_llt_monpulse - peak_time_llt_ptb);
        else peak_llt_diff.push_back(101);
      }
    }

    //char title[100];
    //sprintf(title, "first_llt_diff_event%d.png", int(i));
    //makeHist(first_llt_diff_this_event, title, 151, -500, 500, "Time difference between first LLT 21 for each MonPulse");

    // Value of monpulse when LLT tag
    // loop through LLT Times from MonPulse, get value of monpulse, and push_back
    for (int llt_time : Sim_LLTs21_this_event) mon_vals_monpulse_21.push_back(MonPulses_this_event[llt_time]);
    for (int llt_time : Sim_LLTs20_this_event) mon_vals_monpulse_20.push_back(MonPulses_this_event[llt_time]);
    for (int llt_time : Sim_LLTs18_this_event) mon_vals_monpulse_18.push_back(MonPulses_this_event[llt_time]);
    // loop through LLT Times from PTB, get value of monpulse, and push_back 
    //for (int llt_time : LLT21_Times) mon_vals_ptb_21.push_back(MonPulses_this_event[llt_time]);
    //for (int llt_time : LLT20_Times) mon_vals_ptb_20.push_back(MonPulses_this_event[llt_time]);
    //for (int llt_time : LLT18_Times) mon_vals_ptb_18.push_back(MonPulses_this_event[llt_time]);
    const int tolerance = PTB_TOLERANCE;
    for (int llt_time : LLT21_Times) {

        int lo = std::max(0, llt_time - tolerance);
        int hi = std::min((int)MonPulses_this_event.size() - 1, llt_time + tolerance);

        if (llt_time >= 0 && llt_time < (int)MonPulses_this_event.size()) {
            mon_vals_ptb_21.push_back(MonPulses_this_event[llt_time]);

            int max_val = MonPulses_this_event[lo];
            for (int i = lo + 1; i <= hi; ++i) {
                max_val = std::max(max_val, MonPulses_this_event[i]);
            }
            mon_vals_ptb_21_tolerance.push_back(max_val);

        }
    }
    for (int llt_time : LLT20_Times) {

        int lo = std::max(0, llt_time - tolerance);
        int hi = std::min((int)MonPulses_this_event.size() - 1, llt_time + tolerance);

        if (llt_time >= 0 && llt_time < (int)MonPulses_this_event.size()) {
            mon_vals_ptb_20.push_back(MonPulses_this_event[llt_time]);

            int max_val = MonPulses_this_event[lo];
            for (int i = lo + 1; i <= hi; ++i) {
                max_val = std::max(max_val, MonPulses_this_event[i]);
            }
            mon_vals_ptb_20_tolerance.push_back(max_val);

        }
    }
    for (int llt_time : LLT18_Times) {

        int lo = std::max(0, llt_time - tolerance);
        int hi = std::min((int)MonPulses_this_event.size() - 1, llt_time + tolerance);

        if (llt_time >= 0 && llt_time < (int)MonPulses_this_event.size()) {
            mon_vals_ptb_18.push_back(MonPulses_this_event[llt_time]);

            int max_val = MonPulses_this_event[lo];
            for (int i = lo + 1; i <= hi; ++i) {
                max_val = std::max(max_val, MonPulses_this_event[i]);
            }
            mon_vals_ptb_18_tolerance.push_back(max_val);

        }
    }


    LLTs_per_event.push_back(PTB_LLTs.size());
    SimLLTs_per_event.push_back(All_LLTs_this_event);

    if (PTB_LLTs.size() == 0) events_0.push_back(1);
    else events_0.push_back(0);
    if (All_LLTs_this_event == 0) Simevents_0.push_back(1);
    else Simevents_0.push_back(0);

    LLTs_18_per_event.push_back(LLTs_18_this_event);
    SimLLTs_18_per_event.push_back(SimLLTs_18_this_event);

    LLTs_20_per_event.push_back(LLTs_20_this_event);
    SimLLTs_20_per_event.push_back(SimLLTs_20_this_event);

    LLTs_21_per_event.push_back(LLTs_21_this_event);
    SimLLTs_21_per_event.push_back(SimLLTs_21_this_event);

    if (!PTB_LLTTimes.empty()) {
      max_time = *std::max_element(PTB_LLTTimes.begin(), PTB_LLTTimes.end());
      min_time = *std::min_element(PTB_LLTTimes.begin(), PTB_LLTTimes.end());
      pull_window_lengths.push_back(
          static_cast<float>((max_time - min_time) / 1000000.0));
    }

  }

  f->Close();
  
  auto h16 = makeHistSave(peak_height_missing, "Peak_Height_Missing.png", 65, 0, 65, "Height of MonPulse when LLTs Missing from PTB");
  auto h17 = makeHistSave(peak_height_not_missing, "Peak_Height_Not_Missing.png", 65, 0, 65, "Height of MonPulse when LLTs Not Missing from PTB");
  makeAreaNormComparisonHist(h16, h17, "Missing LLTs from PTB", "Not Missing LLTs from PTB", "Peak_Height_Not_Missing_vs_Missing.png", "", "Height of Peak (# of pairs above threshold)", "MonPulses");

  makeHist(missing_llts_ptb, "Missing_LLTs_from_PTB_by_MonPulse.png", 2, 0, 2, "MonPulse Missing LLTs Entirely from PTB");
  makeHist(missing_llt18_ptb, "Missing_LLT_18_from_PTB_by_MonPulse.png", 2, 0, 2, "MonPulse Missing LLT 18 Entirely from PTB");
  makeHist(missing_llt20_ptb, "Missing_LLT_20_from_PTB_by_MonPulse.png", 2, 0, 2, "MonPulse Missing LLT 20 Entirely from PTB");
  makeHist(missing_llt21_ptb, "Missing_LLT_21_from_PTB_by_MonPulse.png", 2, 0, 2, "MonPulse Missing LLT 21 Entirely from PTB");

  makeHist(missing_flash_this_event_monpulses, "Missing_Flash_Triggers_MonPulses_1.png", 2, 0, 2);
  makeHist(missing_flash_this_event_ptb, "Missing_Flash_Triggers_PTB_1.png", 2, 0, 2);
  makeHist(missing_flash_this_monpulse, "Missing_Flash_Triggers_each_MonPulse.png", 2, 0, 2);
  makeHist(missing_event_this_monpulse, "Missing_Event_Triggers_each_MonPulse.png", 2, 0, 2);
  makeHist(missing_18_this_monpulse, "Missing_18_Triggers_each_MonPulse.png", 2, 0, 2);

  makeHist(time_bn_LLT18, "Time_bn_LLT18_1.png", 25, 0, 500, "Time difference between LLT18 for MonPulses");
  makeHist(time_bn_LLT20, "Time_bn_LLT20_1.png", 25, 0, 500, "Time difference between LLT20 for MonPulses");
  makeHist(time_bn_LLT21, "Time_bn_LLT21_1.png", 25, 0, 500, "Time difference between LLT21 for MonPulses");
  makeHist(time_bn_LLT18_PTB, "Time_bn_LLT18_PTB_1.png", 25, 0, 500, "Time difference between LLT18 for PTB");
  makeHist(time_bn_LLT20_PTB, "Time_bn_LLT20_PTB_1.png", 25, 0, 500, "Time difference between LLT20 for PTB");
  makeHist(time_bn_LLT21_PTB, "Time_bn_LLT21_PTB_1.png", 25, 0, 500, "Time difference between LLT21 for PTB");
  
  makeHist(first_llt_diff, "first_llt_diff.png", 151, -500, 500, "Time difference between first LLT 21 for each MonPulse");
  makeHist(first_llt_diff, "first_llt_diff_with_overflow.png", 151, -500, 500, "Time difference between first LLT 21 for each MonPulse");
  
  makeHist(peak_llt_diff, "peak_llt_diff_with_overflow.png", 205, -102, 102, "Time difference between LLT 21 near peak for each MonPulse");
  makeHist(peak_llt_diff, "peak_llt_diff.png", 501, -2000, 2000, "Time difference between LLT 21 near peak for each MonPulse");
  
  
  makeHist(mon_vals_monpulse_21_direct, "mon_vals_monpulse_21_direct.png", 20, 0, 20, "Value of MonPulse at LLT 21 for MonPulse-extracted LLTs (individual)");
  makeHist(mon_vals_monpulse_20_direct, "mon_vals_monpulse_20_direct.png", 20, 0, 20, "Value of MonPulse at LLT 20 for MonPulse-extracted LLTs (individual)");
  makeHist(mon_vals_monpulse_18_direct, "mon_vals_monpulse_18_direct.png", 20, 0, 20, "Value of MonPulse at LLT 18 for MonPulse-extracted LLTs (individual)");
  

  auto h9 = makeHistSave(passed_events, "trigger_rate.png", 2, 0, 2, "Passed Trigger", "events");
  auto h10 = makeHistSave(passed_events_PTB, "trigger_rate_PTB.png", 2, 0, 2, "Passed Trigger", "events");
  makeComparisonHist(h9, h10, "From MonPulse", "From PTB", "comparison_passed_Trigger.png", "Passed Trigger", "Passed Trigger", "Events");
  /*
  makeHist(deltaT_A, "peaks_A.png", 500, 500, 1000, "MonPulse Peak time", "number of monpulses");
  makeHist(deltaT_B, "first20_B_with_default.png", 500, 500, 1000, "First LLT20 Time", "number of monpulses");
  makeHist(deltaT_B_21, "first21_B_with_default.png", 500, 500, 1000, "First LLT21 Time", "number of monpulses");
  */
  makeHist(mon_vals_monpulse_21, "mon_vals_monpulse_21.png", 20, 0, 20, "Value of MonPulse at LLT 21 for MonPulse-extracted LLTs");
  makeHist(mon_vals_monpulse_20, "mon_vals_monpulse_20.png", 20, 0, 20, "Value of MonPulse at LLT 20 for MonPulse-extracted LLTs");
  makeHist(mon_vals_monpulse_18, "mon_vals_monpulse_18.png", 20, 0, 20, "Value of MonPulse at LLT 18 for MonPulse-extracted LLTs");
  makeHist(mon_vals_ptb_21, "mon_vals_ptb_21.png", 20, 0, 20, "Value of MonPulse at LLT 21 for PTB-extracted LLTs");
  makeHist(mon_vals_ptb_20, "mon_vals_ptb_20.png", 20, 0, 20, "Value of MonPulse at LLT 20 for PTB-extracted LLTs");
  makeHist(mon_vals_ptb_18, "mon_vals_ptb_18.png", 20, 0, 20, "Value of MonPulse at LLT 18 for PTB-extracted LLTs");
  makeHist(mon_vals_ptb_21_tolerance, "mon_vals_ptb_21_tolerance.png", 20, 0, 20, "Max Value of MonPulse in 120 tick window around PTB-extracted LLT 21s");
  makeHist(mon_vals_ptb_20_tolerance, "mon_vals_ptb_20_tolerance.png", 20, 0, 20, "Max Value of MonPulse in 120 tick window around PTB-extracted LLT 20s");
  makeHist(mon_vals_ptb_18_tolerance, "mon_vals_ptb_18_tolerance.png", 20, 0, 20, "Max Value of MonPulse in 120 tick window around PTB-extracted LLT 18s");
  

  auto *h1 = makeHistSave(LLTs_per_event, "num_LLTs_per_event_ptb.png", 50, 0, 500, "Number of LLTs from PTB", "Events");
  auto *h2 = makeHistSave(SimLLTs_per_event, "num_LLTs_per_event_monpulses.png", 50, 0, 500, "Number of LLTs from MonPulses", "Events");
  makeComparisonHist(h1, h2, "From PTB", "From MonPulse", "comparison_Number_of_LLTs_per_Event.png", "Number of LLTs per Event", "Number of LLTs", "Events");
  
  auto *h3 = makeHistSave(LLTs_18_per_event, "num_LLTs_18_per_event_ptb.png", 40, 0, 200, "Number of 18 LLTs from PTB", "Events");
  auto *h4 = makeHistSave(SimLLTs_18_per_event, "num_LLTs_18_per_event_monpulses.png", 40, 0, 200, "Number of 18 LLTs from MonPulses", "Events");
  makeComparisonHist(h3, h4, "From PTB", "From MonPulse", "comparison_Number_of_LLTs_18_per_Event.png", "Number of 18 LLTs per Event", "Number of Event LLTs", "Events");

  auto *h5 = makeHistSave(LLTs_20_per_event, "num_LLTs_20_per_event_ptb.png", 40, 0, 200, "Number of Event LLTs from PTB", "Events");
  auto *h6 = makeHistSave(SimLLTs_20_per_event, "num_LLTs_20_per_event_monpulses.png", 40, 0, 200, "Number of Event LLTs from MonPulses", "Events");
  makeComparisonHist(h5, h6, "From PTB", "From MonPulse", "comparison_Number_of_LLTs_20_per_Event.png", "Number of Event LLTs per Event", "Number of Event LLTs", "Events");

  auto *h7 = makeHistSave(LLTs_21_per_event, "num_LLTs_21_per_event_ptb.png", 40, 0, 200, "Number of Flash LLTs from PTB", "Events");
  auto *h8 = makeHistSave(SimLLTs_21_per_event, "num_LLTs_21_per_event_monpulses.png", 40, 0, 200, "Number of Flash LLTs from MonPulses", "Events");
  makeComparisonHist(h7, h8, "From PTB", "From MonPulse", "comparison_Number_of_LLTs_21_per_Event.png", "Number of Flash LLTs per Event", "Number of Event LLTs", "Events");

  auto *h11 = makeHistSave(events_0, "num_no_LLTs_per_event_ptb.png", 2, 0, 2, "Has 0 LLTs", "Events");
  auto *h12 = makeHistSave(Simevents_0, "num_no_LLTs_per_event_monpulses.png", 2, 0, 2, "Has 0 LLTs", "Events");
  makeComparisonHist(h11, h12, "From PTB", "From MonPulse", "comparison_Events_with_0_LLTs.png", "Number of Events with 0 LLTs", "Has 0 LLTs", "Events");

  //auto *h13 = makeFloatHist(percent_missing, "Fraction of MonPules Missing Flash Triggers by Event", "Fraction of Missing Flash Triggers", "Events", 20, 0, 1, "missing_Flash_trigger.png");
  //auto *h14 = makeFloatHist(percent_missing_merged, "Fraction of Missing Flash Triggers for Merged MonPulses by Event", "Fraction of Missing Flash Triggers", "Events", 20, 0, 1, "missing_merged_Flash_trigger.png");

  auto *h15 = makeFloatHist(pull_window_lengths, "Pull Windows for all Events", "Pull Window (ms)", "Events", 30, 2000, 5000, "pull_window_lengths.png");
  
}

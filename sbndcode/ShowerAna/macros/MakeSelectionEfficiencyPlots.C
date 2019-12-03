#include "TH1.h"
#include "TTree.h"
#include <vector>
#include "TGraph.h"
#include "TH2.h"
#include "TGraph2D.h"
#include "TMath.h"

void DoCutFinding(TH1D* signalhist, TH1D* backgroundhist, TString units, TString axisTitle);
void MakeEfficiencyPlots(TH1F* postcut, TH1F* precut, TH1F* extra, TString postcutname, TString precutname, TString extraname, TString title, TString eff_name);
void MakeResolutionGraphs(std::vector<float>& vals, std::vector<float>& energyval, float min, float max, float num_steps, float width,TString yTitle, TString Metric, TString name,float nbins);

void MakeSelectionEfficiencyPlots(){

  
  //Config Parameters
  /////////////////////////////////////////////////
  bool fVisibleEnergyCut = true;
  bool fVertexEnergyCut  = true; 
  bool fVertexConversionCut = true;
  bool fdEdxCut = true;
  bool fTrackLengthCut =  true; 
  bool fReconstructionEfficiency = true;
  bool fApplyPreviousCuts = false;
  bool fTrackValidation = true;

  float fECut = 200;
  float fConversionGapCut = 2.45;
  float fVertexCut = -1;
  float fdEdxCutVal = 2.9;
  float fMaxLengthCut = 30;

  double signalscale = 1;
  double backgroundscale = 1;

  /////////////////////////////////////////////////

  //Root Dictionaries cos I'm lazy.
  gROOT->ProcessLine("#include <vector>");
  gInterpreter->GenerateDictionary("vector<float>", "vector");
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gStyle->SetPalette(1);

  float Energies[44] = {0,50,75,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500};

  //Get the histogram for the signal
  TFile *f_eff = new TFile("ShowerValidation_osc.root");
  TDirectory* dir_eff = f_eff->GetDirectory("eff");
  TTree *tree_eff = (TTree*)dir_eff->Get("RecoEffMetricTree");  

  //Get the histogram for the background
  TFile *f_bk = new TFile("new/showervalidationGraphs_test.root");
  TDirectory* dir_bk = f_bk->GetDirectory("eff");
  TTree *tree_bk = (TTree*)dir_bk->Get("RecoEffMetricTree");  

  float trueShower_num_branch_eff;
  tree_eff->SetBranchAddress("trueShower_num",&trueShower_num_branch_eff);

  float trueShower_num_branch_bk;
  tree_bk->SetBranchAddress("trueShower_num",&trueShower_num_branch_bk);
  
  float numtrueVtx_branch_bk;
  tree_bk->SetBranchAddress("numtrueVtx_branch",&numtrueVtx_branch_bk);

  float numtrueVtx_branch_eff;
  tree_eff->SetBranchAddress("numtrueVtx_branch",&numtrueVtx_branch_eff);
  
  TFile* file = new TFile("CutFile.root", "RECREATE");
  
  gDirectory->mkdir("ReconstructionEfficiency");
  gDirectory->mkdir("VisibleEnergy");
  gDirectory->mkdir("VertexConversionCut");
  gDirectory->mkdir("dEdxCut");
  gDirectory->mkdir("TrackLengthCut");
  gDirectory->mkdir("TrackValidation");

  //#################################
  //### Reconstruction Efficiency ###
  //#################################
  
  //Very bad analysis every vertex comes froma neutrino hence it is reconstructed.
  //Define reconstructed if there exist a reco pfp particle. 
  if(fReconstructionEfficiency){

    file->cd("ReconstructionEfficiency");

    //Energy cut-off to see how many vertices with true energy above the cut off  are correctly reconstructed reconstructed 
    //    float Energy[35] = {0,100,200,300,400,500,600,800,900,1100,1200,1300,1400,1600,1700,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600};

    //Ratio of the integral above the energy cut off of the vertex that are correctly reconstructed and the total vertices above the cut off 
    //    float efficiency[34];
    //float background[34];

    //Histograms to show the reco and true energy.
    TH1F* reco_neutrino_bk_hist = new TH1F("PFP particle number bk","PFP particle number",30,0,3000);
    TH1F* true_neutrino_bk_hist = new TH1F("True neutrino number bk","True neutrino number",30,0,3000);
    TH1F* extra_reco_neutrino_bk_hist = new TH1F("PFP particle number bk extra","PFP particle number extra",30,0,3000);
    
    //Get the background branches 
    std::vector<float>* nu_E_branch_bk = 0;
    tree_bk->SetBranchAddress("nu_E",&nu_E_branch_bk);

    std::vector<float>* nu_E_numtrue_branch_bk =0;
    tree_bk->SetBranchAddress("nu_E_numtrue",&nu_E_numtrue_branch_bk);


    //Loop over the events 
    Long64_t nentries_bk = tree_bk->GetEntries();
    for (Long64_t evt = 0; evt < nentries_bk; evt++) {
      tree_bk->GetEntry(evt);

      std::vector<float> nu_E_branch_bk_temp;

      //Loop over the true neutrinos.
      for(auto const& neutrino_E: *nu_E_numtrue_branch_bk){
	true_neutrino_bk_hist->Fill(neutrino_E*1000);
      }
   
      //Loop over the reco neutrinos.
      for(auto const& neutrino_E: *nu_E_branch_bk){
	//See if the neutrino has been matched;
	if(std::find(nu_E_branch_bk_temp.begin(),nu_E_branch_bk_temp.end(),neutrino_E) != nu_E_branch_bk_temp.end()){
	  extra_reco_neutrino_bk_hist->Fill(neutrino_E*1000);
	}
	else{
	  nu_E_branch_bk_temp.push_back(neutrino_E);
	  reco_neutrino_bk_hist->Fill(neutrino_E*1000);
	}
      }
    }
    
    //Histograms to show the reco and true energy.
    TH1F* reco_neutrino_eff_hist = new TH1F("PFP particle number sig","PFP particle number",28,200,3000);
    TH1F* true_neutrino_eff_hist = new TH1F("True neutrino number sig","True neutrino number",28,200,3000);
    TH1F* extra_reco_neutrino_eff_hist = new TH1F("PFP particle number eff extra","PFP particle number extra",28,200,3000);

    //Get the background branches 
    std::vector<float>* nu_E_branch_eff = 0;
    tree_eff->SetBranchAddress("nu_E",&nu_E_branch_eff);

    std::vector<float>* nu_E_numtrue_branch_eff =0;
    tree_eff->SetBranchAddress("nu_E_numtrue",&nu_E_numtrue_branch_eff);

    //Loop over the events 
    Long64_t nentries_eff = tree_eff->GetEntries();
    for (Long64_t evt = 0; evt < nentries_eff; evt++) {
      tree_eff->GetEntry(evt);
      
      std::vector<float> nu_E_branch_eff_temp;

      //Loop over the true neutrinos.
      for(auto const& neutrino_E: *nu_E_numtrue_branch_eff){
	true_neutrino_eff_hist->Fill(neutrino_E*1000);
      }
      //Loop over the reco neutrinos.
      for(auto const& neutrino_E: *nu_E_branch_eff){

	//See if the neutrino has been matched
	if(std::find(nu_E_branch_eff_temp.begin(),nu_E_branch_eff_temp.end(),neutrino_E) != nu_E_branch_eff_temp.end()){
	  extra_reco_neutrino_eff_hist->Fill(neutrino_E*1000);
	}
	else{
	  nu_E_branch_eff_temp.push_back(neutrino_E);
	  reco_neutrino_eff_hist->Fill(neutrino_E*1000);
	}
	  
      }
    }

    float maxSignal = true_neutrino_eff_hist->GetBinContent(true_neutrino_eff_hist->GetMaximumBin());
    float maxBackground = true_neutrino_bk_hist->GetBinContent(true_neutrino_bk_hist->GetMaximumBin());
    
    //Make the pretty pictures 
    if(TEfficiency::CheckConsistency(*reco_neutrino_eff_hist,*true_neutrino_eff_hist)){
      //Make the eff plot.
      TEfficiency* reco_eff = new TEfficiency(*reco_neutrino_eff_hist,*true_neutrino_eff_hist);
      reco_eff->SetTitle("Efficiency;Neutrino Energy (MeV);Efficiency");
      reco_eff->SetName("signal eff");
      reco_eff->Write();

      reco_neutrino_eff_hist->Scale(1./maxSignal);
      true_neutrino_eff_hist->Scale(1./maxSignal);
      extra_reco_neutrino_eff_hist->Scale(1./maxSignal);
      
      reco_neutrino_bk_hist->Scale(1./maxBackground);
      true_neutrino_bk_hist->Scale(1./maxBackground);
      extra_reco_neutrino_bk_hist->Scale(1./maxBackground);


      reco_neutrino_eff_hist->SetLineColor(kBlue);
      reco_neutrino_eff_hist->SetFillStyle(3003);
      reco_neutrino_eff_hist->SetFillColor(6);
      reco_neutrino_eff_hist->SetTitle("Reco Number of neutrinos");
      reco_neutrino_eff_hist->Write();
      true_neutrino_eff_hist->SetLineColor(kRed);
      true_neutrino_eff_hist->SetFillStyle(3003);
      true_neutrino_eff_hist->SetFillColor(42);
      true_neutrino_eff_hist->SetTitle("True Number of neutrinos");
      true_neutrino_eff_hist->Write();

      auto sigbk_cut_canvas = new TCanvas("signal_canvas", "signal_canvas", 600, 400);
      true_neutrino_eff_hist->Draw("HIST SAME");
      sigbk_cut_canvas->cd();
      reco_eff->Draw("same");
      reco_neutrino_eff_hist->Draw("HIST SAME");
      sigbk_cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
      sigbk_cut_canvas->Write();
    
      
      reco_neutrino_eff_hist->Scale(maxSignal);
      true_neutrino_eff_hist->Scale(maxSignal);
      extra_reco_neutrino_eff_hist->Scale(maxSignal);
      
      reco_neutrino_bk_hist->Scale(maxBackground);
      true_neutrino_bk_hist->Scale(maxBackground);
      extra_reco_neutrino_bk_hist->Scale(maxBackground);
    }

    //Make the pretty pictures 
    if(TEfficiency::CheckConsistency(*extra_reco_neutrino_eff_hist,*true_neutrino_eff_hist)){
      //Make the eff plot.
      TEfficiency* reco_eff = new TEfficiency(*extra_reco_neutrino_eff_hist,*true_neutrino_eff_hist);
      reco_eff->SetTitle("Efficiency;Neutrino Energy (MeV);Efficiency");
      reco_eff->SetName("signal extra");
      reco_eff->Write();

      reco_neutrino_eff_hist->Scale(1./maxSignal);
      true_neutrino_eff_hist->Scale(1./maxSignal);
      extra_reco_neutrino_eff_hist->Scale(1./maxSignal);
      
      reco_neutrino_bk_hist->Scale(1./maxBackground);
      true_neutrino_bk_hist->Scale(1./maxBackground);
      extra_reco_neutrino_bk_hist->Scale(1./maxBackground);


      extra_reco_neutrino_eff_hist->SetLineColor(kBlue);
      extra_reco_neutrino_eff_hist->SetFillStyle(3003);
      extra_reco_neutrino_eff_hist->SetFillColor(6);
      extra_reco_neutrino_eff_hist->SetTitle("% of Neutrinos with more that one pfp");
      extra_reco_neutrino_eff_hist->Write();
      true_neutrino_eff_hist->SetLineColor(kRed);
      true_neutrino_eff_hist->SetFillStyle(3003);
      true_neutrino_eff_hist->SetFillColor(42);
      true_neutrino_eff_hist->SetTitle("Ture Number of neutrinos");
      true_neutrino_eff_hist->Write();

      auto sigbk_cut_canvas = new TCanvas("signal_canvas extra", "signal_canvas extra", 600, 400);

      sigbk_cut_canvas->cd();
      true_neutrino_eff_hist->Draw("HIST SAME");
      extra_reco_neutrino_eff_hist->Draw("HIST SAME");
      reco_eff->Draw("same");

      sigbk_cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
      sigbk_cut_canvas->Write();

      reco_neutrino_eff_hist->Scale(maxSignal);
      true_neutrino_eff_hist->Scale(maxSignal);
      extra_reco_neutrino_eff_hist->Scale(maxSignal);
      
      reco_neutrino_bk_hist->Scale(maxBackground);
      true_neutrino_bk_hist->Scale(maxBackground);
      extra_reco_neutrino_bk_hist->Scale(maxBackground);
      


    }


    //Make the pretty pictures 
    if(TEfficiency::CheckConsistency(*reco_neutrino_bk_hist,*true_neutrino_bk_hist)){
    
      //Make the bk plot.
      TEfficiency* reco_bk = new TEfficiency(*reco_neutrino_bk_hist,*true_neutrino_bk_hist);

      reco_bk->SetTitle("Efficiency;Neutrino Energy (MeV);Efficiency");
      reco_bk->SetName("background eff");
      reco_bk->Write();

      reco_neutrino_eff_hist->Scale(1./maxSignal);
      true_neutrino_eff_hist->Scale(1./maxSignal);
      extra_reco_neutrino_eff_hist->Scale(1./maxSignal);
      
      reco_neutrino_bk_hist->Scale(1./maxBackground);
      true_neutrino_bk_hist->Scale(1./maxBackground);
      extra_reco_neutrino_bk_hist->Scale(1./maxBackground);

      reco_neutrino_bk_hist->SetLineColor(kBlue);
      reco_neutrino_bk_hist->SetFillStyle(3003);
      reco_neutrino_bk_hist->SetFillColor(6);
      reco_neutrino_bk_hist->SetTitle("Reco Number of neutrinos");
      reco_neutrino_bk_hist->Write();
      true_neutrino_bk_hist->SetLineColor(kRed);
      true_neutrino_bk_hist->SetFillStyle(3003);
      true_neutrino_bk_hist->SetFillColor(42);
      true_neutrino_bk_hist->SetTitle("True Number of neutrinos");
      true_neutrino_bk_hist->Write();

  
      auto sigbk_cut_canvas = new TCanvas("background_canvas", "background_canvas", 600, 400);
      sigbk_cut_canvas->cd();
      true_neutrino_bk_hist->Draw("HIST SAME");
      reco_bk->Draw("same");
      reco_neutrino_bk_hist->Draw("HIST SAME");
      sigbk_cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
      sigbk_cut_canvas->Write();

      reco_neutrino_eff_hist->Scale(maxSignal);
      true_neutrino_eff_hist->Scale(maxSignal);
      extra_reco_neutrino_eff_hist->Scale(maxSignal);
      
      reco_neutrino_bk_hist->Scale(maxBackground);
      true_neutrino_bk_hist->Scale(maxBackground);
      extra_reco_neutrino_bk_hist->Scale(maxBackground);

    }


    //Make the pretty pictures 
    if(TEfficiency::CheckConsistency(*extra_reco_neutrino_bk_hist,*true_neutrino_bk_hist)){
    
      //Make the bk plot.
      TEfficiency* reco_bk = new TEfficiency(*extra_reco_neutrino_bk_hist,*true_neutrino_bk_hist);
      reco_bk->SetTitle("Efficiency;Neutrino Energy (MeV);% of neutrinos reconstructed more than twice");
      reco_bk->SetName("background extra");
      reco_bk->Write();

      reco_neutrino_eff_hist->Scale(1./maxSignal);
      true_neutrino_eff_hist->Scale(1./maxSignal);
      extra_reco_neutrino_eff_hist->Scale(1./maxSignal);
      
      reco_neutrino_bk_hist->Scale(1./maxBackground);
      true_neutrino_bk_hist->Scale(1./maxBackground);
      extra_reco_neutrino_bk_hist->Scale(1./maxBackground);


      extra_reco_neutrino_bk_hist->SetLineColor(kBlue);
      extra_reco_neutrino_bk_hist->SetFillStyle(3003);
      extra_reco_neutrino_bk_hist->SetFillColor(6);
      extra_reco_neutrino_bk_hist->SetTitle("% of Neutrinos with more that one pfp");
      extra_reco_neutrino_bk_hist->Write();
      true_neutrino_bk_hist->SetLineColor(kRed);
      true_neutrino_bk_hist->SetFillStyle(3003);
      true_neutrino_bk_hist->SetFillColor(42);
      true_neutrino_bk_hist->SetTitle("True Number of neutrinos");
      true_neutrino_bk_hist->Write();
  
      auto sigbk_cut_canvas = new TCanvas("background_canvas extra", "background_canvas extra", 600, 400);
      sigbk_cut_canvas->cd();
      true_neutrino_bk_hist->Draw("HIST SAME");
      reco_bk->Draw("same");
      extra_reco_neutrino_bk_hist->Draw("HIST SAME");
      sigbk_cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
      sigbk_cut_canvas->Write();

      reco_neutrino_eff_hist->Scale(maxSignal);
      true_neutrino_eff_hist->Scale(maxSignal);
      extra_reco_neutrino_eff_hist->Scale(maxSignal);
      
      reco_neutrino_bk_hist->Scale(maxBackground);
      true_neutrino_bk_hist->Scale(maxBackground);
      extra_reco_neutrino_bk_hist->Scale(maxBackground);


    }
    


  }



  //##########################
  //### Visible Energy Cut ###
  //##########################

    if(fVisibleEnergyCut){

    std::cout<< "###############################" << std::endl;
    std::cout << "##### Visible Energy Cut #####" <<std::endl;
    std::cout << "###############################" << std::endl;

    file->cd("VisibleEnergy");

    float Energies[44] = {0,50,75,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500};

    float efficiency[44];
    float signal[44];
    float background[44];
    float backgroundRej[44];
    float purity[44];
    float effbk[44];
    float effpur[44];

    std::vector<std::vector<float> >* shower_energy_eff = 0;
    tree_eff->SetBranchAddress("shower_energy",&shower_energy_eff);
    Long64_t nentries_eff = tree_eff->GetEntries();
    //Loop over the energies and calculate the effiecny 
    for(int i=0; i<44; ++i){
      float n_neutrinos = 0;
      float Energy = Energies[i];
      float eff =0;
      //Loop over the events.
      for (Long64_t evt = 0; evt < nentries_eff; evt++) {
	tree_eff->GetEntry(evt);
	n_neutrinos += numtrueVtx_branch_eff;
	//Loop over the neutrinos in the event 
	for(auto const& neutrino: *shower_energy_eff){ 
	  int numshowers = 0;
	  //Loop over the showers associated to the neutrino
	  for(auto const& shower_e: neutrino){
	    //Only keep showers above a specific energy.
	    if(shower_e > Energy){
	      ++numshowers; 
	    }
	  } 
	  //Only keep neutrinos with shower in them.
	  if(numshowers == 1){++eff;}
	}
      }
      efficiency[i] =  eff/n_neutrinos;
      signal[i]     = eff*signalscale;
    }

    float max_effbk = -999;
    int max_i = -999;

    //Repeat for the background events
    std::vector<std::vector<float> >* shower_energy_bk =0;
    tree_bk->SetBranchAddress("shower_energy",&shower_energy_bk);
    Long64_t nentries_bk = tree_bk->GetEntries();
    //Loop over the energy cuts.
    for(int i=0; i<44; ++i){
      float n_neutrinos = 0;
      float Energy = Energies[i];
      float bk=0;
      //Loop over the events 
      for (Long64_t evt = 0; evt < nentries_bk; evt++) {
	tree_bk->GetEntry(evt);
	n_neutrinos += numtrueVtx_branch_bk;
	//Loop over the neutrinos in the event
	for(auto const& neutrino: *shower_energy_bk){
	  int numshowers = 0;
	  //Loop over the showers.
	  for(auto const& shower_e: neutrino){
	    //Only keep showers above the energy in question
	    if(shower_e > Energy){
	      ++numshowers; 
	    }
	  } 
	  //Only keep the neutrinos with one shower. 
	  if(numshowers == 1){++bk;}
	}
      }
      background[i] = bk*backgroundscale;
      backgroundRej[i] = 1-bk/n_neutrinos;
      effbk[i] = backgroundRej[i]*efficiency[i];
      
      //Identify the max values.
      if(effbk[i] > max_effbk){
	max_effbk = effbk[i];
	max_i = i;
      }
    }

    //Calculate the purity:
    int max_i_pure   = -999;
    float max_effpur = -999;
    for(int i=0; i<44; ++i){
      purity[i] = signal[i]/(signal[i]+background[i]);
      effpur[i] = purity[i]*efficiency[i];
      //Get the max value;
      if(max_effpur < effpur[i]){
	max_effpur = effpur[i];
	max_i_pure = i;
      }
    }
    
    //Make eff plot.
    TGraph* efficiencyPlot = new TGraph(44,Energies,efficiency);  
    efficiencyPlot->Draw("a4");
    efficiencyPlot->SetMarkerStyle(8);
    efficiencyPlot->SetMarkerColor(kBlue);
    efficiencyPlot->SetTitle("Efficiency");

    //Make bk plot.
    TGraph* backgroundPlot = new TGraph(44,Energies,backgroundRej); 
    backgroundPlot->Draw("a4");
    backgroundPlot->SetMarkerStyle(8);
    backgroundPlot->SetMarkerColor(kRed);
    backgroundPlot->SetTitle("Bakcground Rejection");

    //Make purity plot.
    TGraph* purityPlot = new TGraph(44,Energies,purity); 
    purityPlot->Draw("a4");
    purityPlot->SetMarkerStyle(8);
    purityPlot->SetMarkerColor(kRed);
    purityPlot->SetTitle("Purity");

    //Make eff*bk plot.
    TGraph* effbackgroundPlot = new TGraph(44,Energies,effbk); 
    effbackgroundPlot->Draw("a4");
    effbackgroundPlot->SetMarkerStyle(8);
    effbackgroundPlot->SetMarkerColor(kGreen);
    effbackgroundPlot->SetTitle("Eff*BkgdRej");

    //Make eff*purity plot.
    TGraph* effpurPlot = new TGraph(44,Energies,effpur); 
    effpurPlot->Draw("a4");
    effpurPlot->SetMarkerStyle(8);
    effpurPlot->SetMarkerColor(kGreen);
    effpurPlot->SetTitle("Eff*BkgdRej");


    //Make the background rej plot.
    auto sigbk_canvas = new TCanvas("sigbk_energycanvas", "sigbk_energy_canvas", 600, 400);
    auto mg2 = new TMultiGraph("energy cut bkrej", "energy cut bkrej");
    mg2->Add(efficiencyPlot);
    mg2->Add(backgroundPlot);
    mg2->Add(effbackgroundPlot);
    mg2->GetYaxis()->SetRangeUser(0, 1);
    mg2->Draw("AP");
    mg2->Write();
    mg2->GetHistogram()->SetTitle(Form("Cut at %0.1f MeV with efficiency %0.2f and background rejection %0.2f", Energies[max_i], efficiency[max_i], backgroundRej[max_i]));
    mg2->GetXaxis()->SetTitle("Min Shower Energy Cut (MeV)");

    sigbk_canvas->BuildLegend();
    TLine* l1 = new TLine(Energies[max_i], 0, Energies[max_i], 1);
    l1->SetLineColor(kBlack);
    l1->Draw();
    sigbk_canvas->Write();


    //Make the purity rej plot.
    auto effpur_canvas = new TCanvas("purity_energy_canvas", "purity_energy_canvas", 600, 400);
    auto mg1 = new TMultiGraph("energy cut pur", "energy cut pur");
    mg1->Add(efficiencyPlot);
    mg1->Add(purityPlot);
    mg1->Add(effpurPlot);
    mg1->GetYaxis()->SetRangeUser(0, 1);
    mg1->Draw("AP");
    mg1->Write();
    mg1->GetHistogram()->SetTitle(Form("Cut at %0.1f MeV with efficiency %0.2f and purity %0.2f", Energies[max_i_pure], efficiency[max_i_pure], purity[max_i_pure]));
    mg1->GetXaxis()->SetTitle("Min Shower Energy Cut (MeV)");
    
    effpur_canvas->BuildLegend();
    TLine* l2 = new TLine(Energies[max_i], 0, Energies[max_i], 1);
    l2->SetLineColor(kBlack);
    l2->Draw();
    effpur_canvas->Write();
    

    std::cout << "max value is: " << max_effbk << " at energy: " << Energies[max_i] << " efficiency: " << efficiency[max_i] << " bkgrnd Rej: " << backgroundRej[max_i] << std::endl;
    std::cout << "max value is: " << max_effpur << " at energy: " << Energies[max_i_pure] << " efficiency: " << efficiency[max_i_pure] << " purity: " << purity[max_i_pure] << std::endl;

  }

  //#########################
  //### Vertex Energy Cut ###
  //#########################

  if(fVertexEnergyCut){

    file->cd("VertexConversionCut");

    //Energy cut-off to see how many vertices with true energy above the cut off  are correctly reconstructed reconstructed 
    float VertexEnergy[34] = {0,10,20,30,40,50,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600};

    //Ratio of the integral above the energy cut off of the vertex that are correctly reconstructed and the total vertices above the cut off 
    float efficiency[34];
    float background[34];

    //Histograms to show the reco and true energy.
    TH1F* vertex_trueKE_bk_hist = new TH1F("True Vertex Energy bk","True Vertex Energy ",30,0,600);
    TH1F* vertex_recoKE_bk_hist = new TH1F("Reco Verte Energy bk","Reco Vertex Energy ",30,0,600);
    
    //Get the background branches 
    std::vector<float>* vertex_recoKE_bk = 0;
    tree_bk->SetBranchAddress("vertex_recoK",&vertex_recoKE_bk);

    std::vector<float>* vertex_trueKE_bk = 0;
    tree_bk->SetBranchAddress("vertex_trueK",&vertex_trueKE_bk);

    std::vector<float>* vertex_reco_bk = 0;
    tree_bk->SetBranchAddress("vertex_reco",&vertex_reco_bk);
    
    std::vector<std::vector<float> >* shower_energy_bk = 0;
    tree_bk->SetBranchAddress("shower_energy",&shower_energy_bk);
	

    //Loop over the events 
    Long64_t nentries_bk = tree_bk->GetEntries();
    for (Long64_t evt = 0; evt < nentries_bk; evt++) {
      tree_bk->GetEntry(evt);

      //Loop over the neutrinos.
      for(int vertex=0; vertex<vertex_trueKE_bk->size(); ++vertex){

	if(fApplyPreviousCuts){
	  //The energy cut will be applied before this so we can apply this cut.
	  int numshowers = 0;
	  for(auto const& shower_e: shower_energy_bk->at(vertex)){
	    //See how many showers are above the cut.
	    if(shower_e > fECut){
	      ++numshowers;
	    }
	  }
	  //If the neutrino only has one shower then keep it.
	  if(numshowers != 1){continue;}
	}
	//Get the true energy.
	float const vertex_trueE = vertex_trueKE_bk->at(vertex);
	vertex_trueKE_bk_hist->Fill(vertex_trueE);
      }

      //Do the same for the reconstructed energy.
      for(int vertex=0; vertex<vertex_recoKE_bk->size(); ++vertex){

	if(fApplyPreviousCuts){
	  int numshowers = 0;
	  for(auto const& shower_e: shower_energy_bk->at(vertex)){
	    if(shower_e > fECut){
	      ++numshowers;
	    }
	  }
	  if(numshowers != 1){continue;}
	}

	float const vertex_recoE = vertex_recoKE_bk->at(vertex);
	vertex_recoKE_bk_hist->Fill(vertex_recoE);
      }


    }

    //Loop over the vertex energies to evaulate the energy.
    for(int i=0; i<34; ++i){
      float n_neutrinos = 0;
      float Energy = VertexEnergy[i];
      float bk=0;
      int reco_num=0;
      //Loop over the events.
      for (Long64_t evt = 0; evt < nentries_bk; evt++) {
	tree_bk->GetEntry(evt);

	//Loop over the neutrinos in the events.
	for(int vertex=0; vertex<vertex_trueKE_bk->size(); ++vertex){
	  
	  //Only keep events with one shower.
	  if(fApplyPreviousCuts){
	    int numshowers = 0;
	    for(auto const& shower_e: shower_energy_bk->at(vertex)){
	      if(shower_e > fECut){
		++numshowers;
	      }
	    }
	    if(numshowers != 1){continue;}
	  }

	  float const vertex_trueE = vertex_trueKE_bk->at(vertex);
	  float const vertex_reco = vertex_reco_bk->at(vertex);
	  if(vertex_trueE < Energy){continue;}
	  ++n_neutrinos;
	  //Check if the vertex was correctly reconstructed.
	  if(vertex_reco > -1){
	    ++bk;
	  }
	}
      }
      //Calculate the effiency 
      background[i] = bk/n_neutrinos;
    }

    //Signal Time 
    TH1F* vertex_trueKE_eff_hist = new TH1F("True Vertex Energy sig","True Vertex Energy",30,0,600);
    TH1F* vertex_recoKE_eff_hist = new TH1F("Reco Vertex Energy sig","Reco Vertex Energy",30,0,600);

    //Get the signal branches.
    std::vector<float>* vertex_recoKE_eff = 0;
    tree_eff->SetBranchAddress("vertex_recoK",&vertex_recoKE_eff);

    std::vector<float>* vertex_trueKE_eff = 0;
    tree_eff->SetBranchAddress("vertex_trueK",&vertex_trueKE_eff);

    std::vector<float>* vertex_reco_eff = 0;
    tree_eff->SetBranchAddress("vertex_reco",&vertex_reco_eff);
    
    std::vector<std::vector<float> >* shower_energy_eff = 0;
    tree_eff->SetBranchAddress("shower_energy",&shower_energy_eff);
	

    Long64_t nentries_eff = tree_eff->GetEntries();
    //Loop over the events
    for (Long64_t evt = 0; evt < nentries_eff; evt++) {
      tree_eff->GetEntry(evt);

      //Loop over the neutrinos
      for(int vertex=0; vertex<vertex_trueKE_eff->size(); ++vertex){

	//Only conider events with one shower above the ECut.
	if(fApplyPreviousCuts){
	  int numshowers = 0;
	  for(auto const& shower_e: shower_energy_eff->at(vertex)){
	    if(shower_e > fECut){
	      ++numshowers;
	    }
	  }
	  if(numshowers != 1){continue;}
	}

	//Get the true Hadronic activity.
	float const vertex_trueE = vertex_trueKE_eff->at(vertex);
	vertex_trueKE_eff_hist->Fill(vertex_trueE);
      }

      //Get the reco vertex.
      for(int vertex=0; vertex<vertex_recoKE_eff->size(); ++vertex){

	if(fApplyPreviousCuts){
	  int numshowers = 0;
	  for(auto const& shower_e: shower_energy_eff->at(vertex)){
	    if(shower_e > fECut){
	      ++numshowers;
	    }
	  }
	  if(numshowers != 1){continue;}
	}

	float const vertex_recoE = vertex_recoKE_eff->at(vertex);
	vertex_recoKE_eff_hist->Fill(vertex_recoE);
      }


    }

    //Loop over the vertex energies to evaulate at.
    for(int i=0; i<31; ++i){
      float n_neutrinos = 0;
      float Energy = VertexEnergy[i];
      float eff=0;
      int reco_num=0;
      //Loop over the events.
      for (Long64_t evt = 0; evt < nentries_eff; evt++) {
	tree_eff->GetEntry(evt);
	
	//Loop over the neutrinos
	for(int vertex=0; vertex<vertex_trueKE_eff->size(); ++vertex){

	  //Only Keep events with one showrr abov the energy cut.
	  if(fApplyPreviousCuts){
	    int numshowers = 0;
	    for(auto const& shower_e: shower_energy_eff->at(vertex)){
	      if(shower_e > fECut){
		++numshowers;
	      }
	    }
	    if(numshowers != 1){continue;}
	  }

	  float const vertex_trueE = vertex_trueKE_eff->at(vertex);
	  float const vertex_reco = vertex_reco_eff->at(vertex);
	  if(vertex_trueE < Energy){continue;}
	  ++n_neutrinos;
	  //Check the vertex was reconstructed correctly.
	  if(vertex_reco > -1){
	    ++eff;
	  }
	}
      }
      //calculate the efficiency.
      efficiency[i] = eff/n_neutrinos;
    }

    //Normalise the histograms so the max bin is at 1.
    double maxScale_eff = TMath::Max(vertex_trueKE_eff_hist->GetBinContent(vertex_trueKE_eff_hist->GetMaximumBin()),vertex_recoKE_eff_hist->GetBinContent(vertex_recoKE_eff_hist->GetMaximumBin()));
    vertex_trueKE_eff_hist->Scale(1/maxScale_eff);
    vertex_trueKE_eff_hist->Write();
    vertex_recoKE_eff_hist->Scale(1/maxScale_eff);
    vertex_recoKE_eff_hist->Write();

    vertex_trueKE_eff_hist->SetLineColor(kBlue);
    vertex_trueKE_eff_hist->SetFillColor(6);
    vertex_trueKE_eff_hist->SetFillStyle(3003);

    vertex_recoKE_eff_hist->SetLineColor(kGreen);
    vertex_recoKE_eff_hist->SetFillColor(42);
    vertex_recoKE_eff_hist->SetFillStyle(3003);

    //Draw the vertex efficiency canvas for the signal.
    auto eff_canvas = new TCanvas("eff_canvas", "eff_canvas", 600, 400);
    TGraph* effnalPlot = new TGraph(31,VertexEnergy,background); 
    effnalPlot->SetTitle("Reco Efficiency");
    effnalPlot->GetXaxis()->SetTitle("True Vertex Energy (MeV)");
    effnalPlot->SetMarkerStyle(8);
    effnalPlot->SetMarkerColor(kRed);
    effnalPlot->GetYaxis()->SetRangeUser(0, 1);
    effnalPlot->Draw();
    vertex_trueKE_eff_hist->Draw("HIST SAME E");
    vertex_recoKE_eff_hist->Draw("HIST SAME E");
    eff_canvas->BuildLegend();
    eff_canvas->Write();


    //Normalise the histograms so the max bin is at 1.
    double maxScale_bk = TMath::Max(vertex_trueKE_bk_hist->GetBinContent(vertex_trueKE_bk_hist->GetMaximumBin()),vertex_recoKE_bk_hist->GetBinContent(vertex_recoKE_bk_hist->GetMaximumBin()));
    vertex_trueKE_bk_hist->Scale(1/maxScale_bk);
    vertex_trueKE_bk_hist->Write();
    vertex_recoKE_bk_hist->Scale(1/maxScale_bk);
    vertex_recoKE_bk_hist->Write();

    vertex_trueKE_bk_hist->SetLineColor(kBlue);
    vertex_trueKE_bk_hist->SetFillColor(6);
    vertex_trueKE_bk_hist->SetFillStyle(3003);

    vertex_recoKE_bk_hist->SetLineColor(kGreen);
    vertex_recoKE_bk_hist->SetFillColor(42);
    vertex_recoKE_bk_hist->SetFillStyle(3003);

    //Draw the vertex efficiency canvas for the background.
    auto bk_canvas = new TCanvas("bk_canvas", "bk_canvas", 600, 400);
    TGraph* backgroundPlot = new TGraph(31,VertexEnergy,background); 
    backgroundPlot->SetMarkerStyle(8);
    backgroundPlot->SetMarkerColor(kRed);
    backgroundPlot->SetTitle("Reco Efficiency");
    backgroundPlot->GetXaxis()->SetTitle("True Vertex Energy (MeV)");
    backgroundPlot->GetYaxis()->SetRangeUser(0, 1);
    backgroundPlot->SetTitle("Reco Efficiency");
    backgroundPlot->GetXaxis()->SetTitle("True Vertex Energy Cut off (MeV)");
    backgroundPlot->Draw();
    vertex_trueKE_bk_hist->Draw("HIST SAME E");
    vertex_recoKE_bk_hist->Draw("HIST SAME E");
    bk_canvas->BuildLegend();
    bk_canvas->Write();
  }

    //#################################
    //### Vertex Conversion Gap Cut ###
    //#################################
    if(fVertexConversionCut){

      std::cout << "###############################" << std::endl;
      std::cout << "####### Conversion Cut ########" <<std::endl;
      std::cout << "###############################" << std::endl;


      file->cd("VertexConversionCut");

      float efficiency[21];
      float background[21];
      float effbk[21];

      //Define the cut position for the conversion gap and vertex
      float VertexEnergy[21] = {0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100.};
      float conversionGap[21]  = {0.,0.5,1.,1.5 ,2.,2.5 ,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.};


      //Get the branches for the signal
      std::vector<float>* vertex_recoKE_eff =0;
      tree_eff->SetBranchAddress("vertex_recoK",&vertex_recoKE_eff);
      
      std::vector<std::vector<float> >* shower_coversion_gap_eff=0;
      tree_eff->SetBranchAddress("shower_coversion_gap",&shower_coversion_gap_eff);

      std::vector<std::vector<float> >* shower_energy_eff = 0;
      tree_eff->SetBranchAddress("shower_energy",&shower_energy_eff);
      
      
      //2D Graph for conversion gap vs vertex Energy.
      TGraph2D *efficency_graph  = new TGraph2D();
      TGraph2D *background_graph = new TGraph2D();
      TGraph2D *effbk_graph      = new TGraph2D();

      //2D Hists for conversion gap vs vertex energy.
      TH2F* efficency_hist  = new TH2F("conversion efficiency","conversion efficiency",20,VertexEnergy,20,conversionGap);
      TH2F* background_hist = new TH2F("conversion bk","conversion bk ",20,VertexEnergy,20,conversionGap);
      TH2F* effbk_hist;

      //Histograms for the conversion gap cut at vertex energy cut at 0. MeV.
      TH1D* signalhist = new TH1D("conversion at 0 MeV signal","conversionsignal",100,0,10);
      TH1D* backgroundhist = new TH1D("conversion at 0 MeV background","conversionbackground",100,0,10);

      //Loop over the conversion gaps to analyze 
      Long64_t nentries_eff = tree_eff->GetEntries();
      for(int j=0; j<21; ++j){
	float convdist = conversionGap[j];
	//Loop over the vertex exnergies to analyze.
	for(int i=0; i<21; ++i){
	  float n_neutrinos = 0;
	  float Energy = VertexEnergy[i];
	  float eff=0;
	  //Loop over the events
	  for (Long64_t evt = 0; evt < nentries_eff; evt++) {
	    tree_eff->GetEntry(evt);

	    //Loop over the neutrinos 
	    for(int vertex=0; vertex<vertex_recoKE_eff->size(); ++ vertex){
	     
	      ++n_neutrinos;

	      //We need a shower.
	      if(shower_energy_eff->at(vertex).size() < 1){continue;}
 
	      //Only analyze events with one shower above fECut.
	      if(fApplyPreviousCuts){
		int numshowers = 0;
		for(auto const& shower_e: shower_energy_eff->at(vertex)){
		  if(shower_e > fECut){
		    ++numshowers;
		  }
		}
		if(numshowers != 1){continue;}
	      }

	      float vertex_recoE = vertex_recoKE_eff->at(vertex);
	      //Remove event if it has vertex lower than cut off.
	      if(vertex_recoE < Energy){++eff;}//continue;}

	      //Get the biggest shower iter
	      float max_energy=-999;
	      int max_energy_iter=-999;
	      for(int shower=0; shower<shower_energy_eff->at(vertex).size(); ++shower){
		if(shower_energy_eff->at(vertex).at(shower) > max_energy){
		  max_energy = shower_energy_eff->at(vertex).at(shower);
		  max_energy_iter = shower;
		}
	      }
	      
	      if(shower_coversion_gap_eff->at(vertex).size() -1 < max_energy_iter){continue;}

	      if(i==0){signalhist->Fill(shower_coversion_gap_eff->at(vertex).at(max_energy_iter));}
	      
	      //Apply the conversion gap to the biggest shower.
	      if(shower_coversion_gap_eff->at(vertex).at(max_energy_iter) < convdist){++eff;}
	      
	    }
	  }
	  if(i==0){efficiency[j] = eff/n_neutrinos;}
	  efficency_graph->SetPoint(efficency_graph->GetN(),Energy,convdist,eff/n_neutrinos);
	  efficency_hist->Fill(Energy,convdist,eff/n_neutrinos);
	}
      }



      //Get the background branches.
      std::vector<float>* vertex_recoKE_bk =0;
      tree_bk->SetBranchAddress("vertex_recoK",&vertex_recoKE_bk);
      
      std::vector<std::vector<float> >* shower_coversion_gap_bk=0;
      tree_bk->SetBranchAddress("shower_coversion_gap",&shower_coversion_gap_bk);

      std::vector<std::vector<float> >* shower_energy_bk = 0;
      tree_bk->SetBranchAddress("shower_energy",&shower_energy_bk);

      float max_effbk = -999;
      float max_eff   = -999;
      float max_bk    = -999;
      int max_i       = -999;

      //Loop over the conversion gaps to be analaysed.
      Long64_t nentries_bk = tree_bk->GetEntries();
      for(int j=0; j<21; ++j){
	float convdist = conversionGap[j];
	//Loop over the vertex energies to be analysed.
	for(int i=0; i<21; ++i){
	  float n_neutrinos = 0;
	  float Energy = VertexEnergy[i];
	  float bk=0;
	  //Loop over the events.
	  for (Long64_t evt = 0; evt < nentries_bk; evt++) {
	    tree_bk->GetEntry(evt);
	    for(int vertex=0; vertex<vertex_recoKE_bk->size(); ++ vertex){
	      
	      ++n_neutrinos;
	      //We need a shower.
	      if(shower_energy_bk->at(vertex).size() < 1){continue;}
	      
	      //Remove the event if there is more than one shower above fECut.
	      if(fApplyPreviousCuts){
		int numshowers = 0;
		for(auto const& shower_e: shower_energy_bk->at(vertex)){
		  if(shower_e > fECut){
		    ++numshowers;
		  }
		}
		if(numshowers != 1){continue;}
	      }

	      //Remove events which do not pass the vertex energy cut 
	      float vertex_recoE = vertex_recoKE_bk->at(vertex);
	      if(vertex_recoE < Energy){++bk;}//continue;}

	      //Get the biggest shower iter
	      float max_energy=-999;
	      int max_energy_iter=-999;
	      for(int shower=0; shower<shower_energy_bk->at(vertex).size(); ++shower){
		if(shower_energy_bk->at(vertex).at(shower) > max_energy){
		  max_energy = shower_energy_bk->at(vertex).at(shower);
		  max_energy_iter = shower;
		}
	      }
	      
	      if((int)(shower_coversion_gap_bk->at(vertex).size() -1) < max_energy_iter){continue;}
	      if(shower_coversion_gap_bk->at(vertex).at(max_energy_iter) < convdist){++bk;}
	      if(i==0){backgroundhist->Fill(shower_coversion_gap_bk->at(vertex).at(max_energy_iter));}
	    }
	  }

	  if(i==0){
	    background[j] = 1- bk/n_neutrinos;
	    effbk[j] = background[j]*efficiency[j];
	  };

	  background_graph->SetPoint(background_graph->GetN(),Energy,convdist,1-bk/n_neutrinos);
	  background_hist->Fill(Energy,convdist, 1- bk/n_neutrinos);
	  
	  double Eff_E=0;
	  double Eff_con=0;
	  double Eff=0;
	  efficency_graph->GetPoint(background_graph->GetN()-1,Eff_E,Eff_con,Eff);
	  effbk_graph->SetPoint(effbk_graph->GetN(),Energy,convdist,Eff*(1-bk/n_neutrinos));

	  //Get the biggest eff*bk background.
	  if(Eff*background[i] > max_effbk){
	    max_effbk = Eff*(1- bk/n_neutrinos); 
	    max_eff = Eff;
	    max_bk =  1- bk/n_neutrinos;
	    max_i = i; 
	  }

	}
      }


      std::cout << "Best Conversion Cut: " << max_effbk << " has conversion gap cut at: " << conversionGap[max_i] << " and vertex energy: " << VertexEnergy[max_i] << " eff: " << max_eff << " bk: " << max_bk << std::endl; 
      
      file->cd("VertexConversionCut");

      //Apply the cutfinder procedure.
      DoCutFinding(signalhist,backgroundhist,"cm","Conversion Gap (cm)");

      file->cd("VertexConversionCut");
      
      //Make the eff*bk hist 
      effbk_hist = (TH2F*) efficency_hist->Clone("conversion eff*bkgr rek hist");
      effbk_hist->Multiply(background_hist);
      effbk_hist->GetXaxis()->SetTitle("Vertex Energy (MeV)");
      effbk_hist->GetYaxis()->SetTitle("Conversion distance (MeV)");

      efficency_hist->GetXaxis()->SetTitle("Vertex Energy (MeV)");
      efficency_hist->GetYaxis()->SetTitle("Conversion distance (MeV)");

      background_hist->GetXaxis()->SetTitle("Vertex Energy (MeV)");
      background_hist->GetYaxis()->SetTitle("Conversion distance (MeV)");
      
      efficency_graph->SetTitle("conversion Efficiency");
      efficency_graph->GetXaxis()->SetTitle("Vertex Energy (MeV)");
      efficency_graph->GetYaxis()->SetTitle("Conversion distance (MeV)");
      efficency_graph->SetName("conversion Efficiency");

      background_graph->SetTitle("conversion Background Rejection");
      background_graph->GetXaxis()->SetTitle("Vertex Energy (MeV)");
      background_graph->GetYaxis()->SetTitle("Conversion distance (MeV)");
      background_graph->SetName("conversion background Rejection");

      effbk_graph->SetTitle("conversion eff*bkgr rej");
      effbk_graph->GetXaxis()->SetTitle("Vertex Energy (MeV)");
      effbk_graph->GetYaxis()->SetTitle("Conversion distance (MeV)");
      effbk_graph->SetName("conversion eff*bkgr re");

      efficency_graph->Write();
      background_graph->Write();
      effbk_graph->Write();
      efficency_hist->Write();
      background_hist->Write();
      effbk_hist->Write();

      TH2D* eff_hist = efficency_graph->GetHistogram();
      eff_hist->SetTitle("Efficiency from TGraph");
      eff_hist->SetName("Efficiency from TGraph");
      eff_hist->GetXaxis()->SetTitle("Vertex Energy (MeV)");
      eff_hist->GetYaxis()->SetTitle("Conversion distance (MeV)");

      TH2D* bk_hist  = background_graph->GetHistogram();
      bk_hist->SetTitle("Background Rejection from TGraph");
      bk_hist->SetName("Background Rejection from TGraph");
      bk_hist->GetXaxis()->SetTitle("Vertex Energy (MeV)");
      bk_hist->GetYaxis()->SetTitle("Conversion distance (MeV)");


      TH2D* effbg_hist = effbk_graph->GetHistogram();
      effbg_hist->SetTitle("Efficiency*Background Rejection from TGraph");
      effbg_hist->SetName("Efficiency*Background Rejection from TGraph");
      effbg_hist->GetXaxis()->SetTitle("Vertex Energy (MeV)");
      effbg_hist->GetYaxis()->SetTitle("Conversion distance (MeV)");


      eff_hist->Write();
      bk_hist->Write();
      effbg_hist->Write();

      
    }

    if(fdEdxCut){

      file->cd("dEdxCut");

      TH1D* signalhist = new TH1D("dEdx signal","dEdx signal",50,0,10);
      TH1D* backgroundhist = new TH1D("dEdx background","dEdx background",50,0,10);

      //Get the signal branches
      std::vector<float>* vertex_recoKE_eff =0;
      tree_eff->SetBranchAddress("vertex_recoK",&vertex_recoKE_eff);
      
      std::vector<std::vector<float> >* shower_coversion_gap_eff=0;
      tree_eff->SetBranchAddress("shower_coversion_gap",&shower_coversion_gap_eff);
      
      std::vector<std::vector<float> >* shower_energy_eff = 0;
      tree_eff->SetBranchAddress("shower_energy",&shower_energy_eff);
      
      std::vector<std::vector<float> >* shower_dEdx_eff=0;
	tree_eff->SetBranchAddress("shower_dEdx",&shower_dEdx_eff);

      Long64_t nentries_eff = tree_eff->GetEntries();
      //Loop over the events
      for (Long64_t evt = 0; evt < nentries_eff; evt++) {
	tree_eff->GetEntry(evt);

	//Loop over the neutrinos.
	for(int vertex=0; vertex<vertex_recoKE_eff->size(); ++ vertex){
	  
	  if(fApplyPreviousCuts){
	    int numshowers = 0;
	    for(int shower_iter=0; shower_iter<shower_energy_eff->at(vertex).size(); ++shower_iter){
	      float shower_e = shower_energy_eff->at(vertex).at(shower_iter);
	      if(shower_e > fECut){
		++numshowers;
	      }
	    }
	    //Check there is only one shower above the fECut.
	    if(numshowers != 1){continue;}

	    //Get the largest shower.
	    int max_i = -999;
	    float max_e = -999;
	    for(int shower_iter=0; shower_iter<shower_energy_eff->at(vertex).size(); ++shower_iter){
	      
	      if(max_e < shower_energy_eff->at(vertex).at(shower_iter)){
		max_i = shower_iter;
		max_e = shower_energy_eff->at(vertex).at(shower_iter);
	      }
	    }
	    

	    //Check the event has a vertex energy greater than the cut off and the conversion gap is less than the cut off.
	    //	    if(shower_coversion_gap_eff->at(vertex).at(max_i) > fConversionGapCut &&
	    //  vertex_recoKE_eff->at(vertex) > fVertexCut){
	    // continue;
	    //}



	  }

	  //Get the largest shower.
	  int max_i = -999;
	  float max_e = -999;
	  for(int shower_iter=0; shower_iter<shower_energy_eff->at(vertex).size(); ++shower_iter){

	    if(max_e < shower_energy_eff->at(vertex).at(shower_iter)){
	      max_i = shower_iter;
	      max_e = shower_energy_eff->at(vertex).at(shower_iter);
	    }
	  }

	  if(max_i < -1){continue;}
	  //Fill the dEdx hist
	  signalhist->Fill(shower_dEdx_eff->at(vertex).at(max_i));
	  
	}
      }


      
      //Get the background branches  
      std::vector<float>* vertex_recoKE_bk =0;
      tree_bk->SetBranchAddress("vertex_recoK",&vertex_recoKE_bk);
      
      std::vector<std::vector<float> >* shower_coversion_gap_bk=0;
      tree_bk->SetBranchAddress("shower_coversion_gap",&shower_coversion_gap_bk);

      std::vector<std::vector<float> >* shower_energy_bk = 0;
      tree_bk->SetBranchAddress("shower_energy",&shower_energy_bk);
      
      std::vector<std::vector<float> >* shower_dEdx_bk=0;
      tree_bk->SetBranchAddress("shower_dEdx",&shower_dEdx_bk);
	
      //Loop over branches.
      Long64_t nentries_bk = tree_bk->GetEntries();
      for (Long64_t evt = 0; evt < nentries_bk; evt++) {
	tree_bk->GetEntry(evt);

	//Loop over neutrinos 
	for(int vertex=0; vertex<vertex_recoKE_bk->size(); ++ vertex){
	  
	  if(fApplyPreviousCuts){
	    int numshowers = 0;
	    for(int shower_iter=0; shower_iter<shower_energy_bk->at(vertex).size(); ++shower_iter){
	      float shower_e = shower_energy_bk->at(vertex).at(shower_iter);
	      if(shower_e > fECut){

		++numshowers;
	      }
	    }
	    if(numshowers != 1){continue;}

	    //Get the largest shower.
	    int max_i = -999;
	    float max_e = -999;
	    for(int shower_iter=0; shower_iter<shower_energy_bk->at(vertex).size(); ++shower_iter){
	      if(max_e < shower_energy_bk->at(vertex).at(shower_iter)){
		max_i = shower_iter;
		max_e = shower_energy_bk->at(vertex).at(shower_iter);
	      }
	    }
	    //Check the event has a vertex energy greater than the cut off and the conversion gap is less than the cut off.
	    if(shower_coversion_gap_bk->at(vertex).at(max_i) > fConversionGapCut &&
	       vertex_recoKE_bk->at(vertex) > fVertexCut){
	      continue;
	    }
	    
	  }

	  //Get the largest shower.
	  int max_i = -999;
	  float max_e = -999;
	  for(int shower_iter=0; shower_iter<shower_energy_bk->at(vertex).size(); ++shower_iter){
	    if(max_e < shower_energy_bk->at(vertex).at(shower_iter)){
	      max_i = shower_iter;
	      max_e = shower_energy_bk->at(vertex).at(shower_iter);
	    }
	  }

	  //Check there is only one shower above the fECut.
	  if(max_i < 0){continue;}  
	  backgroundhist->Fill(shower_dEdx_bk->at(vertex).at(max_i));
	}
      }

      std::cout << "###############################" << std::endl;
      std::cout << "########## dEdx value #########" << std::endl;
      std::cout << "###############################" << std::endl;
      
      //Get the dEdx cut
      DoCutFinding(signalhist,backgroundhist,"MeV/cm","dEdx (MeV/cm)");


    }

    if(fTrackLengthCut){

      std::vector<float> int_type = {1001.,1003.,1004.,1005.,1010.,1011.,1012.,1017.,1021.,1028.,1032.,1039.,1041.,1046.,1048.,1053.,1055.,1062.,1060.,1067.,1070.,1073.,1079.,1080.,1085.,1086.,1090.,1091.,1095.,1097.};
      
      file->cd("TrackLengthCut");

      TH1D* signalhist = new TH1D("TrackLengthCut signal","TrackLengthCut signal",50,0,200);
      TH1D* backgroundhist = new TH1D("TrackLengthCut background","TrackLengthCut background",50,0,200);

      //Get the signal branches
      std::vector<float>* vertex_recoKE_eff =0;
      tree_eff->SetBranchAddress("vertex_recoK",&vertex_recoKE_eff);
      
      std::vector<std::vector<float> >* shower_coversion_gap_eff=0;
      tree_eff->SetBranchAddress("shower_coversion_gap",&shower_coversion_gap_eff);
      
      std::vector<std::vector<float> >* shower_energy_eff = 0;
      tree_eff->SetBranchAddress("shower_energy",&shower_energy_eff);
      
      std::vector<std::vector<float> >* shower_dEdx_eff=0;
      tree_eff->SetBranchAddress("shower_dEdx",&shower_dEdx_eff);

      std::vector<std::vector<float> >* track_lengths_eff =0;
      tree_eff->SetBranchAddress("track_lengths",&track_lengths_eff);

      std::vector<float>* inttype_length_eff=0;
      tree_eff->SetBranchAddress("nu_interaction_type",&inttype_length_eff);

      Long64_t nentries_eff = tree_eff->GetEntries();
      //Loop over the events
      for (Long64_t evt = 0; evt < nentries_eff; evt++) {
	tree_eff->GetEntry(evt);

	//Loop over the neutrinos.
	for(int vertex=0; vertex<vertex_recoKE_eff->size(); ++ vertex){
	  
	  //	  std::cout << " inttype_length_eff->at(vertex): " << inttype_length_eff->at(vertex) << std::endl;
	  //if(std::find(int_type.begin(),int_type.end(),inttype_length_eff->at(vertex)) == int_type.end()){continue;}
	  //std::cout << " inttype_length_eff->at(vertex): " << inttype_length_eff->at(vertex) << std::endl;

	  if(fApplyPreviousCuts){
	    int numshowers = 0;
	    for(int shower_iter=0; shower_iter<shower_energy_eff->at(vertex).size(); ++shower_iter){
	      float shower_e = shower_energy_eff->at(vertex).at(shower_iter);
	      if(shower_e > fECut){
		++numshowers;
	      }
	    }
	    //Check there is only one shower above the fECut.
	    if(numshowers != 1){continue;}
	  
	    //Get the largest shower.
	    int max_i = -999;
	    float max_e = -999;
	    for(int shower_iter=0; shower_iter<shower_energy_eff->at(vertex).size(); ++shower_iter){
	      
	      if(max_e < shower_energy_eff->at(vertex).at(shower_iter)){
		max_i = shower_iter;
		max_e = shower_energy_eff->at(vertex).at(shower_iter);
	      }
	    }

	    //Check the event has a vertex energy greater than the cut off and the conversion gap is less than the cut off.
	    //	    if(shower_coversion_gap_eff->at(vertex).at(max_i) > fConversionGapCut &&
	    //   vertex_recoKE_eff->at(vertex) > fVertexCut){
	    //  continue;
	    // }



	    //See if we pass the dEdx cut.
	    //	    if(shower_dEdx_eff->at(vertex).at(max_i) > fdEdxCutVal){continue;}
	  }

	  //get the biggest track.
	  float max_length = -999;
	  for(int track_iter=0; track_iter<track_lengths_eff->at(vertex).size(); ++track_iter){
	    if(max_length < track_lengths_eff->at(vertex).at(track_iter)){
	      max_length = track_lengths_eff->at(vertex).at(track_iter);
	    }
	  }

	  if(max_length == -999){
	    signalhist->Fill(0);
	  }
	  else{
	    signalhist->Fill(max_length);
	  }

	}
      }


      
      //Get the background branches  
      std::vector<float>* vertex_recoKE_bk =0;
      tree_bk->SetBranchAddress("vertex_recoK",&vertex_recoKE_bk);
      
      std::vector<std::vector<float> >* shower_coversion_gap_bk=0;
      tree_bk->SetBranchAddress("shower_coversion_gap",&shower_coversion_gap_bk);

      std::vector<std::vector<float> >* shower_energy_bk = 0;
      tree_bk->SetBranchAddress("shower_energy",&shower_energy_bk);
      
      std::vector<std::vector<float> >* shower_dEdx_bk=0;
      tree_bk->SetBranchAddress("shower_dEdx",&shower_dEdx_bk);

      std::vector<std::vector<float> >* track_lengths_bk =0;
      tree_bk->SetBranchAddress("track_lengths",&track_lengths_bk);
	
      std::vector<float>* inttype_length_bk=0;
      tree_bk->SetBranchAddress("nu_interaction_type",&inttype_length_bk);

      //Loop over branches.
      Long64_t nentries_bk = tree_bk->GetEntries();
      for (Long64_t evt = 0; evt < nentries_bk; evt++) {
	tree_bk->GetEntry(evt);

	//Loop over neutrinos 
	for(int vertex=0; vertex<vertex_recoKE_bk->size(); ++ vertex){
	  
	  if(std::find(int_type.begin(),int_type.end(),inttype_length_bk->at(vertex)) == int_type.end()){continue;}

	  if(fApplyPreviousCuts){
	    int numshowers = 0;
	    for(int shower_iter=0; shower_iter<shower_energy_bk->at(vertex).size(); ++shower_iter){
	      float shower_e = shower_energy_bk->at(vertex).at(shower_iter);
	      if(shower_e > fECut){
		++numshowers;
	      }
	    }
	    if(numshowers != 1){continue;}
	  

	    //Get the largest shower.
	    int max_i = -999;
	    float max_e = -999;
	    for(int shower_iter=0; shower_iter<shower_energy_bk->at(vertex).size(); ++shower_iter){
	      if(max_e < shower_energy_bk->at(vertex).at(shower_iter)){
		max_i = shower_iter;
		max_e = shower_energy_bk->at(vertex).at(shower_iter);
	      }
	    }

	    //Check the event has a vertex energy greater than the cut off and the conversion gap is less than the cut off.
	    if(shower_coversion_gap_bk->at(vertex).at(max_i) > fConversionGapCut &&
	       vertex_recoKE_bk->at(vertex) > fVertexCut){
	      continue;
	    }
	    

	    //See if we pass the dEdx cut.
	    if(shower_dEdx_bk->at(vertex).at(max_i) > fdEdxCutVal){continue;}
	  }

	  //get the biggest track.
	  float max_length = -999;
	  for(int track_iter=0; track_iter<track_lengths_bk->at(vertex).size(); ++track_iter){
	    if(max_length < track_lengths_bk->at(vertex).at(track_iter)){
	      max_length = track_lengths_bk->at(vertex).at(track_iter);
	    }
	  }
	  
	  if(max_length == -999){
	    backgroundhist->Fill(0);
	  }
	  else{
	    backgroundhist->Fill(max_length);
	  }
	}
      }

      std::cout << "###############################" << std::endl;
      std::cout << "####### Track Length Cut ######" << std::endl;
      std::cout << "###############################" << std::endl;
      
      //Get the dEdx cut
      DoCutFinding(signalhist,backgroundhist,"cm","Max Track Length (cm)");


    }

    //########################################
    //###### Track and Shower Validation #####
    //########################################

    if(fTrackValidation){
   
      file->cd("TrackValidation");

      //Get the branches 
      std::vector<std::vector<float> >* track_E_eff=0;
      tree_eff->SetBranchAddress("track_E",&track_E_eff);

      std::vector<std::vector<float> >* track_trueE_eff=0;
      tree_eff->SetBranchAddress("track_trueE",&track_trueE_eff);

      std::vector<std::vector<float> >* track_track_pdg_eff=0;
      tree_eff->SetBranchAddress("track_pdg",&track_track_pdg_eff);

      std::vector<std::vector<float> >* truepionE_eff=0;
      tree_eff->SetBranchAddress("true_pionE",&truepionE_eff);

      std::vector<std::vector<float> >* trueprotonE_eff=0;
      tree_eff->SetBranchAddress("true_protonE",&trueprotonE_eff);

      std::vector<std::vector<float> >* truekaonE_eff=0;
      tree_eff->SetBranchAddress("true_kaonE",&truekaonE_eff);

      std::vector<std::vector<float> >*truetrackE_eff=0;
      tree_eff->SetBranchAddress("true_trackE",&truetrackE_eff);
      
      std::vector<std::vector<float> >* track_resE_eff=0;
      tree_eff->SetBranchAddress("track_resE",&track_resE_eff);


      std::vector<std::vector<float> >* shower_true_E_eff =0;
      tree_eff->SetBranchAddress("true_energy",&shower_true_E_eff);

      std::vector<std::vector<float> >* truth_pid_eff =0;
      tree_eff->SetBranchAddress("truth_pid",&truth_pid_eff);

      std::vector<std::vector<float> >* shower_energy_eff=0;
      tree_eff->SetBranchAddress("shower_energy",&shower_energy_eff);

      //Get the signal branches
      std::vector<float>* vertex_recoKE_eff =0;
      tree_eff->SetBranchAddress("vertex_recoK",&vertex_recoKE_eff);
      
      std::vector<std::vector<float> >* shower_coversion_gap_eff=0;
      tree_eff->SetBranchAddress("shower_coversion_gap",&shower_coversion_gap_eff);
      
      std::vector<std::vector<float> >* shower_dEdx_eff=0;
      tree_eff->SetBranchAddress("shower_dEdx",&shower_dEdx_eff);

      std::vector<float>* inttype_length_eff=0;
      tree_eff->SetBranchAddress("nu_interaction_type",&inttype_length_eff);

      std::vector<std::vector<float> >* track_lengths_eff=0;
      tree_eff->SetBranchAddress("track_lengths",&track_lengths_eff);

      std::vector<float>* nu_reco_energy_eff=0;
      tree_eff->SetBranchAddress("nu_reco_energy",&nu_reco_energy_eff);

      std::vector<float> * nu_truth_energy_eff=0;
      tree_eff->SetBranchAddress("nu_truth_energy",&nu_truth_energy_eff);
      

      //make the histgorams 
      TH1F* track_res_hist_eff_all = new TH1F("track_res_hist_eff_all","track_res_hist_eff_all",100,-1.5,1.5);
      TH1F* proton_res_hist_eff = new TH1F("proton_res_hist_eff","proton_res_hist_eff",100,-1.5,1.5);
      TH1F* pion_res_hist_eff  = new TH1F("pion_res_hist_eff","pion_res_hist_eff",100,-1.5,1.5);
      TH1F* kaon_res_hist_eff  = new TH1F("kaon_res_hist_eff","kaon_res_hist_eff",100,-1.5,1.5);
      TH1F* muon_res_hist_eff  = new TH1F("muon_res_hist_eff","muon_res_hist_eff",100,-1.5,1.5);

      TH1F* protontrueE_hist_eff = new TH1F("proton_trueE_hist_eff","proton_trueE_hist_eff",100,0,1000);
      TH1F* piontrueE_hist_eff = new TH1F("pion_trueE_hist_eff","pion_trueE_hist_eff",100,0,1000);
      TH1F* kaontrueE_hist_eff = new TH1F("kaon_trueE_hist_eff","kaon_trueE_hist_eff",100,0,1000);
      TH1F* alltrueE_hist_eff = new TH1F("all_trueE_hist_eff","all_trueE_hist_eff",100,0,1000);
      TH1F* proton_reco_hist_eff = new TH1F("proton_reco_hist_eff","proton_reco_hist_eff",100,0,1000);
      TH1F* pion_reco_hist_eff = new TH1F("pion_reco_hist_eff","pion_reco_hist_eff",100,0,1000);
      TH1F* kaon_reco_hist_eff = new TH1F("kaon_reco_hist_eff","kaon_reco_hist_eff",100,0,1000);
      TH1F* all_reco_hist_eff = new TH1F("all_reco_hist_eff","all_reco_hist_eff",100,0,1000);
      TH1F* proton_reco_hist_eff_extra = new TH1F("proton_reco_hist_eff_extra","proton_reco_hist_eff_extra",100,0,1000);
      TH1F* pion_reco_hist_eff_extra = new TH1F("pion_reco_hist_eff_extra","pion_reco_hist_eff_extra",100,0,1000);
      TH1F* kaon_reco_hist_eff_extra = new TH1F("kaon_reco_hist_eff_extra","kaon_reco_hist_eff_extra",100,0,1000);
      TH1F* all_reco_hist_eff_extra = new TH1F("all_reco_hist_eff_extra","all_reco_hist_eff_extra",100,0,1000);

      TH1F* neutrino_true_energy_before_cut_eff = new TH1F("neutrino_true_energy_before_cut_eff","neutrino_true_energy_before_cut_eff",100,0,4000);
      TH1F* neutrino_true_energy_after_cut_eff = new TH1F("neutrino_true_energy_after_cut_eff","neutrino_true_energy_after_cut_eff",100,0,4000);

      TH1F* neutrino_reco_energy_before_cut_eff = new TH1F("neutrino_reco_energy_before_cut_eff","neutrino_reco_energy_before_cut_eff",100,0,4000);
      TH1F* neutrino_reco_energy_after_cut_eff = new TH1F("neutrino_reco_energy_after_cut_eff","neutrino_reco_energy_after_cut_eff",100,0,4000);
      TH1F* neutrino_reco_energy_extra_cut_eff = new TH1F("neutrino_reco_energy_after_extra_cut_eff","neutrino_reco_energy_after_cut_extra_eff",100,0,4000);

      TH1F* neutrino_res_eff = new TH1F("neutrino_res_eff","neutrino_res_eff",100,-5,5);


      
      std::vector<float> all_energy_vec;
      std::vector<float> all_res_vec;
      std::vector<float> proton_energy_vec;
      std::vector<float> proton_res_vec;
      std::vector<float> pion_energy_vec;
      std::vector<float> pion_res_vec;
      std::vector<float> kaon_energy_vec;
      std::vector<float> kaon_res_vec;
      std::vector<float> muon_energy_vec;
      std::vector<float> muon_res_vec;
      

      TH1F* shower_res_hist_eff = new TH1F("shower_res_hist_eff","shower_res_hist_eff",100,-1.5,1.5);
      std::vector<float> shower_res_vec_eff;
      std::vector<float> shower_true_e_vec_eff;

      std::vector<float> neutrino_res_eff_vec;
      std::vector<float> neutrino_true_energy_vec_eff;

      float num_neutrinos_before_cuts_eff = 0;
      float num_neutrinos_after_cuts_eff = 0;
      
      Long64_t nentries_eff = tree_eff->GetEntries();
      //Loop over the events
      for (Long64_t evt = 0; evt < nentries_eff; evt++) {
	tree_eff->GetEntry(evt);

	//Loop over the neutrinos.
	for(int neutrino=0; neutrino<track_E_eff->size(); ++neutrino){

	  ++num_neutrinos_before_cuts_eff;

	  neutrino_true_energy_before_cut_eff->Fill(nu_truth_energy_eff->at(neutrino)*1000);
	  neutrino_reco_energy_before_cut_eff->Fill(nu_reco_energy_eff->at(neutrino));

	  if(fApplyPreviousCuts){
	    int numshowers = 0;
	    for(int shower_iter=0; shower_iter<shower_energy_eff->at(neutrino).size(); ++shower_iter){

	      float shower_e = shower_energy_eff->at(neutrino).at(shower_iter);
	      if(shower_e > fECut){
		++numshowers;
	      }
	    }
	    if(numshowers != 1){continue;}
	  
	    
	    //Get the largest shower.
	    float max_e = -999;
	    int max_i = -999;
	    for(int shower_iter=0; shower_iter<shower_energy_eff->at(neutrino).size(); ++shower_iter){
	      if(max_e < shower_energy_eff->at(neutrino).at(shower_iter)){
	    	max_i = shower_iter;
	    	max_e = shower_energy_eff->at(neutrino).at(shower_iter);
	      }
	    }
	    
	    //Check the event has a vertex energy greater than the cut off and the conversion gap is less than the cut off.
	    if(shower_coversion_gap_eff->at(neutrino).at(max_i) > fConversionGapCut &&
	       vertex_recoKE_eff->at(neutrino) > fVertexCut){
	      continue;
	    }
	    
	    //See if we pass the dEdx cut.
	    if(shower_dEdx_eff->at(neutrino).at(max_i) > fdEdxCutVal){continue;}
	    
	    //get the biggest track.
	    float max_length = -999;
	    for(int track_iter=0; track_iter<track_lengths_eff->at(neutrino).size(); ++track_iter){
	      if(max_length < track_lengths_eff->at(neutrino).at(track_iter)){
	    	max_length = track_lengths_eff->at(neutrino).at(track_iter);
	      }
	    }
	    if(max_length > fMaxLengthCut){continue;}

	  }

	  ++num_neutrinos_after_cuts_eff;
	  neutrino_true_energy_after_cut_eff->Fill(nu_truth_energy_eff->at(neutrino)*1000);
	  neutrino_reco_energy_after_cut_eff->Fill(nu_reco_energy_eff->at(neutrino));
	  neutrino_res_eff->Fill(nu_reco_energy_eff->at(neutrino)/1000*nu_truth_energy_eff->at(neutrino));
	  neutrino_res_eff_vec.push_back(nu_reco_energy_eff->at(neutrino)/1000*nu_truth_energy_eff->at(neutrino));
	  neutrino_true_energy_vec_eff.push_back(nu_truth_energy_eff->at(neutrino));

	  //Loop over the showers. 
	  int max_i = -9999;
	  float max_e = -9999;
	  for(int shower=0; shower<shower_energy_eff->at(neutrino).size(); ++shower){
	    if(shower_energy_eff->at(neutrino).at(shower) > max_e){
	      max_e = shower_energy_eff->at(neutrino).at(shower); 
	      max_i = shower;
	    }
	  }

	  
	  if(max_i == -9999){continue;}

	  float shower_res = (shower_true_E_eff->at(neutrino).at(max_i) - shower_energy_eff->at(neutrino).at(max_i))/shower_true_E_eff->at(neutrino).at(max_i);
	  shower_res_hist_eff->Fill(shower_res);
	  shower_res_vec_eff.push_back(shower_res);
	  shower_true_e_vec_eff.push_back(shower_true_E_eff->at(neutrino).at(max_i)/1000);

	  //Loop over the true tracks 
	  for(int track=0; track<truekaonE_eff->at(neutrino).size(); ++track){
	    //Fill true Hist_Eff
	    kaontrueE_hist_eff->Fill(truekaonE_eff->at(neutrino).at(track));
	  }

	  for(int track=0; track<trueprotonE_eff->at(neutrino).size(); ++track){
	    protontrueE_hist_eff->Fill(trueprotonE_eff->at(neutrino).at(track));
	  }

	  for(int track=0; track<truepionE_eff->at(neutrino).size(); ++track){
	    piontrueE_hist_eff->Fill(truepionE_eff->at(neutrino).at(track));
	  }

	  for(int track=0; track<truetrackE_eff->at(neutrino).size(); ++track){
	    alltrueE_hist_eff->Fill(truetrackE_eff->at(neutrino).at(track));
	  }

	  std::vector<float> truetrackE_eff_temp;
	  //Loop of the reco tracks.
	  for(int track=0; track<track_E_eff->at(neutrino).size(); ++track){

	    float true_trackE  = track_trueE_eff->at(neutrino).at(track);
	    float track_pdg    = track_track_pdg_eff->at(neutrino).at(track);
	    float reco_track_E = track_E_eff->at(neutrino).at(track);
	 
	    track_res_hist_eff_all->Fill((true_trackE-reco_track_E)/true_trackE);
	    all_energy_vec.push_back(true_trackE);
	    all_res_vec.push_back((true_trackE-reco_track_E)/true_trackE);


	    if(track_pdg == 2212){
	      proton_energy_vec.push_back(true_trackE);
	      proton_res_vec.push_back((true_trackE-reco_track_E)/true_trackE);
	      proton_res_hist_eff->Fill((true_trackE-reco_track_E)/true_trackE);
	    }
	    if(track_pdg == 211){
	      pion_energy_vec.push_back(true_trackE);
	      pion_res_vec.push_back((true_trackE-reco_track_E)/true_trackE);
	      pion_res_hist_eff->Fill((true_trackE-reco_track_E)/true_trackE);
	    }
	    if(track_pdg == 321){
	      kaon_energy_vec.push_back(true_trackE);
	      kaon_res_vec.push_back((true_trackE-reco_track_E)/true_trackE);
	      kaon_res_hist_eff->Fill((true_trackE-reco_track_E)/true_trackE);
	    }
	    if(track_pdg == 13){
	      muon_energy_vec.push_back(true_trackE);
	      muon_res_vec.push_back((true_trackE-reco_track_E)/true_trackE);
	      muon_res_hist_eff->Fill((true_trackE-reco_track_E)/true_trackE);
	    }

	    //Only considering "primaries' 
	    if(std::find(truetrackE_eff->at(neutrino).begin(),truetrackE_eff->at(neutrino).end(),true_trackE) == truetrackE_eff->at(neutrino).end()){continue;}

	    //See if the track has been matched yet
	    if(std::find(truetrackE_eff_temp.begin(),truetrackE_eff_temp.end(),true_trackE) != truetrackE_eff_temp.end()){
	      if(track_pdg == 2212){proton_reco_hist_eff_extra->Fill(true_trackE);}
	      if(track_pdg == 211){pion_reco_hist_eff_extra->Fill(true_trackE);}
	      if(track_pdg == 321){kaon_reco_hist_eff_extra->Fill(true_trackE);}
	      all_reco_hist_eff_extra->Fill(true_trackE);
	    }
	    else{
	      truetrackE_eff_temp.push_back(true_trackE);
	      if(track_pdg == 2212){proton_reco_hist_eff->Fill(true_trackE);}
	      if(track_pdg == 211){pion_reco_hist_eff->Fill(true_trackE);}
	      if(track_pdg == 321){kaon_reco_hist_eff->Fill(true_trackE);}
	      all_reco_hist_eff->Fill(true_trackE);
	    }
	  }
	}
      }

      
      //Make the efficency plots 
      MakeEfficiencyPlots(proton_reco_hist_eff,protontrueE_hist_eff,proton_reco_hist_eff_extra,"number of reco protons","number of true protons","number of multiple reconstructions","Efficiency;Track Energy (MeV);","protons_signal");
      MakeEfficiencyPlots(pion_reco_hist_eff,piontrueE_hist_eff,pion_reco_hist_eff_extra,"number of reco pions","number of true pions","number of multiple reconstructions","Efficiency;Track Energy (MeV);","pions_signal");
      MakeEfficiencyPlots(kaon_reco_hist_eff,kaontrueE_hist_eff,kaon_reco_hist_eff_extra,"number of reco kaons","number of true kaons","number of multiple reconstructions","Efficiency;Track Energy (MeV);","kaons_signal");
      MakeEfficiencyPlots(all_reco_hist_eff,alltrueE_hist_eff,all_reco_hist_eff_extra,"number of reco tracks","number of true tracks","number of multiple reconstructions","Efficiency;Track Energy (MeV);","all_signal");
    
      //Make the resolution graphs.
      MakeResolutionGraphs(proton_res_vec,proton_energy_vec,0,500,10,50,"True Energy - Reco Energy/TrueEnergy","track_resE","protons signal",50);
      MakeResolutionGraphs(pion_res_vec,pion_energy_vec,0,500,10,50,"True Energy - Reco Energy/TrueEnergy","track_resE","pions signal",50);
      MakeResolutionGraphs(kaon_res_vec,kaon_energy_vec,0,500,10,50,"True Energy - Reco Energy/TrueEnergy","track_resE","kaons signal",50);
      MakeResolutionGraphs(muon_res_vec,muon_energy_vec,0,1000,10,50,"True Energy - Reco Energy/TrueEnergy","track_resE","muons signal",50);
      MakeResolutionGraphs(all_res_vec,all_energy_vec,0,500,10,50,"True Energy - Reco Energy/TrueEnergy","track_resE","all signal",50);
      
      track_res_hist_eff_all->Write();
      proton_reco_hist_eff->Write();
      pion_reco_hist_eff->Write();
      kaon_reco_hist_eff->Write();
      muon_res_hist_eff->Write();

      neutrino_res_eff->Write();
      MakeEfficiencyPlots(neutrino_true_energy_after_cut_eff,neutrino_true_energy_before_cut_eff,neutrino_reco_energy_extra_cut_eff,"number neutrinos before selection true","number neutrinos after selection true","number of neutrinos extra rrue","Efficiency;Neutrino Energy (MeV);","electron neutrino true");
      MakeEfficiencyPlots(neutrino_reco_energy_after_cut_eff,neutrino_reco_energy_before_cut_eff,neutrino_reco_energy_extra_cut_eff,"number neutrinos before selection reco","number neutrinos after selection reco","number of neutrinos after extra reco","Efficiency;Neutrino Energy (MeV);","electron neutrino reco");
      MakeResolutionGraphs(neutrino_res_eff_vec,neutrino_true_energy_vec_eff,0,4000,10,100,"Reco Energy/True Energy","neutrino_resE","electron neutrino res",50);

      std::cout << " shower res" << std::endl;
      shower_res_hist_eff->Write();
      MakeResolutionGraphs(shower_res_vec_eff,shower_true_e_vec_eff,0.3,1.2,15,0.1,"True Energy - Reco Energy/TrueEnergy","shower_resE","Shower signal",300);



      //Get the branches 
      std::vector<std::vector<float> >* track_E_bk=0;
      tree_bk->SetBranchAddress("track_E",&track_E_bk);

      std::vector<std::vector<float> >* track_trueE_bk=0;
      tree_bk->SetBranchAddress("track_trueE",&track_trueE_bk);

      std::vector<std::vector<float> >* track_track_pdg_bk=0;
      tree_bk->SetBranchAddress("track_pdg",&track_track_pdg_bk);

      std::vector<std::vector<float> >* truepionE_bk=0;
      tree_bk->SetBranchAddress("true_pionE",&truepionE_bk);

      std::vector<std::vector<float> >* trueprotonE_bk=0;
      tree_bk->SetBranchAddress("true_protonE",&trueprotonE_bk);

      std::vector<std::vector<float> >* truekaonE_bk=0;
      tree_bk->SetBranchAddress("true_kaonE",&truekaonE_bk);

      std::vector<std::vector<float> >*truetrackE_bk=0;
      tree_bk->SetBranchAddress("true_trackE",&truetrackE_bk);

      
      std::vector<std::vector<float> >* track_resE_bk=0;
      tree_bk->SetBranchAddress("track_resE",&track_resE_bk);

      std::vector<std::vector<float> >* track_lengths_bk=0;
      tree_bk->SetBranchAddress("track_lengths",&track_lengths_bk);

      std::vector<std::vector<float> >* shower_true_E_bk =0;
      tree_bk->SetBranchAddress("true_energy",&shower_true_E_bk);

      std::vector<std::vector<float> >* shower_energy_bk=0;
      tree_bk->SetBranchAddress("shower_energy",&shower_energy_bk);


      //Get the signal branches
      std::vector<float>* vertex_recoKE_bk =0;
      tree_bk->SetBranchAddress("vertex_recoK",&vertex_recoKE_bk);
      
      std::vector<std::vector<float> >* shower_coversion_gap_bk=0;
      tree_bk->SetBranchAddress("shower_coversion_gap",&shower_coversion_gap_bk);
      
      std::vector<std::vector<float> >* shower_dEdx_bk=0;
      tree_bk->SetBranchAddress("shower_dEdx",&shower_dEdx_bk);

      std::vector<float>* inttype_length_bk=0;
      tree_bk->SetBranchAddress("nu_interaction_type",&inttype_length_bk);

      std::vector<float>* nu_reco_energy_bk=0;
      tree_bk->SetBranchAddress("nu_reco_energy",&nu_reco_energy_bk);

      std::vector<float> * nu_truth_energy_bk=0;
      tree_bk->SetBranchAddress("nu_truth_energy",&nu_truth_energy_bk);
      



      //make the histgorams 
      TH1F* track_res_hist_bk_all = new TH1F("track_res_hist_bk_all","track_res_hist_bk_all",100,-1,1);
      TH1F* proton_res_hist_bk = new TH1F("proton_res_hist_bk","proton_res_hist_bk",100,-1,1);
      TH1F* pion_res_hist_bk  = new TH1F("pion_res_hist_bk","pion_res_hist_bk",100,-1,1);
      TH1F* kaon_res_hist_bk  = new TH1F("kaon_res_hist_bk","kaon_res_hist_bk",100,-1,1);
      TH1F* muon_res_hist_bk  = new TH1F("muon_res_hist_bk","muon_res_hist_bk",100,-1,1);

      TH1F* protontrueE_hist_bk = new TH1F("proton_trueE_hist_bk","proton_trueE_hist_bk",100,0,1000);
      TH1F* piontrueE_hist_bk = new TH1F("pion_trueE_hist_bk","pion_trueE_hist_bk",100,0,1000);
      TH1F* kaontrueE_hist_bk = new TH1F("kaon_trueE_hist_bk","kaon_trueE_hist_bk",100,0,1000);
      TH1F* alltrueE_hist_bk = new TH1F("all_trueE_hist_bk","all_trueE_hist_bk",100,0,1000);
      TH1F* proton_reco_hist_bk = new TH1F("proton_reco_hist_bk","proton_reco_hist_bk",100,0,1000);
      TH1F* pion_reco_hist_bk = new TH1F("pion_reco_hist_bk","pion_reco_hist_bk",100,0,1000);
      TH1F* kaon_reco_hist_bk = new TH1F("kaon_reco_hist_bk","kaon_reco_hist_bk",100,0,1000);
      TH1F* all_reco_hist_bk = new TH1F("all_reco_hist_bk","all_reco_hist_bk",100,0,1000);
      TH1F* proton_reco_hist_bk_extra = new TH1F("proton_reco_hist_bk_extra","proton_reco_hist_bk_extra",100,0,1000);
      TH1F* pion_reco_hist_bk_extra = new TH1F("pion_reco_hist_bk_extra","pion_reco_hist_bk_extra",100,0,1000);
      TH1F* kaon_reco_hist_bk_extra = new TH1F("kaon_reco_hist_bk_extra","kaon_reco_hist_bk_extra",100,0,1000);
      TH1F* all_reco_hist_bk_extra = new TH1F("all_reco_hist_bk_extra","all_reco_hist_bk_extra",100,0,1000);
     
      TH1F* neutrino_true_energy_before_cut_bk = new TH1F("neutrino_true_energy_before_cut_bk","neutrino_true_energy_before_cut_bk",100,0,4000);
      TH1F* neutrino_true_energy_after_cut_bk = new TH1F("neutrino_true_energy_after_cut_bk","neutrino_true_energy_after_cut_bk",100,0,4000);

      TH1F* neutrino_reco_energy_before_cut_bk = new TH1F("neutrino_reco_energy_before_cut_bk","neutrino_reco_energy_before_cut_bk",100,0,4000);
      TH1F* neutrino_reco_energy_after_cut_bk = new TH1F("neutrino_reco_energy_after_cut_bk","neutrino_reco_energy_after_cut_bk",100,0,4000);
      TH1F* neutrino_reco_energy_extra_cut_bk = new TH1F("neutrino_reco_energy_extra_cut_bk","neutrino_reco_energy_extra_cut_bk",100,0,4000);

      TH1F* neutrino_res_bk = new TH1F("neutrino_res_bk","neutrino_res_bk",100,-5,5);
 
      all_energy_vec.clear();
      all_res_vec.clear();
      proton_energy_vec.clear();
      proton_res_vec.clear();
      pion_energy_vec.clear();
      pion_res_vec.clear();
      kaon_energy_vec.clear();
      kaon_res_vec.clear();
      muon_energy_vec.clear();
      muon_res_vec.clear();

      TH1F* shower_res_hist_bk = new TH1F("shower_res_hist_bk","shower_res_hist_bk",100,-1.5,1.5);
      std::vector<float> shower_res_vec_bk;
      std::vector<float> shower_true_e_vec_bk;

      std::vector<float> neutrino_res_bk_vec;
      std::vector<float> neutrino_true_energy_vec_bk;
      
      float num_neutrinos_before_cuts_bk = 0;
      float num_neutrinos_after_cuts_bk = 0;
      
      Long64_t nentries_bk = tree_bk->GetEntries();
      //Loop over the events
      for (Long64_t evt = 0; evt < nentries_bk; evt++) {

	
	tree_bk->GetEntry(evt);

	//Loop over the neutrinos.
	for(int neutrino=0; neutrino<track_E_bk->size(); ++neutrino){

	neutrino_true_energy_before_cut_bk->Fill(nu_truth_energy_bk->at(neutrino)*1000);
	neutrino_reco_energy_before_cut_bk->Fill(nu_reco_energy_bk->at(neutrino));
	++num_neutrinos_before_cuts_bk;

	  if(fApplyPreviousCuts){
	    int numshowers = 0;
	    for(int shower_iter=0; shower_iter<shower_energy_bk->at(neutrino).size(); ++shower_iter){
	      float shower_e = shower_energy_bk->at(neutrino).at(shower_iter);
	      if(shower_e > fECut){

		++numshowers;
	      }
	    }
	    if(numshowers != 1){continue;}
	 
	    //Get the largest shower.
	    int max_i = -999;
	    float max_e = -999;
	    for(int shower_iter=0; shower_iter<shower_energy_bk->at(neutrino).size(); ++shower_iter){
	      if(max_e < shower_energy_bk->at(neutrino).at(shower_iter)){
		max_i = shower_iter;
		max_e = shower_energy_bk->at(neutrino).at(shower_iter);
	      }
	    }

	    //Check the event has a vertex energy greater than the cut off and the conversion gap is less than the cut off.
	    if(shower_coversion_gap_bk->at(neutrino).at(max_i) > fConversionGapCut &&
	       vertex_recoKE_bk->at(neutrino) > fVertexCut){
	      continue;
	    }

	    //See if we pass the dEdx cut.
	    if(shower_dEdx_bk->at(neutrino).at(max_i) < fdEdxCutVal){continue;}
	    
	    //get the biggest track.
	    float max_length = -999;
	    for(int track_iter=0; track_iter<track_lengths_bk->at(neutrino).size(); ++track_iter){
	      if(max_length < track_lengths_bk->at(neutrino).at(track_iter)){
		max_length = track_lengths_bk->at(neutrino).at(track_iter);
	      }
	    }
	    if(max_length > fMaxLengthCut){continue;}

	  }

	  //Loop over the showers. 
	  int max_i = -9999;
	  float max_e = -9999;
	  for(int shower=0; shower<shower_energy_bk->at(neutrino).size(); ++shower){
	    if(shower_energy_bk->at(neutrino).at(shower) > max_e){
	      max_e = shower_energy_bk->at(neutrino).at(shower); 
	      max_i = shower;
	    }
	  }
	  if(max_i == -9999){continue;}

	  ++num_neutrinos_after_cuts_bk;
	  neutrino_true_energy_after_cut_bk->Fill(nu_truth_energy_bk->at(neutrino)*1000);
	  neutrino_reco_energy_after_cut_bk->Fill(nu_reco_energy_bk->at(neutrino));
	  neutrino_res_bk->Fill(nu_reco_energy_bk->at(neutrino)/1000*nu_truth_energy_bk->at(neutrino));
	  neutrino_res_bk_vec.push_back(nu_reco_energy_bk->at(neutrino)/1000*nu_truth_energy_bk->at(neutrino));
	  neutrino_true_energy_vec_bk.push_back(nu_truth_energy_bk->at(neutrino));
	     
	  float shower_res = (shower_true_E_bk->at(neutrino).at(max_i) - shower_energy_bk->at(neutrino).at(max_i))/shower_true_E_bk->at(neutrino).at(max_i);
	  shower_res_hist_bk->Fill(shower_res);
	  shower_res_vec_bk.push_back(shower_res);
	  shower_true_e_vec_bk.push_back(shower_true_E_bk->at(neutrino).at(max_i));


	  //Loop over the true tracks 
	  for(int track=0; track<truekaonE_bk->at(neutrino).size(); ++track){
	    //Fill true Hist_Bk
	    kaontrueE_hist_bk->Fill(truekaonE_bk->at(neutrino).at(track));
	  }

	  for(int track=0; track<trueprotonE_bk->at(neutrino).size(); ++track){
	    protontrueE_hist_bk->Fill(trueprotonE_bk->at(neutrino).at(track));
	  }

	  for(int track=0; track<truepionE_bk->at(neutrino).size(); ++track){
	    piontrueE_hist_bk->Fill(truepionE_bk->at(neutrino).at(track));
	  }

	  for(int track=0; track<truetrackE_bk->at(neutrino).size(); ++track){
	    alltrueE_hist_bk->Fill(truetrackE_bk->at(neutrino).at(track));
	  }

	  std::vector<float> truetrackE_bk_temp;
	  //Loop of the reco tracks.
	  for(int track=0; track<track_E_bk->at(neutrino).size(); ++track){

	    float true_trackE  = track_trueE_bk->at(neutrino).at(track);
	    float track_pdg    = track_track_pdg_bk->at(neutrino).at(track);
	    float reco_track_E = track_E_bk->at(neutrino).at(track);
	 
	    track_res_hist_bk_all->Fill((true_trackE-reco_track_E)/true_trackE);
	    all_energy_vec.push_back(true_trackE);
	    all_res_vec.push_back((true_trackE-reco_track_E)/true_trackE);


	    if(track_pdg == 2212){
	      proton_energy_vec.push_back(true_trackE);
	      proton_res_vec.push_back((true_trackE-reco_track_E)/true_trackE);
	      proton_res_hist_bk->Fill((true_trackE-reco_track_E)/true_trackE);
	    }
	    if(track_pdg == 211){
	      pion_energy_vec.push_back(true_trackE);
	      pion_res_vec.push_back((true_trackE-reco_track_E)/true_trackE);
	      pion_res_hist_bk->Fill((true_trackE-reco_track_E)/true_trackE);
	    }
	    if(track_pdg == 321){
	      kaon_energy_vec.push_back(true_trackE);
	      kaon_res_vec.push_back((true_trackE-reco_track_E)/true_trackE);
	      kaon_res_hist_bk->Fill((true_trackE-reco_track_E)/true_trackE);
	    }
	    if(track_pdg == 13){
	      muon_energy_vec.push_back(true_trackE);
	      muon_res_vec.push_back((true_trackE-reco_track_E)/true_trackE);
	      muon_res_hist_bk->Fill((true_trackE-reco_track_E)/true_trackE);
	    }

	    //Only considering "primaries' 
	    if(std::find(truetrackE_bk->at(neutrino).begin(),truetrackE_bk->at(neutrino).end(),true_trackE) == truetrackE_bk->at(neutrino).end()){continue;}

	    //See if the track has been matched yet
	    if(std::find(truetrackE_bk_temp.begin(),truetrackE_bk_temp.end(),true_trackE) != truetrackE_bk_temp.end()){
	      if(track_pdg == 2212){proton_reco_hist_bk_extra->Fill(true_trackE);}
	      if(track_pdg == 211){pion_reco_hist_bk_extra->Fill(true_trackE);}
	      if(track_pdg == 321){kaon_reco_hist_bk_extra->Fill(true_trackE);}
	      all_reco_hist_bk_extra->Fill(true_trackE);
	    }
	    else{
	      truetrackE_bk_temp.push_back(true_trackE);
	      if(track_pdg == 2212){proton_reco_hist_bk->Fill(true_trackE);}
	      if(track_pdg == 211){pion_reco_hist_bk->Fill(true_trackE);}
	      if(track_pdg == 321){kaon_reco_hist_bk->Fill(true_trackE);}
	      all_reco_hist_bk->Fill(true_trackE);
	    }
	  }
	}
      }

      
      //Make the bkicency plots 
      MakeEfficiencyPlots(proton_reco_hist_bk,protontrueE_hist_bk,proton_reco_hist_bk_extra,"number of reco protons","number of true protons","number of multiple reconstructions","Efficiency;Track Energy (MeV);","protons_background");
      MakeEfficiencyPlots(pion_reco_hist_bk,piontrueE_hist_bk,pion_reco_hist_bk_extra,"number of reco pions","number of true pions","number of multiple reconstructions","Efficiency;Track Energy (MeV);","pions_background");
      MakeEfficiencyPlots(kaon_reco_hist_bk,kaontrueE_hist_bk,kaon_reco_hist_bk_extra,"number of reco kaons","number of true kaons","number of multiple reconstructions","Efficiency;Track Energy (MeV);","kaons_background");
      MakeEfficiencyPlots(all_reco_hist_bk,alltrueE_hist_bk,all_reco_hist_bk_extra,"number of reco tracks","number of true tracks","number of multiple reconstructions","Efficiency;Track Energy (MeV);","all_background");
    
      //Make the resolution graphs.
      MakeResolutionGraphs(proton_res_vec,proton_energy_vec,0,500,10,50,"True Energy-Reco Energy/TrueEnergy","track_resE","protons background",50);
      MakeResolutionGraphs(pion_res_vec,pion_energy_vec,0,500,10,50,"True Energy-Reco Energy/TrueEnergy","track_resE","pions background",50);
      MakeResolutionGraphs(kaon_res_vec,kaon_energy_vec,0,500,10,50,"True Energy-Reco Energ/TrueEnergy","track_resE","kaons background",50);
      MakeResolutionGraphs(muon_res_vec,muon_energy_vec,0,1000,10,50,"True Energy-Reco Energ/TrueEnergy","track_resE","muons background",50);
      MakeResolutionGraphs(all_res_vec,all_energy_vec,0,500,10,50,"True Energy-Reco Energ/TrueEnergy","track_resE","all background",50);
      

      shower_res_hist_bk->Write();
      MakeResolutionGraphs(shower_res_vec_bk,shower_true_e_vec_bk,0,1500,10,50,"True Energy- Reco Energ/TrueEnergy","shower_resE","Shower background",50);

      track_res_hist_bk_all->Write();
      proton_reco_hist_bk->Write();
      pion_reco_hist_bk->Write();
      kaon_reco_hist_bk->Write();
      muon_res_hist_bk->Write();

      neutrino_res_eff->Write();
      MakeEfficiencyPlots(neutrino_true_energy_after_cut_bk,neutrino_true_energy_before_cut_bk,neutrino_reco_energy_extra_cut_bk,"number neutrinos before selection","number neutrinos after selection true","number of neutrinos after selection extra true","Efficiency;Neutrino Energy (MeV);","muon neutrino true");
      MakeEfficiencyPlots(neutrino_reco_energy_after_cut_bk,neutrino_reco_energy_before_cut_bk,neutrino_reco_energy_extra_cut_bk,"number neutrinos before selection reco","number neutrinos after selection reco","number of neutrinos after selection extra reco","Efficiency;Neutrino Energy (MeV);","muon neutrino reco");

      MakeResolutionGraphs(neutrino_res_bk_vec,neutrino_true_energy_vec_bk,0,4000,10,100,"Reco Energy/True Energy","neutrino_resE","muon neutrino res",50);


      std::cout << "total signal efficiency: "<< num_neutrinos_after_cuts_eff/num_neutrinos_before_cuts_eff << std::endl;
      std::cout << "total background rejection: " << num_neutrinos_after_cuts_bk/num_neutrinos_before_cuts_bk << std::endl;

    }

}





void DoCutFinding(TH1D* signalhist, TH1D* backgroundhist, TString units, TString axisTitle)
{

  // TH1::AddDirectory(kFALSE);
  signalhist->Scale(1 / signalhist->GetEntries());
  backgroundhist->Scale(1 /backgroundhist->GetEntries());
  
  //Find the min and max start positions
  float minsig = signalhist->GetBinCenter(0);
  float maxsig = signalhist->GetBinCenter(signalhist->GetNbinsX());
  float minbk = backgroundhist->GetBinCenter(0);
  float maxbk = backgroundhist->GetBinCenter(backgroundhist->GetNbinsX());

  if (minsig != minbk) {
    std::cout << "min point does not match" << std::endl;
    return;
  }
  if (maxbk != maxbk) {
    std::cout << "max point does not match" << std::endl;
    return;
  }
  if (signalhist->GetNbinsX() != backgroundhist->GetNbinsX()) {
    std::cout << "bin numbers sig " << signalhist->GetNbinsX() << " and background " << backgroundhist->GetNbinsX() << "do not match" << std::endl;
    return;
  }

  int max_ind = -999;
  int maxsigbk_ind = -999;
  int n = signalhist->GetNbinsX();
  double bktot_err = 0;
  double bktot = backgroundhist->IntegralAndError(1, n+1, bktot_err);
  double sigtot_err = 0;
  double sigtot = signalhist->IntegralAndError(1, n+1, sigtot_err);
  double maxeffpur = -999;
  double maxsigbk = -999;
  double efficiency[n];
  double efficencyErr[n];
  double purity[n];
  double purityErr[n];
  double bkRej[n];
  double bkRejErr[n];
  double effpur[n];
  double effpurErr[n];
  double sigbk[n];
  double sigbkErr[n];
  double xval[n];
  double xvalErr[n];

  //Get the values.
  for (int i = 1; i < n+1; ++i) {

    double presig_err = 0;
    double prebk_err = 0;
    double presig = signalhist->IntegralAndError(1, i, presig_err);
    double prebk = backgroundhist->IntegralAndError(1, i, prebk_err);

    //std::cout << "presig: " << presig << " prebk: " << prebk << std::endl;

    efficiency[i] = presig / sigtot;
    efficencyErr[i] = efficiency[i] * TMath::Sqrt(((presig_err / presig) * (presig_err / presig)) + (sigtot_err / sigtot) * (sigtot_err / sigtot));

    bkRej[i] = 1 - prebk / bktot;
    bkRejErr[i] = bkRej[i] * TMath::Sqrt((prebk_err / prebk) * (prebk_err / prebk) + (bktot_err / bktot) * (bktot_err / bktot));

    if (presig + prebk != 0) {
      purity[i] = presig / (presig + prebk);
      double purityErr_dem = TMath::Sqrt((presig_err * presig_err) + (prebk_err * prebk_err));
      purityErr[i] = purity[i] * TMath::Sqrt(((presig_err / presig) * (presig_err / presig)) + (purityErr_dem / (presig + prebk)) * (purityErr_dem / (presig + prebk)));
    } else {
      purity[i] = 0;
      purityErr[i] = 0;
    }

    effpur[i] = efficiency[i] * purity[i];
    if (efficiency[i] == 0 || purity[i] == 0) {
      effpurErr[i] = 0;
    } else {
      effpurErr[i] = effpur[i] * TMath::Sqrt(((efficencyErr[i] / efficiency[i]) * (efficencyErr[i] / efficiency[i])) + (purityErr[i] / purity[i]) * (purityErr[i] / purity[i]));
    }

    sigbk[i] = efficiency[i] * bkRej[i];
    if (efficiency[i] == 0 || bkRej[i] == 0) {
      sigbkErr[i] = 0;
    } else {
      sigbkErr[i] = sigbk[i] * TMath::Sqrt(((efficencyErr[i] / efficiency[i]) * (efficencyErr[i] / efficiency[i])) + (bkRejErr[i] / bkRej[i]) * (bkRejErr[i] / bkRej[i]));
    }

    xval[i] = signalhist->GetBinCenter(i);
    xvalErr[i] = signalhist->GetBinWidth(i);

    //    std::cout << "xval " << xval[i] << " efficiency: " << efficiency[i] << " purity: " << purity[i] << " eff * pur "
    //          << effpur[i] << "background Rejection: " << bkRej[i] << " and sig*bkrej: " << sigbk[i] << std::endl;

    //Keep the largest value.
    if (effpur[i] > maxeffpur) {
      max_ind = i;
      maxeffpur = effpur[i];
    }
    if (sigbk[i] > maxsigbk) {
      maxsigbk = sigbk[i];
      maxsigbk_ind = i;
    }
  }

  //Draw the curves
  auto effpur_canvas = new TCanvas("effpur_canvas", "effpur_canvas", 600, 400);
  auto sigbk_canvas = new TCanvas("sigbk_canvas", "sigbk_canvas", 600, 400);

  TGraphErrors* efficiencyPlot = new TGraphErrors(n, xval, efficiency, xvalErr, efficencyErr);
  efficiencyPlot->GetYaxis()->SetTitle("Efficency");
  efficiencyPlot->SetTitle("Efficency");
  efficiencyPlot->SetLineColor(kBlue);
  efficiencyPlot->Draw("a4");
  efficiencyPlot->Write();

  TGraphErrors* bkRejPlot = new TGraphErrors(n, xval, bkRej, xvalErr, bkRejErr);
  bkRejPlot->GetYaxis()->SetTitle("Background Rejection");
  bkRejPlot->SetTitle("Background Rejection");
  bkRejPlot->SetLineColor(kRed);
  bkRejPlot->Draw("a4");
  bkRejPlot->Write();

  TGraphErrors* purityPlot = new TGraphErrors(n, xval, purity, xvalErr, purityErr);
  purityPlot->GetYaxis()->SetTitle("Purity");
  purityPlot->SetTitle("Purity");
  purityPlot->SetLineColor(kRed);
  purityPlot->Draw("a4");
  purityPlot->Write();

  TGraphErrors* effxpurPlot = new TGraphErrors(n, xval, effpur, xvalErr, effpurErr);
  effxpurPlot->GetYaxis()->SetTitle("Effxpur");
  effxpurPlot->SetTitle("Eff*pur");
  effxpurPlot->SetLineColor(kGreen);
  effxpurPlot->Draw("a4");
  effxpurPlot->Write();

  TGraphErrors* sigbkPlot = new TGraphErrors(n, xval, sigbk, xvalErr, sigbkErr);
  sigbkPlot->GetYaxis()->SetTitle("Signal*Background Rejection");
  sigbkPlot->SetTitle("Signal*Background Rejection");
  sigbkPlot->SetLineColor(kGreen);
  sigbkPlot->Draw("a4");
  sigbkPlot->Write();

  //Make the graph to present
  effpur_canvas->cd();
  auto mg1 = new TMultiGraph("mg1", "mg1");
  mg1->Add(efficiencyPlot);
  mg1->Add(purityPlot);
  mg1->Add(effxpurPlot);
  mg1->GetYaxis()->SetRangeUser(0, 1);
  mg1->Write();
  mg1->Draw("AP");
  effpur_canvas->BuildLegend();
  effpur_canvas->Write();

  //Make the graph to present
  sigbk_canvas->cd();
  auto mg2 = new TMultiGraph("mg2", "mg2");
  mg2->Add(efficiencyPlot);
  mg2->Add(bkRejPlot);
  mg2->Add(sigbkPlot);
  mg2->GetYaxis()->SetRangeUser(0, 1);
  mg2->Write();

  mg2->Draw("AP");
  sigbk_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
  sigbk_canvas->Write();

  double maxScale = TMath::Max(backgroundhist->GetBinContent(backgroundhist->GetMaximumBin()),
			       signalhist->GetBinContent(signalhist->GetMaximumBin()));

  signalhist->SetLineColor(kBlue);
  signalhist->SetFillStyle(3003);
  signalhist->SetFillColor(6);
  signalhist->Scale(1. / maxScale);
  signalhist->SetTitle("Signal");
  signalhist->Write();
  backgroundhist->SetLineColor(kRed);
  backgroundhist->SetFillStyle(3003);
  backgroundhist->SetFillColor(42);
  backgroundhist->Scale(1. / maxScale);
  backgroundhist->SetTitle("Background");
  backgroundhist->Write();
  
  auto sigbk_cut_canvas = new TCanvas("sigbk_cut_canvas", "sigbk_cut_canvas", 600, 400);
  sigbk_cut_canvas->cd();
  mg2->GetHistogram()->SetTitle(Form("Cut at %0.1f " + units + " with efficiency %0.2f and background rejection %0.2f", xval[maxsigbk_ind], efficiency[maxsigbk_ind], bkRej[maxsigbk_ind]));
  mg2->GetXaxis()->SetTitle(axisTitle);
  mg2->Draw("AP");
  signalhist->Draw("HIST SAME");
  backgroundhist->Draw("HIST SAME");
  sigbk_cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
  TLine* l1 = new TLine(xval[maxsigbk_ind], 0, xval[maxsigbk_ind], 1);
  l1->SetLineColor(kBlack);
  l1->Draw("same");

  sigbk_cut_canvas->Write();

  //Draw the overlayed samples with the cut.
  auto cut_canvas = new TCanvas("cut_canvas", "cut_canvas", 600, 400);
  cut_canvas->cd();
  mg1->GetHistogram()->SetTitle(Form("Cut at %0.1f " + units + " with efficiency %0.2f and purity %0.2f", xval[max_ind], efficiency[max_ind], purity[max_ind]));
  mg1->GetXaxis()->SetTitle(axisTitle);
  mg1->Draw("AP");
  signalhist->Draw("HIST SAME");
  backgroundhist->Draw("HIST SAME");
  cut_canvas->Update();
  cut_canvas->BuildLegend();

  TLine* l = new TLine(xval[max_ind], 0, xval[max_ind], 1);
  l->SetLineColor(kBlack);
  l->Draw("same");

  cut_canvas->Write();

  std::cout << "Best Cut is at: " << xval[max_ind] << std::endl;
  std::cout << "Efficiency: " << efficiency[max_ind] << " +- " << efficencyErr[max_ind] << std::endl;
  std::cout << "Purity    : " << purity[max_ind] << " +- " << purityErr[max_ind] << std::endl;
  std::cout << "Eff*Purity: " << effpur[max_ind] << " +- " << effpurErr[max_ind] << std::endl;

  std::cout << "Best sig*bkrej Cut is at: " << xval[maxsigbk_ind] << std::endl;
  std::cout << "Efficiency:     " << efficiency[maxsigbk_ind] << " +- " << efficencyErr[maxsigbk_ind] << std::endl;
  std::cout << "Background Rej: " << bkRej[maxsigbk_ind] << " +- " << bkRejErr[maxsigbk_ind] << std::endl;
  std::cout << "Eff*Purity:     " << sigbk[maxsigbk_ind] << " +- " << sigbkErr[maxsigbk_ind] << std::endl;

  delete efficiencyPlot;
  delete bkRejPlot;
  delete purityPlot;
  delete effxpurPlot;
  delete sigbkPlot;

  delete mg1;
  delete mg2;

  delete effpur_canvas;
  delete sigbk_canvas;
  delete sigbk_cut_canvas;
  delete cut_canvas;
}


void MakeEfficiencyPlots(TH1F* postcut, TH1F* precut, TH1F* extra, TString postcutname, TString precutname, TString extraname, TString title, TString eff_name){

  float maxSignal = precut->GetBinContent(precut->GetMaximumBin());
  
  if(TEfficiency::CheckConsistency(*postcut,*precut)){
    //Make the eff plot.
    TEfficiency* reco_eff = new TEfficiency(*postcut,*precut);
    reco_eff->SetTitle(title);
    reco_eff->SetName(eff_name);
    reco_eff->Write();

    postcut->Scale(1./maxSignal);
    precut->Scale(1./maxSignal);
      
    postcut->SetLineColor(kBlue);
    postcut->SetFillStyle(3003);
    postcut->SetFillColor(6);
    postcut->SetTitle(postcutname);
    postcut->Write();
    precut->SetLineColor(kRed);
    precut->SetFillStyle(3003);
    precut->SetFillColor(42);
    precut->SetTitle(precutname);
    precut->Write();

    TString canavs_name = eff_name + "_canvas";
    auto cut_canvas = new TCanvas(canavs_name, canavs_name, 600, 400);
    precut->Draw("HIST SAME");
    cut_canvas->cd();
    reco_eff->Draw("same");
    postcut->Draw("HIST SAME");
    cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
    cut_canvas->Write();
    
      
    postcut->Scale(maxSignal);
    precut->Scale(maxSignal);
  }

  return;

  if(TEfficiency::CheckConsistency(*extra,*precut)){
    //Make the eff plot.
    TEfficiency* reco_eff = new TEfficiency(*extra,*precut);
    reco_eff->SetTitle(title);
    reco_eff->SetName(eff_name);
    reco_eff->Write();

    extra->Scale(1./maxSignal);
    precut->Scale(1./maxSignal);
      
    extra->SetLineColor(kBlue);
    extra->SetFillStyle(3003);
    extra->SetFillColor(6);
    extra->SetTitle(extraname);
    extra->Write();
    precut->SetLineColor(kRed);
    precut->SetFillStyle(3003);
    precut->SetFillColor(42);
    precut->SetTitle(precutname);
    precut->Write();

    TString canavs_name = eff_name + "_canvas_extra";
    auto cut_canvas = new TCanvas(canavs_name, canavs_name, 600, 400);
    precut->Draw("HIST SAME");
    cut_canvas->cd();
    reco_eff->Draw("same");
    extra->Draw("HIST SAME");
    cut_canvas->BuildLegend(0.6, 0.4, 0.9, 0.6);
    cut_canvas->Write();
    
      
    extra->Scale(maxSignal);
    precut->Scale(maxSignal);
  }



}

void MakeResolutionGraphs(std::vector<float>& vals, std::vector<float>& energyval, float min, float max, float num_steps, float width,TString yTitle, TString Metric, TString name, float nbins){


  float step_size = (max-min)/(num_steps);

  float x[(int) num_steps];
  float dx[(int)num_steps];
  float y[(int)num_steps];
  float dy[(int)num_steps];
  float ddy[(int)num_steps];
  float dddy[(int)num_steps];

  TString TrueEnergy = "track_trueE";

  //Loop over the energy step.
  for(int i=0; i<num_steps; ++i){
    
    //Get the widths to access 
    float energylow = min + i*step_size - width/2;
    float energymax = min + i*step_size + width/2;

    TString tempHistname = Metric + "tempHist " + std::to_string(energylow);  
    TH1F* tempHist = new TH1F(tempHistname, tempHistname, nbins, -2, 2);
    for(int val=0; val<vals.size(); ++val){
      //      std::cout << "energyval[val]: " << energyval[val] << " energylow: " << energylow << " energymax: " << energymax  << std::endl;
      if(energyval[val] > energylow && energyval[val] < energymax){
	tempHist->Fill(vals[val]);
      }
    }
    x[i] = (energylow + energymax) / 2.0;
    dx[i] = (energymax - energylow) / 2.0;
    y[i] = tempHist->GetMean();
    dy[i] = tempHist->GetRMS();
    ddy[i] = tempHist->GetRMS();
    dddy[i] = dy[i]/TMath::Sqrt(2*(tempHist->Integral()-1));

    // if(tempHist->GetEntries() > 20){
    //   TF1* fit = new TF1("fit","gaus");
    //   tempHist->Fit("fit");
    //   y[i]  = fit->GetParameter(1);
    //   dy[i] = fit->GetParError(1);
    //   ddy[i]  = fit->GetParameter(2);
    //   dddy[i] = fit->GetParError(2);
      
    //   tempHist->Clear();
    //   for(int val=0; val<vals.size(); ++val){
    // 	//      std::cout << "energyval[val]: " << energyval[val] << " energylow: " << energylow << " energymax: " << energymax  << std::endl;
    // 	if(energyval[val] > energylow && energyval[val] < energymax){
    // 	  if(TMath::Abs(vals[val]-y[i]) <  ddy[i]){ 
    // 	    tempHist->Fill(vals[val]);
    // 	  }
    // 	}
    //   }
    //   delete fit;
    // }

    // if(tempHist->GetEntries() > 20){
    //   TF1* fit1 = new TF1("fit1","gaus");
    //   tempHist->Fit("fit1");
    //   y[i]  = fit1->GetParameter(1);
    //   dy[i] = fit1->GetParError(1);
    //   ddy[i]  = fit1->GetParameter(2);
    //   dddy[i] = fit1->GetParError(2);
    //   tempHist->Write();
    //   delete fit1;
    // }
    
    //    y[i] = tempHist->GetMean();
    //dy[i] = tempHist->GetRMS();
    //ddy[i] = dy[i]/TMath::Sqrt(2*(tempHist->Integral()-1));
    delete tempHist;

    std::cout << "x: " << x[i] << " y: " << y[i]  << " dy: " << dy[i] << " ddy: " << ddy[i] << " dddy[i]: " << dddy[i] << std::endl;
  }


  TString yTitle_mean = yTitle + " mean";
  TString mean_name = yTitle_mean + " " + name;
  TGraphErrors* tempGraph = new TGraphErrors(num_steps, x, y, dx, dy);
  //  tempGraphs->SetLineColor();
  //  tempGraphs->SetMarkerColor(pass + 1);
  tempGraph->GetXaxis()->SetTitle("True Track Energy (MeV)");
  tempGraph->GetYaxis()->SetTitle(yTitle_mean);
  tempGraph->SetTitle(mean_name);
  tempGraph->Write();
 
  TString yTitle_rms = yTitle + " rms"; 
  TString rms_name = yTitle_rms + " " + name;
  TGraphErrors* tempRMSGraph = new TGraphErrors(num_steps, x, ddy, dx, dddy);
  //  tempRMSGraphs->SetLineColor(pass + 1);
  //tempRMSGraphs->SetMarkerColor(pass + 1);
  tempRMSGraph->GetXaxis()->SetTitle("True True Energy (MeV)");
  tempRMSGraph->GetYaxis()->SetTitle(yTitle_rms);
  tempRMSGraph->SetTitle(rms_name);
  

  TString yTitle_mean_canvas = yTitle + " mean_canvas" + name;    
  TCanvas* cMeanGraph = new TCanvas(yTitle_mean_canvas);
  tempGraph->Draw();
  cMeanGraph->Write();

  TString yTitle_rms_canvas = yTitle + " rms_canvas" + name;    
  TCanvas* cRMSGraph = new TCanvas(yTitle_rms_canvas);
  tempRMSGraph->Draw();
  cRMSGraph->Write();
}




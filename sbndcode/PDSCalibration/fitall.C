void fitall() {
  cout << "======" << "SPE ANALYSIS" << "======" << endl << "(all channels)" << endl << "abullock and hollyp for SBND, 2023-2024." << endl << "Launching..." << endl;
  

  //FIT FUNCTION
  TF1 *func = new TF1("fit", "[2]*TMath::Gaus(x,[0],[1]) + [3]*TMath::Gaus(x,2*[0],sqrt(2)*[1]) + [4]*TMath::Gaus(x,3*[0],sqrt(3)*[1]) + [5]*TMath::Exp(-[6]*x)"); //spe integral fit function
  func->SetParNames("Mu", "Sigma", "Peak1", "Peak2", "Peak3", "Pedestal", "Tau"); //name the variables
  func->SetParLimits(0, 0, 500); //limits on mu: this is the value we want
  func->SetParLimits(1, 0, 500); //limits on sigma
  func->SetParLimits(2, 0, 200); func->SetParLimits(3, 0, 200); func->SetParLimits(4, 0, 200); func->SetParLimits(5, 0, 300); //heights: 5 is the issue, wants to be 285
  func->SetParLimits(6, 0, 100); //limits on tau

  //TTREE
  TFile *save_the_trees = new TFile("SPE_Fits.root", "UPDATE");
  TTree *t = (TTree*) save_the_trees->Get("SPE_Analysis"); // get the ttree out of the file
  //branches
  Double_t opc, amp, ampsigma, integ0, integ0sigma, integ1, integ1sigma, integ2, integ2sigma, integ3, integ3sigma, integ4, integ4sigma, integ5, integ5sigma; //gain
  t->SetBranchAddress("OpChannel", &opc); 
  t->SetBranchAddress("Amp", &amp); t->SetBranchAddress("AmpSigma", &ampsigma);
  t->SetBranchAddress("IntegZero", &integ0); t->SetBranchAddress("IntegZeroSigma", &integ0sigma);
  t->SetBranchAddress("IntegThresh", &integ1); t->SetBranchAddress("IntegThreshSigma", &integ1sigma);
  t->SetBranchAddress("IntegManu", &integ2); t->SetBranchAddress("IntegManuSigma", &integ2sigma);
  t->SetBranchAddress("IntegZeroB", &integ3); t->SetBranchAddress("IntegZeroBSigma", &integ3sigma);
  t->SetBranchAddress("IntegThreshB", &integ4); t->SetBranchAddress("IntegThreshBSigma", &integ4sigma);
  t->SetBranchAddress("IntegManuB", &integ5); t->SetBranchAddress("IntegManuBSigma", &integ5sigma);

  //ADDITIONAL TTREES
  TTree *er = (TTree*) save_the_trees->Get("Errors"); // get the ttree out of the file
  TTree *ch = (TTree*) save_the_trees->Get("ChiSquared"); // get the ttree out of the file
  Double_t amp_er, integ0_er, integ1_er, integ2_er, integ3_er, integ4_er, integ5_er;
  Double_t amp_ch, integ0_ch, integ1_ch, integ2_ch, integ3_ch, integ4_ch, integ5_ch;
  er->SetBranchAddress("OpChannel", &opc); 
  ch->SetBranchAddress("OpChannel", &opc); 
  er->SetBranchAddress("AmpError", &amp_er); ch->SetBranchAddress("AmpChisq", &amp_ch);
  er->SetBranchAddress("IntegZeroError", &integ0_er); ch->SetBranchAddress("IntegZeroChisq", &integ0_ch);
  er->SetBranchAddress("IntegThreshError", &integ1_er); ch->SetBranchAddress("IntegThreshChisq", &integ1_ch);
  er->SetBranchAddress("IntegManuError", &integ2_er); ch->SetBranchAddress("IntegManuChisq", &integ2_ch);
  er->SetBranchAddress("IntegZeroBError", &integ3_er); ch->SetBranchAddress("IntegZeroBChisq", &integ3_ch);
  er->SetBranchAddress("IntegThreshBError", &integ4_er); ch->SetBranchAddress("IntegThreshBChisq", &integ4_ch);
  er->SetBranchAddress("IntegManuBError", &integ5_er); ch->SetBranchAddress("IntegManuBChisq", &integ5_ch);

  //LOOP OVER SPEANA
  Int_t channels[20] = {
    7, 17, 295, 305,
    89, 91, 221, 223,
    6, 16, 294, 304,
    88, 90, 220, 222,
    117, 195, 116, 194
  };  //channels to loop over: NEED TO EDIT THIS. CAN I SAVE A LIST WHEN I RUN THE FIRST BIT?
  int tests=0, success=0, failed=0;
  TFile *speana = new TFile("pmt_gain.root", "READ"); //open the big file that the LArSoft module creates with SPE histograms in
  for (int i=0; i<1; i++) { //go back to 20: EDIT
    //TFile *speana = new TFile(Form("ser_opchannel_%d.root", channels[i]), "UPDATE"); //get results of speana, which saves the histograms in a TFile and returns it
    opc = channels[i]; //gain=Gain;
    cout << "Fitting histograms for channel " << channels[i] << "..." << endl;

    //HISTOGRAMS
    TH1D *h[7] = {
      (TH1D*) speana->Get(Form("PMTGain/amp_opchannel_%d", channels[i])),
      (TH1D*) speana->Get(Form("PMTGain/integ_opchannel_%d_zeromode", channels[i])),
      (TH1D*) speana->Get(Form("PMTGain/integ_opchannel_%d_threshmode", channels[i])),
      (TH1D*) speana->Get(Form("PMTGain/integ_opchannel_%d_manualmode", channels[i])),
      (TH1D*) speana->Get(Form("PMTGain/integ_opchannel_%d_zeromodeB", channels[i])),
      (TH1D*) speana->Get(Form("PMTGain/integ_opchannel_%d_threshmodeB", channels[i])),
      (TH1D*) speana->Get(Form("PMTGain/integ_opchannel_%d_manualmodeB", channels[i])),
    };

    //FITTING
    Double_t mu, sigma, error, chisq;
    for (int j=0; j<7; j++) { //go back to 7
      bool fail=false;
      func->SetParameters(h[j]->GetBinCenter(h[j]->GetMaximumBin()), h[j]->GetBinCenter(h[j]->GetMaximumBin())/2, 50, 25, 12, h[j]->GetBinContent(1), h[j]->GetBinCenter(1)); //initial guesses
      // should def be looking for mu by finding that first peak height
      cout << h[j]->GetBinCenter(h[j]->GetMaximumBin()) << " = first peak height?" << endl;
      for (int i = 0; i < 7; i++) {
        cout << "Parameter " << i << ": " << func->GetParameter(i) << endl;
    }


      h[j]->Fit(func, "Q",0); //run the fit with no print and no draw; taken off ,0

      TCanvas *canvas = new TCanvas("fitCanvas", "Fit Result");
      h[j]->Draw();         // Draw the original histogram
      func->Draw("same");    // Draw the fit on top of the histogram
      canvas->Update();      // Update the canvas

    //Save the canvas with a unique filename based on the loop index
      TString filename = TString::Format("fitCanvas_%d.png", j);
      canvas->Print(filename);
      delete canvas;  // Delete the canvas after saving
      
      
      TF1 *fit = h[j]->GetFunction(func->GetName()); //get fitted function
      mu = fit->GetParameter(0); sigma = fit->GetParameter(1); //get fit parameters, and assign them to correct TTree variables
      chisq = fit->GetChisquare(); error = fit->GetParError(0); //get chi squared and mu error
      cout << "peak  =   " << mu << " for hist " << j << endl;
      if (mu==1000 || mu==0) {mu=-1; fail=true;} //problem with mu
      if (sigma==500 || sigma==0) {sigma=-1; fail=true;} //problem with sigma
      if (chisq>200) {mu=-1; sigma=-1; fail=true;} //chisq too big
      //assign to ttree variables
      if (j==0) {amp=mu; ampsigma=sigma; amp_er=error; amp_ch=chisq;}
      else if (j==1) {integ0=mu; integ0sigma=sigma; integ0_er=error, integ0_ch=chisq;}
      else if (j==2) {integ1=mu; integ1sigma=sigma; integ1_er=error, integ1_ch=chisq;}
      else if (j==3) {integ2=mu; integ2sigma=sigma; integ2_er=error, integ2_ch=chisq;}
      else if (j==4) {integ3=mu; integ3sigma=sigma; integ3_er=error, integ3_ch=chisq;}
      else if (j==5) {integ4=mu; integ4sigma=sigma; integ4_er=error, integ4_ch=chisq;}
      else if (j==6) {integ5=mu; integ5sigma=sigma; integ5_er=error, integ5_ch=chisq;}
      tests++; if (fail) {failed++;} else {success++;}
    }

    //ADD TO TREE
    cout << "  Fits successful." << endl;
    t->Fill(); er->Fill(); ch->Fill();
  }

  //SAVE THE TTREES
  cout << "======" << endl;
  cout << "Of " << tests << " fits performed, " << success << " were successful and " << failed << " were not." << endl;
  cout << "Saving TTrees..." << endl;
  save_the_trees->WriteObject(t,"SPE_Analysis");
  save_the_trees->WriteObject(er,"Errors");
  save_the_trees->WriteObject(ch,"ChiSquared");
  cout << "TTrees saved!" << endl;
  save_the_trees->Close();
}

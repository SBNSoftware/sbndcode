//This macro was written by abullock for SBND, August 2023.
//Create a TTree to store fit data in.

bool ttreemaker() {

  TFile *save_the_trees = new TFile("SPE_Fits.root", "RECREATE");
  Double_t opc, amp, ampsigma, integ0, integ0sigma, integ1, integ1sigma, integ2, integ2sigma, integ3, integ3sigma, integ4, integ4sigma, integ5, integ5sigma; //gain
  Double_t amp_er, integ0_er, integ1_er, integ2_er, integ3_er, integ4_er, integ5_er;
  Double_t amp_ch, integ0_ch, integ1_ch, integ2_ch, integ3_ch, integ4_ch, integ5_ch; 
  TTree *t = new TTree("SPE_Analysis","SPE_Analysis");
  t->Branch("OpChannel", &opc); 
  //t->Branch("Gain", &gain); 
  t->Branch("Amp", &amp); t->Branch("AmpSigma", &ampsigma);
  t->Branch("IntegZero", &integ0); t->Branch("IntegZeroSigma", &integ0sigma);
  t->Branch("IntegThresh", &integ1); t->Branch("IntegThreshSigma", &integ1sigma);
  t->Branch("IntegManu", &integ2); t->Branch("IntegManuSigma", &integ2sigma);
  t->Branch("IntegZeroB", &integ3); t->Branch("IntegZeroBSigma", &integ3sigma);
  t->Branch("IntegThreshB", &integ4); t->Branch("IntegThreshBSigma", &integ4sigma);
  t->Branch("IntegManuB", &integ5); t->Branch("IntegManuBSigma", &integ5sigma);
  save_the_trees->WriteObject(t, "SPE_Analysis");
  TTree *er = new TTree("Errors","Errors");
  TTree *ch = new TTree("ChiSquared","ChiSquared");
  er->Branch("OpChannel", &opc); //er->Branch("Gain", &gain); 
  ch->Branch("OpChannel", &opc); //ch->Branch("Gain", &gain);
  er->Branch("AmpError", &amp_er); ch->Branch("AmpChisq", &amp_ch);
  er->Branch("IntegZeroError", &integ0_er); ch->Branch("IntegZeroChisq", &integ0_ch);
  er->Branch("IntegThreshError", &integ1_er); ch->Branch("IntegThreshChisq", &integ1_ch);
  er->Branch("IntegManuError", &integ2_er); ch->Branch("IntegManuChisq", &integ2_ch);
  er->Branch("IntegZeroBError", &integ3_er); ch->Branch("IntegZeroBChisq", &integ3_ch);
  er->Branch("IntegThreshBError", &integ4_er); ch->Branch("IntegThreshBChisq", &integ4_ch);
  er->Branch("IntegManuBError", &integ5_er); ch->Branch("IntegManuBChisq", &integ5_ch); 
  save_the_trees->WriteObject(er, "Errors");
  save_the_trees->WriteObject(ch, "ChiSquared");
  return true;

}

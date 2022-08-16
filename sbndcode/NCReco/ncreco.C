#define ncreco_cxx
#include "ncreco.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <cmath>

double fwhm(TH1* h);
double Median(vector<double>* v);
void makeplot(TH1* h);
void makeplot2d(TH2* h);
void makeplot(TH1* h, vector<double> v);

void ncreco::Loop(){

  double Eh; //hadron E
  double ENuPtcl;
  double Ph; //hadron momentum
  double Phx;
  double Phy;
  double Phz;

  TFile* f = new TFile("ncreco_histos.root", "RECREATE");

  TH1* NCSpectrum = new TH1D("NCSpectrum", "True NC Neutrino Spectrum NC1p;Energy (GeV);", 60, 0, 3); //Spectrum for all NC events, regardless of topology

  TH1* ETrue = new TH1D("ETrue", "True Energy NC1p;Energy (GeV);", 60, 0, 3);
  TH1* EKin = new TH1D("EKin", "Kinematic Energy NC1p;Energy (GeV)", 60, 0, 3);
  TH1* ECal = new TH1D("ECal", "Calorimetric Energy NC1p;Energy (GeV);", 60, 0, 3);

  TH2* KinVsTrue = new TH2D("KinVsTrue", "Kinematic Reco E vs True E NC1p;True Energy (GeV); Kin Reco (GeV)", 30, 0, 3, 30, 0, 3);
  TH2* CalVsTrue = new TH2D("CalVsTrue", "Calorimetric Reco E vs True E NC1p;True Energy (GeV); Cal Reco (GeV)", 30, 0, 3, 30, 0, 3);

  TH1* KinRes = new TH1D("KinRes", "Kinematic resolution NC1p;ETrue-EKin;", 40, -1.25, 2.25); 
  TH1* CalRes = new TH1D("CalRes", "Calorimetric resolution NC1p;ETrue-ECal;", 40, -1.25, 2.25);

  TH1* KinFracRes = new TH1D("KinFracRes", "Kinematic Fractional Resolution NC1p;(ETrue-EKin)/ETrue;", 40, -1, 1); 
  TH1* CalFracRes = new TH1D("CalFracRes", "Calorimetric Fractional Resolution NC1p;(ETrue-ECal)/ETrue;", 40, -1, 1);

  TEfficiency* Eff = new TEfficiency("Eff", "Efficiency NC1p;ETrue (GeV);Eff", 50, 0 ,5); //Efficiency of our cuts


  TH1* NpETrue = new TH1D("NpETrue", "True Energy NCNp;Energy (GeV);", 60, 0, 3);
  TH1* NpEKin = new TH1D("NpEKin", "Kinematic Energy NCNp;Energy (GeV)", 60, 0, 3);
  TH1* NpECal = new TH1D("NpECal", "Calorimetric Energy NCNp;Energy (GeV);", 60, 0, 3);

  TH2* NpKinVsTrue = new TH2D("NpKinVsTrue", "Kinematic Reco E vs True E NCNp;True Energy (GeV); Kin Reco (GeV)", 30, 0, 3, 30, 0, 3);
  TH2* NpCalVsTrue = new TH2D("NpCalVsTrue", "Calorimetric Reco E vs True E NCNp;True Energy (GeV); Cal Reco (GeV)", 30, 0, 3, 30, 0, 3);

  TH1* NpKinRes = new TH1D("NpKinRes", "Kinematic resolution NCNp;ETrue-EKin;", 40, -1.25, 2.25); 
  TH1* NpCalRes = new TH1D("NpCalRes", "Calorimetric resolution NCNp;ETrue-ECal;", 40, -1.25, 2.25);

  TH1* NpKinFracRes = new TH1D("NpKinFracRes", "Kinematic Fractional Resolution NCNp;(ETrue-EKin)/ETrue;", 40, -1, 1); 
  TH1* NpCalFracRes = new TH1D("NpCalFracRes", "Calorimetric Fractional Resolution NCNp;(ETrue-ECal)/ETrue;", 40, -1, 1);

  TEfficiency* NpEff = new TEfficiency("NpEff", "Efficiency NCNp;ETrue (GeV);Eff", 50, 0 ,5); //Efficiency of our cuts


  TH1* IncETrue = new TH1D("IncETrue", "True Energy Inclusive;Energy (GeV);", 60, 0, 3);
  TH1* IncEKin = new TH1D("IncEKin", "Kinematic Energy Inclusive;Energy (GeV)", 60, 0, 3);
  TH1* IncECal = new TH1D("IncECal", "Calorimetric Energy Inclusive;Energy (GeV);", 60, 0, 3);

  TH2* IncKinVsTrue = new TH2D("IncKinVsTrue", "Kinematic Reco E vs True E Inclusive;True Energy (GeV); Kin Reco (GeV)", 30, 0, 3, 30, 0, 3);
  TH2* IncCalVsTrue = new TH2D("IncCalVsTrue", "Calorimetric Reco E vs True E Inclusive;True Energy (GeV); Cal Reco (GeV)", 30, 0, 3, 30, 0, 3);

  TH1* IncKinRes = new TH1D("IncKinRes", "Kinematic resolution Inclusive;ETrue-EKin;", 40, -1.25, 2.25);
  TH1* IncCalRes = new TH1D("IncCalRes", "Calorimetric resolution Inclusive;ETrue-ECal;", 40, -1.25, 2.25);

  TH1* IncKinFracRes = new TH1D("IncKinFracRes", "Kinematic Fractional Resolution Inclusive;(ETrue-EKin)/ETrue;", 50, -1, 1);
  TH1* IncCalFracRes = new TH1D("IncCalFracRes", "Calorimetric Fractional Resolution Inclusive;(ETrue-ECal)/ETrue;", 50, -1, 1);

  TEfficiency* IncEff = new TEfficiency("IncEff", "Efficiency Inclusive;ETrue (GeV);Eff", 50, 0 ,5); //Efficiency of our cuts

    
  vector<double> KinResV;
  vector<double> CalResV;

  vector<double> KinFracResV;
  vector<double> CalFracResV;

  vector<double> KinResVNp;
  vector<double> CalResVNp;

  vector<double> KinFracResVNp;
  vector<double> CalFracResVNp;

  vector<double> KinResVInc;
  vector<double> CalResVInc;

  vector<double> KinFracResVInc;
  vector<double> CalFracResVInc;

  TLorentzVector vTrue; //particle 4-vector
  TLorentzVector vh; //Smeared proton 4-vector

  TRandom3 random;

  const double MN = 0.938919; //Average nucleon mass
  const double MProton = 0.938272; //proton mass
  const double MPi = 0.139570; //proton mass

  int NUnseen = 0; //# of NC events which have no protons, pions or gammas above threshold
  int NNotRealPh = 0;
  int NCalTooBig = 0; //# of events where EVisTot>ENu
  int N0Nu = 0;
  int NTooManyNu = 0;
  int NTooBigProtonEvents = 0;

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  cout<< "start loop\n";

  for (Long64_t jentry=0; jentry<nentries;jentry++) {  //Loop over events

    fChain->GetEntry(jentry);

    if(CCNC==0){ //Checks if NC
      continue;
    }

    int NProton = 0; //# of protons in final state topology
    int NNu = 0;
    int NPion = 0; //Whether there are pions above threshold
    int NGamma = 0; //Whether there are gammas above threshold
    double ENuPtcl = 0;
    double ENu = NuE;
    bool HasRealPh = true; // if smearing gives an energy less than the proton mass, total proton momentum won't be a real number
    double EVisTot = 0; //total visible energy of final state particles
    double PhTot = 0;
    double EhTot = 0;
    double PhxTot = 0;
    double PhyTot = 0;
    double PhzTot = 0;

    Double_t T = 0;

    double TrueETot = 0;
    double TruePTot = 0;
    double TruePzTot = 0;

    double BiggestP = 0;


    for(int i=0; i<Pdg->size(); i++){ //Loop over ptcls

      int pdgcode = Pdg->at(i);

      int statuscode = Status->at(i);
 
      if(statuscode!=1){
	continue;
      }

      double Energy = E->at(i); //true ptcl properties
      double PTrue = P->at(i); 
      double PxTrue = Px->at(i);
      double PyTrue = Py->at(i);
      double PzTrue = Pz->at(i);
//      if(i!=0 && T!=Time.at(i)){
//	cout<<"Non-initial ptcl. T_old="<<T<<", T_new="<<Time.at(i)<<"\n";
  //    }
//      double T = Time.at(i);
 
//      double PFromVec = sqrt(pow(PxTrue, 2)+pow(PyTrue, 2)+pow(PzTrue, 2));
//      double PFracDiff = (PTrue-PFromVec)/PTrue;
//      if(pow(PFracDiff, 2) > 0.00001){
//        std::cout<<"PTrue="<<PTrue<<" from vec="<<PFromVec<<" Frac Diff="<<PFracDiff<<" PDG:"<<pdgcode<<"\n";
//      }
      vTrue.SetE(Energy); //true 4-vector
      vTrue.SetPx(PxTrue);
      vTrue.SetPy(PyTrue);
      vTrue.SetPz(PzTrue);

      double KE = Energy - vTrue.M();

      double EVis = 0; //visible 
      double Eh = 0;  
      double Ph = 0;
      double Phx = 0;
      double Phy = 0;
      double Phz = 0;

      if(pdgcode==12 || pdgcode==14){
//	if(NNu>0){
//	  cout<<"Old ENu="<<ENu<<" Energy="<<Energy<<"\n";
//	}

//        if(Energy>ENuPtcl){ //makes sure that the nu energy selected is from the incoming nu, not the outgoing
//	  if(ENu>0){
//	    cout<<"Old ENu="<<ENu<<" Energy="<<Energy;
//	    ENu = Energy;
//	    cout<<"New ENu="<<ENu<<"\n";
//	  }else{
//	    ENuPtcl = Energy;
//	  }
//	}
        NNu += 1;
        
      }else if(pdgcode==2212 && KE>0.025){ //proton
 
	NProton+=1;

//	Eh = random.Gaus(Energy, 0.06);  //Proton E. has E resolution of 60 MeV
	Eh = Energy;

	if(pow(Eh, 2) - pow(MProton, 2) >= 0){
	  Ph = sqrt(pow(Eh, 2) - pow(MProton, 2)); 
	}else{
	  HasRealPh = false;
	  break;
	}

	double thetaTrue = acos(PzTrue/PTrue);
	double theta = random.Gaus(thetaTrue, 5 * 3.1415/180); //5 degree resolution, but theta is in radians
	double phiTrue = asin(PyTrue/(PTrue * sin(thetaTrue)));
	double phi = random.Gaus(phiTrue, 5 * 3.1415/180);

//	Phx = Ph * sin(theta) * cos(phi); //calculates 4-vec components after accounting for resolution
//	Phy = Ph * sin(theta) * sin(phi);
//	Phz = Ph * cos(theta); 

	Phx = PxTrue; //calculates 4-vec components after accounting for resolution
	Phy = PyTrue;
	Phz = PzTrue; 

      TrueETot+=Energy;
 //     TruePTot+=PTrue;
 //     TruePzTot+=PzTrue;

	EVis = Eh - vTrue.M(); //visible energy of this proton
      
        EhTot+= Eh;
        EVisTot += EVis;
        PhxTot += Phx;
        PhyTot += Phy;
        PhzTot += Phz;

	if(Energy>EVis){
//	  cout<<"Energy="<<Energy<<" BiggestP="<<BiggestP<<"\n";
	  BiggestP = EVis;
	} 

//	std:cout << "proton Eh=" << Eh <<" Proton EVis="<<EVis<<"\n";

      }else if(((pdgcode==211 || pdgcode==-211) && KE>0.010) || pdgcode==111){ // pi+,pi-,pi0

	NPion += 1;

//	Eh = random.Gaus(Energy, Energy*0.1);  //Pion E. has E resolution of 10%
	Eh = Energy;
	
	if(pow(Eh, 2) - pow(MPi, 2) >= 0){
	  Ph = sqrt(pow(Eh, 2) - pow(MPi, 2)); 
	}else{
	  HasRealPh = false;
	  break;
	}

	double thetaTrue = acos(PzTrue/PTrue);
	double theta = random.Gaus(thetaTrue, 2 * 3.1415/180); //2 degree resolution for pions, but theta is in radians
	double phiTrue = asin(PyTrue/(PTrue * sin(thetaTrue)));
	double phi = random.Gaus(phiTrue, 2 * 3.1415/180);

//	Phx = Ph * sin(theta) * cos(phi); //calculates 4-vec components after accounting for resolution
//	Phy = Ph * sin(theta) * sin(phi);
//	Phz = Ph * cos(theta); 

	Phx = PxTrue; //calculates 4-vec components after accounting for resolution
	Phy = PyTrue;
	Phz = PzTrue; 

      TrueETot+=Energy;
//      TruePTot+=PTrue;
//      TruePzTot+=PzTrue;

	EVis = Eh; //visible energy of this ptcl
      
        EhTot+= Eh;
        EVisTot += EVis;
        PhxTot += Phx;
        PhyTot += Phy;
        PhzTot += Phz;

      }else if(pdgcode==22 && Energy>0.030){ //photon

	NGamma += 1;

//	Eh = random.Gaus(Energy, Energy*0.1);  //Photon E. has E resolution of 10%
	Eh = Energy;

	Ph = sqrt(pow(Eh, 2)); 

	double thetaTrue = acos(PzTrue/PTrue);
	double theta = random.Gaus(thetaTrue, 5 * 3.1415/180); //5 degree resolution, but theta is in radians
	double phiTrue = asin(PyTrue/(PTrue * sin(thetaTrue)));
	double phi = random.Gaus(phiTrue, 5 * 3.1415/180);

//	Phx = Ph * sin(theta) * cos(phi); //calculates 4-vec components after accounting for resolution
//	Phy = Ph * sin(theta) * sin(phi);
//	Phz = Ph * cos(theta); 

	Phx = PxTrue; //calculates 4-vec components after accounting for resolution
	Phy = PyTrue;
	Phz = PzTrue; 

        TrueETot+=Energy;
        TruePTot+=PTrue;
        TruePzTot+=PzTrue;

	EVis = Eh; //visible energy of this ptcl      

        EhTot+= Eh;
        EVisTot += EVis;
        PhxTot += Phx;
        PhyTot += Phy;
        PhzTot += Phz;

      }
//	std::cout << "Photon Eh=" << Eh <<" Photon EVis="<<EVis<<"\n";

//      }else{
//	std::cout<<"Particle "<<pdgcode<<" present. E="<<Energy<<" P="<<PTrue<<" Pz="<<PzTrue<<"\n";
    }

    PhTot = sqrt(pow(PhxTot, 2) + pow(PhyTot, 2) + pow(PhzTot, 2));

//    if(NNu != 2){
//      cout << "Warning: event has " << NNu << " neutrinos. Should have 2: incoming and outgoing. Selction info: NProton=" << NProton << " NPion=" << NPion << " NGamma="<<NGamma<<"\n";
//    }
    
    Eff->Fill(NProton==1 && NPion==0 && NGamma==0, ENu);
    NpEff->Fill(NProton>0 && NPion==0 && NGamma==0, ENu);
    IncEff->Fill(NProton>0 || NPion>0 || NGamma>0, ENu);
    NCSpectrum->Fill(ENu);

//    if(ENu!=ENuPtcl){
//      cout<<"ENuPtcl="<<ENuPtcl<<" ENu="<<ENu<<"\n";
//    }

    double kinematic = 0.5*(pow(PhTot, 2)-pow(EhTot-MN, 2))/(MN+PhzTot-EhTot); //include mass in the energy we use here.

    if((NProton>0 || NPion>0 || NGamma>0) && HasRealPh == true){  

      if(NNu==0){
        N0Nu+=1;
	cout<<"0 nu event";
        continue;

      }else if(NNu>1){
        NTooManyNu+=1;
	cout<<NNu<<" nu event\n";
      }

      if(BiggestP>ENu){
        cout << "ENu = "<<ENu<<" BiggestP="<<BiggestP<<" PPiG:"<<NProton<<NPion<<NGamma<<"\n";
	NTooBigProtonEvents+=1;
      }

      if(EVisTot>ENu){
        NCalTooBig+=1;  
        std::cout<<"ENu="<<ENu<<" EVisTot="<<EhTot<<" TrueETot="<<TrueETot<<"PPiG:"<<NProton<<NPion<<NGamma<<"\n"; 
      }

      IncETrue->Fill(ENu);
      IncECal->Fill(EVisTot);
      IncCalVsTrue->Fill(ENu, EVisTot);
      IncCalRes->Fill(ENu-EVisTot);
      IncCalFracRes->Fill((ENu-EVisTot)/ENu);

      CalResVInc.push_back(ENu-EVisTot);
      CalFracResVInc.push_back((ENu-EVisTot)/ENu);

      if(kinematic>0){

        IncEKin->Fill(kinematic);
        IncKinVsTrue->Fill(ENu, kinematic);
        IncKinRes->Fill(ENu-kinematic);
        IncKinFracRes->Fill((ENu-kinematic)/ENu); 
        KinResVInc.push_back(ENu-kinematic);
        KinFracResVInc.push_back((ENu-kinematic)/ENu);

      } 
    }else if(HasRealPh==true){
      NUnseen+=1;
    }else if(HasRealPh==false){
      NNotRealPh+=1;
    }
   
    if(NProton==1 && NPion==0 && NGamma==0 && HasRealPh == true){  

      ETrue->Fill(ENu);
      ECal->Fill(EVisTot);
      CalVsTrue->Fill(ENu, EVisTot);
      CalRes->Fill(ENu-EVisTot);
      CalFracRes->Fill((ENu-EVisTot)/ENu);
      CalResV.push_back(ENu-EVisTot);

      KinFracResV.push_back((ENu-kinematic)/ENu);
      CalFracResV.push_back((ENu-EVisTot)/ENu);
      
      if(kinematic>0){

        EKin->Fill(kinematic);
        KinVsTrue->Fill(ENu, kinematic);
        KinRes->Fill(ENu-kinematic);
        KinFracRes->Fill((ENu-kinematic)/ENu); 
        KinResV.push_back(ENu-kinematic);
        KinFracResV.push_back((ENu-kinematic)/ENu);

      }
    }   

    if(NProton>0 && NPion==0 && NGamma==0 && HasRealPh == true){  

      NpETrue->Fill(ENu);
      NpECal->Fill(EVisTot);
      NpCalVsTrue->Fill(ENu, EVisTot);
      NpCalRes->Fill(ENu-EVisTot);
      NpCalFracRes->Fill((ENu-EVisTot)/ENu);

      CalResVNp.push_back(ENu-EVisTot);
      CalFracResVNp.push_back((ENu-EVisTot)/ENu);

      if(kinematic>0){

        NpEKin->Fill(kinematic);
        NpKinVsTrue->Fill(ENu, kinematic);
        NpKinRes->Fill(ENu-kinematic);
        NpKinFracRes->Fill((ENu-kinematic)/ENu); 
        KinResVNp.push_back(ENu-kinematic);
        KinFracResVNp.push_back((ENu-kinematic)/ENu);
      
      } 
    }
  }
  
  cout<< "end loop\n" << "NUnseen=" << NUnseen <<" NNotRealPh=" << NNotRealPh << " NCalTooBig="<<NCalTooBig<<" N0Nu="<<N0Nu<<" NTooManyNu="<<NTooManyNu<<" NTooBigProtontEvents="<<NTooBigProtonEvents<<"\n";

  makeplot(ETrue);
  makeplot(EKin);
  makeplot(ECal);
  makeplot2d(KinVsTrue);
  makeplot2d(CalVsTrue);
  makeplot(KinRes, KinResV);
  makeplot(CalRes, CalResV);
  makeplot(KinFracRes, KinFracResV);
  makeplot(CalFracRes, CalFracResV);

  makeplot(NpETrue);
  makeplot(NpEKin);
  makeplot(NpECal);
  makeplot2d(NpKinVsTrue);
  makeplot2d(NpCalVsTrue);
  makeplot(NpKinRes, KinResVNp);
  makeplot(NpCalRes, CalResVNp);
  makeplot(NpKinFracRes, KinFracResVNp);
  makeplot(NpCalFracRes, CalFracResVNp);

  makeplot(IncETrue);
  makeplot(IncEKin);
  makeplot(IncECal);
  makeplot2d(IncKinVsTrue);
  makeplot2d(IncCalVsTrue);
  makeplot(IncKinRes, KinResVInc);
  makeplot(IncCalRes, CalResVInc);
  makeplot(IncKinFracRes, KinFracResVInc);
  makeplot(IncCalFracRes, CalFracResVInc);


  TCanvas* c = new TCanvas();
  Eff->Draw();
  NCSpectrum->Scale(6.0/NCSpectrum->Integral());
  NCSpectrum->Draw("same");
  c->SaveAs("/sbnd/app/users/rmacfady/workdir/ncreco/Efficiency.png");

  TCanvas* d = new TCanvas();
  NpEff->Draw();
  NCSpectrum->Scale(6.0/NCSpectrum->Integral());
  NCSpectrum->Draw("same");
  d->SaveAs("/sbnd/app/users/rmacfady/workdir/ncreco/NpEfficiency.png");

  TCanvas* e = new TCanvas();
  IncEff->Draw();
  NCSpectrum->Scale(6.0/NCSpectrum->Integral());
  NCSpectrum->Draw("same");
  d->SaveAs("/sbnd/app/users/rmacfady/workdir/ncreco/IncEfficiency.png");


  Eff->Write();
  ETrue->Write();
  EKin->Write();
  ECal->Write();
  NCSpectrum->Write();
  KinVsTrue->Write();
  CalVsTrue->Write();
  KinRes->Write();
  CalRes->Write();
  KinFracRes->Write();
  CalFracRes->Write();
 
  IncEff->Write();
  IncETrue->Write();
  IncEKin->Write();
  IncECal->Write();
  IncKinVsTrue->Write();
  IncCalVsTrue->Write();
  IncKinRes->Write();
  IncCalRes->Write();
  IncKinFracRes->Write();
  IncCalFracRes->Write();
   
  f->Close();

}

double fwhm(TH1* h){
   int bin1 = h->FindFirstBinAbove(h->GetMaximum()/2);
   int bin2 = h->FindLastBinAbove(h->GetMaximum()/2);
   return h->GetBinCenter(bin2) - h->GetBinCenter(bin1);
}


double Median(vector<double>* v){
   std::sort(v->begin(),v->end());
    int index = v->size()/2;
    if(v->size()%2==0) index-=1;
    double med = v->at(index);
    if(v->size()%2==0){
        med+=v->at(index+1);
        med*=0.5;
    }
    return med;
}


void makeplot(TH1* h){
  TCanvas* c = new TCanvas();
  h->Draw();
  std::stringstream ss;
  ss << "/sbnd/app/users/rmacfady/workdir/ncreco/" << h->GetName() << ".png";
  c->SaveAs(ss.str().c_str());
}

void makeplot2d(TH2* h){
  TCanvas* c = new TCanvas();
  h->SetStats(0);
  h->Draw("colz");
  std::stringstream ss;
  ss << "/sbnd/app/users/rmacfady/workdir/ncreco/" << h->GetName() << ".png";
  c->SaveAs(ss.str().c_str());
}

void makeplot(TH1* h, vector<double> v){
  TCanvas* c = new TCanvas();

  TPaveText *pt = new TPaveText(.15,.75,.4,.9, "NDC");
  std::stringstream ss1;
  std::stringstream ss2;
  ss1 << "Median = " << Median(&v);
  pt->AddText(ss1.str().c_str());
  ss2 << "FWHM = "<<fwhm(h);
  pt->AddText(ss2.str().c_str());

  h->Draw();
  pt->Draw();

  std::stringstream ss3;
  ss3 << "/sbnd/app/users/rmacfady/workdir/ncreco/" << h->GetName() << ".png";
  c->SaveAs(ss3.str().c_str());
}



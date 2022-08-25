//-- Author :  Igor Kreslo, 2017

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// CRTEvent                                                            //
//                                                                      //
// Base event objects for uBooNE/SBND CRT data                          //
// Igor Kreslo, LHEP Uni-Bern (Igor.Kreslo@cern.ch)                     //
//////////////////////////////////////////////////////////////////////////
#define INCH 2.54


#include <iostream>
#include <fstream>
#include <stdio.h>

#include <TSystem.h>
#include <TObject.h>
#include <TClass.h>
#include <TRandom.h>
#include <TMath.h>
#include <TH3.h>
#include <TPolyLine3D.h>

#include "CRTEvent.h"

using namespace std;

typedef struct {
		UShort_t mac5; // ==0xFFFF
		UShort_t flags; // ==0xFFFF
		UShort_t lostcpu;
		UShort_t lostfpga;
		UInt_t ts0; // ==MAGICWORD32
		UInt_t ts1; // ==MAGICWORD32
                Int_t nevsinpoll; 
		UInt_t start_s;
		UInt_t d1;
		UShort_t start_ms;
		UShort_t dd2;
		UInt_t d2;
		UInt_t end_s;
		UInt_t d3;
		UShort_t end_ms;
} EOP_EVENT_t;  // end-of-poll special event


// ClassImp(CRTEvent)
// ClassImp(CRTRawhit)
// ClassImp(CRT2Dhit)
// ClassImp(CRTCalibs)

//______________________________________________________________________________
CRTEvent::CRTEvent()
{
  hits=new TClonesArray("CRTRawhit",100);
  mean_t0=0;
  mean_t1=0;
  Nhits=0;
 // hits->SetOwner(kTRUE);
}


//______________________________________________________________________________
void CRTEvent::Print(Int_t Verbosity) 
{
  for(int i=0; i<Nhits; i++)
  ((CRTRawhit*)(hits->At(i)))->Print(Verbosity);
  printf("%d hits in Event at %d s, <t0>=%f  <t1>=%f.\n",Nhits,s,mean_t0,mean_t1);
}


//______________________________________________________________________________
CRTEvent::~CRTEvent()
{
 // hits->Clear();
 // hits->~TClonesArray();
  //  Clear();
  if(hits) {hits->Delete(); delete hits, hits=0;}
//  SafeDelete(hits);
//if(CL) { CL->Delete(); delete CL; CL=0; }
}
//______________________________________________________________________________
void CRTEvent::Clear()
{
  hits->Clear();
  mean_t0=0;
  mean_t1=0;
  Nhits=0;
  s=0; 
}


//______________________________________________________________________________
void CRTEvent::AddRawhit(CRTRawhit * hit)
{
  CRTRawhit * ht;
if(Nhits>0)
  {
  mean_t0=mean_t0*Nhits;
  mean_t1=mean_t1*Nhits;
  } else 
  {
  mean_t0=0;
  mean_t1=0;
  }
 // ht=(CRTRawhit*)(hits->ConstructedAt(Nhits));
  ht=(CRTRawhit*)(new((*hits)[hits->GetLast()+1])  CRTRawhit());
  ht->Copy(hit);
  Nhits++;
// recalculate mean values
  mean_t0=mean_t0+hit->ts0;
  mean_t1=mean_t1+hit->ts1;
  mean_t0=mean_t0/Nhits;
  mean_t1=mean_t1/Nhits;
  s=hit->s;
}
/*
//______________________________________________________________________________
Int_t CRTEvent::CleanDuplicates()
{
  Int_t removed=0;
  Int_t ind=0;
  while(ind<hits->GetEntriesFast())
  {
   ind++;
  }
  return removed;
}
*/

//______________________________________________________________________________
CRTRawhit::CRTRawhit()
{
  cable=0; 
}

//______________________________________________________________________________
CRTRawhit::~CRTRawhit()
{
}

//______________________________________________________________________________
void CRTRawhit::ReadFromBuffer(void *buf)
{
  memcpy(&mac5, buf, EVLEN);
  cable=0;
}

//______________________________________________________________________________
void CRTRawhit::Print(Int_t Verbosity)
{
int msec;
int sec;
int emsec;
int esec;
int nevsinpoll;

  if(mac5<=0x00FF) //Normal event
  {
  if(Verbosity==0) return;

  if(Verbosity>2) 
   {  
   printf("mac5=%d flags=0x%04x t0=%d t1=%d s=%d ms=%d, cable delay %lf ns; ",mac5,flags,ts0,ts1,s,ms,cable);
   if ((flags & kT0RefOK)==0) printf("T0Ref-NOTOK "); 
   if ((flags & kT1RefOK)==0) printf("T1Ref-NOTOK "); 
   if ((flags & kT0RefEventFlag)>0) printf("T0RefEvent "); 
   if ((flags & kT1RefEventFlag)>0) printf("T1RefEvent "); 
   if ((flags & kBeamEventFlag)>0) printf("BeamEvent ");
   printf("\n"); 
   }
  if(Verbosity>3) 
   {
     printf("ADCs: ");
     for(int i=0;i<32;i++) printf("%d ",adc[i]);
     printf("\n"); 
   }
  if(Verbosity>1) { if((flags & kT0RefEventFlag)>0 ) printf("mac5=%d T0RefEvent t0=%d\n", mac5,ts0);  }
  if(Verbosity>1) { if((flags & kT1RefEventFlag)>0 ) printf("mac5=%d T1RefEvent t1=%d\n", mac5,ts1);  }

  }
  else //EOP (special) event
  {
   	sec=((EOP_EVENT_t*)(&mac5))->start_s;
	msec=((EOP_EVENT_t*)(&mac5))->start_ms;
	esec=((EOP_EVENT_t*)(&mac5))->end_s;
	emsec=((EOP_EVENT_t*)(&mac5))->end_ms;
        nevsinpoll=((EOP_EVENT_t*)(&mac5))->nevsinpoll;
  printf("EOP: start %d s %d ms; end %d s %d ms; %d events in poll.\n",sec,msec,esec,emsec,nevsinpoll);
  }
 
}

//_______________________________________________________________________
  Int_t CRTRawhit::GetSec()
{
  if(mac5<=0x00FF) //Normal event
       return s;
  else  //EOP event
       return ((EOP_EVENT_t*)(&mac5))->end_s;
}

//_______________________________________________________________________
  Int_t CRTRawhit::GetMaxADCIndex()
{
  Int_t max_adc=0, max_i=0;
       for(int i=0;i<32;i++) if(adc[i]>max_adc) { max_adc=adc[i]; max_i=i;}
  return max_i;
}

//_______________________________________________________________________
  Int_t CRTRawhit::GetMaxStripIndex()
{
   Int_t max_ach=0, max_i=0, ch;
     for(int i=0;i<16;i++)
     { 
        ch=adc[i*2]+adc[i*2+1];
        if(ch>max_ach) {max_ach=ch; max_i=i;}
     }
  return max_i;
}

//_______________________________________________________________________
  Int_t CRTRawhit::GetMaxADCValue()
{
  Int_t max_adc=0;
       for(int i=0;i<32;i++) if(adc[i]>max_adc) { max_adc=adc[i];}
  return max_adc;
}

//_______________________________________________________________________
  Int_t CRTRawhit::GetTriggeredADCValue()
{
  Int_t max_adc=0;
  Int_t trig_adc=0;
  Int_t max_i=0;
       for(int i=0;i<32;i++) if(adc[i]>max_adc) { max_adc=adc[i]; max_i=i;}
       if(max_i%2==0) trig_adc=adc[max_i+1];
       else trig_adc=adc[max_i-1]; 
  return trig_adc;
}

//_______________________________________________________________________
  Int_t CRTRawhit::GetTWCorrection(Double_t slope, Double_t offset)
{
  return offset+slope*log(GetTriggeredADCValue()/4096.);
}


//_______________________________________________________________________
  Int_t CRTRawhit::GetMaxStripValue()
{
   Int_t max_ach=0, ch;
     for(int i=0;i<16;i++)
     { 
        ch=adc[i*2]+adc[i*2+1];
        if(ch>max_ach) {max_ach=ch; }
     }
  return max_ach;
}

//______________________________________________________________________________
CRT2Dhit::CRT2Dhit()
{
  calibrated=0;
  Edep=0;
}

//______________________________________________________________________________
CRT2Dhit::CRT2Dhit(CRTRawhit *h1, CRTRawhit *h2, CRTCalibs *cal)
{
 /*   //Reconstructed digits
   Double_t t0;
   Double_t t1;
   Double_t dt0;
   Double_t dt1;
   Int_t s;
   Double_t x,y,z;
   Int_t nhits;
   Int_t plane;
   CRTRawhit * rawhit[2];
*/

   int max1_adc=0;
   int max2_adc=0;
   int max1_ach=0;
   int max2_ach=0;
   int max1_nch=0;
   int max2_nch=0;
   int ch[2][16];
   Double_t calts1[2], calts0[2]; 
   nhits1=0;
   nhits2=0;

   rhit[0].Copy(h1); 
   rhit[1].Copy(h2);

     for(int i=0;i<32;i++)
     { 
        if(rhit[0].adc[i]>max1_adc) max1_adc=rhit[0].adc[i];
        if(rhit[1].adc[i]>max2_adc) max2_adc=rhit[1].adc[i];
     }
     for(int i=0;i<16;i++)
     { 
        ch[0][i]=rhit[0].adc[i*2]+rhit[0].adc[i*2+1];
        if(ch[0][i]>max1_ach) {max1_ach=ch[0][i]; max1_nch=i;}
        ch[1][i]=rhit[1].adc[i*2]+rhit[1].adc[i*2+1];
        if(ch[1][i]>max2_ach) {max2_ach=ch[1][i]; max2_nch=i;}
     }
     for(int i=0;i<16;i++)
     {
      if(ch[0][i]> max1_ach/2) nhits1++;
      if(ch[1][i]> max2_ach/2) nhits2++;
     } 
      
     if(cal==0)
     {
       calts0[0]=rhit[0].ts0;
       calts1[0]=rhit[0].ts1; 
       calts0[1]=rhit[1].ts0;
       calts1[1]=rhit[1].ts1;
    x=0;y=0;z=0; plane1=-1; plane2=-1;
     } 
     else 
     {
      calts0[0]=cal->getHitT1(rhit[0].mac5, max1_nch, rhit[0].ts0, rhit[1].mac5, max2_nch, rhit[1].ts0); 
      calts0[1]=cal->getHitT2(rhit[0].mac5, max1_nch, rhit[0].ts0, rhit[1].mac5, max2_nch, rhit[1].ts0); 
      calts1[0]=cal->getHitT1(rhit[0].mac5, max1_nch, rhit[0].ts1, rhit[1].mac5, max2_nch, rhit[1].ts1); 
      calts1[1]=cal->getHitT2(rhit[0].mac5, max1_nch, rhit[0].ts1, rhit[1].mac5, max2_nch, rhit[1].ts1);
      rhit[0].cable=calts0[0]-rhit[0].ts0;    
      rhit[1].cable=calts0[1]-rhit[1].ts0;    
   x=cal->getHitX(rhit[0].mac5, max1_nch,rhit[1].mac5, max2_nch);
   y=cal->getHitY(rhit[0].mac5, max1_nch,rhit[1].mac5, max2_nch);
   z=cal->getHitZ(rhit[0].mac5, max1_nch,rhit[1].mac5, max2_nch);
   // std::cout << "CRT2Dhit, from " << rhit[0].mac5 << " " << max1_nch << " " << rhit[1].mac5 << " " << max2_nch << std::endl;
   // std::cout << "CRT2Dhit, x " << x << "  y " << y << "  z " << z << std::endl;
   plane1=cal->getHitPlane(rhit[0].mac5);
   plane2=cal->getHitPlane(rhit[1].mac5);
     
     }
    dt0=-calts0[0]+calts0[1];
    dt1=-calts1[0]+calts1[1];
    t0=(calts0[1] + calts0[0])/2.;
    t1=(calts1[1] + calts1[0])/2.;
    s=rhit[0].s;
    Edep=2.1*(rhit[0].GetMaxStripValue() + rhit[1].GetMaxStripValue())/2./1000.;


/*     
//   if(cal==0)
//   {
    calibrated=0;
    t0=(rhit[0].ts0+rhit[1].ts0)/2.;
    t1=(rhit[0].ts1+rhit[1].ts1)/2.;
    dt0=-rhit[0].ts0+rhit[1].ts0;
    dt1=-rhit[0].ts1+rhit[1].ts1;
//    s=(rhit[0].s+rhit[1].s)/2.;
    s=rhit[0].s;
    x=0;y=0;z=0; plane1=-1; plane2=-1;
    Edep=2.1*(rhit[0].GetMaxStripValue() + rhit[1].GetMaxStripValue())/8192;
   }  
   else
   {
   rhit[0].ts0=cal->getHitT1(rhit[0].mac5, max1_nch, rhit[0].ts0, rhit[1].mac5, max2_nch, rhit[1].ts0); 
   rhit[1].ts0=cal->getHitT2(rhit[0].mac5, max1_nch, rhit[0].ts0, rhit[1].mac5, max2_nch, rhit[1].ts0); 
   rhit[0].ts1=cal->getHitT1(rhit[0].mac5, max1_nch, rhit[0].ts1, rhit[1].mac5, max2_nch, rhit[1].ts1); 
   rhit[1].ts1=cal->getHitT2(rhit[0].mac5, max1_nch, rhit[0].ts1, rhit[1].mac5, max2_nch, rhit[1].ts1); 
   dt0=rhit[1].ts0 - rhit[0].ts0;
   dt1=rhit[1].ts1 - rhit[0].ts1;
   t0=(rhit[1].ts0 + rhit[0].ts0)/2.;
   t1=(rhit[1].ts1 + rhit[0].ts1)/2.;
 //  s=(rhit[0].s+rhit[1].s)/2.;
   s=rhit[0].s;
   x=cal->getHitX(rhit[0].mac5, max1_nch,rhit[1].mac5, max2_nch);
   y=cal->getHitY(rhit[0].mac5, max1_nch,rhit[1].mac5, max2_nch);
   z=cal->getHitZ(rhit[0].mac5, max1_nch,rhit[1].mac5, max2_nch);
   plane1=cal->getHitPlane(rhit[0].mac5);
   plane2=cal->getHitPlane(rhit[1].mac5);

  //   pulse=(hit->adc[i] - cal->ADCPedestal[hit->mac5][i])/(cal->ADCGain[hit->mac5][i])*4096.;

   Edep=0;
  for(Int_t i=0;i<32;i++)
  {
    Edep=Edep+2.1*(h1->adc[i] - cal->ADCPedestal[h1->mac5][i])/(1.*cal->ADCGain[h1->mac5][i]);
    Edep=Edep+2.1*(h2->adc[i] - cal->ADCPedestal[h2->mac5][i])/(1.*cal->ADCGain[h2->mac5][i]);
  }

    calibrated=1;
   }
  //   Edep=rhit[0].GetMaxStripValue() + rhit[1].GetMaxStripValue();
  */  
}

//______________________________________________________________________________
CRT2Dhit::~CRT2Dhit()
{
  //  Clear();
}

//______________________________________________________________________________
void CRT2Dhit::Print(Int_t Verbosity)
{
  printf("2D xit:\n");
  rhit[0].Print(Verbosity);
  rhit[1].Print(Verbosity);
  printf("Edep=%lf MeV\n",Edep); 
}


//______________________________________________________________________________
CRTCalibs::CRTCalibs()
{
 FiberDelay=6.2;
}

//______________________________________________________________________________
CRTCalibs::~CRTCalibs()
{
  //  Clear();
}

//______________________________________________________________________________
void CRTCalibs::Print(Bool_t PrPos, Bool_t PrTime, Bool_t PrADC)
{
 int iPos=0;
 int iTime=0;
 int iADC=0;
 int iMod=0;
 if(PrPos) {
   printf("Position calibration data:\n");
   for(int m=0;m<NFEBS;m++) for(int c=0;c<32;c++) 
  {
    if(PositionExists[m][c]) { printf("%d %d (%f %f %f): plane %d\n",m,c,Xs[m][c],Ys[m][c],Zs[m][c], Plane[m][c]); iPos++;}
  }
  printf("%d Entries are printed.\n", iPos);
 } 
 if(PrTime) {
   printf("Cable delay calibration data:\n");
   for(int m=0;m<NFEBS;m++) 
  {
    if(CableDelayExists[m]) { printf("%d %f\n",m,Dts[m]); iTime++;}
  }
  printf("%d Entries are printed.\n", iTime);
 }
 if(PrADC) {
   printf("ADC calibration data:\n");
   for(int m=0;m<NFEBS;m++) for(int c=0;c<32;c++) 
  {
    if(ModuleExists[m]) { printf("%d %d %d %d \n",m,c,ADCPedestal[m][c],ADCGain[m][c]); iADC++;}
  }
  printf("%d Entries are printed.\n", iADC);
 }
   for(int m=0;m<NFEBS;m++) 
  {
    if(ModuleExists[m])  iMod++;
  }
  printf("Found data for %d Modules.\n", iMod);

}

//______________________________________________________________________________
CRTCalibs::CRTCalibs(const char *f_cable_delays, const char *f_positions, const char *f_adccalibs, Double_t StripWidth):TObject()
{
 StripW=StripWidth;
 FiberDelay=6.2;
 Double_t x, y, z, dt, ped, gain;
 Int_t id;
 Int_t pl,ch,p,b,c;
  for(int m=0;m<NFEBS;m++) 
  {
        for(int c=0;c<32;c++) 
        {
             Xs[m][c]=0;
             Ys[m][c]=0;
             Zs[m][c]=0;
             PositionExists[m][c]=0;
             ADCCalibsExists[m][c]=0;
             ADCPedestal[m][c]=0;   
             ADCGain[m][c]=1;   
        }
        Dts[m]=0;
        ModuleExists[m]=0;
        CableDelayExists[m]=0;      
  }

  if(strlen(f_positions)>0)
  { 
  std::ifstream in;
  in.open(f_positions);
  while (!in.eof()) {
    // std::cout << "" << std::endl;
    in>>id>>x>>y>>z>>p>>b>>c;
    pl=id/100; ch=(id-100*int(pl));
    // std::cout << "Emplacing " << x << " "<< y << " "<< z << " for pl " << pl << " and ch " << ch << ", given id " << id << std::endl;
    Xs[pl][ch]=x;
    Ys[pl][ch]=y;
    Zs[pl][ch]=z;
    Plane[pl][ch]=p;
    PositionExists[pl][ch]=1;
    ModuleExists[pl]=1;
    // std::cout << "done id " << id << std::endl;

  }      
  in.close();
  }

  if(strlen(f_cable_delays)>0)
  { 
  std::ifstream in1;
  in1.open(f_cable_delays);
  while (!in1.eof()) {
    in1>>pl>>dt;
    Dts[pl]=dt;
    ModuleExists[pl]=1;
    CableDelayExists[pl]=1;
  }      
  in1.close();
  }

  if(strlen(f_adccalibs)>0)
  { 
  std::ifstream in2;
  in2.open(f_adccalibs);
  while (!in2.eof()) {
    in2>>pl>>ch>>ped>>gain;
    ADCCalibsExists[pl][ch]=1;
    ADCPedestal[pl][ch]=ped;   
    ADCGain[pl][ch]=gain;   
  }      
  in2.close();
  }
}

//_______________________________________________________________________
Double_t CRTCalibs::getDistanceToSIPM1(int mac1, int strip1, int mac2, int strip2)
{
  Double_t L=0;
  if(Xs[mac2][0]!=Xs[mac2][2]) L=Xs[mac2][strip2*2]-Xs[mac1][0]; //coord along the strip minus SiPM position
  if(Ys[mac2][0]!=Ys[mac2][2]) L=Ys[mac2][strip2*2]-Ys[mac1][0]; //coord along the strip minus SiPM position
  if(Zs[mac2][0]!=Zs[mac2][2]) L=Zs[mac2][strip2*2]-Zs[mac1][0]; //coord along the strip minus SiPM position
  if(L<0) L=-L;
  return L;
}

//_______________________________________________________________________
Double_t CRTCalibs::getDistanceToSIPM2(int mac1, int strip1, int mac2, int strip2)
{
  Double_t L=0;
  if(Xs[mac1][0]!=Xs[mac1][2]) L=Xs[mac1][strip1*2]-Xs[mac2][0]; //coord along the strip minus SiPM position
  if(Ys[mac1][0]!=Ys[mac1][2]) L=Ys[mac1][strip1*2]-Ys[mac2][0]; //coord along the strip minus SiPM position
  if(Zs[mac1][0]!=Zs[mac1][2]) L=Zs[mac1][strip1*2]-Zs[mac2][0]; //coord along the strip minus SiPM position
  if(L<0) L=-L;
  return L;
}


//_______________________________________________________________________
Double_t CRTCalibs::getHitT1(int mac1, int strip1,  int t1, int mac2, int strip2, int t2)
{ 
  Double_t ret=0;
  Double_t L=0;
  // std::cout << " getHitT1 mac2 " << mac2 << std::endl;
  if(Xs[mac2][0]!=Xs[mac2][2]) L=Xs[mac2][strip2*2]-Xs[mac1][0]; //coord along the strip minus SiPM position
  // std::cout << " done! " << std::endl;
  if(Ys[mac2][0]!=Ys[mac2][2]) L=Ys[mac2][strip2*2]-Ys[mac1][0]; //coord along the strip minus SiPM position
  if(Zs[mac2][0]!=Zs[mac2][2]) L=Zs[mac2][strip2*2]-Zs[mac1][0]; //coord along the strip minus SiPM position
  if(L<0) L=-L;
  ret=t1-L*FiberDelay/100.; //propagation along the fiber correction, 6.2 ns/m
  ret=ret+Dts[mac1]; //ref PPS cable 
  return ret;
}

//_______________________________________________________________________
Double_t CRTCalibs::getHitT2(int mac1, int strip1,  int t1, int mac2, int strip2, int t2)
{ 
  Double_t ret=0; 
  Double_t L=0;
  if(Xs[mac1][0]!=Xs[mac1][2]) L=Xs[mac1][strip1*2]-Xs[mac2][0]; //coord along the strip minus SiPM position
  if(Ys[mac1][0]!=Ys[mac1][2]) L=Ys[mac1][strip1*2]-Ys[mac2][0]; //coord along the strip minus SiPM position
  if(Zs[mac1][0]!=Zs[mac1][2]) L=Zs[mac1][strip1*2]-Zs[mac2][0]; //coord along the strip minus SiPM position
  if(L<0) L=-L;
  ret=t2-L*FiberDelay/100.; //propagation along the fiber correction,6.2 ns/m
  ret=ret+Dts[mac2];//ref PPS cable  
  return ret;
}

//_______________________________________________________________________
Double_t CRTCalibs::getHitDT(int mac1, int strip1,  int t1, int mac2, int strip2,int t2)
{ 
 return getHitT2(mac1,strip1,t1,mac2,strip2,t2)-getHitT1(mac1,strip1,t1,mac2,strip2,t2);
}

//_______________________________________________________________________
Double_t CRTCalibs::getHitT(int mac1, int strip1,  int t1, int mac2, int strip2,  int t2)
{ 
 return (getHitT2(mac1,strip1,t1,mac2,strip2,t2)+getHitT1(mac1,strip1,t1,mac2,strip2,t2))/2.;
}


//_______________________________________________________________________
Double_t CRTCalibs::getHitX(int mac1, int strip1, int mac2, int strip2)
{ 
  // std::cout << "getHitX, mac1 " << mac1 << ", strip1 " << strip1 << ", mac2 " << mac2 << ", strip2 " << strip2 << std::endl;
  Double_t ret=0; 
  if(Xs[mac1][0]!=Xs[mac1][2]) ret=(Xs[mac1][strip1*2]+Xs[mac1][strip1*2+1])/2.;
  else if(Xs[mac2][0]!=Xs[mac2][2]) ret=(Xs[mac2][strip2*2]+Xs[mac2][strip2*2+1])/2.;
  else ret=(Xs[mac2][strip2*2]+Xs[mac1][strip1*2])/2;
  // std::cout << "ret: " << ret << std::endl;
  return ret;
}

//_______________________________________________________________________
Double_t CRTCalibs::getHitY(int mac1, int strip1, int mac2, int strip2)
{ 
  Double_t ret=0; 
  if(Ys[mac1][0]!=Ys[mac1][2]) ret=(Ys[mac1][strip1*2]+Ys[mac1][strip1*2+1])/2.;
  else if(Ys[mac2][0]!=Ys[mac2][2]) ret=(Ys[mac2][strip2*2]+Ys[mac2][strip2*2+1])/2.;
  else ret=(Ys[mac2][strip2*2]+Ys[mac1][strip1*2])/2;
  return ret;
}

//_______________________________________________________________________
Double_t CRTCalibs::getHitZ(int mac1, int strip1, int mac2, int strip2)
{ 
  Double_t ret=0; 
  if(Zs[mac1][0]!=Zs[mac1][2]) ret=(Zs[mac1][strip1*2]+Zs[mac1][strip1*2+1])/2.;
  else if(Zs[mac2][0]!=Zs[mac2][2]) ret=(Zs[mac2][strip2*2]+Zs[mac2][strip2*2+1])/2.;
  else ret=(Zs[mac2][strip2*2]+Zs[mac1][strip1*2])/2;
  return ret;
}

//_______________________________________________________________________
Int_t CRTCalibs::getHitPlane(int mac1)
{ 
  return Plane[mac1][0];
}


//______________________________________________________________________________
CRTTrack::CRTTrack()
{
  tof=0;
  L=0;
}

//______________________________________________________________________________
CRTTrack::CRTTrack(CRT2Dhit *h1, CRT2Dhit *h2)
{
   hit2d[0].Copy(h1); 
   hit2d[1].Copy(h2);
   tof=h2->t0-h1->t0; 
    Double_t dx,dy,dz;
   dx=h2->x - h1->x;   
   dy=h2->y - h1->y;   
   dz=h2->z - h1->z;   
   L=sqrt(dx*dx+dy*dy+dz*dz)/100.;   //calculate track length in meters 
   
}

//______________________________________________________________________________
void CRTTrack::SetHits(CRT2Dhit *h1, CRT2Dhit *h2)
{
   hit2d[0].Copy(h1); 
   hit2d[1].Copy(h2);
   tof=h2->t0-h1->t0; 
   Double_t dx,dy,dz;
   dx=h2->x - h1->x;   
   dy=h2->y - h1->y;   
   dz=h2->z - h1->z;   
   L=sqrt(dx*dx+dy*dy+dz*dz)/100.;   //calculate track length in meters 
}

//______________________________________________________________________________
CRTTrack::~CRTTrack()
{
  //  Clear();
}

//______________________________________________________________________________
void CRTTrack::Print(Int_t Verbosity)
{
  printf("Passing through track:\n");
  hit2d[0].Print(Verbosity);
  hit2d[1].Print(Verbosity);
}







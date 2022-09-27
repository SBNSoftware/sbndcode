#ifndef ROOT_CRTEvent
#define ROOT_CRTEvent

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// CRTEvent                                                            //
//                                                                      //
// Base event objects for uBooNE/SBND CRT data                          //
// Igor Kreslo, LHEP Uni-Bern (Igor.Kreslo@cern.ch)                     //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TList.h"
#include "TPolyMarker3D.h"
#include "TArrayF.h"
#include "TArrayS.h"

#define MAGICWORD8 0xa5 //marker for the buffer start in the file 
#define MAGICWORD16 0xaa55 //marker for the buffer start in the file 
#define MAGICWORD32 0xaa55aa55 //marker for the buffer start in the file 

#define NFEBS 200
#define EVLEN 80

class CRTEvent;
class CRTRawhit;
class CRT2Dhit;
class CRTTrack;
class CRTCalibs;


//______________________________________________________________________________
class CRTEvent : public TObject {

private:
public:
TClonesArray * hits; // Array of CRTRawhits
Double_t mean_t0; // Averaged t0 for all hits
Double_t mean_t1; // Averaged t1 for all hits
Int_t Nhits;
Int_t s; //Timestamp: Second (linux seconds)
  CRTEvent();
  virtual ~CRTEvent();
  void     Copy(CRTEvent *h) { memcpy( (void*)this, (void*)h, sizeof(CRTEvent) ); }
  void Print(Int_t Verbosity=4);
  void AddRawhit(CRTRawhit * hit);
//  Int_t CleanDuplicates(); //loops over hits and remove duplicates
  void Clear();
  void Dump(){Print(4);}
  // ClassDef(CRTEvent,1)  // CRT combined event: hits within some time window from all detector
};


//______________________________________________________________________________
class CRTRawhit : public TObject {

private:

public:
//Raw FEB digits
   UShort_t mac5; // FEB address (MAC5)
   UShort_t flags; // Event flags
   UShort_t lostcpu; // Lost events in CPU for current poll so far
   UShort_t lostfpga; // Lost events in FPGA for current poll so far
   Int_t ts0; // Time stamp from T0 counter, ns
   Int_t ts1; // Time stamp from T1 counter, ns
   UShort_t adc[32]; // Pointed to SiPM pulse height ADC array
   Int_t s; // Event Time stamp - Seconds (not assigned in raw FEB data)
   Int_t ms; // Event Time stamp - Milliseconds (not assigned in raw FEB data)
 //  CRT2Dhit * Hit;
   Double_t cable; // Time correction (Cable delay - light propagation delay) to be added to ts0 or ts1
  CRTRawhit();
  void ReadFromBuffer(void *buf);
  Bool_t IsT0OK(){ return (flags & kT0RefOK);}
  Bool_t IsT1OK(){ return (flags & kT1RefOK);}
  Bool_t IsT0RefEvent(){ return (flags & kT0RefEventFlag);}
  Bool_t IsT1RefEvent(){ return (flags & kT1RefEventFlag);}
  Bool_t IsBeamEvent(){ return (flags & kBeamEventFlag);}
  Bool_t IsEOPEvent(){ return (mac5==0xFFFF);}
  Int_t GetEventLenght() {return EVLEN;}
  Int_t GetMaxADCIndex();
  Int_t GetMaxStripIndex();
  Int_t GetMaxADCValue(); //maximum ADC pulseheight
  Int_t GetTriggeredADCValue(); //Pulse height in the ADC channel, paired to Maximum. Used for TW correction
  Int_t GetTWCorrection(Double_t slope, Double_t offset); //time walk correction
  Int_t GetMaxStripValue();
  Int_t GetSec(); //Get seconds from the hit If EOP event - returns seconds for end of poll. 
//  Int_t GetDisplayable(Double_t &x, Double_t &y, Double_t &z, Double_t &sx, Double_t &sy, Double_t &sz, Double_t &pulse, CRTCalibs *cal=0);
  virtual ~CRTRawhit();
  void     Copy(CRTRawhit *h) { memcpy( (void*)this, (void*)h, sizeof(CRTRawhit) ); }
  void Print(Int_t Verbosity=4);
  enum {
         kT0RefOK = BIT(0), //T0 Reference pulse is present since less than a second 
         kT1RefOK = BIT(1), //T1 Reference pulse is present since less than a second
         kT0RefEventFlag = BIT(2), //T0 Reference event flag
         kT1RefEventFlag = BIT(3), //T1 Reference event flag
         kBeamEventFlag = BIT(4)   // Indicates, that the event is within selected time window around the Beam trigger (T1). Typically 50us. 
       };  
 
   

  // ClassDef(CRTRawhit,1)  // CRT event from FEB
};
//______________________________________________________________________________

class CRT2Dhit : public TObject {

private:

public:
  //Reconstructed digits
   Double_t t0; //Average t0 time stamp, ns
   Double_t t1; //Average t1 time stamp, ns 
   Double_t dt0; //Difference between raw hit's calibrated ts0
   Double_t dt1; //Difference between raw hit's calibrated ts1
   Int_t s; //Timestamp seconds
   Double_t x,y,z; //calibrated 3D coordinates of the hit
   Int_t nhits1; //Number of ADCs above 50% of max for raw hit#1
   Int_t nhits2; //Number of ADCs above 50% of max for raw hit#2
   Int_t plane1; //Plane to which raw hit#1 belongs (0-Bot, 1-? 2-? 3-Top)
   Int_t plane2; //Plane to which raw hit#2 belongs (0-Bot, 1-? 2-? 3-Top)
   CRTRawhit rhit[2]; //Raw hits
   Bool_t calibrated; //Set to 1 if Calibration object was present during creation.
   Double_t Edep; //Approximate deposited energy, MeV. Assumes calibration is made for ADC=1000 for normal MIP.
    
   
  CRT2Dhit();
  CRT2Dhit(CRTRawhit *h1, CRTRawhit *h2, CRTCalibs *cal=0);

  virtual ~CRT2Dhit();
  void     Copy(CRT2Dhit *h) { memcpy( (void*)this, (void*)h, sizeof(CRT2Dhit) ); }
  void Print(Int_t Verbosity=4);
  
  // ClassDef(CRT2Dhit,1)  // CRT reconstructed 2D hit
};

//______________________________________________________________________________
class CRTCalibs : public TObject {

private:
Double_t StripW;
Int_t Version;

public:
Double_t FiberDelay; //Light in fiber delay, ns/m
  //Calibration data
Int_t ModuleExists[NFEBS];
//Coordinate (position)
Bool_t PositionExists[NFEBS][32];
Double_t Xs[NFEBS][32];
Double_t Ys[NFEBS][32];
Double_t Zs[NFEBS][32];
Int_t Plane[NFEBS][32];
//Cable delays
Bool_t CableDelayExists[NFEBS];
Double_t Dts[NFEBS];
//ADC calibration data
Bool_t ADCCalibsExists[NFEBS][32];
Int_t ADCPedestal[NFEBS][32];   
Int_t ADCGain[NFEBS][32];   
  CRTCalibs();
  CRTCalibs(const char *f_cable_delays, const char *f_positions, const char *f_adccalibs, Double_t StripWidth=10.8);
  virtual ~CRTCalibs();
  void Print(Bool_t PrPos=1, Bool_t PrTime=1, Bool_t PrADC=1);
  void SetVersion(Int_t v){Version=v;}
  Int_t GetVersion() {return Version;}
  void SetStripWidth(Double_t StripWidth) {StripW=StripWidth;}

Double_t getDistanceToSIPM1(int mac1, int strip1, int mac2, int strip2);
Double_t getDistanceToSIPM2(int mac1, int strip1, int mac2, int strip2);
Double_t getHitT1(int mac1, int strip1,  int t1, int mac2, int strip2, int t2);
Double_t getHitT2(int mac1, int strip1,  int t1, int mac2, int strip2, int t2);
Double_t getHitDT(int mac1, int strip1,  int t1, int mac2, int strip2,int t2);
Double_t getHitT(int mac1, int strip1,  int t1, int mac2, int strip2,  int t2);
Double_t getHitX(int mac1, int strip1, int mac2, int strip2);
Double_t getHitY(int mac1, int strip1, int mac2, int strip2);
Double_t getHitZ(int mac1, int strip1, int mac2, int strip2);
Int_t getHitPlane(int mac1);
  
  // ClassDef(CRTCalibs,1)  // CRT calibration data
};

//______________________________________________________________________________
class CRTTrack : public TObject {

private:
public:
  CRT2Dhit  hit2d[2];
  Double_t tof; //time of flight in ns
  Double_t L; //track length in m

  CRTTrack();
  CRTTrack(CRT2Dhit *h1, CRT2Dhit *h2);
  virtual ~CRTTrack();
  void     Copy(CRTTrack *h) { memcpy( (void*)this, (void*)h, sizeof(CRTTrack) ); }
  void Print(Int_t Verbosity=4);
 // void AddRawhit(CRTRawhit * hit);
  void Clear();
  void SetHits(CRT2Dhit *h1, CRT2Dhit *h2);
  void Dump(){Print(4);}

  // ClassDef(CRTTrack,1)  // Track composed of a pair of 2D hits
};


#endif /* ROOT_CRTEvent */




#pragma once
#include <iostream>

struct MuTrack {
  int ID;
  //int TruePDG;
  //int TrueID;
  double Length;
  double StartX, StartY, StartZ;
  double EndX, EndY, EndZ;
};

struct MuHit {
  int PlaneID;
  int TpcID;
  int ChannelID;
  double PeakTime; // Already rescaled by 1/4
  double Integral;
  double RMS;
  int TruePdg;
  int PfpID;
  int TrueID;
};

struct MuCalo {
  int PlaneID;
  int TpcID;
  int ChannelID;
  double PeakTime; // Already rescaled by 1/4
  int TrueID;
  int TruePdg;
  //int PfpID;
  double dEdx;
  double ResRange;
  double X, Y, Z;
};

struct MuBlip {
  // Reco blips variables
  int ID;
  int IsValid;
  int Cryostat;
  int TpcID;
  int NPlanes;
  double MaxWireSpan;
  double TimeTick; // Already rescaled by 1/4
  double Charge;
  double Energy;
  double EnergyESTAR;
  double EnergyPSTAR;
  double ProxTrkDist;
  double ProxTrkID;
  double InCylinder;
  double X, Y, Z;
  double SigmaYZ;
  double dX;
  double dYZ;
};


struct Bragg {
  bool Found = false;
  int Idx;
  int ChannelID;
  double PeakTime;
  double ResRange;
  double TpcID;
  double dEdx;
  double X;
  double Y;
  double Z;
};

struct PCAResult // In coordinates channelID - timeTick
{
  bool Valid = false;
  double dirX;
  double dirY;
};

struct RadialParams {
  float MaxRadius;
  float AngleCut;
  int NSubHits;
  int NRadialIntervals;
};

struct GeoParams {
  float MaxResRange;
  int SlidingWindowSize;
  float MaxDEDx;
};

struct MichelParams {
  GeoParams GeoPset;
  RadialParams RadialPset;
  float AngleMaxRadius;
  int BlipMaxWireSpan;
  float BlipMaxDist;
};


struct MichelFeatures {
  bool featuresComputed = false;
  int nMuonTracks = -1;
  double lowerChi2Ratio = -1;
  double allChi2Ratio = -1;
  double angleQ80 = -1;
  double dropRadius = -1;
  int nBlips = -1;
  double score = -999;
  bool muonInFV = false;

  void Print() const
    {
        std::cout 
        << "NMuonTracks = " << nMuonTracks
        << ", lowerChi2Ratio = " << lowerChi2Ratio
        << ", angleQ80 = " << angleQ80
        << ", dropRadius = " << dropRadius
        << ", nBlips = " << nBlips
        << ", score = " << score
        << ", muonInFV = " << muonInFV
        << std::endl;
    }

  void clear()  
    {
      featuresComputed = false;
      nMuonTracks = -1;
      lowerChi2Ratio = -1;
      allChi2Ratio = -1;
      angleQ80 = -1;
      dropRadius = -1;
      nBlips = -1;
      muonInFV = false;
    }
};
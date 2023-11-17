#ifndef SECONDSHOWERFINDERALG_H_SEEN
#define SECONDSHOWERFINDERALG_H_SEEN

///////////////////////////////////////////////
// SecondShowerFinderAlg.h
//
// Alg intended to find the subleading photon
// missed from the reconstruction of neutral
// pion candidates.
//
// Author: Henry Lay (h.lay@lancaster.ac.uk)
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Shower.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TH1F.h"
#include "TGraph.h"

#include "sbndcode/PiZero/NCPiZeroStructs.h"

typedef std::vector<art::Ptr<recob::Hit>> HitVec;

class SecondShowerFinderAlg
{
  struct HitObj {
    art::Ptr<recob::Hit> hit;
    bool used;
    double x;
    double wire_pos;
  };

  typedef std::vector<HitObj*> ClusterObj;

  std::vector<int> fColours = { kGreen+2, kMagenta+2, kCyan+2, kYellow+2, kAzure+1, kSpring+9, kPink+9, kTeal+9, kOrange+2, kViolet-8 };

 public:
  SecondShowerFinderAlg();

  SecondShowerFinderAlg(fhicl::ParameterSet const& p);

  bool FindSecondShower(const art::Event &e, const HitVec &hits, const HitVec &usedHits, const bool draw);

  bool AnalyseViewHits(const art::Event &e, const HitVec &hits, const HitVec &usedHits, const TString &name, const bool draw);

  void DrawView(const HitVec &hits, const HitVec &usedHits, const std::vector<ClusterObj> clusters, const TString &name);

  double YZtoU(const double y, const double z);

  double YZtoV(const double y, const double z);

  double YZtoW(const double y, const double z);

  void MergeClusters(std::vector<ClusterObj> &clusters);

 private:
  const double wireAngleU = 1.04719758034;
  const double wireAngleV = -1.04719758034;
  const double wireAngleW = 0;

  const double cosU = TMath::Cos(wireAngleU);
  const double sinU = TMath::Sin(wireAngleU);
  const double cosV = TMath::Cos(wireAngleV);
  const double sinV = TMath::Sin(wireAngleV);
  const double cosW = TMath::Cos(wireAngleW);
  const double sinW = TMath::Sin(wireAngleW);
};

#endif

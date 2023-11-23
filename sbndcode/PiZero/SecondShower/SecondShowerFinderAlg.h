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
#include "TLine.h"

#include "sbndcode/PiZero/NCPiZeroStructs.h"

typedef std::vector<art::Ptr<recob::Hit>> HitVec;

class SecondShowerFinderAlg
{
  const double wireAngleU = 1.04719758034;
  const double wireAngleV = -1.04719758034;
  const double wireAngleW = 0;

  const double cosU = TMath::Cos(wireAngleU);
  const double sinU = TMath::Sin(wireAngleU);
  const double cosV = TMath::Cos(wireAngleV);
  const double sinV = TMath::Sin(wireAngleV);
  const double cosW = TMath::Cos(wireAngleW);
  const double sinW = TMath::Sin(wireAngleW);

  const double sinVminusU = TMath::Sin(wireAngleV - wireAngleU);
  const double sinWminusV = TMath::Sin(wireAngleW - wireAngleV);
  const double sinUminusW = TMath::Sin(wireAngleU - wireAngleW);

  struct HitObj {
    art::Ptr<recob::Hit> hit;
    bool used;
    double x;
    double wire_pos;
  };

  typedef std::vector<HitObj*> HitObjVec;

  class ClusterObj {

    HitObjVec hits;
    bool used;
    double min_x;
    double max_x;
    double min_wire_pos;
    double max_wire_pos;
    int tpc;

  public:
    ClusterObj(HitObjVec hitObjVec);

    HitObjVec& Hits() { return hits; }

    const HitObjVec& ConstHits() const { return hits; }

    size_t Size() const { return hits.size(); }

    double MinX() const { return min_x; }

    double MaxX() const { return max_x; }

    double MinWirePos() const { return min_wire_pos; }

    double MaxWirePos() const { return max_wire_pos; }

    int TPC() const { return tpc; }
  };

  typedef std::vector<ClusterObj*> ClusterObjVec;

  std::vector<int> fColours = { kGreen+2, kMagenta+2, kCyan+2, kYellow+2, kAzure+1, kSpring+9, kPink+9, kTeal+9, kOrange+2, kViolet-8 };

  ClusterObjVec ClusterInView(const HitObjVec &hits, const HitObjVec &usedHits, const TString &name, const bool draw);

  void SeparateViews(const art::Event &e, const HitVec &hits, HitObjVec &u_hits, HitObjVec &v_hits, HitObjVec &w_hits);

  void InitialPairings(const HitObjVec &hits, std::vector<HitObjVec> &hitCollections);

  void AddSingleHits(const HitObjVec &hits, std::vector<HitObjVec> &hitCollections);

  void MergeHitCollections(std::vector<HitObjVec> &hitCollections);

  void RemoveClusterBelowLimit(ClusterObjVec &clusters, const size_t limit);

  void TwoDToThreeDMatching(std::vector<ClusterObjVec> &clusters, const bool draw);

  void DrawView(const HitObjVec &hits, const HitObjVec &usedHits, const ClusterObjVec clusters, const TString &name);

  void DrawClusterMatching(const ClusterObj *clusterU, const ClusterObj *clusterV, const ClusterObj *clusterW, const double startX, const double endX);

  void DrawClusterMatchingView(const ClusterObj *cluster, const double maxRangeWirePos, const double rangeWirePos,
                               const double minX, const double maxX, const double startX, const double endX, const TString &name);

  double YZtoU(const double y, const double z, const int tpc);

  double YZtoV(const double y, const double z, const int tpc);

  double YZtoW(const double y, const double z, const int tpc);

  double UVtoW(const double u, const double v, const int tpc);

  double VWtoU(const double v, const double w, const int tpc);

  double WUtoV(const double w, const double u, const int tpc);

  double Dist(const HitObj *hitObjA, const HitObj *hitObjB);

  double GetInterpolatedHitWirePos(const ClusterObj *cluster, const double &x);

  size_t fMinClusterHits;
  double fMaxHitSeparation;

 public:
  SecondShowerFinderAlg();

  SecondShowerFinderAlg(fhicl::ParameterSet const& p);

  std::vector<std::vector<size_t>> FindSecondShower(const art::Event &e, const HitVec &hits, const HitVec &usedHits, const bool draw);
};

#endif

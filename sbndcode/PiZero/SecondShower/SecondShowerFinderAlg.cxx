#include "SecondShowerFinderAlg.h"

SecondShowerFinderAlg::ClusterObj::ClusterObj(HitObjVec hitObjVec)
{
  hits = hitObjVec;
  used = false;

  min_x = std::numeric_limits<double>::max();
  max_x = std::numeric_limits<double>::lowest();

  min_wire_pos = std::numeric_limits<double>::max();
  max_wire_pos = std::numeric_limits<double>::lowest();

  std::set<int> tpcs;

  for(auto const& hitObj : hitObjVec)
    {
      if(hitObj->x < min_x)
        min_x = hitObj->x;

      if(hitObj->x > max_x)
        max_x = hitObj->x;

      if(hitObj->wire_pos < min_wire_pos)
        min_wire_pos = hitObj->wire_pos;

      if(hitObj->wire_pos > max_wire_pos)
        max_wire_pos = hitObj->wire_pos;

      tpcs.insert(hitObj->hit->WireID().asTPCID().deepestIndex());
    }

  if(tpcs.size() == 1)
    tpc = *tpcs.begin();
  else if(tpcs.size() == 2)
    tpc = 2;
  else
    {
      std::cout << "Uh oh... we have: " << tpcs.size() << " TPCs" << std::endl;
      tpc = -1;
    }
}

SecondShowerFinderAlg::SecondShowerFinderAlg()
{
  gROOT->SetBatch(false);

  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameLineWidth(4);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(62);
  gStyle->SetLegendTextSize(0.09);

  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.02);

  gStyle->SetTextFont(62);
  gStyle->SetTextSize(0.07);
  gStyle->SetLabelFont(62,"x");
  gStyle->SetLabelFont(62,"y");
  gStyle->SetLabelFont(62,"z");
  gStyle->SetLabelSize(0,"x");
  gStyle->SetTitleSize(0,"x");
  gStyle->SetLabelSize(0,"y");
  gStyle->SetTitleSize(0,"y");
  gStyle->SetLabelSize(0,"z");
  gStyle->SetTitleSize(0,"z");
  gStyle->SetLabelFont(62,"t");
  gStyle->SetTitleFont(62,"x");
  gStyle->SetTitleFont(62,"y");
  gStyle->SetTitleFont(62,"z");
  gStyle->SetTitleFont(62,"t");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleAlign(23);
  gStyle->SetTitleOffset(1,"y");
  gStyle->SetTitleOffset(1,"z");
  gStyle->SetTitleFontSize(0.07);
  gStyle->SetTitleFont(62,"pad");
  gStyle->SetTitleBorderSize(0);

  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(.25);
  gStyle->SetHistLineWidth(0);
  gStyle->SetLineStyleString(2,"[12 12]");

  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetNdivisions(5000, "X");
  gStyle->SetNdivisions(1000, "Y");

  gStyle->SetPalette(kRainBow);
}

SecondShowerFinderAlg::SecondShowerFinderAlg(fhicl::ParameterSet const& p)
{
  SecondShowerFinderAlg();

  fMinClusterHits   = p.get<size_t>("MinClusterHits");
  fMaxHitSeparation = p.get<double>("MaxHitSeparation");
}

std::vector<std::vector<size_t>> SecondShowerFinderAlg::FindSecondShower(const art::Event &e, const HitVec &hits, const HitVec &usedHits, const bool draw)
{
  HitObjVec u_hits, v_hits, w_hits;
  HitObjVec u_usedHits, v_usedHits, w_usedHits;

  SeparateViews(e, hits, u_hits, v_hits, w_hits);
  SeparateViews(e, usedHits, u_usedHits, v_usedHits, w_usedHits);

  if(draw)
    std::cout << "Extra Hits: " << hits.size()
              << "\tU: " << u_hits.size()
              << "\tV: " << v_hits.size()
              << "\tW: " << w_hits.size() << std::endl;

  std::vector<ClusterObjVec> clusters(3, ClusterObjVec());

  clusters[0] = ClusterInView(u_hits, u_usedHits, "U view", draw);
  clusters[1] = ClusterInView(v_hits, v_usedHits, "V view", draw);
  clusters[2] = ClusterInView(w_hits, w_usedHits, "W view", draw);

  TwoDToThreeDMatching(clusters, draw);

  std::vector<std::vector<size_t>> clusterSizes(3, std::vector<size_t>());

  return clusterSizes;
}

SecondShowerFinderAlg::ClusterObjVec SecondShowerFinderAlg::ClusterInView(const HitObjVec &hits, const HitObjVec &usedHits, const TString &name, const bool draw)
{
  if(draw)
    {
      ClusterObjVec tmpClusters;
      std::cout << "Drawing " << name << " pre-clustering... (" << hits.size() << " unused hits)" << std::endl;
      DrawView(hits, usedHits, tmpClusters, name);
    }

  std::vector<HitObjVec> hitCollections;

  InitialPairings(hits, hitCollections);
  AddSingleHits(hits, hitCollections);
  MergeHitCollections(hitCollections);

  ClusterObjVec clusters;

  for(auto const& hitCollection : hitCollections)
    clusters.push_back(new ClusterObj(hitCollection));

  RemoveClusterBelowLimit(clusters, fMinClusterHits);

  if(draw)
    {
      std::cout << "Drawing " << name << " clustered... (" << clusters.size() << " extra clusters)" << std::endl;
      DrawView(hits, usedHits, clusters, name);
    }

  return clusters;
}

void SecondShowerFinderAlg::SeparateViews(const art::Event &e, const HitVec &hits, HitObjVec &u_hits, HitObjVec &v_hits, HitObjVec &w_hits)
{
  art::ServiceHandle<geo::Geometry const> geom;
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

  for(auto const& hit : hits)
    {
      const geo::WireID wireID = hit->WireID();
      const double x = detProp.ConvertTicksToX(hit->PeakTime(), wireID.Plane, wireID.TPC, wireID.Cryostat);

      const geo::Point_t wireCentre = geom->Wire(wireID).GetCenter();
      const int tpc = wireID.asTPCID().deepestIndex();

      switch(hit->View())
        {
        case geo::kU:
          {
            double wire_pos = YZtoU(wireCentre.Y(), wireCentre.Z(), tpc);
            HitObj *hitObj = new HitObj({hit, false, x, wire_pos});
            u_hits.push_back(hitObj);
            break;
          }
        case geo::kV:
          {
            double wire_pos = YZtoV(wireCentre.Y(), wireCentre.Z(), tpc);
            HitObj *hitObj = new HitObj({hit, false, x, wire_pos});
            v_hits.push_back(hitObj);
            break;
          }
        case geo::kW:
          {
            double wire_pos = YZtoW(wireCentre.Y(), wireCentre.Z(), tpc);
            HitObj *hitObj = new HitObj({hit, false, x, wire_pos});
            w_hits.push_back(hitObj);
            break;
          }
        default:
          std::cout << "Holy shit SBND, you got new planes" << std::endl;
          break;
        }
    }
}

void SecondShowerFinderAlg::InitialPairings(const HitObjVec &hits, std::vector<HitObjVec> &hitCollections)
{
  for(auto const& hitObjA : hits)
    {
      if(hitObjA->used)
        continue;

      for(auto const& hitObjB : hits)
        {
          if(hitObjA->hit == hitObjB->hit)
            continue;

          if(hitObjB->used)
            continue;

          double dist = Dist(hitObjA, hitObjB);

          if(dist < fMaxHitSeparation)
            {
              hitCollections.push_back(HitObjVec());
              hitCollections.back().push_back(hitObjA);
              hitCollections.back().push_back(hitObjB);

              hitObjA->used = true;
              hitObjB->used = true;
              break;
            }
        }
    }
}

void SecondShowerFinderAlg::AddSingleHits(const HitObjVec &hits, std::vector<HitObjVec> &hitCollections)
{
  for(auto const& hitObjA : hits)
    {
      if(hitObjA->used)
        continue;

      std::vector<HitObjVec>::iterator it = hitCollections.begin();

      while(it != hitCollections.end())
        {
          bool add = false;

          for(auto const& hitObjB : *it)
            {
              double dist = Dist(hitObjA, hitObjB);

              if(dist < fMaxHitSeparation)
                add = true;
            }

          if(add)
            {
              it->push_back(hitObjA);
              hitObjA->used = true;
            }

          ++it;
        }
    }
}

void SecondShowerFinderAlg::MergeHitCollections(std::vector<HitObjVec> &hitCollections)
{
  std::vector<HitObjVec>::iterator itA = hitCollections.begin();

  while(itA != hitCollections.end())
    {
      std::vector<HitObjVec>::iterator itB = std::next(itA);

      while(itB != hitCollections.end())
        {
          bool merge = false;

          for(auto const& hitObjA : *itA)
            {
              for(auto const& hitObjB : *itB)
                {
                  double dist = Dist(hitObjA, hitObjB);

                  if(dist < fMaxHitSeparation)
                    merge = true;
                }
            }

          if(merge)
            {
              itA->insert(itA->end(), itB->begin(), itB->end());
              itB = hitCollections.erase(itB);
            }
          else
            ++itB;
        }
      ++itA;
    }
}

void SecondShowerFinderAlg::RemoveClusterBelowLimit(ClusterObjVec &clusters, const size_t limit)
{
  ClusterObjVec::iterator it = clusters.begin();
  while(it != clusters.end())
    {
      if((*it)->Size() < limit)
        it = clusters.erase(it);
      else
        {
          ++it;
        }
    }
}

void SecondShowerFinderAlg::TwoDToThreeDMatching(std::vector<ClusterObjVec> &clusters, const bool draw)
{
  for(auto const& clusterU : clusters[0])
    {
      for(auto const& clusterV : clusters[1])
        {
          for(auto const& clusterW : clusters[2])
            {
              double startX = std::max({clusterU->MinX(), clusterV->MinX(), clusterW->MinX()});
              double endX   = std::min({clusterU->MaxX(), clusterV->MaxX(), clusterW->MaxX()});

              if(startX > endX)
                continue;

              if(draw)
                {
                  std::cout << "Drawing attempt to match cluster set" << std::endl;
                  DrawClusterMatching(clusterU, clusterV, clusterW, startX, endX);
                }
            }
        }
    }
}

void SecondShowerFinderAlg::DrawView(const HitObjVec &hits, const HitObjVec &usedHits, const ClusterObjVec clusters, const TString &name)
{
  TCanvas *c = new TCanvas("c", "", 2100, 1400);
  c->cd();

  TPaveText *t = new TPaveText(.9, .95, .98, 1, "NB NDC");
  t->SetTextAlign(32);
  t->SetTextSize(0.085);
  t->SetTextColor(kBlack);
  t->SetFillStyle(4000);
  t->AddText(name);

  TPaveText *t0 = new TPaveText(.03, .92, .2, .95, "NB NDC");
  t0->SetTextAlign(12);
  t0->SetTextSize(0.06);
  t0->SetTextColor(kBlack);
  t0->SetFillStyle(4000);
  t0->AddText("TPC 0");

  TPaveText *t1 = new TPaveText(.03, .82, .2, .85, "NB NDC");
  t1->SetTextAlign(12);
  t1->SetTextSize(0.06);
  t1->SetTextColor(kBlack);
  t1->SetFillStyle(4000);
  t1->AddText("TPC 1");

  TPad *top    = new TPad("top", "", 0.005, 0.505, 0.995, 0.995);
  TPad *bottom = new TPad("bottom", "", 0.005, 0.005, 0.995, 0.495);
  top->SetBottomMargin(0.002);
  bottom->SetTopMargin(0.002);
  top->Draw();
  bottom->Draw();

  top->cd();

  const double n_wires = name == "W view" ? 1664 : 1984;

  TH1F* base_hist = new TH1F("base_hist", "", 10, 0.5 - 0.5*(1984 - n_wires), 1984.5 - 0.5*(1984 - n_wires));
  base_hist->SetMaximum(3400.0);
  base_hist->Draw();
  t->Draw();
  t1->Draw();

  TGraph *g0 = new TGraph();
  g0->SetMarkerColor(kBlue+2);
  TGraph *g1 = new TGraph();
  g1->SetMarkerColor(kBlue+2);
  TGraph *g0Used = new TGraph();
  g0Used->SetMarkerColor(kRed+2);
  TGraph *g1Used = new TGraph();
  g1Used->SetMarkerColor(kRed+2);

  std::vector<TGraph*> g0Clusters;
  std::vector<TGraph*> g1Clusters;

  for(auto const& hitObj : hits)
    {
      const int tpc = hitObj->hit->WireID().asTPCID().deepestIndex();

      if(tpc == 0)
        g0->SetPoint(g0->GetN(), hitObj->hit->WireID().deepestIndex(), hitObj->hit->PeakTime());
      else if(tpc == 1)
        g1->SetPoint(g1->GetN(), hitObj->hit->WireID().deepestIndex(), 3400 - hitObj->hit->PeakTime());
    }

  for(auto const& hitObj : usedHits)
    {
      const int tpc = hitObj->hit->WireID().asTPCID().deepestIndex();

      if(tpc == 0)
        g0Used->SetPoint(g0Used->GetN(), hitObj->hit->WireID().deepestIndex(), hitObj->hit->PeakTime());
      else if(tpc == 1)
        g1Used->SetPoint(g1Used->GetN(), hitObj->hit->WireID().deepestIndex(), 3400 - hitObj->hit->PeakTime());
    }

  if(g1->GetN())
    g1->Draw("Psame");
  if(g1Used->GetN())
    g1Used->Draw("Psame");

  bottom->cd();
  base_hist->Draw();
  t0->Draw();
  if(g0->GetN())
    g0->Draw("Psame");
  if(g0Used->GetN())
    g0Used->Draw("Psame");

  while(clusters.size() > fColours.size())
    fColours.insert(fColours.end(), fColours.begin(), fColours.end());

  for(auto&& [i, cluster] : enumerate(clusters))
    {
      std::cout << "\tCluster " << i << " of " << cluster->Size() << " hits" << std::endl;
      g0Clusters.push_back(new TGraph());
      g1Clusters.push_back(new TGraph());

      const int colour = fColours[i];
      g0Clusters[i]->SetMarkerColor(colour);
      g1Clusters[i]->SetMarkerColor(colour);

      for(auto const& hitObj : cluster->Hits())
        {
          const int tpc = hitObj->hit->WireID().asTPCID().deepestIndex();

          if(tpc == 0)
            g0Clusters[i]->SetPoint(g0Clusters[i]->GetN(), hitObj->hit->WireID().deepestIndex(), hitObj->hit->PeakTime());
          else if(tpc == 1)
            g1Clusters[i]->SetPoint(g1Clusters[i]->GetN(), hitObj->hit->WireID().deepestIndex(), 3400 - hitObj->hit->PeakTime());
        }

      top->cd();
      if(g1Clusters[i]->GetN())
        g1Clusters[i]->Draw("Psame");
      bottom->cd();
      if(g0Clusters[i]->GetN())
        g0Clusters[i]->Draw("Psame");
    }

  c->Update();
  gSystem->ProcessEvents();

  std::cin.get();

  delete c;
}

void SecondShowerFinderAlg::DrawClusterMatching(const ClusterObj *clusterU, const ClusterObj *clusterV, const ClusterObj *clusterW, const double startX, const double endX)
{
  TCanvas *c = new TCanvas("c", "", 2100, 2100);
  c->cd();

  const double padPadding = 0.05;
  const double padWidth   = (1 - 4 * padPadding) / 3.;

  TPad *padU = new TPad("padU", "", padPadding, padPadding, padPadding + padWidth, 1 - padPadding);
  TPad *padV = new TPad("padV", "", 2 * padPadding + padWidth, padPadding, 2 * padPadding + 2 * padWidth, 1 - padPadding);
  TPad *padW = new TPad("padW", "", 3 * padPadding + 2 * padWidth, padPadding, 3 * padPadding + 3 * padWidth, 1 - padPadding);

  padU->SetBottomMargin(0.002);
  padU->SetTopMargin(0.002);
  padV->SetBottomMargin(0.002);
  padV->SetTopMargin(0.002);
  padW->SetBottomMargin(0.002);
  padW->SetTopMargin(0.002);

  padU->Draw();
  padV->Draw();
  padW->Draw();

  double minX = std::min({clusterU->MinX(), clusterV->MinX(), clusterW->MinX()});
  double maxX = std::max({clusterU->MaxX(), clusterV->MaxX(), clusterW->MaxX()});

  if(minX > maxX)
    std::swap(minX, maxX);

  const double rangeU = clusterU->MaxWirePos() - clusterU->MinWirePos();
  const double rangeV = clusterV->MaxWirePos() - clusterV->MinWirePos();
  const double rangeW = clusterW->MaxWirePos() - clusterW->MinWirePos();

  const double rangeWirePos = 1.2 * std::max({rangeU, rangeV, rangeW});

  padU->cd();
  DrawClusterMatchingView(clusterU, rangeWirePos, rangeU, minX, maxX, startX, endX, "U View");
  padV->cd();
  DrawClusterMatchingView(clusterV, rangeWirePos, rangeV, minX, maxX, startX, endX, "V view");
  padW->cd();
  DrawClusterMatchingView(clusterW, rangeWirePos, rangeW, minX, maxX, startX, endX, "W View");

  c->Update();
  gSystem->ProcessEvents();

  std::cin.get();

  delete c;
}

void SecondShowerFinderAlg::DrawClusterMatchingView(const ClusterObj *cluster, const double maxRangeWirePos, const double rangeWirePos,
                                                    const double minX, const double maxX, const double startX, const double endX, const TString &name)
{
  TPaveText *t = new TPaveText(.9, .92, .98, .97, "NB NDC");
  t->SetTextAlign(32);
  t->SetTextSize(0.065);
  t->SetTextColor(kBlack);
  t->SetFillStyle(4000);
  t->AddText(name);

  const double rangeX = maxX - minX;

  TH1F* base_hist = new TH1F("base_hist", "", 10, cluster->MinWirePos() - 0.5 * (maxRangeWirePos - rangeWirePos), cluster->MaxWirePos() + 0.5 * (maxRangeWirePos - rangeWirePos));
  base_hist->SetMinimum(minX - 0.1 * rangeX);
  base_hist->SetMaximum(maxX + 0.1 * rangeX);
  base_hist->Draw();

  TLine *startLine = new TLine();
  startLine->SetLineColor(kGreen+2);
  startLine->SetLineWidth(5);
  startLine->DrawLine(cluster->MinWirePos() - 0.5 * (maxRangeWirePos - rangeWirePos), startX, cluster->MaxWirePos() + 0.5 * (maxRangeWirePos - rangeWirePos), startX);

  TLine *endLine = new TLine();
  endLine->SetLineColor(kRed+2);
  endLine->SetLineWidth(5);
  endLine->DrawLine(cluster->MinWirePos() - 0.5 * (maxRangeWirePos - rangeWirePos), endX, cluster->MaxWirePos() + 0.5 * (maxRangeWirePos - rangeWirePos), endX);

  TGraph *g = new TGraph();
  g->SetMarkerColor(kMagenta+2);
  g->SetMarkerSize(1.);

  for(auto const& hitObj : cluster->ConstHits())
    g->SetPoint(g->GetN(), hitObj->wire_pos, hitObj->x);

  g->Draw("Psame");
  t->Draw();
}

// Note the use of 'V' angles in TPC0 is due to the fact that view is defined w.r.t. drift direction
// whereas this conversion is working in detector coordinates in which x is pointing the same direction
// regardless of which TPC you're in. The same is true for YZtoV, W doesn't care.
double SecondShowerFinderAlg::YZtoU(const double y, const double z, const int tpc)
{
  if(tpc == 0)
    return z * cosU - y * sinU;
  else if(tpc == 1)
    return z * cosV - y * sinV;

  return std::numeric_limits<double>::lowest();
}

double SecondShowerFinderAlg::YZtoV(const double y, const double z, const int tpc)
{
  if(tpc == 0)
    return z * cosV - y * sinV;
  else if(tpc == 1)
    return z * cosU - y * sinU;

  return std::numeric_limits<double>::lowest();
}

double SecondShowerFinderAlg::YZtoW(const double y, const double z, const int /*tpc*/)
{
  return z * cosW - y * sinW;
}

double SecondShowerFinderAlg::UVtoW(const double u, const double v, const int tpc)
{
  if(tpc == 0)
    return -1. * (u * sinWminusV + v * sinUminusW) / sinVminusU;
  else if(tpc == 1)
    return -1. * (v * sinWminusV + u * sinUminusW) / sinVminusU;

  return std::numeric_limits<double>::lowest();
}

double SecondShowerFinderAlg::VWtoU(const double v, const double w, const int tpc)
{
  if(tpc == 0)
    return -1. * (v * sinUminusW + w * sinVminusU) / sinWminusV;
  else if(tpc == 1)
    return -1. * (v * sinWminusV + w * sinVminusU) / sinUminusW;

  return std::numeric_limits<double>::lowest();
}

double SecondShowerFinderAlg::WUtoV(const double w, const double u, const int tpc)
{
  if(tpc == 0)
    return -1. * (u * sinWminusV + w * sinVminusU) / sinUminusW;
  else if(tpc == 1)
    return -1. * (u * sinUminusW + w * sinVminusU) / sinWminusV;

  return std::numeric_limits<double>::lowest();
}

double SecondShowerFinderAlg::Dist(const HitObj *hitObjA, const HitObj *hitObjB)
{
  return sqrt((hitObjA->x - hitObjB->x) * (hitObjA->x - hitObjB->x) +
              (hitObjA->wire_pos - hitObjB->wire_pos) * (hitObjA->wire_pos - hitObjB->wire_pos));
}

double SecondShowerFinderAlg::GetInterpolatedHitWirePos(const ClusterObj *cluster, const double &x)
{
  if(x < cluster->MinX() || x > cluster->MaxX())
    return std::numeric_limits<double>::lowest();

  double closestLowerXDiff     = std::numeric_limits<double>::max();
  double closestLowerXWirePos  = std::numeric_limits<double>::lowest();
  double closestHigherXDiff    = std::numeric_limits<double>::max();
  double closestHigherXWirePos = std::numeric_limits<double>::lowest();

  bool foundLower = false, foundHigher = false;

  for(auto const& hitObj : cluster->ConstHits())
    {
      if(hitObj->x <= x && std::abs(x - hitObj->x) < closestLowerXDiff)
        {
          closestLowerXDiff    = std::abs(x - hitObj->x);
          closestLowerXWirePos = hitObj->wire_pos;
          foundLower = true;
        }

      if(hitObj->x >= x && std::abs(x - hitObj->x) < closestHigherXDiff)
        {
          closestHigherXDiff    = std::abs(x - hitObj->x);
          closestHigherXWirePos = hitObj->wire_pos;
          foundHigher = true;
        }
    }

  if(foundLower && foundHigher)
    {
      double fractionalDistance = closestLowerXDiff + closestHigherXDiff < std::numeric_limits<double>::epsilon() ? 0 :
        closestLowerXDiff / (closestLowerXDiff + closestHigherXDiff);

      return closestLowerXWirePos + fractionalDistance * (closestHigherXWirePos - closestLowerXWirePos);
    }
  else if(foundLower)
    std::cout << "Just lower....!" << std::endl;
  else if(foundHigher)
    std::cout << "Just higher....!" << std::endl;
  else
    std::cout << "Found neither...!" << std::endl;

  return std::numeric_limits<double>::lowest();
}

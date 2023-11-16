#include "SecondShowerFinderAlg.h"

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
}

bool SecondShowerFinderAlg::FindSecondShower(const art::Event &e, const HitVec &hits, const HitVec &usedHits, const bool draw)
{
  HitVec u_hits, v_hits, w_hits;
  HitVec u_usedHits, v_usedHits, w_usedHits;

  for(auto const& hit : hits)
    {
      switch(hit->View())
        {
        case geo::kU:
          u_hits.push_back(hit);
          break;
        case geo::kV:
          v_hits.push_back(hit);
          break;
        case geo::kW:
          w_hits.push_back(hit);
          break;
        default:
          std::cout << "Holy shit SBND, you got new planes" << std::endl;
          break;
        }
    }

  for(auto const& hit : usedHits)
    {
      switch(hit->View())
        {
        case geo::kU:
          u_usedHits.push_back(hit);
          break;
        case geo::kV:
          v_usedHits.push_back(hit);
          break;
        case geo::kW:
          w_usedHits.push_back(hit);
          break;
        default:
          std::cout << "Holy shit SBND, you got new planes" << std::endl;
          break;
        }
    }

  std::cout << "Extra Hits: " << hits.size()
            << "\tU: " << u_hits.size()
            << "\tV: " << v_hits.size()
            << "\tW: " << w_hits.size() << std::endl;

  bool u_yes = AnalyseViewHits(e, u_hits, u_usedHits, "U view", draw);
  bool v_yes = AnalyseViewHits(e, v_hits, v_usedHits, "V view", draw);
  bool w_yes = AnalyseViewHits(e, w_hits, w_usedHits, "W view", draw);

  if(u_yes && v_yes && w_yes)
    std::cout << "YES!- GOOD JOB" << std::endl;

  return false;
}

bool SecondShowerFinderAlg::AnalyseViewHits(const art::Event &e, const HitVec &hits, const HitVec &usedHits, const TString &name, const bool draw)
{
  art::ServiceHandle<geo::Geometry const> geom;
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(e);

  ClusterObj allHits;
  std::vector<ClusterObj> clusters;

  if(draw)
    DrawView(hits, usedHits, clusters, name);

  for(auto const& hit : hits)
    {
      const geo::WireID wireID = hit->WireID();
      const double x = detProp.ConvertTicksToX(hit->PeakTime(), wireID.Plane, wireID.TPC, wireID.Cryostat);

      auto const xyz = geom->Wire(wireID).GetCenter();
      double wire_pos = std::numeric_limits<double>::lowest();

      switch(hit->View())
        {
        case geo::kU:
          wire_pos = YZtoU(xyz.Y(), xyz.Z());
          break;
        case geo::kV:
          wire_pos = YZtoV(xyz.Y(), xyz.Z());
          break;
        case geo::kW:
          wire_pos = YZtoW(xyz.Y(), xyz.Z());
          break;
        default:
          std::cout << "Holy shit SBND, you got new planes" << std::endl;
          break;
        }

      allHits.push_back(new HitObj({hit, false, x, wire_pos}));
    }

  for(auto const& hitObjA : allHits)
    {
      if(hitObjA->used)
        continue;

      for(auto const& hitObjB : allHits)
        {
          if(hitObjA->hit == hitObjB->hit)
            continue;

          if(hitObjB->used)
            continue;

          double dist = sqrt( (hitObjA->x - hitObjB->x) * (hitObjA->x - hitObjB->x) +
                              (hitObjA->wire_pos - hitObjB->wire_pos) * (hitObjA->wire_pos - hitObjB->wire_pos));

          if(dist < 0.9)
            {
              clusters.push_back(ClusterObj());
              clusters.back().push_back(hitObjA);
              clusters.back().push_back(hitObjB);

              hitObjA->used = true;
              hitObjB->used = true;
            }
        }
    }

  int total = 0;

  if(draw)
    {
      for(auto const& cluster : clusters)
        total += cluster.size();

      std::cout << allHits.size() << " " << clusters.size() << " (" << total << ")\n\t";
      for(auto const& cluster : clusters)
        std::cout << cluster.size() << " ";
      std::cout << std::endl;
    }

  MergeClusters(clusters);

  if(draw)
    {
      total = 0;
      for(auto const& cluster : clusters)
        total += cluster.size();

      std::cout << allHits.size() << " " << clusters.size() << " (" << total << ")\n\t";
      for(auto const& cluster : clusters)
        std::cout << cluster.size() << " ";
      std::cout << std::endl;
    }

  std::vector<ClusterObj>::iterator it = clusters.begin();
  while(it != clusters.end())
    {
      if(it->size() < 10)
        it = clusters.erase(it);
      else
        {
          ++it;
        }
    }

  if(draw)
    {
      total = 0;
      for(auto const& cluster : clusters)
        total += cluster.size();

      std::cout << allHits.size() << " " << clusters.size() << " (" << total << ")\n\t";
      for(auto const& cluster : clusters)
        std::cout << cluster.size() << " ";
      std::cout << std::endl;

      DrawView(hits, usedHits, clusters, name);
    }

  return clusters.size();
}

void SecondShowerFinderAlg::DrawView(const HitVec &hits, const HitVec &usedHits, const std::vector<ClusterObj> clusters, const TString &name)
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

  std::vector<TGraph*> g0Clusters(clusters.size(), new TGraph());
  std::vector<TGraph*> g1Clusters(clusters.size(), new TGraph());

  for(auto const& hit : hits)
    {
      const int tpc = hit->WireID().asTPCID().deepestIndex();

      if(tpc == 0)
        g0->SetPoint(g0->GetN(), hit->WireID().deepestIndex(), hit->PeakTime());
      else if(tpc == 1)
        g1->SetPoint(g1->GetN(), hit->WireID().deepestIndex(), 3400 - hit->PeakTime());
    }

  for(auto const& hit : usedHits)
    {
      const int tpc = hit->WireID().asTPCID().deepestIndex();

      if(tpc == 0)
        g0Used->SetPoint(g0Used->GetN(), hit->WireID().deepestIndex(), hit->PeakTime());
      else if(tpc == 1)
        g1Used->SetPoint(g1Used->GetN(), hit->WireID().deepestIndex(), 3400 - hit->PeakTime());
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

  int color = 6;

  for(auto&& [i, cluster] : enumerate(clusters))
    {
      g0Clusters[i]->SetMarkerColor(color);
      g1Clusters[i]->SetMarkerColor(color);

      for(auto const& hitObj : cluster)
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

      std::cout << i << " - color: " << color << std::endl;

      ++color;
    }

  c->Update();
  gSystem->ProcessEvents();

  std::cin.get();

  delete c;
}

double SecondShowerFinderAlg::YZtoU(const double y, const double z)
{
  return z * cosU - y * sinU;
}

double SecondShowerFinderAlg::YZtoV(const double y, const double z)
{
  return z * cosV - y * sinV;
}

double SecondShowerFinderAlg::YZtoW(const double y, const double z)
{
  return z * cosW - y * sinW;
}

void SecondShowerFinderAlg::MergeClusters(std::vector<ClusterObj> &clusters)
{
  std::vector<ClusterObj>::iterator itA = clusters.begin();

  while(itA != clusters.end())
    {
      std::vector<ClusterObj>::iterator itB = clusters.begin();

      while(itB != clusters.end())
        {
          if(itA == itB)
            {
              ++itB;
              continue;
            }

          bool merge = false;

          for(auto const& hitObjA : *itA)
            {
              for(auto const& hitObjB : *itB)
                {
                  double dist = sqrt( (hitObjA->x - hitObjB->x) * (hitObjA->x - hitObjB->x) +
                                      (hitObjA->wire_pos - hitObjB->wire_pos) * (hitObjA->wire_pos - hitObjB->wire_pos));

                  if(dist < 0.9)
                    merge = true;
                }
            }

          if(merge)
            {
              itA->insert(itA->end(), itB->begin(), itB->end());
              itB = clusters.erase(itB);
            }
          else
            ++itB;
        }
      ++itA;
    }
}

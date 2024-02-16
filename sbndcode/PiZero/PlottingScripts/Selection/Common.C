#pragma once

#include "/exp/sbnd/app/users/hlay/plotting_utils/Structs.h"
#include "Enumerate.h"

const TString baseSaveDir = "/exp/sbnd/data/users/hlay/ncpizero/plots";
const TString baseFileDir = "/pnfs/sbnd/persistent/users/hlay/ncpizero";

const double goalPOT     = 10e20;
const double potPerSpill = 5e12;
const double goalSpills  = goalPOT / potPerSpill;

const double nTargets = 4.36e31;
const double intFlux  = 1.73e13;

double GetPOT(TChain *subruns)
{
  double sum = 0., pot = 0;

  subruns->SetBranchAddress("pot", &pot);

  for(size_t i = 0; i < subruns->GetEntries(); ++i)
    {
      subruns->GetEntry(i);
      sum += pot;
    }

  return sum;
}

int GetGenEvents(TChain *subruns)
{
  int sum = 0., ngenevts = 0;

  subruns->SetBranchAddress("ngenevts", &ngenevts);

  for(size_t i = 0; i < subruns->GetEntries(); ++i)
    {
      subruns->GetEntry(i);
      sum += ngenevts;
    }

  return sum;
}

void GetScaling(TChain *rockboxSubruns, TChain *intimeSubruns, double &rockboxScaling, double &intimeScaling)
{
  const double rockboxPOT = GetPOT(rockboxSubruns);
  const int rockboxSpills = GetGenEvents(rockboxSubruns);
  const int intimeSpills  = GetGenEvents(intimeSubruns);

  rockboxScaling = goalPOT / rockboxPOT;

  const double scaledRockboxSpills = rockboxScaling * rockboxSpills;
  
  intimeScaling = (goalSpills - scaledRockboxSpills) / intimeSpills;
}

TString POTString()
{
  TString potString = Form(" (%g POT)", goalPOT);
  potString.ReplaceAll("e+", "x10^{");
  potString.ReplaceAll(" POT", "} POT");

  return potString;
}

Cut TotalCut(const std::vector<Cut> &cuts)
{
  TCut totalCut = "";

  for(auto const& cut : cuts)
    totalCut += cut.cut;

  Cut cut = { "full_selection", totalCut, "Full Selection" };

  return cut;
}

Cut TotalCutExcluding(const std::vector<Cut> &cuts, const int excludingCut)
{
  TCut totalCut = "";

  for(auto&& [ cut_i, cut ] : enumerate(cuts))
    {
      if(cut_i != excludingCut)
        totalCut += cut.cut;
    }

  Cut cut = { "full_selection", totalCut, "Full Selection" };

  return cut;
}

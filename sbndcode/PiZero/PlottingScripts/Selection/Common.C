const double goalPOT     = 10e20;
const double potPerSpill = 5e12;
const double goalSpills  = goalPOT / potPerSpill;

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

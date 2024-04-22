#include "Common.C"

double totalpot = 0;
int totalevts = 0;

float PerFile(const TString fileName);

void PerFile()
{
  //  std::fstream file("/pnfs/sbnd/scratch/users/hlay/ncpizero/NCPiZeroBv3/ncpizero/ana/filesana.list");
  std::fstream file("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv17/ncpizero/ana/filesana.list");

  TH1F *hist = new TH1F("hist", ";NC1#pi^{0} Events at 1#times10^{21} POT;Files", 40, 200000, 350000);

  if(file.is_open())
    {
      std::string a;
      while(getline(file, a))
        hist->Fill(PerFile(a));
    }

  hist->Draw("histe");

  std::cout << totalevts * (1e21/totalpot) << std::endl;
}

float PerFile(const TString fileName)
{
  TChain *neutrinos = new TChain("ncpizeroxsectrees/neutrinos");
  neutrinos->AddFile(fileName);

  TChain *subruns = new TChain("ncpizeroxsectrees/subruns");
  subruns->AddFile(fileName);

  const double pot = GetPOT(subruns);

  const int signal = neutrinos->GetEntries("event_type_incl==0 || event_type_incl==1");

  std::cout << signal << " " << pot << " " << signal/pot << " " << signal*(1e21/pot) << std::endl;
  totalpot += pot;
  totalevts += signal;

  return signal*(1e21/pot);
}

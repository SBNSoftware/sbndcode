void ObservablesResolution(const TString productionVersion)
{
  const TString saveDir = "/sbnd/data/users/hlay/ncpizero/plots/" + productionVersion + "/observables/" + saveDirExt;
  gSystem->Exec("mkdir -p " + saveDir);

  const TString rockboxFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_rockbox.root";
  const TString intimeFile = "/pnfs/sbnd/persistent/users/hlay/ncpizero/" + productionVersion + "/" + productionVersion + "_intime.root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *rockboxEvents = new TChain("ncpizeroana/events");
  rockboxEvents->Add(rockboxFile);

  std::vector<std::vector<dou

  for(int i = 0; i < rockboxEvents->GetEntries(); ++i)
    {
      rockboxEvents->GetEntry(i);


    }
}

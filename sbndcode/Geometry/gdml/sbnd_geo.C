typedef struct _drawopt 
{
  const char* volume;
  int         color;
} drawopt;

sbnd_geo(TString volName="")
{
  gSystem->Load("libGeom");
  gSystem->Load("libGdml");

  TGeoManager::Import("sbnd_v00_08.gdml");

  drawopt optuboone[] = {
    {"volWorld",        0},
    {"volDetEnclosure", kWhite},
    {"volCryostat",     kOrange},
    {"volTPC",          kOrange-5},
    {"volTPCBackWall",  kRed},
    {"volTPCVertWall",  kCyan-5},
    {"volTPCHorizWall", kOrange},
    {0, 0}
  };

  // for (int i=0;; ++i) {
  //   if (optuboone[i].volume==0) break;
  //     gGeoManager->FindVolumeFast(optuboone[i].volume)->SetLineColor(optuboone[i].color);
  // }
  TList* mat = gGeoManager->GetListOfMaterials();
  TIter next(mat);
  TObject *obj;
  while (obj = next()) {
    obj->Print();
  }

  gGeoManager->GetTopNode();
  gGeoManager->CheckOverlaps(0.01);
  gGeoManager->PrintOverlaps();
  gGeoManager->SetMaxVisNodes(70000);

  //gGeoManager->GetTopVolume()->Draw();
  if ( ! volName.IsNull() ) gGeoManager->FindVolumeFast(volName)->Draw("ogl");

  TFile *tf = new TFile("sbnd.root", "RECREATE");
 
  gGeoManager->Write();

  tf->Close();
}

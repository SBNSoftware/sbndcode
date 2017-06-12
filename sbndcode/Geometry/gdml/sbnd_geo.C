#include "TGeoManager.h"

typedef struct {
  const char* volume;
  int color;
} drawopt;

int sbnd_geo(TString volName="") {
  gSystem->Load("libGeom");
  gSystem->Load("libGdml");

  TGeoManager::Import("sbnd_v00_08.gdml");

  drawopt opts[] = {
    {"volTPCActive", kGreen+2},
    {"volCryostat", kOrange+4},
    {0, 0}
  };

  for (int i=0;; i++) {
    if (opts[i].volume == 0) break;
    std::cout << "Painting '" << opts[i].volume << "'" << std::endl;
      gGeoManager->FindVolumeFast(opts[i].volume)->SetLineColor(opts[i].color);
  }

  TList* mat = gGeoManager->GetListOfMaterials();
  TIter next(mat);
  TObject* obj;
  while ((obj = next())) {
    obj->Print();
  }

  gGeoManager->GetTopNode();
//  gGeoManager->CheckOverlaps(10e-11);
  gGeoManager->PrintOverlaps();
  gGeoManager->SetMaxVisNodes(70000);

  if (!volName.IsNull())
    gGeoManager->FindVolumeFast(volName)->Draw("ogl");
  return 0;
  
  TGeoVolume* TPC = gGeoManager->FindVolumeFast("volTPCActive");
  float m_tpc = TPC->Weight();
  TGeoVolume* Cathode = gGeoManager->FindVolumeFast("volCathodePlate");
  float m_cathode = Cathode->Weight();

  float m_tpc_argon = m_tpc - m_cathode ;
  cout << "LAr weight in TPC = " << m_tpc_argon << " kg\n" <<endl;

  TFile* tf = new TFile("sbnd.root", "RECREATE");
  gGeoManager->Write();
  tf->Close();
  return 0;
}


typedef struct _drawopt {
  const char* volume;
  int         color;
} drawopt;

void lariat_geo_wdisk(TString volName="volCryostat"){

gSystem->Load("libGeom");
gSystem->Load("libGdml");

TGeoManager::Import("lariat_wdisk.gdml");

drawopt optArgoNeuT[] = {
//   {"volWorld",                 0},
//   {"volDetEnclosure",          kWhite},
//   {"volCryostat",              kOrange},
//   {"volTPCWirePlaneLengthSide", kCyan+3},
//   {"volTPCWirePlaneWidthSide", kRed},
  {"volTPCWire45", kRed},
  {"volTPCWire0", kBlue},
  {"volTPCWidthFace", kRed},
  {"volTPCLengthFace", kCyan+5},
  {"volTPCBottomFace", kOrange},
{"volTubBottom", kOrange+7},
{"voltheX", kOrange+7},
{"volTPCShieldPlane", kBlue},
{"volArgon_solid_L", kRed},
{"volArgon_cap_L", kOrange},
{"volArgon_cap_front", kOrange},
{"volTPC", kOrange},
{"volDetEnclosure", kBlue},
{"volMND", kBlue},
{"volTPCActive",kGreen},
  {0, 0}
};

for (int i=0;; ++i) {
  if (optArgoNeuT[i].volume==0) break;
    gGeoManager->FindVolumeFast(optArgoNeuT[i].volume)->SetLineColor(optArgoNeuT[i].color);
}

TList* mat = gGeoManager->GetListOfMaterials();
TIter next(mat);
TObject *obj;
// while (obj = next()) {
//  obj->Print();
// }

 gGeoManager->CheckOverlaps(0.01);
 gGeoManager->PrintOverlaps();
 gGeoManager->SetMaxVisNodes(70000);

 //gGeoManager->GetTopVolume()->Draw();
 //gGeoManager->FindVolumeFast(volName)->Draw();

 TFile *tf = new TFile("lariat_wdisk.root", "RECREATE");
 gGeoManager->Write();
 tf->Close();
}

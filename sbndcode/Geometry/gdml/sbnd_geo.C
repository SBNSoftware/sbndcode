#include "TGeoManager.h"

typedef struct {
  const char* volume;
  int color;
} drawopt;

int sbnd_geo(TString volName="")
{
	gSystem->Load("libGeom");
	gSystem->Load("libGdml");

	TGeoManager::Import("sbnd_v01_01.gdml");

drawopt optsbnd[] = {
	{"volTheBuilding", kYellow-8},
	{"volBuildingBase",kYellow-8},
	{"volGlassWindow", kBlue+4},
	{"volMezzanine",kYellow-8},
	{"volDetectorHall",kYellow-5},
	{"volGroundLevel0",kOrange+3},
	{"volGroundLevel1",kOrange+4},
	{"volGroundLevel2",kOrange+5},
	{"volOpDetSensitive", kYellow},
	{"volTPCActive",  kCyan-9},
	{"volTPCPlane_U", kBlue-4},
	{"volTPCPlane_V", kGreen-6},
	{"volTPCPlane_Y", kRed-7},
	{"volOneAPA", 	  kMagenta},	
	{"volInsulation", kMagenta},
	{"volStructureLid", kBlue},
	{"volAuxDetSensitiveCRT_X", kBlue},
	{"volAuxDetSensitiveCRT_Z", kBlue},
	{"volShieldingLid", kYellow-5},
	{"volShieldingTop", kWhite},
	{"volMezzanineLid", kWhite},
	{0, 0}
};


TGeoVolume *vol;
for (int i=0;; ++i) {
	if (optsbnd[i].volume==0) break;
	vol=gGeoManager->FindVolumeFast(optsbnd[i].volume);
	if(vol){
		vol->SetLineColor(optsbnd[i].color);
	//	vol->SetTransparency(optsbnd[i].transparency);
	}
	}
  TList* mat = gGeoManager->GetListOfMaterials();
  TIter next(mat);
  TObject* obj;
  while ((obj = next())) {
    obj->Print();
  }

  gGeoManager->GetTopNode();
  gGeoManager->CheckOverlaps(10e-11);
  gGeoManager->PrintOverlaps();
  gGeoManager->SetMaxVisNodes(70000);

  if (!volName.IsNull())
    gGeoManager->FindVolumeFast(volName)->Draw("ogl");
  
  TGeoVolume* TPC = gGeoManager->FindVolumeFast("volTPCActive");
  float m_tpc = TPC->Weight();
//  TGeoVolume* Cathode = gGeoManager->FindVolumeFast("volWorld");
//  float m_cathode = Cathode->Weight();

//  float m_tpc_argon = m_tpc - m_cathode ;
//  cout << "LAr weight in TPC = " << m_tpc_argon << " kg\n" <<endl;

  TFile* tf = new TFile("sbnd.root", "RECREATE");
  gGeoManager->Write();
  tf->Close();
  return 0;
}


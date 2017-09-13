typedef struct _drawopt 
{
  const char* volume;
  int         color;
  int	      transparency;
} drawopt;

void geoVis_sbnd(TString volName="volWorld")
{
	gSystem->Load("libGeom");
	gSystem->Load("libGdml");

	TGeoManager::Import("sbnd_v01_01_nowires.gdml");

drawopt optsbnd[] = {
	{"volTheBuilding", kYellow-8,0},
	{"volBuildingBase",kYellow-8,0},
	{"volGlassWindow", kBlue+4,50},
	{"volMezzanine",kYellow-8,0},
	{"volDetectorHall",kYellow-5,0},
	{"volGroundLevel0",kOrange+3,80},
	{"volGroundLevel1",kOrange+4,80},
	{"volGroundLevel2",kOrange+5,80},
	{"volOpDetSensitive", kYellow, 10},
	{"volTPCActive",  kCyan-9,50},
	{"volTPCPlane_U", kBlue-4,   80},
	{"volTPCPlane_V", kGreen-6,  80},
	{"volTPCPlane_Y", kRed-7,    80},
	{"volOneAPA", 	  kMagenta,    80},	
	{"volInsulation", kMagenta,    80},
	{"volStructureLid", kBlue-2,    0},
	{"volAuxDetSensitiveCRT_X", kBlue,0},
	{"volAuxDetSensitiveCRT_Z", kBlue,0},
	{"volShieldingLid", kYellow-5, 0},
	{"volShieldingTop", kWhite, 0},
	{"volMezzanineLid", kWhite, 0},
	{0, 0}
};


TGeoVolume *vol;
for (int i=0;; ++i) {
	if (optsbnd[i].volume==0) break;
	vol=gGeoManager->FindVolumeFast(optsbnd[i].volume);
	if(vol){
		vol->SetLineColor(optsbnd[i].color);
		vol->SetTransparency(optsbnd[i].transparency);
	}
}

gGeoManager->GetTopNode();
gGeoManager->CheckOverlaps(10e-11);
gGeoManager->PrintOverlaps();
gGeoManager->SetMaxVisNodes(70000);

if ( ! volName.IsNull() ) {gGeoManager->FindVolumeFast(volName)->Draw("ogl");} 
else  {gGeoManager->FindVolumeFast(volName)->Draw("ogl");}

float mLArTPC = gGeoManager->FindVolumeFast("volTPCActive")->Weight();
float mCPA = gGeoManager->FindVolumeFast("volCPA")->Weight();

//cout << mLArTPC << endl;
//cout << mCPA << endl;

}

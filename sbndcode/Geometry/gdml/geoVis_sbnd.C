#include "TString.h"
#include "TSystem.h"
#include "TGeoManager.h"
#include "TGeoTube.h"

#include <iostream>

typedef struct _drawopt 
{
  const char* volume;
  int         color;
  int	      transparency;
} drawopt;

void geoVis_sbnd(const TString &inputFileName,const TString &volName="volWorld")
{
	gSystem->Load("libGeom");
	gSystem->Load("libGdml");

	TGeoManager::Import(inputFileName);

drawopt optsbnd[] = {
	{"volBuildingBase", kGray+2,0},	
	{"volBuldingTube_H_1", kRed+1,0},
	{"volBuldingTube_H_2", kRed+1,0},
	{"volBuldingTube_V_1", kRed+1,0},
	{"volBuldingTube_V_2", kRed+1,0},
	{"volBuildingIBeam_V", kRed+1,0},
	{"volBuildingIBeam_H", kRed+1,0}, 
	{"volBuildingTranlucentWall", kBlue-10,40}, 
	{"volBuildingWindow", kBlue-10,60}, 
	{"volCeilingInsulation_External",kCyan+4,0},
	{"volCeilingInsulation",kWhite,0},	
	{"volCeilingSkylightCover",kYellow,80},
	{"volWallMetalPlate_NorthSouth",kGray+2,0},
	{"volWallMetalPlate_West",kGray+2,0},
	{"volWallMetalPlate_East",kGray+2,0}, //
	{"volWallInsulation_NorthSouth",kGray+2,0},
	{"volWallInsulation_West",kGray+2,0},
	{"volWallInsulation_East",kGray+2,0}, 
	{"volGlassWindow", kBlue+4,50},
	{"volMezzanine",kYellow-8,0},
	{"volDetectorHall",kYellow-5,0},
	{"volDummyTPCActive", kCyan-9, 10},
	{"volGroundLevel0",kOrange+3,80},
	{"volGroundLevel1",kOrange+4,80},
	{"volGroundLevel2",kOrange+5,80},
	{"volOpDetSensitive", kYellow, 10},
	{"volTPCActive_East",  kCyan-9,30},
	{"volTPCActive_West",  kCyan-9,30},
	{"volOneAPAFrame", 	  kMagenta-2,    0},	
	{"volCryoExt", kOrange, 0}, 
	{"volCryoInsulation", kOrange-10, 90},
	{"volCryoMembrane", kBlue-5,    0},
	{"volAuxDetSensitiveCRT_X", kBlue,0},
	{"volAuxDetSensitiveCRT_Z", kBlue,0},
	{"volShieldingLid", kYellow-5, 0},
	{"volShieldingTop", kWhite, 0},
	{"volMezzanineLid", kWhite, 0},
	{"volCryostatLidInsulation_Large", kBlue-9,    0},
//	{"volCryostatLid_Larger", kBlue-5,    0},
//	{"volDummyTPC", kCyan-9, 90},
//	{"volDummyCryostat", kOrange-3, 40},
//	{"volDummyDetEnclosure", kViolet+10, 80},
	{"volTPCPlane_U", kBlue-4,   80},
	{"volTPCPlane_V", kGreen-6,  80},
	{"volTPCPlaneVert", kRed-7,    80},
	{"volTPCWireVert", kRed-7, 0},
	{"volFC_IStructure", kPink-7, 0},
	{"volFC_strip_hor", kPink-6, 0},
	{"volFC_strip_ver", kPink-6, 0},
	{"volCPA_Around", kOrange-3, 0},
	{"volCPA_East", kOrange-3, 0},
	{"volCPA_West", kOrange-3, 0},
	{"volCPAMesh", kOrange-3, 90},
	{"volCPAFoil", kWhite, 90},
	{"volCPATPB", kCyan, 90},
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


struct wire{
	TString Name;
	float X,Y,Z;
	wire(const TString &name, const float &x,const float &y,const float &z) : Name(name), X(x), Y(y), Z(z) {}
};

std::vector<wire> wires;
TGeoNode *aNode;
bool checkWDist = true;
vol = gGeoManager->FindVolumeFast("volTPCPlaneVert");
if( vol ) { // check Y plane
	wires.clear();
	if(vol->GetNodes())for(int i=0; i<vol->GetNodes()->GetSize(); i++) {
		aNode = vol->GetNode(i);
		if(aNode){
			auto matrix = aNode->GetMatrix();
			wires.emplace_back( "volTPCWire_"+TString::Itoa(i+1,10), matrix->GetTranslation()[0], matrix->GetTranslation()[1], matrix->GetTranslation()[2] );
		} else {
			std::cerr << "Wire " << i+1 << " not found.\n";
		}
	}

	float dist, c1, c2, m = tan(-30.0f*M_PI/180.0f);
	TString wireName;
	for(int i=0; i<wires.size()-1; i++) {
		dist = abs(wires[i+1].Z-wires[i].Z);

		if( abs(dist*10.0f-3.0f) >= 1.0e-3 ) {
			std::cout  << "Distance between "<< wires[i].Name << " and " << wires[i+1].Name << " is " << dist*10.0f << "mm\n";
			checkWDist &= false;
		}
		
	}
} else {
	std::cout << "No Y wireplane was found.\n";
}

if(vol && checkWDist) std::cout << "Y plane: All wire distances equal to 3mm (1e-3 precision). \n";

checkWDist = true;
vol = gGeoManager->FindVolumeFast("volTPCPlane_U");
if( vol ) { // check U plane

	wires.clear();
	if(vol->GetNodes())for(int i=0; i<vol->GetNodes()->GetSize(); i++) {
		aNode = vol->GetNode(i);
		aNode->GetVolume()->SetLineColor(kBlue-4);
		if(aNode){
			auto matrix = aNode->GetMatrix();
			wires.emplace_back( "volTPCWire_"+TString::Itoa(i+1,10), matrix->GetTranslation()[0], matrix->GetTranslation()[1], matrix->GetTranslation()[2] );
		} else {
			std::cerr << "Wire " << i+1 << " not found.\n";
		}
	}

	float dist, c1, c2, m = tan(-30.0f*M_PI/180.0f);
	TString wireName;
	for(int i=0; i<wires.size()-1; i++) {
		c1 = wires[i+1].Y - m* wires[i+1].Z;
		c2 = wires[i].Y - m* wires[i].Z;
		dist = abs(c1-c2)/sqrt(m*m+1);

		if( abs(dist*10.0f-3.0f) >= 1.0e-3 ) {
			std::cout  << "Distance between "<< wires[i].Name << " and " << wires[i+1].Name << " is " << dist*10.0f << "mm\n";
			checkWDist &= false;
		}
		
	}
} else {
	std::cout << "No U wireplane was found.\n";
}

if(vol && checkWDist) std::cout << "U plane: All wire distances equal to 3mm (1e-3 precision). \n";

checkWDist = true;
vol = gGeoManager->FindVolumeFast("volTPCPlane_V");
if( vol ) { // check U plane
	wires.clear();
	if(vol->GetNodes())for(int i=0; i<vol->GetNodes()->GetSize(); i++) {
		aNode = vol->GetNode(i);
		aNode->GetVolume()->SetLineColor(kGreen-6);
		if(aNode){
			auto matrix = aNode->GetMatrix();
			wires.emplace_back( "volTPCWire_"+TString::Itoa(i+1,10), matrix->GetTranslation()[0], matrix->GetTranslation()[1], matrix->GetTranslation()[2] );
		} else {
			std::cerr << "Wire " << i+1 << " not found.\n";
		}
	}

	float dist, c1, c2, m = tan(+30.0f*M_PI/180.0f);
	TString wireName;
	for(int i=0; i<wires.size()-1; i++) {
		c1 = wires[i+1].Y - m* wires[i+1].Z;
		c2 = wires[i].Y - m* wires[i].Z;
		dist = abs(c1-c2)/sqrt(m*m+1);

		if( abs(dist*10.0f-3.0f) >= 1.0e-3 ) {
			std::cout  << "Distance between "<< wires[i].Name << " and " << wires[i+1].Name << " is " << dist*10.0f << "mm\n";
			checkWDist &= false;
		}
		
	}
} else {
	std::cout << "No V wireplane was found.\n";
}

if(vol && checkWDist) std::cout << "V plane: All wire distances equal to 3mm (1e-3 precision). \n";


gGeoManager->GetTopNode();
gGeoManager->CheckOverlaps(10e-11);
gGeoManager->PrintOverlaps();
gGeoManager->SetMaxVisNodes(70000);

if ( ! volName.IsNull() ) {gGeoManager->FindVolumeFast(volName)->Draw("ogl");} 
else  {gGeoManager->FindVolumeFast(volName)->Draw("ogl");}

//float mLArTPC = gGeoManager->FindVolumeFast("volTPCActive")->Weight();
//std::cout <<  gGeoManager->FindVolumeFast("volBuilding")->Weight() << "\n";

//cout << mLArTPC << endl;
//cout << mCPA << endl;

}

typedef struct _drawopt 
{
  const char* volume;
  int         color;
} drawopt;

lar1_geo(TString volName="")
{
  gSystem->Load("libGeom");
  gSystem->Load("libGdml");

  TGeoManager::Import("lar1nd_nowires.gdml");

  drawopt optuboone[] = {
	{"volHorizontalBeam",	  kGreen+2},
	{"volGround",				kOrange+4},
	//{"volTPCFrame",			  kGray+1},
    //"volCathodePlate",          kGray+1}, 
    {0, 0}
  };

  for (int i=0;; ++i) {
 		if (optuboone[i].volume==0) break;
        gGeoManager->FindVolumeFast(optuboone[i].volume)->SetLineColor(optuboone[i].color);
  }
  TList* mat = gGeoManager->GetListOfMaterials();
  TIter next(mat);
  TObject *obj;
  while (obj = next()) {
    obj->Print();
  }

  gGeoManager->GetTopNode();
  //gGeoManager->CheckOverlaps(0.0000001);
  //need to fix some overlaps, for now save some time
  gGeoManager->CheckOverlaps(10e-25);
  gGeoManager->PrintOverlaps();
  gGeoManager->SetMaxVisNodes(70000);

  //gGeoManager->GetTopVolume()->Draw();
  if ( ! volName.IsNull() ) gGeoManager->FindVolumeFast(volName)->Draw("ogl");
  //gGeoManager->FindVolumeFast("volWorld")->Draw("ogl");

  TGeoVolume *TPC = gGeoManager->FindVolumeFast("volTPC");
  float m_tpc = TPC->Weight();
//  TGeoVolume *TPC2 = gGeoManager->FindVolumeFast("volTPC2");
//  float m_tpc = TPC2->Weight();
  TGeoVolume *Cathode = gGeoManager->FindVolumeFast("volCathodePlate");
  float m_cathode = Cathode->Weight();

/*  TGeoVolume *Ground = gGeoManager->FindVolumeFast("volGroundPlate");
  float m_ground = Ground->Weight();
  TGeoVolume *UVPlane = gGeoManager->FindVolumeFast("volTPCPlane");
  float m_uvplane = UVPlane->Weight();
  TGeoVolume *YPlane = gGeoManager->FindVolumeFast("volTPCPlaneVert");
  float m_yplane = YPlane->Weight();
  TGeoVolume *FieldCageH = gGeoManager->FindVolumeFast("volFieldCageTubeTop");
  float m_fchoriz = FieldCageH->Weight();
  TGeoVolume *FieldCageV = gGeoManager->FindVolumeFast("volFieldCageTubeFront");
  float m_fcvert = FieldCageV->Weight();
*/
  float m_tpc_argon = m_tpc - m_cathode ;
//( m_cathode -m_yplane + m_ground + 2*m_uvplane + m_yplane + 50*(m_fchoriz + m_fcvert));
  cout << "LAr weight in TPC = " << m_tpc_argon << " kg\n" <<endl;

  TFile *tf = new TFile("lar1.root", "RECREATE");
 
  gGeoManager->Write();

  tf->Close();
}

{
  const TString saveDirectory = "PUT YOUR SAVE DIRECTORY HERE";
  using namespace std;
 
  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TFile *file = new TFile("/sbnd/app/users/hlay/rachel_build/work/basictree_ana_sbnd.root","READ");
  TTree *events = (TTree*) file->Get("basictreeana/EventTree");

  int N_slices, N_neutrinos;
  vector<float> *tX = 0, *tY = 0, *tZ = 0,
    *sl_VX = 0, *sl_VY = 0, *sl_VZ = 0, *xCorrection = 0;

  events->SetBranchAddress("nNeutrinos",&N_neutrinos);
  events->SetBranchAddress("trueVX",&tX);
  events->SetBranchAddress("trueVY",&tY);
  events->SetBranchAddress("trueVZ",&tZ);
  events->SetBranchAddress("xCorrection",&xCorrection);

  events->SetBranchAddress("nUsedSlices",&N_slices);
  events->SetBranchAddress("sl_VX",&sl_VX);
  events->SetBranchAddress("sl_VY",&sl_VY);
  events->SetBranchAddress("sl_VZ",&sl_VZ);
  
  int N_events = events->GetEntries();
  cout << N_events << " events..." << endl;

  TH1F *error_hist = new TH1F("error_hist",";Vertex Error (cm); Entries",50,0,10);

  float xOutEdge = 185, xInEdge = 1.5, yEdge = 185, zFrontEdge = 30, zBackEdge = 435;

  for(int event_i = 0; event_i < N_events; ++event_i){
    events->GetEntry(event_i);

    if(N_neutrinos != 1) continue;

    TVector3 true_vertex(tX->at(0),tY->at(0),tZ->at(0));
    if(TMath::Abs(tX->at(0)) > xOutEdge || TMath::Abs(tX->at(0)) < xInEdge || 
       TMath::Abs(tY->at(0)) > yEdge || tZ->at(0) < zFrontEdge || 
       tZ->at(0) > zBackEdge) continue;

    float selected_error = -std::numeric_limits<float>::max();

    for(int slice_i = 0; slice_i < N_slices; ++slice_i){
      float x = sl_VX->at(slice_i);
      if(tX->at(0) < 0) x -= xCorrection->at(0);
      else x += xCorrection->at(0);
      TVector3 reco_vertex(x,sl_VY->at(slice_i),sl_VZ->at(slice_i));

      float error = (reco_vertex - true_vertex).Mag();
      if(error < selected_error)
	selected_error = error;
    }

    if(N_slices > 0)
      error_hist->Fill(selected_error);
  }

  TCanvas *errorCanvas = new TCanvas("errorCanvas","errorCanvas");
  errorCanvas->cd();

  error_hist->GetYaxis()->SetTitleOffset(1.4);
  error_hist->Draw();
}

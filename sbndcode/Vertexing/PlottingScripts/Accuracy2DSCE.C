void Accuracy2DSCE(const TString sample, const bool cc_only, const int tpc = -1);

void Accuracy2DSCE()
{
  Accuracy2DSCE("rockbox", true);
  Accuracy2DSCE("rockbox", false);

  Accuracy2DSCE("rockbox", true, 0);
  Accuracy2DSCE("rockbox", false, 0);

  Accuracy2DSCE("rockbox", true, 1);
  Accuracy2DSCE("rockbox", false, 1);

  Accuracy2DSCE("intrnue", true);
  Accuracy2DSCE("intrnue", false);

  Accuracy2DSCE("intrnue", true, 0);
  Accuracy2DSCE("intrnue", false, 0);

  Accuracy2DSCE("intrnue", true, 1);
  Accuracy2DSCE("intrnue", false, 1);
}

void Accuracy2DSCE(const TString sample, const bool cc_only, const int tpc)
{
  const TString saveDir = "/exp/sbnd/data/users/hlay/vertexing_fix/validation/" + sample;
  gSystem->Exec("mkdir -p " + saveDir);

  const TString file = "/pnfs/sbnd/persistent/users/hlay/vertexing_fix/vertexing_fix_" + sample + ".root";

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *tree = new TChain("vertexana/SliceTree");
  tree->Add(file);

  TCanvas *canvas = new TCanvas("canvas", "canvas");
  canvas->cd();
  canvas->SetRightMargin(.25);

  TH2F *hist = new TH2F("hist", ";True x (cm);#splitline{Distance between true and}{reconstructed vertex (cm)};Events",
                        55, -220, 220, 40, 0, 10);

  TCut evSel = "is_fv";

  if(cc_only)
    {
      evSel += "cc";
      if(sample == "rockbox")
        evSel += "numu";
      else if(sample == "intrnue")
        evSel += "nue";
    }
  
  if(tpc == 0)
    evSel += "true_vtx_x<0";
  else if(tpc == 1)
    evSel += "true_vtx_x>0";
  
  tree->Draw("dr_x_corr:true_vtx_x>>hist", evSel);

  hist->GetYaxis()->SetTitleOffset(1.25);
  hist->GetXaxis()->SetTitleOffset(1.25);
  hist->GetYaxis()->SetTitleSize(.05);
  hist->GetXaxis()->SetTitleSize(.06);
  hist->Draw("colz");

  TString addedName = "";

  if(cc_only)
    {
      addedName += "_cc";
      if(sample == "rockbox")
        addedName += "_numu";
      else if(sample == "intrnue")
        addedName += "_nue";
    }

  if(tpc != -1)
    addedName += Form("_tpc_%i", tpc);

  canvas->SaveAs(saveDir + "/vertex_accuracy_true_vtx_x_" + sample + addedName + ".png");
  canvas->SaveAs(saveDir + "/vertex_accuracy_true_vtx_x_" + sample + addedName + ".pdf");
}

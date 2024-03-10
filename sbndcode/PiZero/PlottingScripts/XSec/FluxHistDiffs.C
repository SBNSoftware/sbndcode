#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"
#include "Common.C"

void FluxHistDiffs(const TString productionVersion)
{
  const TString baseDir  = baseSaveDir + "/" + productionVersion + "/flux_hists";
  const TString frontDir = baseDir;
  const TString backDir  = baseDir + "_back_face";
  const TString effZDir  = baseDir + "_eff_z";

  const TString saveDir = baseDir + "_comparisons";
  gSystem->Exec("mkdir -p " + saveDir);

  TFile *frontFile = new TFile(frontDir + "/sbnd_flux.root", "READ");
  TFile *backFile = new TFile(backDir + "/sbnd_flux.root", "READ");
  TFile *effZFile = new TFile(effZDir + "/sbnd_flux.root", "READ");

  const std::vector<TString> flavours = { "numu", "anumu", "nue", "anue" };

  for(auto const &flavour : flavours)
    {
      TCanvas *c = new TCanvas(Form("c%s", flavour.Data()), Form("c%s", flavour.Data()));
      c->cd();
      c->Draw();

      TPad *p1 = new TPad(Form("p1%s", flavour.Data()), Form("p1%s", flavour.Data()),
                          0., .3, 1., 1.);
      p1->Draw();
      p1->SetLeftMargin(.15);
      p1->SetBottomMargin(.04);
      p1->cd();

      TH1F* front = (TH1F*) frontFile->Get("flux_sbnd_" + flavour);
      TH1F* back  = (TH1F*) backFile->Get("flux_sbnd_" + flavour);
      TH1F* effZ  = (TH1F*) effZFile->Get("flux_sbnd_" + flavour);

      front->SetLineColor(kAzure-8);
      back->SetLineColor(kViolet-6);
      effZ->SetLineColor(kOrange+7);

      front->GetXaxis()->SetLabelSize(0);
      front->GetXaxis()->SetLabelOffset(999);

      front->Draw("hist][");
      back->Draw("histsame][");
      effZ->Draw("histsame][");

      TLegend *leg1 = new TLegend(.55, .65, .8, .8);
      leg1->AddEntry(front, "Front Face (z = 0cm)", "l");
      leg1->AddEntry(back, "Back Face (z = 500cm)", "l");
      leg1->AddEntry(effZ, "Effective Z (z = 227.8cm)", "l");
      leg1->Draw();

      c->cd();

      TPad *p2 = new TPad(Form("p2%s", flavour.Data()), Form("p2%s", flavour.Data()),
                          0., 0., 1., .31);
      p2->Draw();
      p2->SetLeftMargin(.15);
      p2->SetBottomMargin(.35);
      p2->SetTopMargin(.04);
      p2->cd();

      TH1F* frontBack = (TH1F*) back->Clone();
      frontBack->Divide(front);
      frontBack->SetLineColor(kGreen+3);
      frontBack->GetXaxis()->SetTitleSize(0.16);
      frontBack->GetYaxis()->SetTitleSize(0.13);
      frontBack->GetXaxis()->SetTitleOffset(0.65);
      frontBack->GetYaxis()->SetTitleOffset(0.4);
      frontBack->GetXaxis()->SetLabelSize(0.14);
      frontBack->GetYaxis()->SetLabelSize(0.12);
      frontBack->GetXaxis()->SetLabelOffset(0);
      frontBack->GetYaxis()->SetLabelOffset(0);
      frontBack->SetTitle(";E_{#nu} (GeV);Ratio");
      frontBack->SetMinimum(0.82);
      frontBack->SetMaximum(1.08);
      frontBack->Draw("hist][");

      TH1F* frontEffZ = (TH1F*) effZ->Clone();
      frontEffZ->Divide(front);
      frontEffZ->SetLineColor(kPink-6);
      frontEffZ->Draw("histsame][");

      TLegend *leg2 = new TLegend(.73, .5, .89, .8);
      leg2->AddEntry(frontBack, "Back / Front", "l");
      leg2->AddEntry(frontEffZ, "Eff Z / Front", "l");
      leg2->Draw();

      c->SaveAs(saveDir + "/flux_z_dependence_" + flavour + ".pdf");
      c->SaveAs(saveDir + "/flux_z_dependence_" + flavour + ".png");
    }
}

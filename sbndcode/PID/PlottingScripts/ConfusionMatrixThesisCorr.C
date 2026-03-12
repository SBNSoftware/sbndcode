#include "/exp/sbnd/app/users/hlay/plotting_utils/Plotting.C"
#include "/exp/sbnd/app/users/hlay/plotting_utils/HistUtils.C"

#include "RazzledHeaders.h"
#pragma link C++ class vector<vector<bool> >+;

void ConfusionMatrix(const PIDTraining &training, const bool efficiency_mode = true,
                     const bool purity_mode = false, const double comp_thresh = 0.5,
                     const double pur_thresh = 0.5, const int ke_range = -1);

void ConfusionMatrixThesisCorr()
{
  ConfusionMatrix(razzled_v12, true, false, 0.5, 0.5, 0);
  ConfusionMatrix(razzled_v12, false, true, 0.5, 0.5, 0);
  ConfusionMatrix(razzled_v12, true, true, 0.5, 0.5, 0);

  ConfusionMatrix(razzled_v12, true, false, 0.8, 0.8, 0);
  ConfusionMatrix(razzled_v12, false, true, 0.8, 0.8, 0);
  ConfusionMatrix(razzled_v12, true, true, 0.8, 0.8, 0);

  ConfusionMatrix(razzled_v12, true, false, 0.5, 0.5, 1);
  ConfusionMatrix(razzled_v12, false, true, 0.5, 0.5, 1);
  ConfusionMatrix(razzled_v12, true, true, 0.5, 0.5, 1);

  ConfusionMatrix(razzled_v12, true, false, 0.8, 0.8, 1);
  ConfusionMatrix(razzled_v12, false, true, 0.8, 0.8, 1);
  ConfusionMatrix(razzled_v12, true, true, 0.8, 0.8, 1);

  ConfusionMatrix(razzled_v12, true, false, 0.5, 0.5, 2);
  ConfusionMatrix(razzled_v12, false, true, 0.5, 0.5, 2);
  ConfusionMatrix(razzled_v12, true, true, 0.5, 0.5, 2);

  ConfusionMatrix(razzled_v12, true, false, 0.8, 0.8, 2);
  ConfusionMatrix(razzled_v12, false, true, 0.8, 0.8, 2);
  ConfusionMatrix(razzled_v12, true, true, 0.8, 0.8, 2);

  return;

  ConfusionMatrix(razzled_v12, true, false);
  ConfusionMatrix(razzled_v12, false, true);
  ConfusionMatrix(razzled_v12, true, true);

  ConfusionMatrix(razzled_v12, true, false, 0.8, 0.8);
  ConfusionMatrix(razzled_v12, false, true, 0.8, 0.8);
  ConfusionMatrix(razzled_v12, true, true, 0.8, 0.8);
}

void ConfusionMatrix(const PIDTraining &training, const bool efficiency_mode,
                     const bool purity_mode, const double comp_thresh,
                     const double pur_thresh, const int ke_range)
{
  TString save_dir = "/exp/sbnd/data/users/hlay/thesis/ncpizero/plots/NCPiZeroAv12/razzled/confusionmatrices";
  gSystem->Exec("mkdir -p " + save_dir);

  gROOT->SetStyle("henrySBND");
  gROOT->ForceStyle();

  TChain *events = new TChain("ncpizeroana/events");
  events->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_rockbox.root");
  events->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_intrnue.root");
  events->Add("/pnfs/sbnd/persistent/users/hlay/ncpizero/NCPiZeroAv12/NCPiZeroAv12_intime.root");

  std::vector<std::vector<int>> *truePDGV = 0, *recoPDGV = 0;
  std::vector<std::vector<float>> *compV = 0, *purityV = 0;
  std::vector<std::vector<double>> *slc_pfp_track_start_x = 0, *slc_pfp_track_start_y = 0, *slc_pfp_track_start_z = 0,
    *slc_pfp_shower_start_x = 0, *slc_pfp_shower_start_y = 0, *slc_pfp_shower_start_z = 0, *slc_pfp_track_length = 0,
    *slc_pfp_shower_energy = 0, *slc_pfp_true_ke = 0;
  std::vector<bool> *unambiguousSlice = 0;
  std::vector<std::vector<bool>> *recoPrimary = 0;

  events->SetBranchAddress("slc_pfp_razzled_pdg", &recoPDGV);
  events->SetBranchAddress("slc_pfp_true_pdg", &truePDGV);
  events->SetBranchAddress("slc_pfp_comp", &compV);
  events->SetBranchAddress("slc_pfp_pur", &purityV);
  events->SetBranchAddress("slc_pfp_track_start_x", &slc_pfp_track_start_x);
  events->SetBranchAddress("slc_pfp_track_start_y", &slc_pfp_track_start_y);
  events->SetBranchAddress("slc_pfp_track_start_z", &slc_pfp_track_start_z);
  events->SetBranchAddress("slc_pfp_shower_start_x", &slc_pfp_shower_start_x);
  events->SetBranchAddress("slc_pfp_shower_start_y", &slc_pfp_shower_start_y);
  events->SetBranchAddress("slc_pfp_shower_start_z", &slc_pfp_shower_start_z);
  events->SetBranchAddress("slc_pfp_track_length", &slc_pfp_track_length);
  events->SetBranchAddress("slc_pfp_shower_energy", &slc_pfp_shower_energy);
  events->SetBranchAddress("slc_is_clear_cosmic", &unambiguousSlice);
  events->SetBranchAddress("slc_pfp_primary_child", &recoPrimary);
  events->SetBranchAddress("slc_pfp_true_ke", &slc_pfp_true_ke);

  const unsigned N = events->GetEntries();

  unsigned n_categories_true = 5, n_categories_reco = 5;

  TH2F *hConfusionMatrix = new TH2F("hConfusionMatrix",
                                    ";True Class;Razzled Class;",
                                    n_categories_true, 0, n_categories_true,
                                    n_categories_reco, 0, n_categories_reco);

  for(unsigned i = 0; i < N; ++i)
    {
      events->GetEntry(i);

      for(int slc_i = 0; slc_i < unambiguousSlice->size(); ++slc_i)
	{
	  if(unambiguousSlice->at(slc_i))
	    continue;

	  for(int pfp_i = 0; pfp_i < recoPrimary->at(slc_i).size(); ++pfp_i)
	    {
	      if(!(recoPrimary->at(slc_i).at(pfp_i)))
		continue;

	      int truePDG = truePDGV->at(slc_i).at(pfp_i);
	      int recoPDG = recoPDGV->at(slc_i).at(pfp_i);

	      if(abs(truePDG) != 11 && abs(truePDG) != 13 && abs(truePDG) != 22 
		 && abs(truePDG) != 211 && abs(truePDG) != 2212)
		continue;

	      if(abs(recoPDG) != 11 && abs(recoPDG) != 13 && abs(recoPDG) != 22 
		 && abs(recoPDG) != 211 && abs(recoPDG) != 2212)
		continue;

	      float comp = compV->at(slc_i).at(pfp_i);
	      float purity = purityV->at(slc_i).at(pfp_i);

	      if(comp < comp_thresh || purity < pur_thresh)
		continue;

	      double trackStartX = slc_pfp_track_start_x->at(slc_i).at(pfp_i);
	      double trackStartY = slc_pfp_track_start_y->at(slc_i).at(pfp_i);
	      double trackStartZ = slc_pfp_track_start_z->at(slc_i).at(pfp_i);

	      double showerStartX = slc_pfp_shower_start_x->at(slc_i).at(pfp_i);
	      double showerStartY = slc_pfp_shower_start_y->at(slc_i).at(pfp_i);
	      double showerStartZ = slc_pfp_shower_start_z->at(slc_i).at(pfp_i);

	      if(abs(trackStartX) > 180 || abs(trackStartY) > 180 || trackStartZ < 10 || trackStartZ > 450 ||
		 abs(showerStartX) > 180 || abs(showerStartY) > 180 || showerStartZ < 10 || showerStartZ > 450)
		continue;

	      double trk_length = slc_pfp_track_length->at(slc_i).at(pfp_i);
	      double showerEnergy = slc_pfp_shower_energy->at(slc_i).at(pfp_i);

	      if(trk_length < 3 && showerEnergy<10)
		continue;

	      if(ke_range > -1)
		{
		  double ke = slc_pfp_true_ke->at(slc_i).at(pfp_i);

		  if(ke_range == 0 && ke > 0.1)
		    continue;

		  if(ke_range == 1 && (ke <= 0.1 || ke > 1))
		    continue;

		  if(ke_range == 2 && ke <= 1)
		    continue;
		}

	      int trueClass = razzledMap.at(abs(truePDG));
	      int recoClass = razzledMap.at(abs(recoPDG));

	      hConfusionMatrix->Fill(trueClass, recoClass);
	    }
	}
    }
  
  TCanvas *c = new TCanvas("c", "c");
  c->cd();
  
  if(efficiency_mode)
    NormaliseEntriesByXTotal(hConfusionMatrix);
  if(purity_mode)
    NormaliseEntriesByYTotal(hConfusionMatrix);

  std::vector<TString> axisLabels = razzledAxisLabels;

  for(int i = 1; i <= hConfusionMatrix->GetNbinsX(); ++i)
    hConfusionMatrix->GetXaxis()->SetBinLabel(i, axisLabels[i]);
  for(int i = 1; i <= hConfusionMatrix->GetNbinsY(); ++i)
    hConfusionMatrix->GetYaxis()->SetBinLabel(i, axisLabels[i]);

  gStyle->SetPaintTextFormat("1.2g");
  hConfusionMatrix->SetMarkerSize(3);
  hConfusionMatrix->Draw("col text");

  TString file_name = "razzled";

  if(efficiency_mode)
    file_name += "_efficiency";
  if(purity_mode)
    file_name += "_purity";
  if(comp_thresh > 0.5 && pur_thresh > 0.5)
    file_name += "_high_quality";
  if(comp_thresh < 0.5 && pur_thresh < 0.5)
    file_name += "_low_quality";
  if(ke_range == 0)
    file_name += "_low_ke";
  if(ke_range == 1)
    file_name += "_mod_ke";
  if(ke_range == 2)
    file_name += "_high_ke";

  c->SaveAs(save_dir + "/" + file_name + ".pdf");
  c->SaveAs(save_dir + "/" + file_name + ".png");
}

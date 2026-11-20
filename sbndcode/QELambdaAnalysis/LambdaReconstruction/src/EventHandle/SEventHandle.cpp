////////////////////////////////////////////////////////////////////////////
//
// \file TPCLineParameters.h
//
// \brief Definition of TPCLinesParameters
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "SEventHandle.h"

std::vector<TPad*> buildpadcanvas(int nx, int ny){
    std::vector<TPad*> Tp;
    double x=0, y=1, dx=1./nx, dy=1./ny;
    TPad *pad = new TPad("padFrame","padFrame", 0, 0, 1, 1, -1, -1, -1);
    Tp.push_back(pad);
    for(int i=1; i<=nx; i++){
      y=1;
      for(int j=1; j<=ny; j++){
        TPad *pad = new TPad( ("pad"+std::to_string(i)+"_"+std::to_string(j)).c_str(), ("pad"+std::to_string(i)+"_"+std::to_string(j)).c_str(), x, y-dy, x+dx, y, -1, -1, -1);
        Tp.push_back(pad);
        y-=dy;
      }
      x+=dx;
    }
    for(int i=0; i<=nx*ny; i++){
      Tp.at(i)->Draw();
      Tp.at(i)->SetBottomMargin(0.15);
    }
    return Tp;
}


std::vector<TString> GetInputFileList(std::string file_name,  std::string ext, std::string dirPath){
    // Get the candidate Files
    std::vector<TString> fFilePaths;
    std::cout<<" Looking in directory: "<<dirPath<<std::endl;
    TSystemDirectory dir(dirPath.c_str(), dirPath.c_str());
    TList *files = dir.GetListOfFiles();
    TString targetFileName(file_name);
    TString targetExtension(ext);
    std::cout<<" Target file name "<<targetFileName<<std::endl;
    gSystem->Exec("ls");
    if (files){
        TSystemFile *file;
        TString fname;
        TIter next(files);
        while ((file=(TSystemFile*)next())) {
            fname = file->GetName();
            
            if (!file->IsDirectory() && fname.EndsWith(targetExtension)){
                if(fname.Contains(targetFileName)){
                if(dirPath!=".")
                    fFilePaths.push_back(dirPath+fname);
                else
                    fFilePaths.push_back(fname);
                }
            }
        }
    }

    // sort in alphabetical order
    std::sort(fFilePaths.begin(), fFilePaths.end());

    // cout
    std::cout<<" Found "<<fFilePaths.size()<<" files"<<std::endl;
    for(auto &f:fFilePaths){
        std::cout<<"  "<<f<<std::endl;
    }

    return fFilePaths;
}


EfficiencyCalculator::EfficiencyCalculator():
      fEventList({}),
      nEvents(0),
      nEventsSkipped(0),
      nProcessedEvents(0),
      nEventsSelected(0)
      {
          //LABELS SIZE AND FONT
          gStyle->SetLabelFont(132, "XYZ");
          gStyle->SetTitleSize(0.06, "TXYZ");
          gStyle->SetLabelSize(0.06, "XYZ");

          hNOrigins = new TH1F("hNOrigins", "NOrigins;# origins; # entries", 10, 0, 10);
          hNOriginsMult1 = new TH1F("hNOriginsMult1", "NOrigins mult=1;# origins; # entries", 10, 0, 10);
          hNOriginsMult2 = new TH1F("hNOriginsMult2", "NOrigins mult=2;# origins; # entries", 10, 0, 10);
          hNOriginsMultGt3 = new TH1F("hNOriginsMultGt3", "NOrigins mult>=3;# origins; # entries", 10, 0, 10);
          hProf2DMatrix = new TH2F("hProf2DMatrix", "# origins; origin 1 multiplicity; origin 2 multiplicity", 10, 0, 10, 10, 0, 10);

          hHitDensity = new TH1F("hHitDensity", ";HitDensity [#hits/wire]; # entries", 30, 0, 10);
      }


void EfficiencyCalculator::UpdateHistograms(SEvent recoEv){
    hNOrigins->Fill(recoEv.GetNOrigins());
    hNOriginsMult1->Fill(recoEv.GetNOriginsMult(1));
    hNOriginsMult2->Fill(recoEv.GetNOriginsMult(2));
    hNOriginsMultGt3->Fill(recoEv.GetNOriginsMultGt(3));
    
    /*for(int m1=0; m1<=8; m1++){
        for(int m2=0; m2<=8; m2++){
            if(recoEv.GetNOriginsMult(m1)>0 && recoEv.GetNOriginsMult(m2)>0){
                hProf2DMatrix->Fill(m1, m2, 1);
            }
            
        }
    }*/

    // use only the two main origins
    std::vector<SOrigin> origins = recoEv.GetOrigins();

    if(origins.size()==1){
        hProf2DMatrix->Fill(origins[0].Multiplicity(), 0 );
    }
    else{
        int maxMult = std::max(origins[0].Multiplicity(), origins[1].Multiplicity());
        int minMult = std::min(origins[0].Multiplicity(), origins[1].Multiplicity());
        hProf2DMatrix->Fill(maxMult, minMult);
    }
    
    hHitDensity->Fill(recoEv.HitDensity());
}

void EfficiencyCalculator::UpdateValues(){
    nEvents = fEventList.size();
    nEventsSkipped=0;
    nProcessedEvents=0;
    nEventsSelected=0;
    for(auto &pair:fEventList){
        if(pair.second.NSkipped()==pair.second.NSlices()){
            nEventsSkipped++;
        }
        else{
            nProcessedEvents++;
            if(pair.second.NSelected()>0){
                nEventsSelected++;
            }
        }
    }
}

void EfficiencyCalculator::DrawHistograms(TCanvas *c){
    
    gStyle->SetOptStat(1);
    
    c->cd();

    std::vector<TPad*> padV = buildpadcanvas(3, 2);
    double fLeftMargin=0.15;

    padV[1]->cd();
    padV[1]->SetLeftMargin(fLeftMargin);
    hNOrigins->Draw();
    
    padV[2]->cd();
    padV[2]->SetLeftMargin(fLeftMargin);
    hNOriginsMult1->Draw();

    padV[3]->cd();
    padV[3]->SetLeftMargin(fLeftMargin);
    hNOriginsMult2->Draw();

    padV[4]->cd();
    padV[4]->SetLeftMargin(fLeftMargin);
    hNOriginsMultGt3->Draw();

    padV[5]->cd();
    padV[5]->SetLeftMargin(fLeftMargin);
    hHitDensity->Draw();

    padV[6]->cd();
    padV[6]->SetLeftMargin(fLeftMargin);
    hProf2DMatrix->SetStats(0);

    padV[6]->SetRightMargin(0.12);
    // Normalize the histogram to the number of entries
    if (hProf2DMatrix->GetEntries() > 0) {
        hProf2DMatrix->Scale(1.0 / hProf2DMatrix->GetEntries());
    }

    hProf2DMatrix->Draw("colz");

    int maxMult=4;
    // Add percentage values to each bin as labels
    for (Int_t i = 1; i <= maxMult; ++i) {
        for (Int_t j = 1; j <maxMult; ++j) {
            Double_t binContent = hProf2DMatrix->GetBinContent(i, j);
            // two to replace the symmetric matrix, double counting
            Double_t percentage =  binContent * 100.0;

            // Create a TText object for the label
            TText* label = new TText(hProf2DMatrix->GetXaxis()->GetBinCenter(i), hProf2DMatrix->GetYaxis()->GetBinCenter(j), Form("%.2f%%", percentage));
            label->SetTextSize(0.05); // Adjust the text size as needed
            label->SetTextAlign(22);  // Centered alignment
            label->SetTextColor(kRed);
            label->Draw();
        }
    }

    hProf2DMatrix->GetXaxis()->SetRangeUser(0, 4);
    hProf2DMatrix->GetYaxis()->SetRangeUser(0, 4);
    hProf2DMatrix->GetXaxis()->SetNdivisions(4);
    hProf2DMatrix->GetYaxis()->SetNdivisions(4);

    c->Update();
    c->cd();
    c->WaitPrimitive();
}
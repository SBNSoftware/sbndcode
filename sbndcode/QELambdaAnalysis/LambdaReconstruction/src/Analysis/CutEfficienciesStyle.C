
#ifndef CUTEFFICIENCIES_STYLE_H
#define CUTEFFICIENCIES_STYLE_H

#include "TString.h"

// CutStyler class
class CutStyler{
    public:
        CutStyler(bool batchMode):
            fBatchMode(batchMode),
            fColors({kBlack, kRed-3,  kAzure-5, kGreen+3, kOrange-3, kMagenta+1, kCyan- 3, kYellow+2, kViolet-1, kTeal-1,}),
            //fColors({kRed+1, kBlue-4, kGreen+4, kOrange+5, kMagenta+2, kCyan+3, kYellow+2, kViolet-1, kTeal-1, kAzure-5}),
            fLineStyle({0, 1, 2, 9, 3, 6, 10, 7, 4, 5, 8}),
            fMarkerStyle({0, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30})
        {
            SetGStyle();
        }
        void SetGStyle();
        size_t GetColor(size_t i){return fColors[i];}
        size_t GetLineStyle(size_t i){return fLineStyle[i];}
        size_t GetMarkerStyle(size_t i){return fMarkerStyle[i];}

    private:
        bool fBatchMode;
        std::vector<int> fColors;
        std::vector<int> fLineStyle;
        std::vector<int> fMarkerStyle;
};


void CutStyler::SetGStyle(){
    
  gStyle->SetTitleFont(62, "TXYZ");
  gStyle->SetTitleSize(0.05,"TXYZ");

  // LABELS SIZE AND FONT
  gStyle->SetLabelFont(62, "TXYZ");
  gStyle->SetLabelSize(0.05,"TXYZ");

  // AXIS OFFSETS AND SIZES
  gStyle->SetTitleXOffset(0.80);
  gStyle->SetTitleXSize(0.06);

  gStyle->SetTitleYOffset(1.2);
  gStyle->SetTitleYSize(0.06);

  gStyle->SetMarkerStyle(33);
  gStyle->SetMarkerSize(1.5);
  gStyle->SetLineColor(46);
  gStyle->SetLineWidth(1);

  gStyle->SetHistLineWidth(2);
  gStyle->SetMarkerSize(1);

  // Draw horizontal and vertical grids
  //gStyle->SetPadGridX(kTRUE);                     
  //gStyle->SetPadGridY(kTRUE);

  // set left margin to 10%
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.1);

  // N divisions
  gStyle->SetNdivisions(505);

  //gStyle->SetOptLogy(1);
}

vector<TPad*> buildpadcanvas(int nx, int ny){
    vector<TPad*> Tp;
    double x=0, y=1, dx=1./nx, dy=1./ny;
    TPad *pad = new TPad("","", 0, 0, 1, 1, -1, -1, -1);
    Tp.push_back(pad);
    for(int i=1; i<=nx; i++){
      y=1;
      for(int j=1; j<=ny; j++){
        //TPad *pad = new TPad("a", "a",x, y, x+dx,y+dy);
        //TPad pad("a", "a",x, y, x+dx,y+dy,1, 1, 2);
        //cout<<x<<" "<<y<<endl;
        TPad *pad = new TPad("","", x, y-dy, x+dx, y, -1, -1, -1);
        Tp.push_back(pad);
        y-=dy;
      }
      x+=dx;
    }
    for(int i=0; i<=nx*ny; i++){
      Tp.at(i)->Draw();
    }
    return Tp;
}

#endif
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
        // Bottom margin
        pad->SetBottomMargin(0.15);
        // Left margin
        pad->SetLeftMargin(0.15);
      }
      x+=dx;
      
    }
    for(int i=0; i<=nx*ny; i++){
      Tp.at(i)->Draw();
    }
    return Tp;
}

void SetFRANSStyle(){
  // Axis offsets
  gStyle->SetTitleOffset(1.3, "Y");
  // Set divisions
  gStyle->SetNdivisions(505, "X");
  // Line widths
  gStyle->SetLineWidth(2);
  // Graphs line widths
  gStyle->SetHistLineWidth(5);
}
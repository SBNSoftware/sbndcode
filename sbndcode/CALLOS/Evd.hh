// Some functions needed to debug the code
// #include "TTimer.h"
#include "AverageWaveform.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLine.h"
#include "nuevdb/EventDisplayBase/Canvas.h"
#include <iostream>

// #include "lareventdisplay/EventDisplay/DrawingPad.h"

void lets_pause()
{
  // TTimer *timer = new TTimer("gSystem->ProcessEvents();", 50, kFALSE);
  // timer->TurnOn();
  // timer->Reset();
  std::cout << "q/Q to quit, other to continuee: ";
  char kkey;
  std::cin.get(kkey);
  if (kkey == 'q' || kkey == 'Q')
    throw std::exception(); // std::exit(0); //gSystem->Exit(0); //
  // timer->TurnOff();
  // delete timer;
}

void event_display(const std::vector<float>& inputArray, int starting_tick){
    // Get the number of samples in the waveform
    int size = inputArray.size();
    
    auto c1 = new TCanvas("c1","c1",600,400);

    // draw the waveform using a TH1F
    auto h = new TH1F("h1f","Test random numbers",size,0,size);
    for (int i = 0; i < size; i++) {
        h->SetBinContent(i+1, inputArray[i]);
    };
    h->Draw();

    // draw a vertical line at the starting tick
    TLine *v_line= new TLine(starting_tick, h->GetMinimum(), starting_tick, h->GetMaximum());
    v_line->SetLineColor(kRed);
    v_line->SetLineWidth(2);
    v_line->SetLineStyle(kDashed);
    //draw the line

    v_line->Draw("same");
    c1->Update();
    c1->Draw();
    c1->SaveAs("my_evd.C");
    
    
    c1->Close();
    delete c1;
    delete h;
    delete v_line;

    bool display=false;

    if (display) system("root -l my_evd.C"); // open file with root, closes after execution of line arguments
    
    lets_pause();
}
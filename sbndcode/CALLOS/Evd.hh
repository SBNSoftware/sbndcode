// Some functions needed to debug the code
// #include "TTimer.h"
#include "AverageWaveform.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLine.h"
#include "nuevdb/EventDisplayBase/Canvas.h"
#include <iostream>
#include <vector>
#include <TMarker.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include "TLatex.h"  

// #include "lareventdisplay/EventDisplay/DrawingPad.h"

void lets_pause()
{
  // TTimer *timer = new TTimer("gSystem->ProcessEvents();", 50, kFALSE);
  // timer->TurnOn();
  // timer->Reset();
  std::cout << "q/Q to quit, other to continue: ";
  char kkey;
  std::cin.get(kkey);
  if (kkey == 'q' || kkey == 'Q')
    throw std::exception(); // std::exit(0); //gSystem->Exit(0); //
  // timer->TurnOff();
  // delete timer;
}

void event_display(const std::vector<float>& inputArray, int starting_tick, int ending_tick, double baseline, double baseline_std, double peak_height, double peak_pos, int wfChannel, int FirstBin, int LastBin) {
    int size = inputArray.size();
    //baseline = 0;
    std::string title = "Waveform for Channel " + std::to_string(wfChannel);
    
    auto c1 = new TCanvas("c1",title.c_str(),600,400);
    
    // draw the waveform using a TH1F
    auto h = new TH1F("h1f",title.c_str(),size,0,size);
    for (int i = 0; i < size; i++) {
        h->SetBinContent(i+1, inputArray[i]);
    };
    
    h->GetXaxis()->SetTitle("Ticks [16 ns/tick]");
    h->GetYaxis()->SetTitle("ADC counts");
    h->Draw();

    // // Draw the smoothed waveform using a TH1F
    // auto h_smoothed = new TH1F("h_smoothed", "Smoothed Waveform", size, 0, size);
    // for (int i = 0; i < size; i++) {
    //     h_smoothed->SetBinContent(i + 1+starting_tick, Smoothed[i]);
    // }
    // h_smoothed->SetLineColor(kRed);
    // h_smoothed->Draw("same");

    // // Draw the derivative waveform using a TH1F
    // auto h_derivative = new TH1F("h_derivative", "Derivative Waveform", size, 0, size);
    // for (int i = 0; i < size; i++) {
    //     h_derivative->SetBinContent(i + 1 +starting_tick, Derivate[i]);
    // }
    // h_derivative->SetLineColor(kGreen);
    // h_derivative->Draw("same");
    // draw a vertical line at the starting tick
    TLine *v_line= new TLine(starting_tick, h->GetMinimum(), starting_tick, h->GetMaximum());
    v_line->SetLineColor(kRed);
    v_line->SetLineWidth(2);
    v_line->SetLineStyle(kDashed);
    //draw the line

    v_line->Draw("same");

    // draw a vertical line at the ending tick
    TLine *v_line2= new TLine(ending_tick, h->GetMinimum(), ending_tick, h->GetMaximum());
    v_line2->SetLineColor(kRed);
    v_line2->SetLineWidth(2);
    v_line2->SetLineStyle(kDashed);
    //draw the line

    v_line2->Draw("same");

    TLine *v_line3= new TLine(FirstBin, h->GetMinimum(), FirstBin, h->GetMaximum());
    v_line3->SetLineColor(kBlack);
    v_line3->SetLineWidth(2);
    v_line3->SetLineStyle(kDashed);

    v_line3->Draw("same");

    TLine *v_line4= new TLine(LastBin, h->GetMinimum(), LastBin, h->GetMaximum());
    v_line4->SetLineColor(kBlue);
    v_line4->SetLineWidth(2);
    v_line4->SetLineStyle(kDashed);

    v_line4->Draw("same");

    // draw a horizontal line at the baseline
    TLine *h_line= new TLine(0, baseline, size, baseline);
    h_line->SetLineColor(kMagenta);
    h_line->SetLineWidth(2);
    h_line->SetLineStyle(kDashed);
    //draw the line

    h_line->Draw("same");

    // draw a horizontal line at the baseline + baseline_std
    TLine *h_line2= new TLine(0, baseline + baseline_std, size, baseline + baseline_std);
    h_line2->SetLineColor(kGreen);
    h_line2->SetLineWidth(2);
    h_line2->SetLineStyle(kDashed);
    //draw the line

    h_line2->Draw("same");

    // draw a horizontal line at the baseline - baseline_std
    TLine *h_line3= new TLine(0, baseline - baseline_std, size, baseline - baseline_std);
    h_line3->SetLineColor(kGreen);
    h_line3->SetLineWidth(2);
    h_line3->SetLineStyle(kDashed);
    //draw the line

    h_line3->Draw("same");

    // draw a horizontal star at the peak
    TLine *h_line4= new TLine(0, peak_height, size, peak_height);
    h_line4->SetLineColor(kBlue);
    h_line4->SetLineWidth(2);
    h_line4->SetLineStyle(kDashed);
    //draw the line

    h_line4->Draw("same");

    // Add a dot at coordinates (x, y)
    TMarker *marker = new TMarker(peak_pos, peak_height, kFullCircle);
    marker->SetMarkerColor(kRed);
    marker->SetMarkerSize(1);
    marker->Draw("same");


        // Crear la leyenda
    auto legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Coordenadas (x1, y1, x2, y2) en el canvas
    // legend->SetHeader("Legend", "C"); // Título opcional de la leyenda
    legend->SetBorderSize(1); // Tamaño del borde
    legend->SetFillColor(0); // Fondo transparente

    // Añadir elementos a la leyenda
    // legend->AddEntry(h, "Original Waveform", "l"); // "l" indica línea
    legend->AddEntry(v_line, "ROI", "l");
    legend->AddEntry(v_line2, "End Tick", "l");
    legend->AddEntry(h_line, "Baseline", "l");
    legend->AddEntry(h_line2, "Baseline STD", "l");
    legend->AddEntry(h_line4, "Signal amplitude", "l");
    legend->AddEntry(marker, "Alignment point", "p"); // "p" indica punto
    legend->AddEntry(v_line3, "First Bin", "l");
    legend->AddEntry(v_line4, "Last Bin", "l");
    // Dibujar la leyenda
    legend->Draw();

   
    c1->Update();
    c1->Draw();
    c1->SaveAs("my_evd.C");


    c1->Close();
    delete c1;
    delete h;
    delete v_line;
    delete v_line2;
    delete v_line3;
    delete v_line4;
    delete h_line;
    delete h_line2;
    delete h_line3;
    delete h_line4;
    delete marker;
    delete legend;

    bool display=false;

    if (display) system("root -l my_evd.C"); // open file with root, closes after execution of line arguments

    lets_pause();
}

void event_display_fft(const std::vector<float>& inputArray, int wfChannel) {
    // Get the number of samples in the waveform
    int size = inputArray.size();
    std::string title = "Waveform for Channel " + std::to_string(wfChannel);
    auto c1 = new TCanvas("c1","c1",600,400);

    // draw the waveform using a TH1F
    auto h = new TH1F("h1f",title.c_str(),size,0,size);
    for (int i = 0; i < size; i++) {
        h->SetBinContent(i+1, inputArray[i]);
    };
    h->GetXaxis()->SetTitle("Ticks [16 ns/tick]");
    h->GetYaxis()->SetTitle("ADC counts");
    h->Draw();


    c1->Update();
    c1->Draw();
    c1->SaveAs("my_evd_fft.C");


    c1->Close();
    delete c1;
    delete h;

    lets_pause();
}

 

void event_display_paper(const std::vector<float>& inputArray, int starting_tick, int ending_tick, double baseline, double baseline_std) {
    int size = inputArray.size();
    //baseline = 0;

    // Conversion factor: 16 ns/tick -> microseconds
    const double tick_to_us = 16.0 / 1000.0;
    double x_max = size * tick_to_us;
    double start_us = starting_tick * tick_to_us;
    double end_us   = ending_tick * tick_to_us;

    

    auto c1 = new TCanvas("c1", "c1", 600, 400);

    // draw the waveform using a TH1F, X axis in microseconds, no title
    auto h = new TH1F("h1f", "", size, 0, x_max);
    for (int i = 0; i < size; i++) {
        h->SetBinContent(i + 1, inputArray[i]);
    }
    h->SetTitle("");
    h->SetStats(0);   // also remove the stats box
    h->GetXaxis()->SetTitle("Time [#mus]");
    h->GetYaxis()->SetTitle("ADC counts");
    h->Draw();

    // vertical line at the starting tick (in us)
    TLine *v_line = new TLine(start_us, h->GetMinimum(), start_us, h->GetMaximum());
    v_line->SetLineColor(kRed);
    v_line->SetLineWidth(2);
    v_line->SetLineStyle(kDashed);
    v_line->Draw("same");

    // vertical line at the ending tick (in us)
    TLine *v_line2 = new TLine(end_us, h->GetMinimum(), end_us, h->GetMaximum());
    v_line2->SetLineColor(kRed);
    v_line2->SetLineWidth(2);
    v_line2->SetLineStyle(kDashed);
    v_line2->Draw("same");

    // horizontal line at the baseline
    TLine *h_line = new TLine(0, baseline, x_max, baseline);
    h_line->SetLineColor(kMagenta+2);
    h_line->SetLineWidth(2);
    h_line->SetLineStyle(kDashed);
    h_line->Draw("same");

    // horizontal line at baseline + baseline_std
    TLine *h_line2 = new TLine(0, baseline + baseline_std, x_max, baseline + baseline_std);
    h_line2->SetLineColor(kGreen+2);
    h_line2->SetLineWidth(2);
    h_line2->SetLineStyle(kDashed);
    h_line2->Draw("same");

    // horizontal line at baseline - baseline_std
    TLine *h_line3 = new TLine(0, baseline - baseline_std, x_max, baseline - baseline_std);
    h_line3->SetLineColor(kGreen+2);
    h_line3->SetLineWidth(2);
    h_line3->SetLineStyle(kDashed);
    h_line3->Draw("same");

    // Crear la leyenda
    auto legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    // Pass a null object so the entry does not try to inherit attributes
    // from the TLine (that inheritance is what was making every entry
    // show up red); set color/style/width on the entry directly instead.
    auto* entry_roi      = legend->AddEntry((TObject*)nullptr, "ROI", "l");
    entry_roi->SetLineColor(kRed);
    entry_roi->SetLineStyle(kDashed);
    entry_roi->SetLineWidth(2);
    auto* entry_baseline = legend->AddEntry((TObject*)nullptr, "Baseline", "l");
    entry_baseline->SetLineColor(kMagenta+2);
    entry_baseline->SetLineStyle(kDashed);
    entry_baseline->SetLineWidth(2);
    auto* entry_std      = legend->AddEntry((TObject*)nullptr, "Baseline STD", "l");
    entry_std->SetLineColor(kGreen+2);
    entry_std->SetLineStyle(kDashed);
    entry_std->SetLineWidth(2);
    legend->Draw();

    // "SBND Data" label (top-right corner, NDC coordinates)
    TLatex *sbnd_label = new TLatex();
    sbnd_label->SetNDC();
    sbnd_label->SetTextFont(62);
    sbnd_label->SetTextSize(0.045);
    sbnd_label->SetTextAlign(33);
    sbnd_label->DrawLatex(0.90, 0.96, "SBND Data");

    c1->Update();
    c1->Draw();
    c1->SaveAs("my_evd.C");
    c1->Close();
    delete c1;
    delete h;
    delete v_line;
    delete v_line2;
    delete h_line;
    delete h_line2;
    delete h_line3;
    delete legend;
    delete sbnd_label;
    bool display = false;
    if (display) system("root -l my_evd.C");
    lets_pause();
}
////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesAlgo.h
//
// \brief Definition of TPCLinesDisplay
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCLinesDisplay.h"


TPCLinesDisplay::TPCLinesDisplay(bool show, std::string outputPath):
            fOutputPath(outputPath)
{
    fColors = {kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5, kOrange+7, kGreen+2, kRed, kViolet+2, kOrange+3, kPink+9, kGray+2, kBlue+10, kBlue-4, kYellow+4, kRed-2, kGreen-5};

    fColorsOrigins = {40, 42, 46, 30, 35, 40, 42, 46, 30, 35, 40, 42, 46, 30, 35, 40, 42, 46, 30, 35, 40, 42, 46, 30, 35, 40, 42, 46, 30, 35};

    fHitsMarkers = {20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};

    std::iota(std::begin(fAlbet), std::end(fAlbet), 'a');

    if(show==false){
        gROOT->SetBatch(true);
    }
    
}

void TPCLinesDisplay::SetStyle(){
    //TITLES SIZE AND FONT
    gStyle->SetPalette(112,0);
    gStyle->SetTitleFont(132, "TXYZ");
    gStyle->SetTitleSize(0.05, "TXYZ");

    gStyle->SetTitleFont(132, "titleFont"); 
    
    // Off stats
    gStyle->SetOptStat(0);     // Turn off statistics box

    
    //LABELS SIZE AND FONT
    gStyle->SetLabelFont(132, "XYZ");
    gStyle->SetLabelSize(0.05, "XYZ");

    gStyle->SetTitleYOffset (1.4);

}

TH2F* TPCLinesDisplay::GetFrame(std::vector<SHit> hitsV, std::string label){

    std::vector<double> x, y, err;
    for(size_t ix=0; ix<hitsV.size(); ix++){
        x.push_back(hitsV[ix].X());
        y.push_back(hitsV[ix].Y());
        err.push_back(hitsV[ix].Width());
    }
    double minX = *std::min_element(x.begin(), x.end());
    double maxX = *std::max_element(x.begin(), x.end());
    double minY = *std::min_element(y.begin(), y.end());
    double maxY = *std::max_element(y.begin(), y.end());
    double overX = 0.1*(maxX-minX);
    double overY = 0.1*(maxY-minY);

    TH2F* hFrame = new TH2F("frame", (label+";WireId;TimeTick [0.5 #mu s]").c_str(), 200, minX-overX, maxX+overX, 200, minY-overY, maxY+overY);
    hFrame->SetStats(0);
    return hFrame;
}

void TPCLinesDisplay::DrawLine(LineEquation line, double xmin, double xmax, TLegend * leg, std::string label, int color, int style){
    // Draw a horizontal line from x = 1 to x = 4 at y = 2
    double y1 = line.EvaluateX(xmin);
    double y2 = line.EvaluateX(xmax);
    TLine* horizontalLine = new TLine(xmin, y1, xmax, y2);
    horizontalLine->SetLineColor(color); // Set line color
    horizontalLine->SetLineStyle(style); // Set line color
    horizontalLine->SetLineWidth(2);     // Set line width
    horizontalLine->Draw("same");

    if(label!="")
        leg->AddEntry(horizontalLine, label.c_str(), "l");
}

void TPCLinesDisplay::DrawHitScatter(std::vector<SHit> hitsV, TLegend * leg, std::string label, int color, int style, double size, double errorAlpha, bool markerByCluster, bool addClusterIDLabels){
    
    if(hitsV.size()>0){
        // map with vector of hits by ClusterID
        std::map<int, std::vector<SHit>> hitMap;
        for(size_t ix=0; ix<hitsV.size(); ix++){
            if( hitMap.find(hitsV[ix].ClusterId())==hitMap.end())
                hitMap[hitsV[ix].ClusterId()] = {hitsV[ix]};
            else
                hitMap[hitsV[ix].ClusterId()].push_back(hitsV[ix]);
        }

        int clusterCounter = 0;
        for(auto & cluster:hitMap){
            std::vector<SHit> hits = cluster.second;
            std::vector<double> x, y, err;
            for(size_t ix=0; ix<hits.size(); ix++){
                x.push_back(hits[ix].X());
                y.push_back(hits[ix].Y());
                err.push_back(hits[ix].Width());
            }

            TGraphErrors *g = new TGraphErrors(x.size(),&x[0],&y[0], nullptr, &err[0]); 
            g->SetMarkerColorAlpha(color, 0.5);
            if(markerByCluster) g->SetMarkerStyle(fHitsMarkers[cluster.first]);
            else g->SetMarkerStyle(style);
            g->SetMarkerSize(size);
            g->SetLineColorAlpha(kGray, errorAlpha);
            g->Draw("p");
            
            if(!addClusterIDLabels){
                if(label!="" && clusterCounter==0)
                    leg->AddEntry(g, label.c_str(), "p");
            }
            else{
                leg->AddEntry(g, ("PFP "+std::to_string(cluster.first)).c_str(), "p");
            }

            clusterCounter++;
        }
    }

    return;
}

void TPCLinesDisplay::DrawVertex(SVertex vertex, TLegend * leg, std::string label, int color, int marker, double alpha){
    TGraph *pointGraph = new TGraph();
    pointGraph->SetPoint(0, vertex.Point().X(), vertex.Point().Y()); // Set the point's coordinates
    
    // Set marker style and size for the point
    pointGraph->SetMarkerStyle(marker);
    pointGraph->SetMarkerSize(2.);
    //pointGraph->SetMarkerColorAlpha(color, alpha);
    pointGraph->SetMarkerColor(color);

    if(label!="" && leg!=nullptr)
        leg->AddEntry(pointGraph, label.c_str(), "p");

    // Draw the point
    pointGraph->Draw("P same");
}

void TPCLinesDisplay::DrawTriangle(STriangle tri, TLegend * leg, std::string label, int colorP, int color, double alpha){ 
    // Define the triangle's vertices
    Double_t x[3] = {tri.GetMainVertex().X(), tri.GetVertexB().X(), tri.GetVertexC().X()};
    Double_t y[3] = {tri.GetMainVertex().Y(), tri.GetVertexB().Y(), tri.GetVertexC().Y()};

    // Create a polyline representing the triangle
    TPolyLine *triangle = new TPolyLine(3, x, y, "F");
    triangle->SetFillColorAlpha(color, alpha);  // Set fill color

    TGraph *pointGraph = new TGraph();
    pointGraph->SetPoint(0, tri.GetMainVertex().X(), tri.GetMainVertex().Y()); // Set the point's coordinates
    // Set marker style and size for the point
    pointGraph->SetMarkerStyle(29);
    pointGraph->SetMarkerSize(2.2);
    pointGraph->SetMarkerColor(colorP);

    if(label!="")
        leg->AddEntry(pointGraph, label.c_str(), "p");

    // Draw the point
    pointGraph->Draw("P same");

    // Draw the triangle
    triangle->Draw("F same");
}

void TPCLinesDisplay::DrawLinearCluster(SLinearCluster cluster, int cID, TLegend * leg, std::string label, int color, double size, int style){
    if(cluster.NHits()>0){
        std::vector<SHit> hits = cluster.GetHits();
        std::vector<double> x, y;
        for(size_t ix=0; ix<hits.size(); ix++){
            x.push_back(hits[ix].X());
            y.push_back(hits[ix].Y());
        }

        TGraph *g = new TGraph(x.size(),&x[0],&y[0]); 

        //g->SetMarkerColorAlpha(color, 0.5);
        g->SetMarkerColor(color);
        g->SetMarkerStyle(style);;
        g->SetMarkerSize(size);
        g->SetLineWidth(20);
        //g->SetLineColorAlpha(kGray, errorAlpha);
        g->Draw("p same");



        DrawLine(cluster.GetTrackEquation(), cluster.GetMinX(), cluster.GetMaxX(), leg, "", color, kSolid);
        DrawLine(cluster.GetTrackEquationStart(), cluster.GetMinX(), cluster.GetMaxX(), leg, "", color, kDashed);
        DrawLine(cluster.GetTrackEquationEnd(), cluster.GetMinX(), cluster.GetMaxX(), leg, "", color, kDashed);

        int preClusterCounter=0;
        if(cluster.GetId()!=-1)
            leg->AddEntry(g, ( "Cluster "+std::to_string(cID)).c_str(), "p");
        else{
            leg->AddEntry(g, ("Pre-cluster "+std::to_string(cID)).c_str(), "p");
            preClusterCounter++;
        }
    }

    
    return;
}



void TPCLinesDisplay::Show(
    bool wait,
    std::string eventLabel,
    std::vector<SHit> allHitsV,
    LineEquation houghLine,
    std::vector<SHit> selectedHitsV,
    std::vector<SLinearCluster> clustersV,
    SLinearCluster mainDirection,
    std::vector<STriangle> originAngles,
    SVertex recoVertex,
    SVertex trueVertex,
    std::vector<SOrigin> origins,
    TCanvas *canvas)
{

    SetStyle();

    if(allHitsV.size()==0){return;}
    
    if(canvas==nullptr)
        canvas = new TCanvas("c", "c",  0, 0, 1000, 800);

    canvas->cd();
    
    TPad *pad1 = new TPad("pad1", "Graph Pad", 0.0, 0., .8, 1.0);
    TPad *pad2 = new TPad("pad2", "Legend Pad", 0.75, 0., 1.0, 1.);
    pad1->SetLeftMargin(0.15);
    pad2->SetLeftMargin(-0.2);
    pad1->Draw();
    pad2->Draw();

    TLegend * legend = new TLegend(0., 0.1, 0.99, 0.9);

    pad1->cd();

    // general frame
    TH2F * hFrame=GetFrame(allHitsV, eventLabel);
    hFrame->Draw();

    // Main direction
    if(mainDirection.NHits()>0){
        // selected hits scatter
        std::vector<SHit> mainDirHits = mainDirection.GetHits();
        DrawHitScatter(mainDirHits, legend, "MainLine", kBlack, 8, 1., 0);   
    }

    // all hits scatter
    bool addClusterIDLabels = clustersV.size()>0? false:true;
    DrawHitScatter(allHitsV, legend, "AllHits", 65, 8, 1.5, 0.6, true, addClusterIDLabels);
    
    // selected hits
    if(selectedHitsV.size()!=0){
        // selected hits scatter
        DrawHitScatter(selectedHitsV, legend, "Selected hits", kRed, 24, 1.1, 0);

        // hough line
        DrawLine(houghLine, hFrame->GetXaxis()->GetXmin(), hFrame->GetXaxis()->GetXmax(), legend, "Track line", kBlack, kDashed);
    }
   
    // linear clusters
    if(clustersV.size()>0){

        for(size_t cIx=0; cIx<clustersV.size(); cIx++){
            int cIDLabel = ( clustersV[cIx].GetId()!=-1)?  clustersV[cIx].GetId():cIx;
            DrawLinearCluster(clustersV[cIx], cIDLabel, legend, "Cluster", fColors[cIx], .75, 8);
        }
    }

    // triangles
    if(originAngles.size()>0){
        for(size_t oIx=0; oIx<originAngles.size(); oIx++){
            std::string VLabel = "V_{"+std::to_string(oIx)+"} ";
            VLabel = VLabel + originAngles[oIx].GetTrack1().GetId()+"-"+originAngles[oIx].GetTrack2().GetId()+"#rightarrow";
            VLabel = VLabel + originAngles[oIx].GetMainTrack().GetId();

            DrawTriangle(originAngles[oIx], legend, VLabel, fColorsOrigins[oIx], 90, 0.5);
        }
    }


    if(origins.size()>0){
        for(size_t oIx=0; oIx<origins.size(); oIx++){
            // get letter abs for the origin

            std::string originLabel = "#upsilon_{"+std::string(1, fAlbet[oIx])+"}^{"+std::to_string(origins[oIx].Multiplicity())+"} (";
            bool first = true;
             int sign = (origins[oIx].IsEdgeOrigin())? 1:-1;
            for(SLinearCluster & trk:origins[oIx].GetTracks()){  
                if(!first) originLabel += ", ";
                originLabel += std::to_string(sign*trk.GetId());
                first = false;
            }
            originLabel += ")";
            if(sign==1)
                DrawVertex(origins[oIx].GetPoint(), legend, originLabel.c_str(), fColorsOrigins[oIx], 90, 0.5);
            else
                DrawVertex(origins[oIx].GetPoint(), nullptr, originLabel.c_str(), fColorsOrigins[oIx], 90, 0.5);
        }
    }


    // draw PANDORA vertex
    if(recoVertex.Point().X()!=-1 && recoVertex.Point().Y()!=-1){
        DrawVertex(recoVertex, legend, "PANDORA vtx", kPink-9, 33, 1);
    }

    // draw PANDORA vertex
    if(trueVertex.Point().X()!=-1 && trueVertex.Point().Y()!=-1){
        DrawVertex(trueVertex, legend, "True vtx", kPink-12, 110, 1);
    }
    
    // Create a legend
    pad2->cd();
    legend->SetBorderSize(0);
    legend->SetTextSize(fLegendFontSize); 
    legend->Draw("same");

    canvas->cd();
    canvas->Update();
    if(wait)
        canvas->WaitPrimitive();


    return;
}

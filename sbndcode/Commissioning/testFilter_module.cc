// where you stopped: 

// Framework includes 
#include "art/Framework/Core/EDFilter.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes 
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larcoreobj/SummaryData/POTSummary.h"

// ROOT includes 
#include "TRandom3.h"

// C++ includes
#include <map>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>
#include <bitset>

const int DEFAULT_VALUE = -9999;
#define PI 3.14159265

using std::vector;

class testFilter : public art::EDFilter {
public:
   explicit testFilter(fhicl::ParameterSet const& p);
   virtual bool filter(art::Event& e) override;
   void    reconfigure(fhicl::ParameterSet const& p);
   virtual ~testFilter() { }
   
private:
   // declare functions 
   /// Resets variables that are saved to the tree 
   void ResetWireHitsVars(int n);
   /// function that calculates distance 
   float Distance(int x1, int y1, int x2, int y2);
   /// function that performs the Hough transform 
   void Hough(vector<vector<int>> coords, int threshold, int max_gap, int range, float min_length, int muon_length, int nentry, vector<vector<int>>& lines);

   // Wire hits variables
   int                 nhits;   

   art::ServiceHandle<art::TFileService> tfs;
   
   std::string fHitsModuleLabel; 
   int max_hits;

   int fHoughThreshold;
   int fHoughMaxGap;
   int fHoughRange;
   int fHoughMinLength;
   int fHoughMuonLength;
};  

testFilter::testFilter(fhicl::ParameterSet const& p): EDFilter{p} 
{
   this->reconfigure(p);
}

void testFilter::reconfigure(fhicl::ParameterSet const& p)
{
   fHitsModuleLabel     = p.get<std::string>("HitsModuleLabel");
   max_hits             = p.get<int>("MaxHits", 50000);

   // Hough Transform parameters 
   fHoughThreshold      = p.get<int>("HoughThreshold",10);
   fHoughMaxGap         = p.get<int>("HoughMaxGap",30);
   fHoughRange          = p.get<int>("HoughRange",100);
   fHoughMinLength      = p.get<int>("HoughMinLength",500);
   fHoughMuonLength     = p.get<int>("HoughMuonLength",2500);
}

bool testFilter::filter(art::Event& evt)
{ 
   bool pass = false;
   int event = evt.id().event();

   // get nhits 
   art::Handle<vector<recob::Hit>> hitListHandle;
   vector<art::Ptr<recob::Hit>> hitlist;
   if (evt.getByLabel(fHitsModuleLabel,hitListHandle)) {
      art::fill_ptr_vector(hitlist, hitListHandle);
      nhits = hitlist.size();
   }
   else {
      std::cout << "Failed to get recob::Hit data product." << std::endl;
      nhits = 0;
    }
   if (nhits > max_hits) {
      std::cout << "Available hits are " << nhits 
                << ", which is above the maximum number allowed to store." << std::endl;
      std::cout << "Will only store " << max_hits << "hits." << std::endl;
      nhits = max_hits;
    }

   //ResetWireHitsVars(nhits);
   vector<vector<int>> collection_hits0, collection_hits1;
   collection_hits0.reserve(3000); // set reserved length somewhere else in parameters? 
   collection_hits1.reserve(3000); 

   for (int i = 0; i < nhits; ++i) {
      geo::WireID wireid = hitlist[i]->WireID();
      int hit_wire = int(wireid.Wire), hit_peakT = int(hitlist[i]->PeakTime()), hit_plane = wireid.Plane, hit_tpc = wireid.TPC;
      if (hit_plane==2 && hit_peakT>0){ // if collection plane and only positive peakT 
         vector<int> v{hit_wire,hit_peakT};
         if (hit_tpc==0)
            collection_hits0.push_back(v);  
         else
            collection_hits1.push_back(v); 
      }
   } // end of nhit loop
   collection_hits0.shrink_to_fit(); collection_hits1.shrink_to_fit(); 
   vector<vector<int>> lines0, lines1; 
   lines0.reserve(10); lines1.reserve(10); 
    
   //transform to find lines 
   Hough(collection_hits0,fHoughThreshold,fHoughMaxGap,fHoughRange,fHoughMinLength,fHoughMuonLength,event,lines0);
   Hough(collection_hits1,fHoughThreshold,fHoughMaxGap,fHoughRange,fHoughMinLength,fHoughMuonLength,event,lines1);

   lines0.shrink_to_fit(); lines1.shrink_to_fit();

   if (lines0.empty() == false || lines1.empty() == true){
      // if an AC muon is found, perform hough transform on hits in plane==0
      pass = true;
      if (lines0.empty == false){ // if a line was found in TPC0

      }
   }
   return pass; 
} // end of event loop 

float testFilter::Distance(int x1, int y1, int x2, int y2){
   return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) * 1.0);
}

void testFilter::Hough(vector<vector<int>> coords, int threshold, int max_gap, int range, float min_length, int muon_length, int nentry, vector<vector<int>>& lines){
   //set global variables 
   TRandom3 rndgen;
   const int h = 3500; const int w = 2000; //range of hit_wire
   constexpr int accu_h = h + w + 1 ; const int accu_w = 180; 
   const int x_c = (w/2); const int y_c = (h/2); 

   //create accumulator/pointer to accumulator 
   int accu[accu_h][accu_w] = {{0}};
   int (*adata)[accu_w];
   adata = &accu[(accu_h-1)/2]; // have pointer point to the middle of the accumulator array (then we can use negative indices)

   //declare necessary vectors 
   vector<vector<int>> data = coords; // points will not be removed 
   vector<vector<int>> deaccu; //deaccumulator
   vector<vector<int>> outlines; 

   // loop over points and perform transform 
   int count = coords.size(); 
   for ( ; count>0; count--){ 
      int idx = rndgen.Uniform(count);
      int max_val = threshold-1;
      if ((coords.at(idx)).empty())
         continue; 
      int x = coords[idx][0], y = coords[idx][1], rho = 0, theta = 0;
      vector<int> v{x,y}; 
      deaccu.push_back(v);
      //loop over all angles and fill the accumulator 
      for (int j=0; j<accu_w; j++){ 
         int r = int(round((x-x_c)*cos(j*PI/accu_w) + (y-y_c)*sin(j*PI/accu_w)));
         int val = ++(adata[r][j]);       
         if (max_val < val){
            max_val = val;
            rho = r;
            theta = j*180/accu_w;
         }
      }
      if (max_val < threshold){
         (coords.at(idx)).clear(); 
         continue;
      }
      //start at point and walk the corridor on both sides 
      vector<vector<int>> endpoint(2, vector<int>(2));
      for (int k=0; k<2;k++){ 
         int i=0, gap=0;
         while (gap < max_gap){ 
            (k==0)? i++ : i--; 
            if ( (idx+i) == int(data.size()) || (idx+i) <0) // if we reach the edges of the data set 
               break;
            if ((data.at(idx+i)).empty()) // if the point has already been removed 
               continue;
            int x1 = data[idx+i][0], y1 = int(data[idx+i][1]); 
            int last_x, diffx; 
            if (endpoint[k][0]!= 0){ // ensure we don't jump large x-values 
               last_x = endpoint[k][0];
               diffx = abs(last_x - x1);
               if (diffx > 30){
                  break;
               }
            }
            int y_val = int(round((rho - (x1 - x_c)*cos(theta*PI/180.0))/sin(theta*PI/180.0) + y_c));
            if (y1 >= (y_val-range) && y1 <= (y_val+range)){
               gap = 0;
               endpoint[k] = {x1, y1};
               (coords.at(idx+i)).clear();
               (data.at(idx+i)).clear();
            }
            else
               gap++;
         } // end of while loop 
      } // end of k loop 

      // unvote from the accumulator 
      for (int n = (deaccu.size()-1); n>=0; n--){ 
         int x1 = deaccu[n][0], y1 = int(deaccu[n][1]);
         int y_val = int(round((rho - (x1 - x_c)*cos(theta*PI/180.0))/sin(theta*PI/180.0) + y_c));
         if (y1 >= (y_val-range) && y1 <= (y_val+range)){
            for (int m=0; m<accu_w; m++){
               int r = int(round((x1-x_c)*cos(m*PI/accu_w) + (y1-y_c)*sin(m*PI/accu_w)));
               (adata[r][m])--;
            }
            deaccu.erase(deaccu.begin() + n);
         }
      } // end of deaccumulator loop

      int x0_end = endpoint[0][0], y0_end = endpoint[0][1], x1_end = endpoint[1][0], y1_end = endpoint[1][1];
      if ((x0_end==0 && y0_end==0) || (x1_end==0 && y1_end==0)) // don't add the (0,0) points 
         continue;
      vector<int> line = {x0_end, y0_end, x1_end, y1_end, rho, theta};
      outlines.push_back(line);

   } // end of point loop 
   // combine lines that are split 
   for (int i=0; i<int(outlines.size()); i++){
      bool same = false;
      for (int j=i+1; j<int(outlines.size()) && same == false;j++){ 
         int xi_coords[2] = {outlines[i][0], outlines[i][2]}; int xj_coords[2] = {outlines[j][0], outlines[j][2]};
         int yi_coords[2] = {outlines[i][1], outlines[i][3]}; int yj_coords[2] = {outlines[j][1], outlines[j][3]};
         int rhoi = outlines[i][4], rhoj = outlines[j][4];
         int thetai = outlines[i][5], thetaj = outlines[j][5]; 

         int var = 100;
         int rho_var = 30;
         int theta_var = 20; 
         for (int k=0; k<2 && same == false; k++){
            for (int l=0; l<2 && same == false; l++){
               int counter = 0; 
               if ((xi_coords[k] < (xj_coords[l] + var)) && (xi_coords[k] > (xj_coords[l] - var)))
                  counter++;
               if ((yi_coords[k] < (yj_coords[l] + var)) && (yi_coords[k] > (yj_coords[l] - var)))
                  counter++ ;
               if ((rhoi < (rhoj + rho_var)) && (rhoi > (rhoj - rho_var)))
                  counter++; 
               if ((thetai < (thetaj + theta_var)) && (thetai > (thetaj - theta_var)))
                  counter++;
               if (counter >= 3){ // if at least three of the conditions are fulfilled 
                  if(k==0){
                     if(l==0){
                        outlines[j][2] = outlines[i][0];
                        outlines[j][3] = outlines[i][1];
                     }
                     else{
                        outlines[j][0] = outlines[i][0];
                        outlines[j][1] = outlines[i][1];
                     }
                  }
                  else{
                     if(l==0){
                        outlines[j][2] = outlines[i][2]; 
                        outlines[j][3] = outlines[i][3];                        
                     }
                     else{
                        outlines[j][0] = outlines[i][2];
                        outlines[j][1] = outlines[i][3]; 
                     }  
                  }
                  same = true;
                  (outlines.at(i)).clear();
                  //std::fill ((outlines.at(i)).begin(),(outlines.at(i)).end(),0);
               } 
            }
         }
      } // end of j loop 
   } // end of i loop 

   for (int i=0; i<int(outlines.size()); i++){
      if ((outlines.at(i)).empty())
         continue;
      int x0_end = outlines[i][0], y0_end = outlines[i][1], x1_end = outlines[i][2], y1_end = outlines[i][3];
      if (muon_length!=0){
         int y_diff = abs(y0_end-y1_end);
         if (y_diff > muon_length){
            lines.push_back(outlines.at(i));
         }
      }
      else{
         float length = Distance(x0_end,y0_end,x1_end,y1_end);
         if (length > min_length){
            lines.push_back(outlines.at(i));
         }
      }
   }
   //free memory 
   data.clear(); deaccu.clear(); outlines.clear();
} // end of hough 
  // A macro required for a JobControl module.
  DEFINE_ART_MODULE(testFilter)
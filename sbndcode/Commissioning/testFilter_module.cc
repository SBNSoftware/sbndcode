// where you stopped: 
// added wire information into the lines 

// for each AC muon, want to store:
// anode endpoint, cathode endpoint, t0, thetaxz, thetayz 

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
   void ResetCollectionHitVectors();
   void ResetInductionHitVectors();
   void FindWireIntersection(vector<vector<int>> lines_col, vector<vector<int>> lines_ind, vector<art::Ptr<recob::Hit>> hitlist); 
   float Distance(int x1, int y1, int x2, int y2);
   void Hough(vector<vector<int>> coords, int threshold, int max_gap, int range, float min_length, int muon_length, int nentry, vector<vector<int>>& lines);

   // Wire hits variables
   int   nhits;
   vector<vector<int>> hit_00; // tpc0, plane0
   vector<vector<int>> hit_01; // tpc0, plane1
   vector<vector<int>> hit_02; // tpc0, plane2

   vector<vector<int>> hit_10; // tpc1, plane0
   vector<vector<int>> hit_11; // tpc1, plane1
   vector<vector<int>> hit_12; // tpc1, plane2 

   vector<vector<int>> lines_00; 
   vector<vector<int>> lines_01;
   vector<vector<int>> lines_02; 

   vector<vector<int>> lines_10;
   vector<vector<int>> lines_11;
   vector<vector<int>> lines_12; 

   art::ServiceHandle<art::TFileService> tfs;
   
   std::string fHitsModuleLabel; 
   int max_hits;

   int fHoughThreshold;
   int fHoughMaxGap;
   int fHoughRange;
   int fHoughMinLength;
   int fHoughMuonLength;

   geo::GeometryCore const* fGeometryService;
};

testFilter::testFilter(fhicl::ParameterSet const& p): EDFilter{p} 
{
   fGeometryService = lar::providerFrom<geo::Geometry>();

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

   ResetCollectionHitVectors();

   for (int i = 0; i < nhits; ++i) {
      geo::WireID wireid = hitlist[i]->WireID();
      int hit_wire = int(wireid.Wire), hit_peakT = int(hitlist[i]->PeakTime()), hit_plane = wireid.Plane, hit_tpc = wireid.TPC;
      if (hit_plane==2 && hit_peakT>0){ // if collection plane and only positive peakT 
         vector<int> v{hit_wire,hit_peakT,i};
         if (hit_tpc==0)
            hit_02.push_back(v);  
         else
            hit_12.push_back(v); 
      }
   } // end of nhit loop
   hit_02.shrink_to_fit(); hit_12.shrink_to_fit(); 
   Hough(hit_02,fHoughThreshold,fHoughMaxGap,fHoughRange,fHoughMinLength,fHoughMuonLength,event,lines_02);
   Hough(hit_12,fHoughThreshold,fHoughMaxGap,fHoughRange,fHoughMinLength,fHoughMuonLength,event,lines_12);

   //flags 
   bool ac_tpc0 = !(lines_02.empty()); // will be true if an ac muon was detected in tpc0 
   bool ac_tpc1 = !(lines_12.empty()); // will be true if an ac muon was detected in tpc1

   ResetInductionHitVectors();
   //find induction plane hits
   if (ac_tpc0 == true || ac_tpc1 == true){
      std::cout << "ac muon found in tpc0: " << ac_tpc0 << std::endl;
      std::cout << "ac muon found in tpc1: " << ac_tpc1 << std::endl;
      pass = true;
      for (int i = 0; i < nhits; ++i) {
         geo::WireID wireid = hitlist[i]->WireID();
         int hit_wire = int(wireid.Wire), hit_peakT = int(hitlist[i]->PeakTime()), hit_tpc = wireid.TPC, hit_plane = wireid.Plane;
         vector<int> v{hit_wire,hit_peakT,i};
         if (ac_tpc0 == true){ //if ac muon was found in tpc0 
            if (hit_plane==0 && hit_tpc==0 && hit_peakT>0) 
               hit_00.push_back(v);
            if (hit_plane==1 && hit_tpc==0 && hit_peakT>0)
               hit_01.push_back(v);
         }
         if (ac_tpc1 == true){ // if ac muon was found in tpc 1
            if (hit_plane != 2){
               std::cout << hit_plane << std::endl;
            }
            if (hit_plane==0){ //&& hit_tpc==1 && hit_peakT>0){
               std::cout << "here(1)!" << std::endl;
               //hit_10.push_back(v);
            }
            if (hit_plane==1){ //&& hit_tpc==1 && hit_peakT>0){
               std::cout << "here(2)!" << std::endl;
               //hit_11.push_back(v);
            }
         } 
      }
      std::cout << "sizes of tpc1 hit coords: " << hit_10.size() << ", " << hit_11.size() << std::endl;
      if (ac_tpc0){
         Hough(hit_00,fHoughThreshold,fHoughMaxGap,fHoughRange,fHoughMinLength,fHoughMuonLength,event,lines_00);
         if (lines_00.empty() == false){
            FindWireIntersection(lines_02,lines_00,hitlist);
         }
         else{
            Hough(hit_01,fHoughThreshold,fHoughMaxGap,fHoughRange,fHoughMinLength,fHoughMuonLength,event,lines_01);
            if (lines_01.empty()==false){
               FindWireIntersection(lines_02,lines_01,hitlist);
            }
            else{
               std::cout << "tpc0: induction lines not found" << std::endl;
            }
         }
      }
      if (ac_tpc1){
         Hough(hit_10,fHoughThreshold,fHoughMaxGap,fHoughRange,fHoughMinLength,fHoughMuonLength,event,lines_10);
         if (lines_10.empty() == false){
            std::cout << "tpc1: line found on plane0 "  << std::endl;
            FindWireIntersection(lines_12,lines_10,hitlist);
         }
         else{
            Hough(hit_11,fHoughThreshold,fHoughMaxGap,fHoughRange,fHoughMinLength,fHoughMuonLength,event,lines_11);
            if (lines_11.empty() == false){
               std::cout << "tpc1: line found on plane 0 "  << std::endl;
               FindWireIntersection(lines_12,lines_11,hitlist);
            }
            else{
               std::cout << "tpc1: induction lines not found" << std::endl;
            }
         }
      }
   }
   
   // structure of vector in lines = (x0_end, y0_end, x1_end, y1_end, wire0_end, wire1_end, rho, theta);
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
      vector<vector<int>> endpoint(2, vector<int>(3));
      for (int k=0; k<2;k++){ 
         int i=0, gap=0;
         while (gap < max_gap){ 
            (k==0)? i++ : i--; 
            if ( (idx+i) == int(data.size()) || (idx+i) <0) // if we reach the edges of the data set 
               break;
            if ((data.at(idx+i)).empty()) // if the point has already been removed 
               continue;
            int x1 = data[idx+i][0], y1 = int(data[idx+i][1]), wire_idx = data[idx+i][2]; 
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
               endpoint[k] = {x1, y1, wire_idx};
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
      int wire0_end = endpoint[0][2], wire1_end = endpoint[1][2]; 
      if ((x0_end==0 && y0_end==0) || (x1_end==0 && y1_end==0)) // don't add the (0,0) points 
         continue;
      vector<int> outline = {x0_end, y0_end, x1_end, y1_end, wire0_end, wire1_end, rho, theta};
      outlines.push_back(outline);

   } // end of point loop 
   // combine lines that are split 
   for (int i=0; i<int(outlines.size()); i++){
      bool same = false;
      for (int j=i+1; j<int(outlines.size()) && same == false;j++){ 
         int xi_coords[2] = {outlines[i][0], outlines[i][2]}; int xj_coords[2] = {outlines[j][0], outlines[j][2]};
         int yi_coords[2] = {outlines[i][1], outlines[i][3]}; int yj_coords[2] = {outlines[j][1], outlines[j][3]};
         int rhoi = outlines[i][6], rhoj = outlines[j][6];
         int thetai = outlines[i][7], thetaj = outlines[j][7]; 

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

void testFilter::ResetCollectionHitVectors() {

   hit_02.clear(); 
   hit_12.clear(); 
   lines_02.clear(); 
   lines_12.clear(); 
  
   hit_02.reserve(3000); 
   hit_12.reserve(3000); 
   lines_02.reserve(10); 
   lines_12.reserve(10); 
}

void testFilter::ResetInductionHitVectors(){
   hit_00.clear();
   hit_01.clear();
   hit_10.clear();
   hit_11.clear();

   lines_00.clear(); 
   lines_01.clear();
   lines_10.clear(); 
   lines_11.clear();

   hit_00.reserve(5000);
   hit_01.reserve(5000);
   hit_10.reserve(5000);
   hit_11.reserve(5000);

   lines_00.reserve(1); 
   lines_01.reserve(1);
   lines_10.reserve(1); 
   lines_11.reserve(1);
}

void testFilter::FindWireIntersection( vector<vector<int>> lines_col, vector<vector<int>> lines_ind, vector<art::Ptr<recob::Hit>> hitlist){
   for (int i=0; i<int(lines_col.size()); i++){
      for (int j=0; j<int(lines_ind.size()); j++){
         int peakT0_col, peakT1_col, peakT0_ind, peakT1_ind; 
         int wire0_col, wire1_col, wire0_ind, wire1_ind;
         if (lines_col[i][1] < lines_col[i][3]){
            peakT0_col = lines_col[i][1];
            peakT1_col = lines_col[i][3];
            wire0_col = lines_col[i][4]; 
            wire1_col = lines_col[i][5]; 
         }
         else{
            peakT0_col = lines_col[i][3];
            peakT1_col = lines_col[i][1];
            wire0_col = lines_col[i][5]; 
            wire1_col = lines_col[i][4]; 
         }

         if (lines_ind[j][1] < lines_ind[j][3]){
            peakT0_ind = lines_ind[j][1];
            peakT1_ind = lines_ind[j][3];
            wire0_ind = lines_ind[j][4]; 
            wire1_ind = lines_ind[j][5]; 
         }
         else{
            peakT0_ind = lines_ind[j][3];
            peakT1_ind = lines_ind[j][1];
            wire0_ind = lines_ind[j][5]; 
            wire1_ind = lines_ind[j][4]; 
         }

         int peakT_range = 10; 
         if ( (peakT0_col + peakT_range > peakT0_ind) && (peakT0_col - peakT_range < peakT0_ind) && 
              (peakT1_col + peakT_range > peakT1_ind) && (peakT1_col + peakT_range > peakT1_ind)){
               
            std::cout << "wire intersection found!" << std::endl;
            geo::WireID awire_col = hitlist[wire0_col]->WireID(); 
            geo::WireID awire_ind = hitlist[wire0_ind]->WireID();
            geo::WireID cwire_col = hitlist[wire1_col]->WireID();
            geo::WireID cwire_ind = hitlist[wire1_ind]->WireID();
            geo::Point_t apoint, cpoint; 

            if (fGeometryService->WireIDsIntersect(awire_col, awire_ind, apoint) == true){
               std::cout << "intersection coord of awire:" << apoint.X() << ", " << apoint.Y() << ", " << apoint.Z() << std::endl;
            }
            else{
               std::cout << "intersection of awire not found by WireIDsIntersect" << std::endl;
            }

            if (fGeometryService->WireIDsIntersect(cwire_col, cwire_ind, cpoint) == true){
               std::cout << "intersection coord of cwire:" << cpoint.X() << ", " << cpoint.Y() << ", " << cpoint.Z() << std::endl;
            }
            else{
               std::cout << "intersection of cwire not found by WireIDsIntersect"<< std::endl;
            }
         }
      }
   }
}
// A macro required for a JobControl module.
DEFINE_ART_MODULE(testFilter)
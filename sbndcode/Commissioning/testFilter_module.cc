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
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"

// SBN/SBND includes
#include "sbnobj/SBND/CRT/CRTData.hh"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/CRTHitRecoAlg.h"

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
   // Define functions 
   // Resets hit information for collection plane
   void ResetCollectionHitVectors();
   // Resets hit information for induction planes
   void ResetInductionHitVectors();
   // Resets variables for AC Crossing muons 
   void ResetACVariables(); 
   // Finds distance between two points 
   float Distance(int x1, int y1, int x2, int y2);
   // Performs Hough Transform for collection planes 
   void Hough_col(vector<vector<int>> coords, int threshold, int max_gap, int range, float min_length, int muon_length, int nentry, 
                  vector<vector<int>>& lines, vector<vector<int>>& colwire, vector<vector<int>>& colpeakT);
   // Performs Hough Transform for induction planes 
   void Hough_ind(vector<vector<int>> coords, int threshold, int max_gap, int range, float min_length, int muon_length, int nentry, 
                  vector<vector<int>>& lines);
   // Finds collection plane hit info for AC muon, stores them in vectors 
   // Finds t0, stores them in a vector<vector<double>> 
   void Findt0(vector<vector<int>> lines, vector<double>& ac_t0); 
   // Finds endpoints, stores them in a vector<vector<geo::Point_t>>    
   void FindEndpoints(vector<vector<int>> lines_col, vector<vector<int>> lines_ind, int range, vector<art::Ptr<recob::Hit>> hitlist, vector<vector<geo::Point_t>>& ac_endpoints); 
   // Finds trajectories, stores them in a vector<vector<double>> 
   void FindTrajectories(vector<vector<geo::Point_t>> ac_endpoints,vector<vector<double>>& ac_trajectories); 

   // Define variables 
   // Wire hits variables
   int   nhits;
   // will store hit_wire and hit_peakT for hit_<tpc><plane> 
   vector<vector<int>> hit_00;
   vector<vector<int>> hit_01;
   vector<vector<int>> hit_02;

   vector<vector<int>> hit_10;
   vector<vector<int>> hit_11;
   vector<vector<int>> hit_12; 

   // will store hough transform output for lines_<tpc><plane>
   vector<vector<int>> lines_00; 
   vector<vector<int>> lines_01;
   vector<vector<int>> lines_02; 
   vector<vector<int>> colwire_02; 
   vector<vector<int>> colpeakT_02;

   vector<vector<int>> lines_10;
   vector<vector<int>> lines_11;
   vector<vector<int>> lines_12; 
   vector<vector<int>> colwire_12; 
   vector<vector<int>> colpeakT_12;

   // AC crossing muon variables 
   vector<int> ac_tpc; 
   vector<double> ac_t0;
   vector<vector<geo::Point_t>> ac_endpoints;
   vector<vector<double>> ac_trajectories; 

   // AC variables stored in CRTTrack 
   vector<sbn::crt::CRTTrack> ac_tracks;

   // parameters from the fcl file 
   std::string fHitsModuleLabel; 
   int fMax_Hits;
   int fHoughThreshold;
   int fHoughMaxGap;
   int fHoughRange;
   int fHoughMinLength;
   int fHoughMuonLength;
   int fEndpointRange; 

   // services 
   art::ServiceHandle<art::TFileService> tfs;
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
   fMax_Hits             = p.get<int>("MaxHits", 50000);

   // Hough parameters 
   fHoughThreshold      = p.get<int>("HoughThreshold",10);
   fHoughMaxGap         = p.get<int>("HoughMaxGap",30);
   fHoughRange          = p.get<int>("HoughRange",100);
   fHoughMinLength      = p.get<int>("HoughMinLength",500);
   fHoughMuonLength     = p.get<int>("HoughMuonLength",2500);

   // FindEndpoints parameters 
   fEndpointRange       = p.get<int>("EndpointRange",30);
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
   if (nhits > fMax_Hits){
      std::cout << "Available hits are " << nhits 
                << ", which is above the maximum number allowed to store." << std::endl;
      std::cout << "Will only store " << fMax_Hits << "hits." << std::endl;
      nhits = fMax_Hits;
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
   Hough_col(hit_02,fHoughThreshold,fHoughMaxGap,fHoughRange,fHoughMinLength,fHoughMuonLength,event,lines_02, colwire_02,colpeakT_02);
   Hough_col(hit_12,fHoughThreshold,fHoughMaxGap,fHoughRange,fHoughMinLength,fHoughMuonLength,event,lines_12, colwire_12,colpeakT_12);

   //flags 
   bool ac_in_tpc0 = !(lines_02.empty()); // will be true if an ac muon was detected in tpc0 
   bool ac_in_tpc1 = !(lines_12.empty()); // will be true if an ac muon was detected in tpc1

   ResetInductionHitVectors();
   ResetACVariables();
   //find induction plane hits
   if (ac_in_tpc0 == true || ac_in_tpc1 == true){
      std::cout << "ac muon found in tpc0: " << ac_in_tpc0 << std::endl;
      std::cout << "ac muon found in tpc1: " << ac_in_tpc1 << std::endl;
      pass = true;
      for (int i = 0; i < nhits; ++i) {
         geo::WireID wireid = hitlist[i]->WireID();
         int hit_wire = int(wireid.Wire), hit_peakT = int(hitlist[i]->PeakTime()), hit_tpc = wireid.TPC, hit_plane = wireid.Plane;
         vector<int> v{hit_wire,hit_peakT,i};
         if (ac_in_tpc0 == true){ //if ac muon was found in tpc0 
            if (hit_plane==0 && hit_tpc==0 && hit_peakT>0) 
               hit_00.push_back(v);
            if (hit_plane==1 && hit_tpc==0 && hit_peakT>0)
               hit_01.push_back(v);
         }
         if (ac_in_tpc1 == true){ // if ac muon was found in tpc 1
            if (hit_plane==0){ //&& hit_tpc==1 && hit_peakT>0){
               hit_10.push_back(v);
            }
            if (hit_plane==1){ //&& hit_tpc==1 && hit_peakT>0){
               hit_11.push_back(v);
            }
         } 
      }
      if (ac_in_tpc0){
         ac_tpc.insert(ac_tpc.end(),0,lines_02.size());
         Findt0(lines_02, ac_t0);
         Hough_ind(hit_00,fHoughThreshold,fHoughMaxGap,fHoughRange,fHoughMinLength,fHoughMuonLength,event,lines_00);
         if (lines_00.empty() == false){
            FindEndpoints(lines_02,lines_00,fEndpointRange,hitlist,ac_endpoints);
         }
         else{
            Hough_ind(hit_01,fHoughThreshold,fHoughMaxGap,fHoughRange,fHoughMinLength,fHoughMuonLength,event,lines_01);
            if (lines_01.empty()==false){
               FindEndpoints(lines_02,lines_01,fEndpointRange,hitlist,ac_endpoints);
            }
            else{
               std::cout << "tpc0: induction lines not found" << std::endl;
            }
         }
      }
      if (ac_in_tpc1){
         ac_tpc.insert(ac_tpc.end(),0,lines_12.size());
         Findt0(lines_12, ac_t0);
         Hough_ind(hit_10,fHoughThreshold,fHoughMaxGap,fHoughRange,fHoughMinLength,fHoughMuonLength,event,lines_10);
         if (lines_10.empty() == false){
            std::cout << "tpc1: line found on plane0 "  << std::endl;
            FindEndpoints(lines_12,lines_10,fEndpointRange,hitlist, ac_endpoints);
         }
         else{
            Hough_ind(hit_11,fHoughThreshold,fHoughMaxGap,fHoughRange,fHoughMinLength,fHoughMuonLength,event,lines_11);
            if (lines_11.empty() == false){
               std::cout << "tpc1: line found on plane1 "  << std::endl;
               FindEndpoints(lines_12,lines_11,fEndpointRange,hitlist, ac_endpoints);
            }
            else{
               std::cout << "tpc1: induction lines not found" << std::endl;
            }
         }
      }
      FindTrajectories(ac_endpoints, ac_trajectories);
      for (int i=0; i<int(ac_t0.size()); i++){
         double t0 = ac_t0.at(i); 
         geo::Point_t apoint = (ac_endpoints.at(i)).at(0), cpoint = (ac_endpoints.at(i)).at(1); 
         float apoint_x = float(apoint.X()), apoint_y = float(apoint.Y()), apoint_z = float(apoint.Z()); 
         float cpoint_x = float(cpoint.X()), cpoint_y = float(cpoint.Y()), cpoint_z = float(cpoint.Z()); 

         sbn::crt::CRTTrack mytrack;

         mytrack.ts0_ns = t0*1e-3; //t0 is currently in us
         mytrack.x1_pos = apoint_x; 
         mytrack.y1_pos = apoint_y; 
         mytrack.z1_pos = apoint_z; 
         mytrack.x2_pos = cpoint_x; 
         mytrack.y2_pos = cpoint_y; 
         mytrack.z2_pos = cpoint_z; 

         ac_tracks.push_back(mytrack);
      }
   } // end of finding ind hits
   return pass;
} // end of event loop 

float testFilter::Distance(int x1, int y1, int x2, int y2){
   return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) * 1.0);
}

void testFilter::Hough_col(vector<vector<int>> coords, int threshold, int max_gap, int range, float min_length, int muon_length, int nentry, 
                       vector<vector<int>>& lines, vector<vector<int>>& colwire, vector<vector<int>>& colpeakT){
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
   vector<vector<int>> outcolwire; 
   vector<vector<int>> outcolpeakT; 

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
      vector<vector<int>> endpoint(2, vector<int>(4));
      vector<int> wires; 
      vector<int> peakTs; 
      wires.clear(); peakTs.clear(); 
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
               endpoint[k] = {x1, y1, wire_idx, idx+i};
               (coords.at(idx+i)).clear();
               (data.at(idx+i)).clear();
               wires.push_back(x1);
               peakTs.push_back(y1); 
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
      int idx0_end = endpoint[0][3], idx1_end = endpoint[1][3];
      if ((x0_end==0 && y0_end==0) || (x1_end==0 && y1_end==0)) // don't add the (0,0) points 
         continue;
      vector<int> outline = {x0_end, y0_end, x1_end, y1_end, wire0_end, wire1_end, idx0_end, idx1_end, rho, theta};

      outlines.push_back(outline);
      outcolwire.push_back(wires);
      outcolpeakT.push_back(peakTs); 

   } // end of point loop 
   // combine lines that are split 
   for (int i=0; i<int(outlines.size()); i++){
      bool same = false;
      for (int j=i+1; j<int(outlines.size()) && same == false; j++){ 
         int xi_coords[2] = {outlines[i][0], outlines[i][2]}; int xj_coords[2] = {outlines[j][0], outlines[j][2]};
         int yi_coords[2] = {outlines[i][1], outlines[i][3]}; int yj_coords[2] = {outlines[j][1], outlines[j][3]};
         int rhoi = outlines[i][8], rhoj = outlines[j][8];
         int thetai = outlines[i][9], thetaj = outlines[j][9]; 

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
                  (outcolwire.at(j)).insert( (outcolwire.at(j)).end(),  (outcolwire.at(i)).begin(),  (outcolwire.at(i)).end()); 
                  (outcolpeakT.at(j)).insert( (outcolpeakT.at(j)).end(),  (outcolpeakT.at(i)).begin(),  (outcolpeakT.at(i)).end());
                  
                  (outlines.at(i)).clear();
                  (outcolwire.at(i)).clear();
                  (outcolpeakT.at(i)).clear();
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
            colwire.push_back(outcolwire.at(i)); 
            colpeakT.push_back(outcolpeakT.at(i));
         }
      }
      else{
         float length = Distance(x0_end,y0_end,x1_end,y1_end);
         if (length > min_length){
            lines.push_back(outlines.at(i));
            colwire.push_back(outcolwire.at(i)); 
            colpeakT.push_back(outcolpeakT.at(i));
         }
      }
   }
   //free memory 
   data.clear(); deaccu.clear(); outlines.clear(); outcolwire.clear(); outcolpeakT.clear(); 
} // end of hough 

void testFilter::Hough_ind(vector<vector<int>> coords, int threshold, int max_gap, int range, float min_length, int muon_length, int nentry, vector<vector<int>>& lines){
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
   colwire_02.clear(); 
   colpeakT_02.clear(); 
   colwire_12.clear(); 
   colpeakT_12.clear();   
  
   hit_02.reserve(3000); 
   hit_12.reserve(3000); 
   lines_02.reserve(10); 
   lines_12.reserve(10); 
   colwire_02.reserve(3000); 
   colpeakT_02.reserve(3000); 
   colwire_12.reserve(3000); 
   colpeakT_12.reserve(3000); 
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

   lines_00.reserve(5); 
   lines_01.reserve(5);
   lines_10.reserve(5); 
   lines_11.reserve(5);
}

void testFilter::ResetACVariables(){
   ac_t0.clear();
   ac_endpoints.clear();
   ac_trajectories.clear(); 
   
   ac_tracks.clear();

   ac_t0.reserve(5); 
   ac_endpoints.reserve(5);
   ac_trajectories.reserve(5);

   ac_tracks.reserve(5);
}

// void testFilter::FindCollectionPlaneHits(vector<vector<int>> colwire, vector<vector<int>> colpeakT){
//    for (int i=0; i<int(colwire.size()); i++){
//       vector<int> cwire = (colwire.at(i)); 
//       vector<int> cpeakT = (colpeakT.at(i)); 
//    }
// }

void testFilter::Findt0(vector<vector<int>> lines,
                        vector<double>& ac_t0){
   for (int i=0; i<int(lines.size()); i++){
      double t_i, t_j; 
      t_i = (lines[i][1]-500)*0.5;
      t_j = (lines[i][3]-500)*0.5;
      double t0 = (t_i<t_j)? t_i:t_j;
      std::cout << "interaction times (in us): " << t0 << std::endl;
      ac_t0.push_back(t0);
   }           
}

void testFilter::FindEndpoints( vector<vector<int>> lines_col, vector<vector<int>> lines_ind, 
                                int range, vector<art::Ptr<recob::Hit>> hitlist, 
                                vector<vector<geo::Point_t>>& ac_endpoints){
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
         // std:: cout << "peakT0_col, peakT1_col, peakT0_ind, peakT1_ind: " << peakT0_col << ", " << peakT1_col << ", "<< peakT0_ind << ", "<< peakT1_ind << std::endl;
         // std:: cout << "wire0_col, wire1_col, wire0_ind, wire1_ind: " << wire0_col << ", " << wire1_col << ", "<< wire0_ind << ", "<< wire1_ind << std::endl;

         int peakT_range = range; 
         if ( (peakT0_col + peakT_range > peakT0_ind) && (peakT0_col - peakT_range < peakT0_ind) && 
              (peakT1_col + peakT_range > peakT1_ind) && (peakT1_col + peakT_range > peakT1_ind)){
               
            geo::WireID awire_col = hitlist[wire0_col]->WireID(); 
            geo::WireID awire_ind = hitlist[wire0_ind]->WireID();
            geo::WireID cwire_col = hitlist[wire1_col]->WireID();
            geo::WireID cwire_ind = hitlist[wire1_ind]->WireID();
            geo::Point_t apoint, cpoint;

            bool aintersect = fGeometryService->WireIDsIntersect(awire_col, awire_ind, apoint);
            bool cintersect = fGeometryService->WireIDsIntersect(cwire_col, cwire_ind, cpoint); 

            if (aintersect)
               std::cout << "anode endpoint: " << apoint.X() << ", " << apoint.Y() << ", " << apoint.Z() << std::endl;
            else
               std::cout << "intersection of awire not found by WireIDsIntersect" << std::endl;
            
            if (cintersect){
               cpoint.SetX(0);
               std::cout << "cathode endpoint: " << cpoint.X() << ", " << cpoint.Y() << ", " << cpoint.Z() << std::endl;
            }
            else
               std::cout << "intersection of cwire not found by WireIDsIntersect"<< std::endl;
            
            vector<geo::Point_t> pair{apoint,cpoint};
            ac_endpoints.push_back(pair); 
         }
      } // end of j loop 
   } // end of i loop 
} // end of function 

void testFilter::FindTrajectories(vector<vector<geo::Point_t>> ac_endpoints, vector<vector<double>>& ac_trajectories){
   if (ac_endpoints.empty() == false){
      for (int i=0; i<int(ac_endpoints.size()); i++){
         vector<geo::Point_t> pair = ac_endpoints[i];
         float x_a = float(pair[0].X()), y_a = float(pair[0].Y()), z_a = float(pair[0].Z()); 
         float x_c = float(pair[1].X()), y_c = float(pair[1].Y()), z_c = float(pair[1].Z()); 
         float dx = x_a - x_c, dy = y_a - y_c, dz = z_a - z_c; 
         double theta_xz = atan2(dx,dz) * 180/PI;
         double theta_yz = atan2(dy,dz) * 180/PI;
         std::cout << "theta_xz: " << theta_xz << std::endl;
         std::cout << "theta_yz: " << theta_yz << std::endl;
         vector<double> trajectories{theta_xz, theta_yz};
         ac_trajectories.push_back(trajectories);
      }
   }
}
// A macro required for a JobControl module.
DEFINE_ART_MODULE(testFilter)
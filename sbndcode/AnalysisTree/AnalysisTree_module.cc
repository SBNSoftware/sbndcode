////////////////////////////////////////////////////////////////////////
// $Id: DBSCANfinderAna.cxx,v 1.36 2010/09/15  bpage Exp $
//
// module to create a TTree for analysis
//
// \author tjyang@fnal.gov, sowjanyag@phys.ksu.edu
//
////////////////////////////////////////////////////////////////////////
// To reduce memory usage:
// [x] create the data structure connected to the tree only when needed
// [x] reduce the size of the elemental items (Double_t => Float_t could damage precision)
// [x] create a different structure for each tracker, allocate only what needed
// [x] use variable size array buffers for each tracker datum instead of [kMaxTrack]
// [x] turn the truth/GEANT information into vectors
// [ ] move hit_trkid into the track information, remove kMaxTrackers
// [ ] turn the hit information into vectors (~1 MB worth), remove kMaxHits
// [ ] fill the tree branch by branch
// 
// Current implementation:
// There is one tree only, with one set of branches for each tracking algorithm.
// The data structure which hosts the addresses of the tree branches is
// dynamically allocated on demand, and it can be optionally destroyed at the
// end of each event.
// The data structure (AnalysisTreeDataStruct) firectly contains the truth and
// simulation information as C arrays. The data from tracking algorithms is the
// largest, and it is contained in a C++ vector of structures (TrackDataStruct),
// one per algorithm. These structures can also be allocated on demand.
// Each of these structures is connected to a set of branches, one branch per
// data member. Data members are vectors of numbers or vectors of fixed-size
// C arrays. The vector index represents the tracks reconstructed by the
// algorithm, and each has a fixed size pool for hits (do ROOT t3rees support
// branches with more than one dimension with variable size?).
// The data structures can assign default values to their data, connect to a
// ROOT tree (creating the branches they need) and resize.
// The AnalysisTreeDataStruct is constructed with as many tracking algorithms as
// there are named in the module configuration (even if they are not backed by
// any available tracking data).
// By default construction, TrackDataStruct is initialized in a state which does
// not allow any track (maximum tracks number is zero), and in such state trying
// to connect to a tree has no effect. This is done so that the
// AnalysisTreeDataStruct can be initialized first (and with unusable track data
// structures), and then the TrackDataStruct instances are initialized one by
// one when the number of tracks needed is known.
// A similar mechanism is implemented for the truth information.
// 
// The "UseBuffers: false" mode assumes that on each event a new
// AnalysisTreeDataStruct is created with unusable tracker data, connected to
// the ROOT tree (the addresses of the available branches are assigned), then
// each of the tracking algorithm data is resized to host the correct number
// of reconstructed tracks and connected to the tree. Then the normal process of
// filling the event data and then the tree take place. Finally, the whole
// data structure is freed and the tree is left in a invalid state (branch
// addresses are invalid). It could be possible to make the tree in a valid
// state by resetting the addresses, but there is no advantage in that.
// 
// The "UseBuffers: true" mode assumes that on the first event a new
// AnalysisTreeDataStruct is created and used just as in the other mode
// described above. At the end of the first event, the data structure is left
// around (and the tree is in a valid state). On the next event, all the
// addresses are checked, then for each tracker the data is resized to
// accomodate the right number of tracks for tis event. If the memory is
// increased, the address will be changed. All the branches are reconnected to
// the data structure, and the procedure goes on as normal.
// 
// Note that reducing the maximum number of tracks in a TrackDataStruct does not
// necessarily make memory available, because of how std::vector::resize()
// works; that feature can be implemented, but it currently has not been.
// 
// The BoxedArray<> class is a wrapper around a normal C array; it is needed
// to be able to include such structure in a std::vector. This container
// requires its objects to be default-constructable and copy-constructable,
// and a C array is neither. BoxedArray<> is: the default construction leaves it
// uninitialized (for speed reasons) while the copy construction is performed
// as in a Plain Old Data structure (memcpy; really!).
// 
////////////////////////////////////////////////////////////////////////

#ifndef ANALYSISTREE_H
#define ANALYSISTREE_H

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/BeamInfo.h"
// #include "larcore/DetectorInfoServices/LArPropertiesService.h"
// #include "larcore/Utilities/AssociationUtil.h"
// #include "larcore/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larreco/RecoAlg/PMAlg/PmaTrack3D.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "sbndcode/RecoUtils/RecoUtils.h"

#include <cstring> // std::memcpy()
#include <vector>
#include <map>
#include <utility>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <functional> // std::mem_fn()
#include <typeinfo>
#include <cmath>

#include "TTree.h"
#include "TTimeStamp.h"

constexpr int kNplanes       = 3;     //number of wire planes
constexpr int kMaxHits       = 25000; //maximum number of hits;
constexpr int kMaxTrackHits  = 2000;  //maximum number of hits on a track
constexpr int kMaxTrackers   = 15;    //number of trackers passed into fTrackModuleLabel
constexpr unsigned short kMaxVertices = 100;    //max number of 3D vertices
constexpr unsigned short kMaxShowers = 100;    //max number of 3D showers
constexpr unsigned short kMaxAuxDets  = 4; ///< max number of auxiliary detector cells per MC particle

/// total_extent\<T\>::value has the total number of elements of an array
template <typename T>
struct total_extent {
  using value_type = size_t;
  static constexpr value_type value
    = sizeof(T) / sizeof(typename std::remove_all_extents<T>::type);
}; // total_extent<>


namespace sbnd {

  /// Data structure with all the tree information.
  /// 
  /// Can connect to a tree, clear its fields and resize its data.
  class AnalysisTreeDataStruct {
      public:
    
    /// A wrapper to a C array (needed to embed an array into a vector)
    template <typename Array_t>
    class BoxedArray {
        protected:
      Array_t array; // actual data
      
        public:
      using This_t = BoxedArray<Array_t>;
      typedef typename std::remove_all_extents<Array_t>::type Data_t;
      
      BoxedArray() {} // no initialization
      BoxedArray(const This_t& from)
        { std::memcpy((char*) &(data()), (char*) &(from.data()), sizeof(Array_t)); }
      
      Array_t& data() { return array; }
      const Array_t& data() const { return array; }
      
      //@{
      /// begin/end interface
      static constexpr size_t size() { return total_extent<Array_t>::value; }
      Data_t* begin() { return reinterpret_cast<Data_t*>(&array); }
      const Data_t* begin() const { return reinterpret_cast<const Data_t*>(&array); }
      Data_t* end() { return begin() + size(); }
      const Data_t* end() const { return begin() + size(); }
      //@}
      
      //@{
      /// Array interface
      auto operator[] (size_t index) -> decltype(*array) { return array[index]; }
      auto operator[] (size_t index) const -> decltype(*array) { return array[index]; }
      auto operator+ (ptrdiff_t index) -> decltype(&*array) { return array + index; }
      auto operator+ (ptrdiff_t index) const -> decltype(&*array) { return array + index; }
      auto operator- (ptrdiff_t index) -> decltype(&*array) { return array - index; }
      auto operator- (ptrdiff_t index) const -> decltype(&*array) { return array - index; }
      auto operator* () -> decltype(*array) { return *array; }
      auto operator* () const -> decltype(*array) { return *array; }
      
      operator decltype(&array[0]) () { return &array[0]; }
      operator decltype(&array[0]) () const { return &array[0]; }
      //@}
      
    }; // BoxedArray
    
    /// Tracker algorithm result
    /// 
    /// Can connect to a tree, clear its fields and resize its data.
    class TrackDataStruct {
        public:
      /* Data structure size:
       *
       * TrackData_t<Short_t>                    :  2  bytes/track
       * TrackData_t<Float_t>                    :  4  bytes/track
       * PlaneData_t<Float_t>, PlaneData_t<Int_t>: 12  bytes/track
       * HitData_t<Float_t>                      : 24k bytes/track
       * HitCoordData_t<Float_t>                 : 72k bytes/track
       */
      template <typename T>
      using TrackData_t = std::vector<T>;
      template <typename T>
      using PlaneData_t = std::vector<BoxedArray<T[kNplanes]>>;
      template <typename T>
      using HitData_t = std::vector<BoxedArray<T[kNplanes][kMaxTrackHits]>>;
      template <typename T>
      using HitCoordData_t = std::vector<BoxedArray<T[kNplanes][kMaxTrackHits][3]>>;
      
      size_t MaxTracks; ///< maximum number of storable tracks
      
      Short_t  ntracks;             //number of reconstructed tracks
      PlaneData_t<Float_t>    trkke;
      PlaneData_t<Float_t>    trkrange;
      PlaneData_t<Int_t>      trkidtruth_recoutils_totaltrueenergy;  //true geant trackid from TrueParticleIDFromTotalTrueEnergy in RecoUtils 
      PlaneData_t<Int_t>      trkidtruth_recoutils_totalrecocharge;  //true geant trackid from TrueParticleIDFromTotalRecoCharge in RecoUtils 
      PlaneData_t<Int_t>      trkidtruth_recoutils_totalrecohits;  //true geant trackid from TrueParticleIDFromTotalTrueHits in RecoUtils 
      PlaneData_t<Int_t>      trkidtruth;  //true geant trackid
      PlaneData_t<Short_t>    trkorigin;   //_ev_origin 0: unknown, 1: cosmic, 2: neutrino, 3: supernova, 4: singles
      PlaneData_t<Int_t>      trkpdgtruth; //true pdg code
      PlaneData_t<Float_t>    trkefftruth; //completeness
      PlaneData_t<Float_t>    trksimIDEenergytruth;
      PlaneData_t<Float_t>    trksimIDExtruth;
      PlaneData_t<Float_t>    trksimIDEytruth;
      PlaneData_t<Float_t>    trksimIDEztruth;
      PlaneData_t<Float_t>    trkpurtruth; //purity of track
      PlaneData_t<Float_t>    trkpitchc;
      PlaneData_t<Short_t>    ntrkhits;
      HitData_t<Float_t>      trkdedx;
      HitData_t<Float_t>      trkdqdx;
      HitData_t<Float_t>      trkresrg;
      HitCoordData_t<Float_t> trkxyz;

      // more track info
      TrackData_t<Short_t> trkId;
      TrackData_t<Short_t> trkncosmictags_tagger;
      TrackData_t<Float_t> trkcosmicscore_tagger;
      TrackData_t<Short_t> trkcosmictype_tagger;
      TrackData_t<Short_t> trkncosmictags_flashmatch;
      TrackData_t<Float_t> trkcosmicscore_flashmatch;
      TrackData_t<Short_t> trkcosmictype_flashmatch;
      TrackData_t<Float_t> trkstartx;     // starting x position.
      TrackData_t<Float_t> trkstarty;     // starting y position.
      TrackData_t<Float_t> trkstartz;     // starting z position.
      TrackData_t<Float_t> trkstartd;     // starting distance to boundary.
      TrackData_t<Float_t> trkendx;       // ending x position.
      TrackData_t<Float_t> trkendy;       // ending y position.
      TrackData_t<Float_t> trkendz;       // ending z position.
      TrackData_t<Float_t> trkendd;       // ending distance to boundary.
      TrackData_t<Float_t> trktheta;      // theta.
      TrackData_t<Float_t> trkphi;        // phi.
      TrackData_t<Float_t> trkstartdcosx;
      TrackData_t<Float_t> trkstartdcosy;
      TrackData_t<Float_t> trkstartdcosz;
      TrackData_t<Float_t> trkenddcosx;
      TrackData_t<Float_t> trkenddcosy;
      TrackData_t<Float_t> trkenddcosz;
      TrackData_t<Float_t> trkthetaxz;    // theta_xz.
      TrackData_t<Float_t> trkthetayz;    // theta_yz.
      TrackData_t<Float_t> trkmom;        // momentum.
      TrackData_t<Float_t> trklen;        // length.
      TrackData_t<Float_t> trkmomrange;   // track momentum from range using CSDA tables
      TrackData_t<Float_t> trkmommschi2;  // track momentum from multiple scattering Chi2 method
      TrackData_t<Float_t> trkmommsllhd;  // track momentum from multiple scattering LLHD method
      TrackData_t<Short_t> trksvtxid;     // Vertex ID associated with the track start
      TrackData_t<Short_t> trkevtxid;     // Vertex ID associated with the track end
      PlaneData_t<Int_t>   trkpidpdg;     // particle PID pdg code
      PlaneData_t<Float_t> trkpidchi;
      PlaneData_t<Float_t> trkpidchipr;   // particle PID chisq for proton
      PlaneData_t<Float_t> trkpidchika;   // particle PID chisq for kaon
      PlaneData_t<Float_t> trkpidchipi;   // particle PID chisq for pion
      PlaneData_t<Float_t> trkpidchimu;   // particle PID chisq for muon
      PlaneData_t<Float_t> trkpidpida;    // particle PIDA
      TrackData_t<Short_t> trkpidbestplane; // this is defined as the plane with most hits     
      
      /// If SaveHierarchyInfo is true, there are additional variables:
      TrackData_t<Int_t>   trkisprimary;   // If the track is a daughter of the neutrino
      TrackData_t<Int_t>   trkndaughters;  // Number of daughters the track has
      TrackData_t<Int_t>   trkpfpid;       // The track's pfparticle ID
      TrackData_t<Int_t>   trkparentpfpid; // The parent of the track's pfparticle ID

      /// Creates an empty tracker data structure
      TrackDataStruct(): MaxTracks(0) { Clear(); }
      /// Creates a tracker data structure allowing up to maxTracks tracks
      TrackDataStruct(size_t maxTracks): MaxTracks(maxTracks) { Clear(); }
      void Clear();
      void SetMaxTracks(size_t maxTracks)
        { MaxTracks = maxTracks; Resize(MaxTracks); }
      void Resize(size_t nTracks);
      void SetAddresses(TTree* pTree, std::string tracker, bool isCosmics, bool saveHierarchyInfo);
      
      size_t GetMaxTracks() const { return MaxTracks; }
      size_t GetMaxPlanesPerTrack(int /* iTrack */ = 0) const
        { return (size_t) kNplanes; }
      size_t GetMaxHitsPerTrack(int /* iTrack */ = 0, int /* ipl */ = 0) const
        { return (size_t) kMaxTrackHits; }
      
    }; // class TrackDataStruct
    
    /*
     * Added by Rhiannon 8th October 2018
     */
    /// Shower algorithm result
    /// 
    /// Can connect to a tree, clear its fields and resize its data.
    class ShowerDataStruct {
        public:
      
      template <typename T>
      using ShowerData_t = std::vector<T>;
      
      size_t MaxShowers; ///< maximum number of storable showers
      
      Short_t  nshowers;                  // Number of reconstructed showers
      ShowerData_t<Int_t>   shwId;        // shower ID
      ShowerData_t<Int_t>   shwbestplane; // shower best plane
      ShowerData_t<Float_t> shwlength;    // Length of the shower
      ShowerData_t<Float_t> shwopenangle; // Opening angle of the shower
      ShowerData_t<Float_t> shwstartx;    // Shower start position, x 
      ShowerData_t<Float_t> shwstarty;    // Shower start position, y 
      ShowerData_t<Float_t> shwstartz;    // Shower start position, z 
      ShowerData_t<Float_t> shwdirx;      // Shower direction, x 
      ShowerData_t<Float_t> shwdiry;      // Shower direction, y 
      ShowerData_t<Float_t> shwdirz;      // Shower direction, z
      
      /// If SaveHierarchyInfo is true, there are additional variables:
      ShowerData_t<Int_t>   shwisprimary; // If the shower is a primary particle
      ShowerData_t<Int_t>   shwndaughters; // Number of daughters the shower has
      ShowerData_t<Int_t>   shwpfpid;      // PFParticle shower ID
      ShowerData_t<Int_t>   shwparentpfpid;      // PFParticle shower ID
      
      /// Creates an empty showerer data structure
      ShowerDataStruct(): MaxShowers(0) { Clear(); }
      /// Creates a shower data structure allowing up to maxShowers showers
      ShowerDataStruct(size_t maxShowers): MaxShowers(maxShowers) { Clear(); }
      void Clear();
      void SetMaxShowers(size_t maxShowers)
        { MaxShowers = maxShowers; Resize(MaxShowers); }
      void Resize(size_t nShowers);
      void SetShowerAddresses(TTree* pTree, std::string showerLabel, bool saveHierarchyInfo);
      
      size_t GetMaxShowers() const { return MaxShowers; }
      size_t GetMaxPlanesPerShower(int /* iShower */ = 0) const
        { return (size_t) kNplanes; }
      
    }; // class ShowerDataStruct
    
    /*
     * Added by Rhiannon 8th October 2018
     */
    /// Vertices
    /// 
    /// Can connect to a tree, clear its fields and resize its data.
    class VertexDataStruct {
        public:
      
      template <typename T>
      using VertexData_t = std::vector<T>;
      template <typename T>
      using VtxCoordData_t = std::vector<BoxedArray<T[3]> >;
      
      size_t MaxVertices; ///< maximum number of storable vertices
      
      // vertex information
      Short_t  nvtx;                  // Number of reconstructed vertices
      VtxCoordData_t<Float_t> vtx;    // 3D position of the vertex
      VertexData_t<Int_t> primaryvtx; // Whether the vertex is primary or not

      /// Creates an empty tracker data structure
      VertexDataStruct(): MaxVertices(0) { Clear(); }
      /// Creates a vertex data structure allowing up to maxVertices vertices
      VertexDataStruct(size_t maxVertices): MaxVertices(maxVertices) { Clear(); }
      void Clear();
      void SetMaxVertices(size_t maxVertices)
        { MaxVertices = maxVertices; Resize(MaxVertices); }
      void Resize(size_t nVertices);
      void SetAddresses(TTree* pTree, std::string vertexLabel, bool saveHierarchyInfo);
      
      size_t GetMaxVertices() const { return MaxVertices; }
      
    }; // class VertexDataStruct
 
    enum DataBits_t: unsigned int {
      tdAuxDet = 0x01,
      tdCry = 0x02,
      tdGenie = 0x04,
      tdGeant = 0x08,
      tdHit = 0x10,
      tdTrack = 0x20,
      tdShower = 0x20,
      tdVtx = 0x40,
      tdDefault = 0
    }; // DataBits_t
    
/*    /// information from the run
    struct RunData_t {
        public:
      RunData_t() { Clear(); }
      void Clear() {}
    }; // struct RunData_t
*/
    /// information from the subrun
    struct SubRunData_t {
      SubRunData_t() { Clear(); }
      void Clear() { pot = -99999.; }
      Double_t pot; //protons on target
    }; // struct SubRunData_t

//    RunData_t    RunData; ///< run data collected at begin of run
    SubRunData_t SubRunData; ///< subrun data collected at begin of subrun

    //run information
    Int_t      run;                  //run number
    Int_t      subrun;               //subrun number
    Int_t      event;                //event number
    Double_t   evttime;              //event time in sec
    Double_t   beamtime;             //beam time
  //  Double_t   pot;                  //protons on target moved in subrun data
    // Double_t   taulife;              //electron lifetime
    Char_t     isdata;               //flag, 0=MC 1=data

    // hit information (non-resizeable, 45x kMaxHits = 900k bytes worth)
    Int_t    no_hits;                  //number of hits
    Short_t  hit_tpc[kMaxHits];        //tpc number
    Short_t  hit_plane[kMaxHits];      //plane number
    Short_t  hit_wire[kMaxHits];       //wire number
    Short_t  hit_channel[kMaxHits];    //channel ID
    Float_t  hit_peakT[kMaxHits];      //peak time
    Float_t  hit_charge[kMaxHits];     //charge (area)
    Float_t  hit_ph[kMaxHits];         //amplitude
    Float_t  hit_startT[kMaxHits];     //hit start time
    Float_t  hit_endT[kMaxHits];       //hit end time
    Float_t  hit_nelec[kMaxHits];     //hit number of electrons
    Float_t  hit_energy[kMaxHits];       //hit energy
    Short_t  hit_trkid[kMaxHits];      //is this hit associated with a reco track?

    // track information
    Char_t kNTracker;
    std::vector<TrackDataStruct> TrackData;

    // vertex information
    std::vector<VertexDataStruct> VertexData;

    // shower information
    ShowerDataStruct ShowerData;

    //mctruth information
    size_t MaxMCNeutrinos;     ///! The number of MCNeutrinos there is currently room for
    Int_t     mcevts_truth;    //number of neutrino Int_teractions in the spill
    std::vector<Int_t>     nuScatterCode_truth; //Scattering code given by Genie for each neutrino
    std::vector<Int_t>     nuID_truth;          //Unique ID of each true neutrino
    std::vector<Int_t>     nuPDG_truth;         //neutrino PDG code
    std::vector<Int_t>     ccnc_truth;          //0=CC 1=NC
    std::vector<Int_t>     mode_truth;          //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
    std::vector<Float_t>   enu_truth;           //true neutrino energy
    std::vector<Float_t>   Q2_truth;            //Momentum transfer squared
    std::vector<Float_t>   W_truth;             //hadronic invariant mass
    std::vector<Int_t>     hitnuc_truth;        //hit nucleon
    std::vector<Float_t>   nuvtxx_truth;        //neutrino vertex x
    std::vector<Float_t>   nuvtxy_truth;        //neutrino vertex y
    std::vector<Float_t>   nuvtxz_truth;        //neutrino vertex z
    std::vector<Float_t>   nu_dcosx_truth;      //neutrino dcos x
    std::vector<Float_t>   nu_dcosy_truth;      //neutrino dcos y
    std::vector<Float_t>   nu_dcosz_truth;      //neutrino dcos z
    std::vector<Float_t>   lep_mom_truth;       //lepton momentum
    std::vector<Float_t>   lep_dcosx_truth;     //lepton dcos x
    std::vector<Float_t>   lep_dcosy_truth;     //lepton dcos y
    std::vector<Float_t>   lep_dcosz_truth;     //lepton dcos z

    //flux information
    std::vector<Float_t>  tpx_flux;        //Px of parent particle leaving BNB target
    std::vector<Float_t>  tpy_flux;        //Py of parent particle leaving BNB target
    std::vector<Float_t>  tpz_flux;        //Pz of parent particle leaving BNB target
    std::vector<Int_t>     tptype_flux;     //Type of parent particle leaving BNB target

    //genie information
    size_t MaxGeniePrimaries = 0;
    Int_t     genie_no_primaries;
    std::vector<Int_t>     genie_primaries_pdg;
    std::vector<Float_t>  genie_Eng;
    std::vector<Float_t>  genie_Px;
    std::vector<Float_t>  genie_Py;
    std::vector<Float_t>  genie_Pz;
    std::vector<Float_t>  genie_P;
    std::vector<Int_t>     genie_status_code;
    std::vector<Float_t>  genie_mass;
    std::vector<Int_t>     genie_trackID;
    std::vector<Int_t>     genie_ND;
    std::vector<Int_t>     genie_mother;
    
    //cosmic cry information
    Int_t     mcevts_truthcry;    //number of neutrino Int_teractions in the spill
    Int_t     cry_no_primaries;
    std::vector<Int_t>    cry_primaries_pdg;
    std::vector<Float_t>  cry_Eng;
    std::vector<Float_t>  cry_Px;
    std::vector<Float_t>  cry_Py;
    std::vector<Float_t>  cry_Pz;
    std::vector<Float_t>  cry_P;
    std::vector<Float_t>  cry_StartPointx;
    std::vector<Float_t>  cry_StartPointy;
    std::vector<Float_t>  cry_StartPointz;
    std::vector<Int_t>    cry_status_code;
    std::vector<Float_t>  cry_mass;
    std::vector<Int_t>    cry_trackID;
    std::vector<Int_t>    cry_ND;
    std::vector<Int_t>    cry_mother;
    
    //geant information
    size_t MaxGEANTparticles = 0; ///! how many particles there is currently room for
    Int_t     no_primaries;      //number of primary geant particles
    Int_t     geant_list_size;  //number of all geant particles
    Int_t     geant_list_size_in_tpcAV;
    std::vector<Int_t>    pdg;
    std::vector<Int_t>    status;    
    std::vector<Float_t>  Eng;
    std::vector<Float_t>  EndE;
    std::vector<Float_t>  Mass;
    std::vector<Float_t>  Px;
    std::vector<Float_t>  Py;
    std::vector<Float_t>  Pz;
    std::vector<Float_t>  P;
    std::vector<Float_t>  StartPointx;
    std::vector<Float_t>  StartPointy;
    std::vector<Float_t>  StartPointz;
    std::vector<Float_t>  StartT;  
    std::vector<Float_t>  EndT;          
    std::vector<Float_t>  EndPointx;
    std::vector<Float_t>  EndPointy;
    std::vector<Float_t>  EndPointz;
    std::vector<Float_t>  theta;    
    std::vector<Float_t>  phi;    
    std::vector<Float_t>  theta_xz;    
    std::vector<Float_t>  theta_yz;    
    std::vector<Float_t>  pathlen;    
    std::vector<Int_t>    inTPCActive;    
    std::vector<Float_t>  StartPointx_tpcAV;
    std::vector<Float_t>  StartPointy_tpcAV;
    std::vector<Float_t>  StartPointz_tpcAV;
    std::vector<Float_t>  EndPointx_tpcAV;
    std::vector<Float_t>  EndPointy_tpcAV;
    std::vector<Float_t>  EndPointz_tpcAV;
    std::vector<Int_t>    NumberDaughters;
    std::vector<Int_t>    TrackId;
    std::vector<Int_t>    Mother;
    std::vector<Int_t>    process_primary;
    std::vector<std::string> processname;
    std::vector<Int_t>    MergedId; //geant track segments, which belong to the same particle, get the same
    std::vector<Int_t>    MotherNuId; //The unique ID of the mother neutrino
    
    // Auxiliary detector variables saved for each geant track
    // This data is saved as a vector (one item per GEANT particle) of C arrays
    // (wrapped in a BoxedArray for technical reasons), one item for each
    // affected detector cell (which one is saved in AuxDetID
    template <typename T>
    using AuxDetMCData_t = std::vector<BoxedArray<T[kMaxAuxDets]>>;
    
    std::vector<UShort_t> NAuxDets;         ///< Number of AuxDets crossed by this particle
    AuxDetMCData_t<Short_t> AuxDetID;       ///< Which AuxDet this particle went through
    AuxDetMCData_t<Float_t> entryX;         ///< Entry X position of particle into AuxDet
    AuxDetMCData_t<Float_t> entryY;         ///< Entry Y position of particle into AuxDet
    AuxDetMCData_t<Float_t> entryZ;         ///< Entry Z position of particle into AuxDet
    AuxDetMCData_t<Float_t> entryT;         ///< Entry T position of particle into AuxDet
    AuxDetMCData_t<Float_t> exitX;          ///< Exit X position of particle out of AuxDet
    AuxDetMCData_t<Float_t> exitY;          ///< Exit Y position of particle out of AuxDet
    AuxDetMCData_t<Float_t> exitZ;          ///< Exit Z position of particle out of AuxDet
    AuxDetMCData_t<Float_t> exitT;          ///< Exit T position of particle out of AuxDet
    AuxDetMCData_t<Float_t> exitPx;         ///< Exit x momentum of particle out of AuxDet
    AuxDetMCData_t<Float_t> exitPy;         ///< Exit y momentum of particle out of AuxDet
    AuxDetMCData_t<Float_t> exitPz;         ///< Exit z momentum of particle out of AuxDet
    AuxDetMCData_t<Float_t> CombinedEnergyDep; ///< Sum energy of all particles with this trackID (+ID or -ID) in AuxDet
    
    unsigned int bits; ///< complementary information
 
    /// Returns whether we have auxiliary detector data
    bool hasAuxDetector() const { return bits & tdAuxDet; }
    
    /// Returns whether we have Cry data
    bool hasCryInfo() const { return bits & tdCry; }
    
    /// Returns whether we have Genie data
    bool hasGenieInfo() const { return bits & tdGenie; }
    
    /// Returns whether we have Hit data
    bool hasHitInfo() const { return bits & tdHit; }

    /// Returns whether we have Track data
    bool hasTrackInfo() const { return bits & tdTrack; }
    
    /// Returns whether we have Track data
    bool hasShowerInfo() const { return bits & tdShower; }
    
    /// Returns whether we have Vertex data
    bool hasVertexInfo() const { return bits & tdVtx; }
    
    /// Returns whether we have Geant data
    bool hasGeantInfo() const { return bits & tdGeant; }

    /// Sets the specified bits
    void SetBits(unsigned int setbits, bool unset = false)
      { if (unset) bits &= ~setbits; else bits |= setbits; }
      
    /// Constructor; clears all fields
    AnalysisTreeDataStruct(size_t nTrackers = 0): bits(tdDefault) 
      { SetTrackers(nTrackers); SetVertices(nTrackers); Clear(); }

    TrackDataStruct& GetTrackerData(size_t iTracker)
      { return TrackData.at(iTracker); }
    const TrackDataStruct& GetTrackerData(size_t iTracker) const
      { return TrackData.at(iTracker); }
    
    VertexDataStruct& GetVertexData(size_t iTracker)
      { return VertexData.at(iTracker); }
    
    ShowerDataStruct& GetShowerData()
      { return ShowerData; }
    
    /// Clear all fields if this object (not the tracker algorithm data)
    void ClearLocalData();
    
    /// Clear all fields
    void Clear();
    
    /// Allocates data structures for the given number of trackers (no Clear())
    void SetTrackers(size_t nTrackers) { TrackData.resize(nTrackers); }

    /// Allocates data structures for the given number of trackers (no Clear())
    void SetVertices(size_t nTrackers) { VertexData.resize(nTrackers); }

    /// Resize the data structure for MCNeutrino particles
    void ResizeMCNeutrino(int nNeutrinos);
    
    /// Resize the data strutcure for GEANT particles
    void ResizeGEANT(int nParticles);
    
    /// Resize the data strutcure for Genie primaries
    void ResizeGenie(int nPrimaries);
    
    /// Resize the data strutcure for Cry primaries
    void ResizeCry(int nPrimaries);
    
    /// Connect this object with a tree
    void SetAddresses(TTree* pTree, const std::vector<std::string>& trackers, std::string showerLabel, const std::vector< std::string> &vertexLabels, bool isCosmics, const std::vector<bool>& saveHierarchyInfo, bool saveShowerHierarchyInfo);
    
    /// Connect this object with a tree
    void SetShowerAddresses(TTree* pTree, std::string showerLabel, bool saveHierarchyInfo);
    
    /// Returns the number of trackers for which data structures are allocated
    size_t GetNTrackers() const { return TrackData.size(); }
    
    /// Returns the number of hits for which memory is allocated
    size_t GetMaxHits() const { return kMaxHits; }
    
    /// Returns the number of trackers for which memory is allocated
    size_t GetMaxTrackers() const { return TrackData.capacity(); }
    
    /// Returns the number of GEANT particles for which memory is allocated
    size_t GetMaxGEANTparticles() const { return MaxGEANTparticles; }
        
    /// Returns the number of GENIE primaries for which memory is allocated
    size_t GetMaxGeniePrimaries() const { return MaxGeniePrimaries; }
    
    
  private:
    /// Little helper functor class to create or reset branches in a tree
    class BranchCreator {
        public:
      TTree* pTree; ///< the tree to be worked on
      BranchCreator(TTree* tree): pTree(tree) {}

      //@{
      /// Create a branch if it does not exist, and set its address
      void operator()
        (std::string name, void* address, std::string leaflist /*, int bufsize = 32000 */)
        {
          if (!pTree) return;
          TBranch* pBranch = pTree->GetBranch(name.c_str());
          if (!pBranch) {
            pTree->Branch(name.c_str(), address, leaflist.c_str() /*, bufsize */);
            MF_LOG_DEBUG("AnalysisTreeStructure")
              << "Creating branch '" << name << " with leaf '" << leaflist << "'";
          }
          else if (pBranch->GetAddress() != address) {
            pBranch->SetAddress(address);
            MF_LOG_DEBUG("AnalysisTreeStructure")
              << "Reassigning address to branch '" << name << "'";
          }
          else {
            MF_LOG_DEBUG("AnalysisTreeStructure")
              << "Branch '" << name << "' is fine";
          }
        } // operator()
      void operator()
        (std::string name, void* address, const std::stringstream& leaflist /*, int bufsize = 32000 */)
        { return this->operator() (name, address, leaflist.str() /*, int bufsize = 32000 */); }
      template <typename T>
      void operator()
        (std::string name, std::vector<T>& data, std::string leaflist /*, int bufsize = 32000 */)
        { return this->operator() (name, (void*) data.data(), leaflist /*, int bufsize = 32000 */); }

      template <typename T>
      void operator() (std::string name, std::vector<T>& data)
        {
          // overload for a generic object expressed directly by reference
          // (as opposed to a generic object expressed by a pointer or
          // to a simple leaf sequence specification);
          // TTree::Branch(name, T* obj, Int_t bufsize, splitlevel) and
          // TTree::SetObject() are used.
          if (!pTree) return;
          TBranch* pBranch = pTree->GetBranch(name.c_str());
          if (!pBranch) {
            pTree->Branch(name.c_str(), &data);
            // ROOT needs a TClass definition for T in order to create a branch,
            // se we are sure that at this point the TClass exists
            MF_LOG_DEBUG("AnalysisTreeStructure")
              << "Creating object branch '" << name
              << " with " << TClass::GetClass(typeid(T))->ClassName();
          }
          else if
            (*(reinterpret_cast<std::vector<T>**>(pBranch->GetAddress())) != &data)
          {
            // when an object is provided directly, the address of the object
            // is assigned in TBranchElement::fObject (via TObject::SetObject())
            // and the address itself is set to the address of the fObject
            // member. Here we check that the address of the object in fObject
            // is the same as the address of our current data type
            pBranch->SetObject(&data);
            MF_LOG_DEBUG("AnalysisTreeStructure")
              << "Reassigning object to branch '" << name << "'";
          }
          else {
            MF_LOG_DEBUG("AnalysisTreeStructure")
              << "Branch '" << name << "' is fine";
          }
        } // operator()
      //@}
    }; // class BranchCreator

  }; // class AnalysisTreeDataStruct
  
  
  /**
   * @brief Creates a simple ROOT tree with tracking and calorimetry information
   * 
   * <h2>Configuration parameters</h2>
   * - <b>UseBuffers</b> (default: false): if enabled, memory is allocated for
   *   tree data for all the run; otherwise, it's allocated on each event, used
   *   and freed; use "true" for speed, "false" to save memory
   * - <b>SaveAuxDetInfo</b> (default: false): if enabled, auxiliary detector
   *   data will be extracted and included in the tree
   */
  class AnalysisTree : public art::EDAnalyzer {

  public:

    explicit AnalysisTree(fhicl::ParameterSet const& pset);
    virtual ~AnalysisTree();

    /// read access to event
    void analyze(const art::Event& evt);
  //  void beginJob() {}
    void beginSubRun(const art::SubRun& sr);

  private:

    void   HitsPurity(detinfo::DetectorClocksData const& clockData,
                      std::vector< art::Ptr<recob::Hit> > const& hits, Int_t& trackid, Float_t& purity, double& maxe);
    double length(const recob::Track& track);
    double length(const simb::MCParticle& part, TVector3& start, TVector3& end);
    double bdist(const recob::Track::Point_t& pos);

    TTree* fTree;

    // event information is huge and dynamic;
    // run information is much smaller and we still store it statically
    // in the event
    AnalysisTreeDataStruct* fData;
//    AnalysisTreeDataStruct::RunData_t RunData;
    AnalysisTreeDataStruct::SubRunData_t SubRunData;

    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fCalDataModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fCryGenModuleLabel;
    std::string fG4ModuleLabel;
    std::string fPFParticleModuleLabel;
    std::string fShowerModuleLabel;
    std::vector<std::string> fVertexModuleLabel;
    std::vector<std::string> fTrackModuleLabel;
    std::vector<std::string> fCalorimetryModuleLabel;
    std::vector<std::string> fParticleIDModuleLabel;
    std::string fPOTModuleLabel;
    bool fUseBuffer; ///< whether to use a permanent buffer (faster, huge memory)    
    bool fSaveAuxDetInfo; ///< whether to extract and save auxiliary detector data
    bool fSaveCryInfo; ///whether to extract and save CRY particle data
    bool fSaveGenieInfo; ///whether to extract and save Genie information
    bool fSaveGeantInfo; ///whether to extract and save Geant information
    bool fSaveHitInfo; ///whether to extract and save Hit information
    bool fSaveTrackInfo; ///whether to extract and save Track information
    bool fSaveShowerInfo; ///whether to extract and save Shower information
    bool fSaveVertexInfo; ///whether to extract and save Vertex information
    std::vector<bool> fSaveHierarchyInfo;  ///< if the user wants to access the tracks with their hierarchy for each tracker
    bool fSaveShowerHierarchyInfo;  ///< if the user wants to access the showers with their hierarchy
    
    std::vector<std::string> fCosmicTaggerAssocLabel;
    std::vector<std::string> fFlashMatchAssocLabel;

    bool isCosmics;        ///< if it contains cosmics
    bool fSaveCaloCosmics; ///< save calorimetry information for cosmics
    float fG4minE;         ///< Energy threshold to save g4 particle info
    
    /// Returns the number of trackers configured
    size_t GetNTrackers() const { return fTrackModuleLabel.size(); }
    
    /// Creates the structure for the tree data; optionally initializes it
    void CreateData(bool bClearData = false)
      {
        if (!fData) {
          fData = new AnalysisTreeDataStruct(GetNTrackers());
          fData->SetBits(AnalysisTreeDataStruct::tdAuxDet, !fSaveAuxDetInfo);
          fData->SetBits(AnalysisTreeDataStruct::tdCry, !fSaveCryInfo);	  
          fData->SetBits(AnalysisTreeDataStruct::tdGenie, !fSaveGenieInfo);
          fData->SetBits(AnalysisTreeDataStruct::tdGeant, !fSaveGeantInfo); 
        }
        else {
          fData->SetBits(AnalysisTreeDataStruct::tdHit, !fSaveHitInfo);	
          fData->SetBits(AnalysisTreeDataStruct::tdTrack, !fSaveTrackInfo);	
          fData->SetBits(AnalysisTreeDataStruct::tdShower, !fSaveShowerInfo);	
          fData->SetBits(AnalysisTreeDataStruct::tdVtx, !fSaveVertexInfo);	  	  	    	  	    	  	    	  
          fData->SetTrackers(GetNTrackers());
          fData->SetVertices(GetNTrackers());
          if (bClearData) fData->Clear();
        }
      } // CreateData()
   
    /// Sets the addresses of all the tree branches, creating the missing ones
    void SetAddresses()
      {
        CheckData("SetAddress()"); CheckTree("SetAddress()");
        fData->SetAddresses(fTree, fTrackModuleLabel, fShowerModuleLabel, fVertexModuleLabel, isCosmics, fSaveHierarchyInfo, fSaveShowerHierarchyInfo);
      } // SetAddresses()
    
    /// Sets the addresses of all the tree branches of the specified tracking algo,
    /// creating the missing ones
    void SetTrackerAddresses(size_t iTracker)
      {
        CheckData("SetTrackerAddresses()"); CheckTree("SetTrackerAddresses()");
        if (iTracker >= fData->GetNTrackers()) {
          throw art::Exception(art::errors::LogicError)
            << "AnalysisTree::SetTrackerAddresses(): no tracker #" << iTracker
            << " (" << fData->GetNTrackers() << " available)";
        }
        fData->GetTrackerData(iTracker) \
          .SetAddresses(fTree, fTrackModuleLabel[iTracker], isCosmics, fSaveHierarchyInfo[iTracker]);
      } // SetTrackerAddresses()
    
    void SetVerticesAddresses(size_t iTracker)
      {
        CheckData("SetVerticesAddresses()"); CheckTree("SetVerticesAddresses()");
        if (iTracker >= fData->GetNTrackers()) {
          throw art::Exception(art::errors::LogicError)
            << "AnalysisTree::SetVerticesAddresses(): no tracker #" << iTracker
            << " (" << fData->GetNTrackers() << " available)";
        }
        fData->GetVertexData(iTracker) \
          .SetAddresses(fTree, fVertexModuleLabel[iTracker], fSaveHierarchyInfo[iTracker]);
      } // SetVerticesAddresses()
    
    /// Sets the addresses of all the tree branches, creating the missing ones
    void SetShowerAddresses()
      {
        CheckData("SetShowerAddress()"); CheckTree("SetShowerAddress()");
        fData->ShowerData.SetShowerAddresses(fTree, fShowerModuleLabel, fSaveShowerHierarchyInfo);
      } // SetShowerAddresses()
    
    /// Create the output tree and the data structures, if needed
    void CreateTree(bool bClearData = false);
    
    /// Destroy the local buffers (existing branches will point to invalid address!)
    void DestroyData() { if (fData) { delete fData; fData = nullptr; } }
    
    /// Helper function: throws if no data structure is available
    void CheckData(std::string caller) const
      {
        if (fData) return;
        throw art::Exception(art::errors::LogicError)
          << "AnalysisTree::" << caller << ": no data";
      } // CheckData()
    /// Helper function: throws if no tree is available
    void CheckTree(std::string caller) const
      {
        if (fTree) return;
        throw art::Exception(art::errors::LogicError)
          << "AnalysisTree::" << caller << ": no tree";
      } // CheckData()
  }; // class sbnd::AnalysisTree
} // namespace sbnd


namespace { // local namespace
  /// Simple stringstream which empties its buffer on operator() call
  class AutoResettingStringSteam: public std::ostringstream {
      public:
    AutoResettingStringSteam& operator() () { str(""); return *this; }
  }; // class AutoResettingStringSteam

  /// Fills a sequence of TYPE elements
  template <typename ITER, typename TYPE>
  inline void FillWith(ITER from, ITER to, TYPE value)
    { std::fill(from, to, value); }

  /// Fills a sequence of TYPE elements
  template <typename ITER, typename TYPE>
  inline void FillWith(ITER from, size_t n, TYPE value)
    { std::fill(from, from + n, value); }

  /// Fills a container with begin()/end() interface
  template <typename CONT, typename V>
  inline void FillWith(CONT& data, const V& value)
    { FillWith(std::begin(data), std::end(data), value); }

} // local namespace


//------------------------------------------------------------------------------
//---  AnalysisTreeDataStruct::TrackDataStruct
//---

void sbnd::AnalysisTreeDataStruct::TrackDataStruct::Resize(size_t nTracks)
{
  MaxTracks = nTracks;
  
  trkId.resize(MaxTracks);
  trkncosmictags_tagger.resize(MaxTracks);
  trkcosmicscore_tagger.resize(MaxTracks);
  trkcosmictype_tagger.resize(MaxTracks);
  trkncosmictags_flashmatch.resize(MaxTracks);
  trkcosmicscore_flashmatch.resize(MaxTracks);
  trkcosmictype_flashmatch.resize(MaxTracks);
  trkstartx.resize(MaxTracks);
  trkstarty.resize(MaxTracks);
  trkstartz.resize(MaxTracks);
  trkstartd.resize(MaxTracks);
  trkendx.resize(MaxTracks);
  trkendy.resize(MaxTracks);
  trkendz.resize(MaxTracks);
  trkendd.resize(MaxTracks);
  trktheta.resize(MaxTracks);
  trkphi.resize(MaxTracks);
  trkstartdcosx.resize(MaxTracks);
  trkstartdcosy.resize(MaxTracks);
  trkstartdcosz.resize(MaxTracks);
  trkenddcosx.resize(MaxTracks);
  trkenddcosy.resize(MaxTracks);
  trkenddcosz.resize(MaxTracks);
  trkthetaxz.resize(MaxTracks);
  trkthetayz.resize(MaxTracks);
  trkmom.resize(MaxTracks);
  trkmomrange.resize(MaxTracks);
  trkmommschi2.resize(MaxTracks);
  trkmommsllhd.resize(MaxTracks);  
  trklen.resize(MaxTracks);
  trksvtxid.resize(MaxTracks);
  trkevtxid.resize(MaxTracks);
  // PID variables
  trkpidpdg.resize(MaxTracks);
  trkpidchi.resize(MaxTracks);
  trkpidchipr.resize(MaxTracks);
  trkpidchika.resize(MaxTracks);
  trkpidchipi.resize(MaxTracks);
  trkpidchimu.resize(MaxTracks);
  trkpidpida.resize(MaxTracks);
  trkpidbestplane.resize(MaxTracks);
  
  trkke.resize(MaxTracks);
  trkrange.resize(MaxTracks);
  trkidtruth_recoutils_totaltrueenergy.resize(MaxTracks);
  trkidtruth_recoutils_totalrecocharge.resize(MaxTracks);
  trkidtruth_recoutils_totalrecohits.resize(MaxTracks);
  trkidtruth.resize(MaxTracks);
  trkorigin.resize(MaxTracks);
  trkpdgtruth.resize(MaxTracks);
  trkefftruth.resize(MaxTracks);
  trksimIDEenergytruth.resize(MaxTracks);
  trksimIDExtruth.resize(MaxTracks);
  trksimIDEytruth.resize(MaxTracks);
  trksimIDEztruth.resize(MaxTracks);
  trkpurtruth.resize(MaxTracks);
  trkpitchc.resize(MaxTracks);
  ntrkhits.resize(MaxTracks);
  
  trkisprimary.resize(MaxTracks);  
  trkndaughters.resize(MaxTracks); 
  trkpfpid.resize(MaxTracks);      
  trkparentpfpid.resize(MaxTracks);

  trkdedx.resize(MaxTracks);
  trkdqdx.resize(MaxTracks);
  trkresrg.resize(MaxTracks);
  trkxyz.resize(MaxTracks);
  
} // sbnd::AnalysisTreeDataStruct::TrackDataStruct::Resize()

void sbnd::AnalysisTreeDataStruct::TrackDataStruct::Clear() {
  Resize(MaxTracks);
  ntracks = 0;
  
  FillWith(trkId        , -9999  );
  FillWith(trkncosmictags_tagger, -9999  );
  FillWith(trkcosmicscore_tagger, -99999.);
  FillWith(trkcosmictype_tagger, -9999  );
  FillWith(trkncosmictags_flashmatch, -9999  );
  FillWith(trkcosmicscore_flashmatch, -99999.);
  FillWith(trkcosmictype_flashmatch, -9999  );
  FillWith(trkstartx    ,   -99999.);
  FillWith(trkstarty    ,   -99999.);
  FillWith(trkstartz    ,   -99999.);
  FillWith(trkstartd    ,   -99999.);
  FillWith(trkendx      ,   -99999.);
  FillWith(trkendy      ,   -99999.);
  FillWith(trkendz      ,   -99999.);
  FillWith(trkendd      ,   -99999.);
  FillWith(trktheta     ,   -99999.);
  FillWith(trkphi       ,   -99999.);
  FillWith(trkstartdcosx,   -99999.);
  FillWith(trkstartdcosy,   -99999.);
  FillWith(trkstartdcosz,   -99999.);
  FillWith(trkenddcosx  ,   -99999.);
  FillWith(trkenddcosy  ,   -99999.);
  FillWith(trkenddcosz  ,   -99999.);
  FillWith(trkthetaxz   ,   -99999.);
  FillWith(trkthetayz   ,   -99999.);
  FillWith(trkmom       ,   -99999.);
  FillWith(trkmomrange  ,   -99999.);  
  FillWith(trkmommschi2 ,   -99999.);  
  FillWith(trkmommsllhd ,   -99999.);  
  FillWith(trklen       ,   -99999.);
  FillWith(trksvtxid    ,   -1);
  FillWith(trkevtxid    ,   -1);
  FillWith(trkpidbestplane, -1); 
  FillWith(trkisprimary,   -9999);  
  FillWith(trkndaughters,  -9999);  
  FillWith(trkpfpid,       -9999);  
  FillWith(trkparentpfpid, -9999);  

 
  for (size_t iTrk = 0; iTrk < MaxTracks; ++iTrk){
    
    // the following are BoxedArray's;
    // their iterators traverse all the array dimensions
    FillWith(trkke[iTrk]      , -99999.);
    FillWith(trkrange[iTrk]   , -99999.);
    FillWith(trkidtruth_recoutils_totaltrueenergy[iTrk] , -99999 );
    FillWith(trkidtruth_recoutils_totalrecocharge[iTrk] , -99999 );
    FillWith(trkidtruth_recoutils_totalrecohits[iTrk] , -99999 );
    FillWith(trkidtruth[iTrk] , -99999 );
    FillWith(trkorigin[iTrk]  , -1 );
    FillWith(trkpdgtruth[iTrk], -99999 );
    FillWith(trkefftruth[iTrk], -99999.);
    FillWith(trksimIDEenergytruth[iTrk], -99999.);
    FillWith(trksimIDExtruth[iTrk], -99999.);
    FillWith(trksimIDEytruth[iTrk], -99999.);
    FillWith(trksimIDEztruth[iTrk], -99999.);
    FillWith(trkpurtruth[iTrk], -99999.);
    FillWith(trkpitchc[iTrk]  , -99999.);
    FillWith(ntrkhits[iTrk]   ,  -9999 );
    
    FillWith(trkdedx[iTrk], 0.);
    FillWith(trkdqdx[iTrk], 0.);
    FillWith(trkresrg[iTrk], 0.);
    
    FillWith(trkxyz[iTrk], 0.);
 
    FillWith(trkpidpdg[iTrk]    , -1);
    FillWith(trkpidchi[iTrk]    , -99999.);
    FillWith(trkpidchipr[iTrk]  , -99999.);
    FillWith(trkpidchika[iTrk]  , -99999.);
    FillWith(trkpidchipi[iTrk]  , -99999.);
    FillWith(trkpidchimu[iTrk]  , -99999.);
    FillWith(trkpidpida[iTrk]   , -99999.);
    
  } // for track
  
} // sbnd::AnalysisTreeDataStruct::TrackDataStruct::Clear()


void sbnd::AnalysisTreeDataStruct::TrackDataStruct::SetAddresses(
  TTree* pTree, std::string tracker, bool isCosmics, bool saveHierarchyInfo
) {
  if (MaxTracks == 0) return; // no tracks, no tree!
  
  sbnd::AnalysisTreeDataStruct::BranchCreator CreateBranch(pTree);

  AutoResettingStringSteam sstr;
  sstr() << kMaxTrackHits;
  std::string MaxTrackHitsIndexStr("[" + sstr.str() + "]");
  
  std::string TrackLabel = tracker;
  std::string BranchName;

  BranchName = "ntracks_" + TrackLabel;
  CreateBranch(BranchName, &ntracks, BranchName + "/S");
  std::string NTracksIndexStr = "[" + BranchName + "]";
  
  BranchName = "trkId_" + TrackLabel;
  CreateBranch(BranchName, trkId, BranchName + NTracksIndexStr + "/S");
  
  BranchName = "trkncosmictags_tagger_" + TrackLabel;
  CreateBranch(BranchName, trkncosmictags_tagger, BranchName + NTracksIndexStr + "/S");
  
  BranchName = "trkcosmicscore_tagger_" + TrackLabel;
  CreateBranch(BranchName, trkcosmicscore_tagger, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkcosmictype_tagger_" + TrackLabel;
  CreateBranch(BranchName, trkcosmictype_tagger, BranchName + NTracksIndexStr + "/S");

  BranchName = "trkncosmictags_flashmatch_" + TrackLabel;
  CreateBranch(BranchName, trkncosmictags_flashmatch, BranchName + NTracksIndexStr + "/S");
  
  BranchName = "trkcosmicscore_flashmatch_" + TrackLabel;
  CreateBranch(BranchName, trkcosmicscore_flashmatch, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkcosmictype_flashmatch_" + TrackLabel;
  CreateBranch(BranchName, trkcosmictype_flashmatch, BranchName + NTracksIndexStr + "/S");
  
  BranchName = "trkke_" + TrackLabel;
  CreateBranch(BranchName, trkke, BranchName + NTracksIndexStr + "[3]/F");
  
  BranchName = "trkrange_" + TrackLabel;
  CreateBranch(BranchName, trkrange, BranchName + NTracksIndexStr + "[3]/F");
   
  BranchName = "trkidtruth_recoutils_totaltrueenergy_" + TrackLabel;
  CreateBranch(BranchName, trkidtruth_recoutils_totaltrueenergy, BranchName + NTracksIndexStr + "[3]/I");

  BranchName = "trkidtruth_recoutils_totalrecocharge_" + TrackLabel;
  CreateBranch(BranchName, trkidtruth_recoutils_totalrecocharge, BranchName + NTracksIndexStr + "[3]/I");

  BranchName = "trkidtruth_recoutils_totalrecohits_" + TrackLabel;
  CreateBranch(BranchName, trkidtruth_recoutils_totalrecohits, BranchName + NTracksIndexStr + "[3]/I");

  BranchName = "trkidtruth_" + TrackLabel;
  CreateBranch(BranchName, trkidtruth, BranchName + NTracksIndexStr + "[3]/I");

  BranchName = "trkorigin_" + TrackLabel;
  CreateBranch(BranchName, trkorigin, BranchName + NTracksIndexStr + "[3]/S");
  
  BranchName = "trkpdgtruth_" + TrackLabel;
  CreateBranch(BranchName, trkpdgtruth, BranchName + NTracksIndexStr + "[3]/I");
  
  BranchName = "trkefftruth_" + TrackLabel;
  CreateBranch(BranchName, trkefftruth, BranchName + NTracksIndexStr + "[3]/F");
 
  BranchName = "trksimIDEenergytruth_" + TrackLabel;
  CreateBranch(BranchName, trksimIDEenergytruth, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trksimIDExtruth_" + TrackLabel;
  CreateBranch(BranchName, trksimIDExtruth, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trksimIDEytruth_" + TrackLabel;
  CreateBranch(BranchName, trksimIDEytruth, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trksimIDEztruth_" + TrackLabel;
  CreateBranch(BranchName, trksimIDEztruth, BranchName + NTracksIndexStr + "[3]/F");
 
  BranchName = "trkpurtruth_" + TrackLabel;
  CreateBranch(BranchName, trkpurtruth, BranchName + NTracksIndexStr + "[3]/F");
  
  BranchName = "trkpitchc_" + TrackLabel;
  CreateBranch(BranchName, trkpitchc, BranchName + NTracksIndexStr + "[3]/F");
  
  BranchName = "ntrkhits_" + TrackLabel;
  CreateBranch(BranchName, ntrkhits, BranchName + NTracksIndexStr + "[3]/S");
  
  if (!isCosmics){
    BranchName = "trkdedx_" + TrackLabel;
    CreateBranch(BranchName, trkdedx, BranchName + NTracksIndexStr + "[3]" + MaxTrackHitsIndexStr + "/F");
  
    BranchName = "trkdqdx_" + TrackLabel;
    CreateBranch(BranchName, trkdqdx, BranchName + NTracksIndexStr + "[3]" + MaxTrackHitsIndexStr + "/F");
    
    BranchName = "trkresrg_" + TrackLabel;
    CreateBranch(BranchName, trkresrg, BranchName + NTracksIndexStr + "[3]" + MaxTrackHitsIndexStr + "/F");
    
    BranchName = "trkxyz_" + TrackLabel;
    CreateBranch(BranchName, trkxyz, BranchName + NTracksIndexStr + "[3]" + MaxTrackHitsIndexStr + "/F");
  }

  BranchName = "trkstartx_" + TrackLabel;
  CreateBranch(BranchName, trkstartx, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkstarty_" + TrackLabel;
  CreateBranch(BranchName, trkstarty, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkstartz_" + TrackLabel;
  CreateBranch(BranchName, trkstartz, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkstartd_" + TrackLabel;
  CreateBranch(BranchName, trkstartd, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkendx_" + TrackLabel;
  CreateBranch(BranchName, trkendx, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkendy_" + TrackLabel;
  CreateBranch(BranchName, trkendy, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkendz_" + TrackLabel;
  CreateBranch(BranchName, trkendz, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkendd_" + TrackLabel;
  CreateBranch(BranchName, trkendd, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trktheta_" + TrackLabel;
  CreateBranch(BranchName, trktheta, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkphi_" + TrackLabel;
  CreateBranch(BranchName, trkphi, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkstartdcosx_" + TrackLabel;
  CreateBranch(BranchName, trkstartdcosx, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkstartdcosy_" + TrackLabel;
  CreateBranch(BranchName, trkstartdcosy, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkstartdcosz_" + TrackLabel;
  CreateBranch(BranchName, trkstartdcosz, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkenddcosx_" + TrackLabel;
  CreateBranch(BranchName, trkenddcosx, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkenddcosy_" + TrackLabel;
  CreateBranch(BranchName, trkenddcosy, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkenddcosz_" + TrackLabel;
  CreateBranch(BranchName, trkenddcosz, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkthetaxz_" + TrackLabel;
  CreateBranch(BranchName, trkthetaxz, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkthetayz_" + TrackLabel;
  CreateBranch(BranchName, trkthetayz, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkmom_" + TrackLabel;
  CreateBranch(BranchName, trkmom, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trkmomrange_" + TrackLabel;
  CreateBranch(BranchName, trkmomrange, BranchName + NTracksIndexStr + "/F");

  BranchName = "trkmommschi2_" + TrackLabel;
  CreateBranch(BranchName, trkmommschi2, BranchName + NTracksIndexStr + "/F");

  BranchName = "trkmommsllhd_" + TrackLabel;
  CreateBranch(BranchName, trkmommsllhd, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trklen_" + TrackLabel;
  CreateBranch(BranchName, trklen, BranchName + NTracksIndexStr + "/F");
  
  BranchName = "trksvtxid_" + TrackLabel;
  CreateBranch(BranchName, trksvtxid, BranchName + NTracksIndexStr + "/S");
  
  BranchName = "trkevtxid_" + TrackLabel;
  CreateBranch(BranchName, trkevtxid, BranchName + NTracksIndexStr + "/S");

  BranchName = "trkpidpdg_" + TrackLabel;
  CreateBranch(BranchName, trkpidpdg, BranchName + NTracksIndexStr + "[3]/I");

  BranchName = "trkpidchi_" + TrackLabel;
  CreateBranch(BranchName, trkpidchi, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trkpidchipr_" + TrackLabel;
  CreateBranch(BranchName, trkpidchipr, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trkpidchika_" + TrackLabel;
  CreateBranch(BranchName, trkpidchika, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trkpidchipi_" + TrackLabel;
  CreateBranch(BranchName, trkpidchipi, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trkpidchimu_" + TrackLabel;
  CreateBranch(BranchName, trkpidchimu, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trkpidpida_" + TrackLabel;
  CreateBranch(BranchName, trkpidpida, BranchName + NTracksIndexStr + "[3]/F");

  BranchName = "trkpidbestplane_" + TrackLabel;
  CreateBranch(BranchName, trkpidbestplane, BranchName + NTracksIndexStr + "/S");

  if(saveHierarchyInfo){
    BranchName = "trkisprimary_" + TrackLabel;
    CreateBranch(BranchName, trkisprimary, BranchName + NTracksIndexStr + "/O");
  
    BranchName = "trkndaughters_" + TrackLabel;
    CreateBranch(BranchName, trkndaughters, BranchName + NTracksIndexStr + "/I");
  
    BranchName = "trkpfpid_" + TrackLabel;
    CreateBranch(BranchName, trkpfpid, BranchName + NTracksIndexStr + "/I");
  
    BranchName = "trkparentpfpid_" + TrackLabel;
    CreateBranch(BranchName, trkparentpfpid, BranchName + NTracksIndexStr + "/I");
  }
} // sbnd::AnalysisTreeDataStruct::TrackDataStruct::SetAddresses()

//------------------------------------------------------------------------------
//---  AnalysisTreeDataStruct::ShowerDataStruct
//---

void sbnd::AnalysisTreeDataStruct::ShowerDataStruct::Resize(size_t nShowers)
{
  MaxShowers = nShowers;
  
  shwId.resize(MaxShowers);  
  shwbestplane.resize(MaxShowers);
  shwlength.resize(MaxShowers);   
  shwopenangle.resize(MaxShowers);
  shwstartx.resize(MaxShowers);   
  shwstarty.resize(MaxShowers);   
  shwstartz.resize(MaxShowers);   
  shwdirx.resize(MaxShowers);     
  shwdiry.resize(MaxShowers);     
  shwdirz.resize(MaxShowers);     
  shwisprimary.resize(MaxShowers);
  shwndaughters.resize(MaxShowers);
  shwpfpid.resize(MaxShowers);
  shwparentpfpid.resize(MaxShowers);

} // AnalysisTreeDataStruct::ShowerDataStruct.Resize()

void sbnd::AnalysisTreeDataStruct::ShowerDataStruct::Clear() {
  Resize(MaxShowers);
  nshowers = 0;

  // For variables with the ShowerData_t form
  FillWith(shwId,          -9999);
  FillWith(shwbestplane,   -9999);
  FillWith(shwlength,      -9999.);
  FillWith(shwopenangle,   -9999.);
  FillWith(shwstartx,      -9999.);
  FillWith(shwstarty,      -9999.);
  FillWith(shwstartz,      -9999.);
  FillWith(shwdirx,        -9999.);
  FillWith(shwdiry,        -9999.);
  FillWith(shwdirz,        -9999.);
  FillWith(shwisprimary,   -9999);
  FillWith(shwndaughters,  -9999);
  FillWith(shwpfpid,       -9999);
  FillWith(shwparentpfpid, -9999);
  
} // sbnd::AnalysisTreeDataStruct::ShowerDataStruct::Clear()

void sbnd::AnalysisTreeDataStruct::ShowerDataStruct::SetShowerAddresses(
  TTree* pTree, std::string showerLabel, bool saveHierarchyInfo
) {
  if (MaxShowers == 0) return; // no tracks, no tree!
  
  sbnd::AnalysisTreeDataStruct::BranchCreator CreateBranch(pTree);

  std::string ShowerLabel = showerLabel;
  std::string BranchName;

  BranchName = "nshowers_" + ShowerLabel;
  CreateBranch(BranchName, &nshowers, BranchName + "/S");
  std::string NShowersIndexStr = "[" + BranchName + "]";
  
  BranchName = "shwId_" + ShowerLabel;
  CreateBranch(BranchName, shwId, BranchName + NShowersIndexStr + "/I");
  
  BranchName = "shwbestplane_" + ShowerLabel;
  CreateBranch(BranchName, shwbestplane, BranchName + NShowersIndexStr + "/I");
  
  BranchName = "shwlength_" + ShowerLabel;
  CreateBranch(BranchName, shwlength, BranchName + NShowersIndexStr + "/F");
  
  BranchName = "shwlength_" + ShowerLabel;
  CreateBranch(BranchName, shwlength, BranchName + NShowersIndexStr + "/F");
  
  BranchName = "shwopenangle_" + ShowerLabel;
  CreateBranch(BranchName, shwopenangle, BranchName + NShowersIndexStr + "/F");
  
  BranchName = "shwstartx_" + ShowerLabel;
  CreateBranch(BranchName, shwstartx, BranchName + NShowersIndexStr + "/F");
  
  BranchName = "shwstarty_" + ShowerLabel;
  CreateBranch(BranchName, shwstarty, BranchName + NShowersIndexStr + "/F");
  
  BranchName = "shwstartz_" + ShowerLabel;
  CreateBranch(BranchName, shwstartz, BranchName + NShowersIndexStr + "/F");
  
  BranchName = "shwdirx_" + ShowerLabel;
  CreateBranch(BranchName, shwdirx, BranchName + NShowersIndexStr + "/F");
  
  BranchName = "shwdiry_" + ShowerLabel;
  CreateBranch(BranchName, shwdiry, BranchName + NShowersIndexStr + "/F");
  
  BranchName = "shwdirz_" + ShowerLabel;
  CreateBranch(BranchName, shwdirz, BranchName + NShowersIndexStr + "/F");
  
  if(saveHierarchyInfo){
    BranchName = "shwisprimary_" + ShowerLabel;
    CreateBranch(BranchName, shwisprimary, BranchName + NShowersIndexStr + "/O");
  
    BranchName = "shwndaughters_" + ShowerLabel;
    CreateBranch(BranchName, shwndaughters, BranchName + NShowersIndexStr + "/I");
  
    BranchName = "shwpfpid_" + ShowerLabel;
    CreateBranch(BranchName, shwpfpid, BranchName + NShowersIndexStr + "/I");
  
    BranchName = "shwparentpfpid_" + ShowerLabel;
    CreateBranch(BranchName, shwparentpfpid, BranchName + NShowersIndexStr + "/I");
  }
} // sbnd::AnalysisTreeDataStruct::ShowerDataStruct::SetShowerAddresses()

//------------------------------------------------------------------------------
//---  AnalysisTreeDataStruct::VertexDataStruct
//---

void sbnd::AnalysisTreeDataStruct::VertexDataStruct::Resize(size_t nVertices)
{
  MaxVertices = nVertices;
 
  vtx.resize(MaxVertices);
  primaryvtx.resize(MaxVertices);
} // AnalysisTreeDataStruct::VertexDataStruct::Resize()

void sbnd::AnalysisTreeDataStruct::VertexDataStruct::Clear() {
  Resize(MaxVertices);
  nvtx = 0;

  for (size_t iVtx = 0; iVtx < MaxVertices; ++iVtx){
    FillWith(vtx[iVtx], -9999.);
  }
  FillWith(primaryvtx, -9999);  
} // AnalysisTreeDataStruct::VertexDataStruct::Clear()

void sbnd::AnalysisTreeDataStruct::VertexDataStruct::SetAddresses(
  TTree* pTree, std::string vertices, bool saveHierarchyInfo
) {
  if (MaxVertices == 0) return; // no vertices, no tree!
  
  sbnd::AnalysisTreeDataStruct::BranchCreator CreateBranch(pTree);

  std::string VertexLabel = vertices;
  std::string BranchName;

  BranchName = "nvtx_" + VertexLabel;
  CreateBranch(BranchName, &nvtx, BranchName + "/S");
  std::string NVertexIndexStr = "[" + BranchName + "]";
 
  if(saveHierarchyInfo){
    BranchName = "primaryvtx_" + VertexLabel;
    CreateBranch(BranchName, primaryvtx, BranchName + NVertexIndexStr + "/I");
  }

  BranchName = "vtx_" + VertexLabel;
  CreateBranch(BranchName, vtx, BranchName + NVertexIndexStr + "[3]/F");
} // AnalysisTreeDataStruct::VertexDataStruct::SetAddresses()

//------------------------------------------------------------------------------
//---  AnalysisTreeDataStruct
//---

void sbnd::AnalysisTreeDataStruct::ClearLocalData() {

//  RunData.Clear();
  SubRunData.Clear();

  run = -99999;
  subrun = -99999;
  event = -99999;
  evttime = -99999;
  beamtime = -99999;
  isdata = -99;
  //taulife = -99999;

  no_hits = 0;
 
  std::fill(hit_tpc, hit_tpc + sizeof(hit_tpc)/sizeof(hit_tpc[0]), -9999);
  std::fill(hit_plane, hit_plane + sizeof(hit_plane)/sizeof(hit_plane[0]), -9999);
  std::fill(hit_wire, hit_wire + sizeof(hit_wire)/sizeof(hit_wire[0]), -9999);
  std::fill(hit_channel, hit_channel + sizeof(hit_channel)/sizeof(hit_channel[0]), -9999);
  std::fill(hit_peakT, hit_peakT + sizeof(hit_peakT)/sizeof(hit_peakT[0]), -99999.);
  std::fill(hit_charge, hit_charge + sizeof(hit_charge)/sizeof(hit_charge[0]), -99999.);
  std::fill(hit_ph, hit_ph + sizeof(hit_ph)/sizeof(hit_ph[0]), -99999.);
  std::fill(hit_startT, hit_startT + sizeof(hit_startT)/sizeof(hit_startT[0]), -99999.);
  std::fill(hit_endT, hit_endT + sizeof(hit_endT)/sizeof(hit_endT[0]), -99999.);
  std::fill(hit_trkid, hit_trkid + sizeof(hit_trkid)/sizeof(hit_trkid[0]), -9999);
  std::fill(hit_nelec, hit_nelec + sizeof(hit_nelec)/sizeof(hit_nelec[0]), -99999.);
  std::fill(hit_energy, hit_energy + sizeof(hit_energy)/sizeof(hit_energy[0]), -99999.);

  /*
  nvtx = 0;
  for (size_t ivtx = 0; ivtx < kMaxVertices; ++ivtx) {
    std::fill(vtx[ivtx], vtx[ivtx]+3, -99999.);
  }
  FillWith(primaryvtx, -9999);
*/

  mcevts_truth = 0;
  mcevts_truthcry = -99999;
  FillWith(nuScatterCode_truth,-99999);
  FillWith(nuID_truth,-99999);
  FillWith(nuPDG_truth,-99999);
  FillWith(ccnc_truth,-99999);
  FillWith(mode_truth,-99999);
  FillWith(enu_truth,-99999);
  FillWith(Q2_truth,-99999);
  FillWith(W_truth,-99999);
  FillWith(hitnuc_truth,-99999);
  FillWith(nuvtxx_truth,-99999);
  FillWith(nuvtxy_truth,-99999);
  FillWith(nuvtxz_truth,-99999);
  FillWith(nu_dcosx_truth,-99999);
  FillWith(nu_dcosy_truth,-99999);
  FillWith(nu_dcosz_truth,-99999);
  FillWith(lep_mom_truth,-99999);
  FillWith(lep_dcosx_truth,-99999);
  FillWith(lep_dcosy_truth,-99999);
  FillWith(lep_dcosz_truth,-99999);
  FillWith(tpx_flux,-99999);
  FillWith(tpy_flux,-99999);
  FillWith(tpz_flux,-99999);
  FillWith(tptype_flux,-99999);

  genie_no_primaries = 0;
  cry_no_primaries = 0;
  no_primaries = 0;
  geant_list_size=0;
  geant_list_size_in_tpcAV = 0;
  
  FillWith(pdg, -99999);
  FillWith(status, -99999);
  FillWith(Mass, -99999.);
  FillWith(Eng, -99999.);
  FillWith(EndE, -99999.);
  FillWith(Px, -99999.);
  FillWith(Py, -99999.);
  FillWith(Pz, -99999.);
  FillWith(P, -99999.);
  FillWith(StartPointx, -99999.);
  FillWith(StartPointy, -99999.);
  FillWith(StartPointz, -99999.);
  FillWith(StartT, -99999.);
  FillWith(EndT, -99999.);    
  FillWith(EndPointx, -99999.);
  FillWith(EndPointy, -99999.);
  FillWith(EndPointz, -99999.);
  FillWith(EndT, -99999.);
  FillWith(theta, -99999.);
  FillWith(phi, -99999.);
  FillWith(theta_xz, -99999.);
  FillWith(theta_yz, -99999.);
  FillWith(pathlen, -99999.);
  FillWith(inTPCActive, -99999);
  FillWith(StartPointx_tpcAV, -99999.);
  FillWith(StartPointy_tpcAV, -99999.);
  FillWith(StartPointz_tpcAV, -99999.);
  FillWith(EndPointx_tpcAV, -99999.);
  FillWith(EndPointy_tpcAV, -99999.);
  FillWith(EndPointz_tpcAV, -99999.);  
  FillWith(NumberDaughters, -99999);
  FillWith(Mother, -99999);
  FillWith(TrackId, -99999);
  FillWith(process_primary, -99999);
  FillWith(processname, "noname");
  FillWith(MergedId, -99999);
  FillWith(MotherNuId, -99999);
  FillWith(genie_primaries_pdg, -99999);
  FillWith(genie_Eng, -99999.);
  FillWith(genie_Px, -99999.);
  FillWith(genie_Py, -99999.);
  FillWith(genie_Pz, -99999.);
  FillWith(genie_P, -99999.);
  FillWith(genie_status_code, -99999);
  FillWith(genie_mass, -99999.);
  FillWith(genie_trackID, -99999);
  FillWith(genie_ND, -99999);
  FillWith(genie_mother, -99999);
  FillWith(cry_primaries_pdg, -99999);
  FillWith(cry_Eng, -99999.);
  FillWith(cry_Px, -99999.);
  FillWith(cry_Py, -99999.);
  FillWith(cry_Pz, -99999.);
  FillWith(cry_P, -99999.);
  FillWith(cry_StartPointx, -99999.);
  FillWith(cry_StartPointy, -99999.);
  FillWith(cry_StartPointz, -99999.);  
  FillWith(cry_status_code, -99999);
  FillWith(cry_mass, -99999.);
  FillWith(cry_trackID, -99999);
  FillWith(cry_ND, -99999);
  FillWith(cry_mother, -99999);
  
  // auxiliary detector information;
  FillWith(NAuxDets, 0);
  // - set to -9999 all the values of each of the arrays in AuxDetID;
  //   this auto is BoxedArray<Short_t>
  for (auto& partInfo: AuxDetID) FillWith(partInfo, -9999);
  // - pythonish C++: as the previous line, for each one in a list of containers
  //   of the same type (C++ is not python yet), using pointers to avoid copy;
  for (AuxDetMCData_t<Float_t>* cont: {
   &entryX, &entryY, &entryZ, &entryT,
   &exitX , &exitY , &exitZ, &exitT, &exitPx, &exitPy, &exitPz,
   &CombinedEnergyDep
   })
  {
    // this auto is BoxedArray<Float_t>
    for (auto& partInfo: *cont) FillWith(partInfo, -99999.);
  } // for container
  
} // sbnd::AnalysisTreeDataStruct::ClearLocalData()


void sbnd::AnalysisTreeDataStruct::Clear() {
  ClearLocalData();
  std::for_each
    (TrackData.begin(), TrackData.end(), std::mem_fn(&TrackDataStruct::Clear));
  std::for_each
    (VertexData.begin(), VertexData.end(), std::mem_fn(&VertexDataStruct::Clear));
  std::mem_fn(&ShowerDataStruct::Clear);
} // sbnd::AnalysisTreeDataStruct::Clear()

void sbnd::AnalysisTreeDataStruct::ResizeMCNeutrino(int nNeutrinos){

  //min size is 1, to guarantee an address
  MaxMCNeutrinos = (size_t) std::max(nNeutrinos, 1);
  nuScatterCode_truth.resize(MaxMCNeutrinos);
  nuID_truth.resize(MaxMCNeutrinos);
  nuPDG_truth.resize(MaxMCNeutrinos);
  ccnc_truth.resize(MaxMCNeutrinos);
  mode_truth.resize(MaxMCNeutrinos);
  enu_truth.resize(MaxMCNeutrinos);
  Q2_truth.resize(MaxMCNeutrinos);
  W_truth.resize(MaxMCNeutrinos);
  hitnuc_truth.resize(MaxMCNeutrinos);
  nuvtxx_truth.resize(MaxMCNeutrinos);
  nuvtxy_truth.resize(MaxMCNeutrinos);
  nuvtxz_truth.resize(MaxMCNeutrinos);
  nu_dcosx_truth.resize(MaxMCNeutrinos);
  nu_dcosy_truth.resize(MaxMCNeutrinos);
  nu_dcosz_truth.resize(MaxMCNeutrinos);
  lep_mom_truth.resize(MaxMCNeutrinos);
  lep_dcosx_truth.resize(MaxMCNeutrinos);
  lep_dcosy_truth.resize(MaxMCNeutrinos);
  lep_dcosz_truth.resize(MaxMCNeutrinos);
  //Also resize the flux information here as it's a 1:1 with the MCNeutrino
  tpx_flux.resize(MaxMCNeutrinos);
  tpy_flux.resize(MaxMCNeutrinos);
  tpz_flux.resize(MaxMCNeutrinos);
  tptype_flux.resize(MaxMCNeutrinos);

  return;
} // sbnd::AnalysisTreeDataStruct::ResizeMCNeutrino()

void sbnd::AnalysisTreeDataStruct::ResizeGEANT(int nParticles) {

  // minimum size is 1, so that we always have an address
  MaxGEANTparticles = (size_t) std::max(nParticles, 1);
  
  pdg.resize(MaxGEANTparticles);
  status.resize(MaxGEANTparticles);
  Mass.resize(MaxGEANTparticles);  
  Eng.resize(MaxGEANTparticles);
  EndE.resize(MaxGEANTparticles);
  Px.resize(MaxGEANTparticles);
  Py.resize(MaxGEANTparticles);
  Pz.resize(MaxGEANTparticles);
  P.resize(MaxGEANTparticles);
  StartPointx.resize(MaxGEANTparticles);
  StartPointy.resize(MaxGEANTparticles);
  StartPointz.resize(MaxGEANTparticles);
  StartT.resize(MaxGEANTparticles); 
  EndT.resize(MaxGEANTparticles);    
  EndPointx.resize(MaxGEANTparticles);
  EndPointy.resize(MaxGEANTparticles);
  EndPointz.resize(MaxGEANTparticles);
  EndT.resize(MaxGEANTparticles);  
  theta.resize(MaxGEANTparticles);
  phi.resize(MaxGEANTparticles);
  theta_xz.resize(MaxGEANTparticles);
  theta_yz.resize(MaxGEANTparticles);
  pathlen.resize(MaxGEANTparticles);
  inTPCActive.resize(MaxGEANTparticles);
  StartPointx_tpcAV.resize(MaxGEANTparticles);
  StartPointy_tpcAV.resize(MaxGEANTparticles);
  StartPointz_tpcAV.resize(MaxGEANTparticles);
  EndPointx_tpcAV.resize(MaxGEANTparticles);
  EndPointy_tpcAV.resize(MaxGEANTparticles);
  EndPointz_tpcAV.resize(MaxGEANTparticles);    
  NumberDaughters.resize(MaxGEANTparticles);
  Mother.resize(MaxGEANTparticles);
  TrackId.resize(MaxGEANTparticles);
  process_primary.resize(MaxGEANTparticles);
  processname.resize(MaxGEANTparticles);
  MergedId.resize(MaxGEANTparticles);
  MotherNuId.resize(MaxGEANTparticles);
  
  // auxiliary detector structure
  NAuxDets.resize(MaxGEANTparticles);
  AuxDetID.resize(MaxGEANTparticles);
  entryX.resize(MaxGEANTparticles);
  entryY.resize(MaxGEANTparticles);
  entryZ.resize(MaxGEANTparticles);
  entryT.resize(MaxGEANTparticles);
  exitX.resize(MaxGEANTparticles);
  exitY.resize(MaxGEANTparticles);
  exitZ.resize(MaxGEANTparticles);
  exitT.resize(MaxGEANTparticles);
  exitPx.resize(MaxGEANTparticles);
  exitPy.resize(MaxGEANTparticles);
  exitPz.resize(MaxGEANTparticles);
  CombinedEnergyDep.resize(MaxGEANTparticles);
  
} // sbnd::AnalysisTreeDataStruct::ResizeGEANT()

void sbnd::AnalysisTreeDataStruct::ResizeGenie(int nPrimaries) {
  
  // minimum size is 1, so that we always have an address
  MaxGeniePrimaries = (size_t) std::max(nPrimaries, 1);
  genie_primaries_pdg.resize(MaxGeniePrimaries);
  genie_Eng.resize(MaxGeniePrimaries);
  genie_Px.resize(MaxGeniePrimaries);
  genie_Py.resize(MaxGeniePrimaries);
  genie_Pz.resize(MaxGeniePrimaries);
  genie_P.resize(MaxGeniePrimaries);
  genie_status_code.resize(MaxGeniePrimaries);
  genie_mass.resize(MaxGeniePrimaries);
  genie_trackID.resize(MaxGeniePrimaries);
  genie_ND.resize(MaxGeniePrimaries);
  genie_mother.resize(MaxGeniePrimaries);

} // sbnd::AnalysisTreeDataStruct::ResizeGenie()

// minimum size is 1, so that we always have an address
void sbnd::AnalysisTreeDataStruct::ResizeCry(int nPrimaries) {

  cry_primaries_pdg.resize(nPrimaries);
  cry_Eng.resize(nPrimaries);
  cry_Px.resize(nPrimaries);
  cry_Py.resize(nPrimaries);
  cry_Pz.resize(nPrimaries);
  cry_P.resize(nPrimaries);
  cry_StartPointx.resize(nPrimaries);
  cry_StartPointy.resize(nPrimaries);
  cry_StartPointz.resize(nPrimaries);  
  cry_status_code.resize(nPrimaries);
  cry_mass.resize(nPrimaries);
  cry_trackID.resize(nPrimaries);
  cry_ND.resize(nPrimaries);
  cry_mother.resize(nPrimaries);

} // sbnd::AnalysisTreeDataStruct::ResizeCry()

void sbnd::AnalysisTreeDataStruct::SetAddresses(
  TTree* pTree,
  const std::vector<std::string>& trackers,
  const std::string showerLabel,
  const std::vector<std::string>& vertexLabels,
  bool isCosmics,
  const std::vector<bool>& saveHierarchyInfo,
  bool saveShowerHierarchyInfo
) {
  BranchCreator CreateBranch(pTree);

  CreateBranch("run",&run,"run/I");
  CreateBranch("subrun",&subrun,"subrun/I");
  CreateBranch("event",&event,"event/I");
  CreateBranch("evttime",&evttime,"evttime/D");
  CreateBranch("beamtime",&beamtime,"beamtime/D");
  CreateBranch("pot",&SubRunData.pot,"pot/D");

  CreateBranch("isdata",&isdata,"isdata/B");
  //CreateBranch("taulife",&taulife,"taulife/D");

  if (hasHitInfo()){
    CreateBranch("no_hits",&no_hits,"no_hits/I");
    CreateBranch("hit_tpc",hit_tpc,"hit_tpc[no_hits]/S");
    CreateBranch("hit_plane",hit_plane,"hit_plane[no_hits]/S");
    CreateBranch("hit_wire",hit_wire,"hit_wire[no_hits]/S");
    CreateBranch("hit_channel",hit_channel,"hit_channel[no_hits]/S");
    CreateBranch("hit_peakT",hit_peakT,"hit_peakT[no_hits]/F");
    CreateBranch("hit_charge",hit_charge,"hit_charge[no_hits]/F");
    CreateBranch("hit_ph",hit_ph,"hit_ph[no_hits]/F");
    CreateBranch("hit_startT",hit_startT,"hit_startT[no_hits]/F");
    CreateBranch("hit_endT",hit_endT,"hit_endT[no_hits]/F");
    CreateBranch("hit_trkid",hit_trkid,"hit_trkid[no_hits]/S");
    CreateBranch("hit_nelec",hit_nelec,"hit_nelec[no_hits]/F");
    CreateBranch("hit_energy",hit_energy,"hit_energy[no_hits]/F");
  }

  if (hasVertexInfo()){
    /*
    CreateBranch("nvtx",&nvtx,"nvtx/S");
    CreateBranch("vtx",vtx,"vtx[nvtx][3]/F");
    if(saveHierarchyInfo)
      CreateBranch("primaryvtx",primaryvtx,"primaryvtx[nvtx]/I");
      */

    kNTracker = trackers.size();
    for(int i=0; i<kNTracker; i++){
      std::string VertexLabel = vertexLabels[i];

      // note that if the tracker data has maximum number of tracks 0,
      // nothing is initialized (branches are not even created)
      VertexData[i].SetAddresses(pTree, VertexLabel, saveHierarchyInfo[i]);   
    }
  }  

  if (hasTrackInfo()){
    AutoResettingStringSteam sstr;
    sstr() << kMaxTrackHits;
    std::string MaxTrackHitsIndexStr("[" + sstr.str() + "]");

    kNTracker = trackers.size();
    CreateBranch("kNTracker",&kNTracker,"kNTracker/B");
    for(int i=0; i<kNTracker; i++){
      std::string TrackLabel = trackers[i];
      std::string BranchName;

      // note that if the tracker data has maximum number of tracks 0,
      // nothing is initialized (branches are not even created)
      TrackData[i].SetAddresses(pTree, TrackLabel, isCosmics, saveHierarchyInfo[i]);    
    } // for trackers
  } 
  
  if (hasShowerInfo()){
    std::string ShowerLabel = showerLabel;
    ShowerData.SetShowerAddresses(pTree, ShowerLabel, saveShowerHierarchyInfo);    
  }

  if (hasGenieInfo()){
    CreateBranch("mcevts_truth",&mcevts_truth,"mcevts_truth/I");
    CreateBranch("nuScatterCode_truth",nuScatterCode_truth,"nuScatterCode_truth[mcevts_truth]/I");
    CreateBranch("nuID_truth",nuID_truth,"nuID_truth[mcevts_truth]/I");
    CreateBranch("nuPDG_truth",nuPDG_truth,"nuPDG_truth[mcevts_truth]/I");
    CreateBranch("ccnc_truth",ccnc_truth,"ccnc_truth[mcevts_truth]/I");
    CreateBranch("mode_truth",mode_truth,"mode_truth[mcevts_truth]/I");
    CreateBranch("enu_truth",enu_truth,"enu_truth[mcevts_truth]/F");
    CreateBranch("Q2_truth",Q2_truth,"Q2_truth[mcevts_truth]/F");
    CreateBranch("W_truth",W_truth,"W_truth[mcevts_truth]/F");
    CreateBranch("hitnuc_truth",hitnuc_truth,"hitnuc_truth[mcevts_truth]/I");
    CreateBranch("nuvtxx_truth",nuvtxx_truth,"nuvtxx_truth[mcevts_truth]/F");
    CreateBranch("nuvtxy_truth",nuvtxy_truth,"nuvtxy_truth[mcevts_truth]/F");
    CreateBranch("nuvtxz_truth",nuvtxz_truth,"nuvtxz_truth[mcevts_truth]/F");
    CreateBranch("nu_dcosx_truth",nu_dcosx_truth,"nu_dcosx_truth[mcevts_truth]/F");
    CreateBranch("nu_dcosy_truth",nu_dcosy_truth,"nu_dcosy_truth[mcevts_truth]/F");
    CreateBranch("nu_dcosz_truth",nu_dcosz_truth,"nu_dcosz_truth[mcevts_truth]/F");
    CreateBranch("lep_mom_truth",lep_mom_truth,"lep_mom_truth[mcevts_truth]/F");
    CreateBranch("lep_dcosx_truth",lep_dcosx_truth,"lep_dcosx_truth[mcevts_truth]/F");
    CreateBranch("lep_dcosy_truth",lep_dcosy_truth,"lep_dcosy_truth[mcevts_truth]/F");
    CreateBranch("lep_dcosz_truth",lep_dcosz_truth,"lep_dcosz_truth[mcevts_truth]/F");

    CreateBranch("tpx_flux",tpx_flux,"tpx_flux[mcevts_truth]/F");
    CreateBranch("tpy_flux",tpy_flux,"tpy_flux[mcevts_truth]/F");
    CreateBranch("tpz_flux",tpz_flux,"tpz_flux[mcevts_truth]/F");
    CreateBranch("tptype_flux",tptype_flux,"tptype_flux[mcevts_truth]/I");

    CreateBranch("genie_no_primaries",&genie_no_primaries,"genie_no_primaries/I");
    CreateBranch("genie_primaries_pdg",genie_primaries_pdg,"genie_primaries_pdg[genie_no_primaries]/I");
    CreateBranch("genie_Eng",genie_Eng,"genie_Eng[genie_no_primaries]/F");
    CreateBranch("genie_Px",genie_Px,"genie_Px[genie_no_primaries]/F");
    CreateBranch("genie_Py",genie_Py,"genie_Py[genie_no_primaries]/F");
    CreateBranch("genie_Pz",genie_Pz,"genie_Pz[genie_no_primaries]/F");
    CreateBranch("genie_P",genie_P,"genie_P[genie_no_primaries]/F");
    CreateBranch("genie_status_code",genie_status_code,"genie_status_code[genie_no_primaries]/I");
    CreateBranch("genie_mass",genie_mass,"genie_mass[genie_no_primaries]/F");
    CreateBranch("genie_trackID",genie_trackID,"genie_trackID[genie_no_primaries]/I");
    CreateBranch("genie_ND",genie_ND,"genie_ND[genie_no_primaries]/I");
    CreateBranch("genie_mother",genie_mother,"genie_mother[genie_no_primaries]/I");
  }

   if (hasCryInfo()){
    CreateBranch("mcevts_truthcry",&mcevts_truthcry,"mcevts_truthcry/I");  
    CreateBranch("cry_no_primaries",&cry_no_primaries,"cry_no_primaries/I");
    CreateBranch("cry_primaries_pdg",cry_primaries_pdg,"cry_primaries_pdg[cry_no_primaries]/I");
    CreateBranch("cry_Eng",cry_Eng,"cry_Eng[cry_no_primaries]/F");
    CreateBranch("cry_Px",cry_Px,"cry_Px[cry_no_primaries]/F");
    CreateBranch("cry_Py",cry_Py,"cry_Py[cry_no_primaries]/F");
    CreateBranch("cry_Pz",cry_Pz,"cry_Pz[cry_no_primaries]/F");
    CreateBranch("cry_P",cry_P,"cry_P[cry_no_primaries]/F");
    CreateBranch("cry_StartPointx",cry_StartPointx,"cry_StartPointx[cry_no_primaries]/F");
    CreateBranch("cry_StartPointy",cry_StartPointy,"cry_StartPointy[cry_no_primaries]/F");
    CreateBranch("cry_StartPointz",cry_StartPointz,"cry_StartPointz[cry_no_primaries]/F");   
    CreateBranch("cry_status_code",cry_status_code,"cry_status_code[cry_no_primaries]/I");
    CreateBranch("cry_mass",cry_mass,"cry_mass[cry_no_primaries]/F");
    CreateBranch("cry_trackID",cry_trackID,"cry_trackID[cry_no_primaries]/I");
    CreateBranch("cry_ND",cry_ND,"cry_ND[cry_no_primaries]/I");
    CreateBranch("cry_mother",cry_mother,"cry_mother[cry_no_primaries]/I");
  }  

  if (hasGeantInfo()){  
    CreateBranch("no_primaries",&no_primaries,"no_primaries/I");
    CreateBranch("geant_list_size",&geant_list_size,"geant_list_size/I");
    CreateBranch("geant_list_size_in_tpcAV",&geant_list_size_in_tpcAV,"geant_list_size_in_tpcAV/I");  
    CreateBranch("pdg",pdg,"pdg[geant_list_size]/I");
    CreateBranch("status",status,"status[geant_list_size]/I");
    CreateBranch("Mass",Mass,"Mass[geant_list_size]/F");
    CreateBranch("Eng",Eng,"Eng[geant_list_size]/F");
    CreateBranch("EndE",EndE,"EndE[geant_list_size]/F");
    CreateBranch("Px",Px,"Px[geant_list_size]/F");
    CreateBranch("Py",Py,"Py[geant_list_size]/F");
    CreateBranch("Pz",Pz,"Pz[geant_list_size]/F");
    CreateBranch("P",P,"P[geant_list_size]/F");
    CreateBranch("StartPointx",StartPointx,"StartPointx[geant_list_size]/F");
    CreateBranch("StartPointy",StartPointy,"StartPointy[geant_list_size]/F");
    CreateBranch("StartPointz",StartPointz,"StartPointz[geant_list_size]/F");
    CreateBranch("StartT",StartT,"StartT[geant_list_size]/F");
    CreateBranch("EndPointx",EndPointx,"EndPointx[geant_list_size]/F");
    CreateBranch("EndPointy",EndPointy,"EndPointy[geant_list_size]/F");
    CreateBranch("EndPointz",EndPointz,"EndPointz[geant_list_size]/F");
    CreateBranch("EndT",EndT,"EndT[geant_list_size]/F");
    CreateBranch("theta",theta,"theta[geant_list_size]/F");
    CreateBranch("phi",phi,"phi[geant_list_size]/F");
    CreateBranch("theta_xz",theta_xz,"theta_xz[geant_list_size]/F");
    CreateBranch("theta_yz",theta_yz,"theta_yz[geant_list_size]/F");
    CreateBranch("pathlen",pathlen,"pathlen[geant_list_size]/F");
    CreateBranch("inTPCActive",inTPCActive,"inTPCActive[geant_list_size]/I");  
    CreateBranch("StartPointx_tpcAV",StartPointx_tpcAV,"StartPointx_tpcAV[geant_list_size]/F");
    CreateBranch("StartPointy_tpcAV",StartPointy_tpcAV,"StartPointy_tpcAV[geant_list_size]/F");
    CreateBranch("StartPointz_tpcAV",StartPointz_tpcAV,"StartPointz_tpcAV[geant_list_size]/F");
    CreateBranch("EndPointx_tpcAV",EndPointx_tpcAV,"EndPointx_tpcAV[geant_list_size]/F");
    CreateBranch("EndPointy_tpcAV",EndPointy_tpcAV,"EndPointy_tpcAV[geant_list_size]/F");
    CreateBranch("EndPointz_tpcAV",EndPointz_tpcAV,"EndPointz_tpcAV[geant_list_size]/F");
    CreateBranch("NumberDaughters",NumberDaughters,"NumberDaughters[geant_list_size]/I");
    CreateBranch("Mother",Mother,"Mother[geant_list_size]/I");
    CreateBranch("TrackId",TrackId,"TrackId[geant_list_size]/I");
    CreateBranch("MergedId", MergedId, "MergedId[geant_list_size]/I");
    CreateBranch("MotherNuId", MotherNuId, "MotherNuId[geant_list_size]/I");
    CreateBranch("process_primary",process_primary,"process_primary[geant_list_size]/I");
    CreateBranch("processname", processname);
  }

  if (hasAuxDetector()) {
    // Geant information is required to fill aux detector information.
    // if fSaveGeantInfo is not set to true, show an error message and quit!
    if (!hasGeantInfo()){
      throw art::Exception(art::errors::Configuration)
      << "Saving Auxiliary detector information requies saving GEANT information, "
      <<"please set fSaveGeantInfo flag to true in your fhicl file and rerun.\n"; 
    }    
    std::ostringstream sstr;
    sstr << "[" << kMaxAuxDets << "]";
    std::string MaxAuxDetIndexStr = sstr.str();
    CreateBranch("NAuxDets",     NAuxDets, "NAuxDets[geant_list_size]/s");
    CreateBranch("AuxDetID",     AuxDetID, "AuxDetID[geant_list_size]" + MaxAuxDetIndexStr + "/S");
    CreateBranch("AuxDetEntryX", entryX,   "AuxDetEntryX[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetEntryY", entryY,   "AuxDetEntryY[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetEntryZ", entryZ,   "AuxDetEntryZ[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetEntryT", entryT,   "AuxDetEntryT[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitX",  exitX,    "AuxDetExitX[geant_list_size]"  + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitY",  exitY,    "AuxDetExitY[geant_list_size]"  + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitZ",  exitZ,    "AuxDetExitZ[geant_list_size]"  + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitT",  exitT,    "AuxDetExitT[geant_list_size]"  + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitPx", exitPx,   "AuxDetExitPx[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitPy", exitPy,   "AuxDetExitPy[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("AuxDetExitPz", exitPz,   "AuxDetExitPz[geant_list_size]" + MaxAuxDetIndexStr + "/F");
    CreateBranch("CombinedEnergyDep", CombinedEnergyDep,
      "CombinedEnergyDep[geant_list_size]" + MaxAuxDetIndexStr + "/F");
  } // if hasAuxDetector
  
} // sbnd::AnalysisTreeDataStruct::SetAddresses()


//------------------------------------------------------------------------------
//---  AnalysisTree
//---

sbnd::AnalysisTree::AnalysisTree(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset),
  fTree(nullptr), fData(nullptr),
  fDigitModuleLabel         (pset.get< std::string >("DigitModuleLabel")                   ),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel")                    ),
  fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel")                ),
  fCalDataModuleLabel       (pset.get< std::string >("CalDataModuleLabel")                 ),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel")                ),
  fCryGenModuleLabel        (pset.get< std::string >("CryGenModuleLabel")                  ), 
  fG4ModuleLabel            (pset.get< std::string >("G4ModuleLabel")                      ),
  fPFParticleModuleLabel    (pset.get< std::string> ("PFParticleModuleLabel")              ),
  fShowerModuleLabel        (pset.get< std::string> ("ShowerModuleLabel")                  ),
  fVertexModuleLabel        (pset.get< std::vector<std::string> >("VertexModuleLabel")     ),
  fTrackModuleLabel         (pset.get< std::vector<std::string> >("TrackModuleLabel")      ),
  fCalorimetryModuleLabel   (pset.get< std::vector<std::string> >("CalorimetryModuleLabel")),
  fParticleIDModuleLabel    (pset.get< std::vector<std::string> >("ParticleIDModuleLabel") ),
  fPOTModuleLabel           (pset.get< std::string >("POTModuleLabel")                     ),

  fUseBuffer                (pset.get< bool >("UseBuffers", false)),
  fSaveAuxDetInfo           (pset.get< bool >("SaveAuxDetInfo", false)),
  fSaveCryInfo              (pset.get< bool >("SaveCryInfo", false)),  
  fSaveGenieInfo	          (pset.get< bool >("SaveGenieInfo", false)), 
  fSaveGeantInfo	          (pset.get< bool >("SaveGeantInfo", false)), 
  fSaveHitInfo	            (pset.get< bool >("SaveHitInfo", false)), 
  fSaveTrackInfo	          (pset.get< bool >("SaveTrackInfo", false)), 
  fSaveShowerInfo	          (pset.get< bool >("SaveShowerInfo", false)), 
  fSaveVertexInfo	          (pset.get< bool >("SaveVertexInfo", false)),
  fSaveHierarchyInfo        (pset.get< std::vector<bool> >("SaveHierarchyInfo", {false})),
  fSaveShowerHierarchyInfo  (pset.get< bool >("SaveShowerHierarchyInfo", false)),
  //fCosmicTaggerAssocLabel  (pset.get<std::vector< std::string > >("CosmicTaggerAssocLabel") ),
  //fFlashMatchAssocLabel (pset.get<std::vector< std::string > >("FlashMatchAssocLabel") ),
  isCosmics(false),
  fSaveCaloCosmics          (pset.get< bool >("SaveCaloCosmics",false)),
  fG4minE                   (pset.get< float>("G4minE",0.01))
{
  if (fSaveAuxDetInfo == true) fSaveGeantInfo = true;
  mf::LogInfo("AnalysisTree") << "Configuration:"
    << "\n  UseBuffers: " << std::boolalpha << fUseBuffer
    ;
  if (GetNTrackers() > kMaxTrackers) {
    throw art::Exception(art::errors::Configuration)
      << "AnalysisTree currently supports only up to " << kMaxTrackers
      << " tracking algorithms, but " << GetNTrackers() << " are specified."
      << "\nYou can increase kMaxTrackers and recompile.";
  } // if too many trackers
  if (fTrackModuleLabel.size() != fVertexModuleLabel.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fVertexModuleLabel.size() = "<<fVertexModuleLabel.size();
  }
  if (fTrackModuleLabel.size() != fCalorimetryModuleLabel.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fCalorimetryModuleLabel.size() = "<<fCalorimetryModuleLabel.size();
  }
  if (fTrackModuleLabel.size() != fParticleIDModuleLabel.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fParticleIDModuleLabel.size() = "<<fParticleIDModuleLabel.size();
  }
  if (fTrackModuleLabel.size() != fSaveHierarchyInfo.size()){
    throw art::Exception(art::errors::Configuration)
      << "fTrackModuleLabel.size() = "<<fTrackModuleLabel.size()<<" does not match "
      << "fSaveHierarchyInfo.size() = "<<fSaveHierarchyInfo.size();
  }
} // sbnd::AnalysisTree::AnalysisTree()

//-------------------------------------------------
sbnd::AnalysisTree::~AnalysisTree()
{
  DestroyData();
}

void sbnd::AnalysisTree::CreateTree(bool bClearData /* = false */) {
  if (!fTree) {
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("anatree","analysis tree");
  }
  CreateData(bClearData);
  SetAddresses();
} // sbnd::AnalysisTree::CreateTree()


void sbnd::AnalysisTree::beginSubRun(const art::SubRun& sr)
{
  art::Handle< sumdata::POTSummary > potListHandle;
  //sr.getByLabel(fPOTModuleLabel,potListHandle);

  if(sr.getByLabel(fPOTModuleLabel,potListHandle))
    SubRunData.pot=potListHandle->totpot;
  else
    SubRunData.pot=0.;
}

void sbnd::AnalysisTree::analyze(const art::Event& evt)
{
  //services
  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  // auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  // auto const* LArProp = lar::providerFrom<detinfo::LArPropertiesService>();

  // collect the sizes which might me needed to resize the tree data structure:
  bool isMC = !evt.isRealData();
  
  // * hits
  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);


  // * MC truth information
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);
    
  // *MC truth cosmic generator information
  art::Handle< std::vector<simb::MCTruth> > mctruthcryListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclistcry;
  if (fSaveCryInfo){
    if (evt.getByLabel(fCryGenModuleLabel,mctruthcryListHandle))
      art::fill_ptr_vector(mclistcry, mctruthcryListHandle);
  }       
  
  art::Ptr<simb::MCTruth> mctruthcry;
  int nCryPrimaries = 0;
   
  if (fSaveCryInfo){
    mctruthcry = mclistcry[0];      
    nCryPrimaries = mctruthcry->NParticles();  
  } 
  
  int nGeniePrimaries = 0, nGEANTparticles = 0, nMCNeutrinos = 0;
  
  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

  art::Ptr<simb::MCTruth> mctruth;
  //Brailsford 2017/10/16
  //Fix for issue 17917
  //With the code change to accept multiple neutrinos per TTree::Entry into the TTree, this int is no longer needed (it makes compilation fail due to a warning)
  //int imc = 0;
  if (isMC) { //is MC
    // GENIE
    if (!mclist.empty()){//at least one mc record
      //if (fSaveGenieInfo){
        //if (mclist[0]->NeutrinoSet()){//is neutrino
        //sometimes there can be multiple neutrino interactions in one spill
        //this is trying to find the primary interaction
        //by looking for the highest energy deposition
        //std::map<art::Ptr<simb::MCTruth>,double> mctruthemap;
      static bool isfirsttime = true;
      if (isfirsttime){
	      for (size_t i = 0; i<hitlist.size(); i++){
	        //if (hitlist[i]->View() == geo::kV){//collection view
          // tbrooks: use TrackIDEs rather than eveTrackIDEs because the eve ID doesn't always seem to correspond to the g4 track FIXME may need further investigation
                std::vector<sim::TrackIDE> eveIDs = bt_serv->HitToTrackIDEs(clockData, hitlist[i]);
	        for (size_t e = 0; e<eveIDs.size(); e++){
	          try {
		    art::Ptr<simb::MCTruth> ev_mctruth = pi_serv->TrackIdToMCTruth_P(eveIDs[e].trackID);
		    //mctruthemap[ev_mctruth]+=eveIDs[e].energy;
		    if (ev_mctruth->Origin() == simb::kCosmicRay) isCosmics = true;
		  } 
		  catch(...){
		    std::cout<<"Exception thrown matching eveTrackId "<<eveIDs[e].trackID<<" to MCTruth\n";
		  }
		}
	        //}
	      }
	      isfirsttime = false;
	      if (fSaveCaloCosmics) isCosmics = false; //override to save calo info
      }

//        double maxenergy = -1;
//        int imc0 = 0;
//        for (std::map<art::Ptr<simb::MCTruth>,double>::iterator ii=mctruthemap.begin(); ii!=mctruthemap.end(); ++ii){
//          if ((ii->second)>maxenergy){
//            maxenergy = ii->second;
//            mctruth = ii->first;
//            imc = imc0;
//          }
//          imc0++;
//        }

	  //imc = 0; //set imc to 0 to solve a confusion for BNB+cosmic files where there are two MCTruth
      mctruth = mclist[0];

      if (mctruth->NeutrinoSet()) nGeniePrimaries = mctruth->NParticles();
	//} //end (fSaveGenieInfo)
      
      const sim::ParticleList& plist = pi_serv->ParticleList();
      nGEANTparticles = plist.size();

      // to know the number of particles in AV would require
      // looking at all of them; so we waste some memory here
    } // if have MC truth
    MF_LOG_DEBUG("AnalysisTree") << "Expected "
      << nGEANTparticles << " GEANT particles, "
      << nGeniePrimaries << " GENIE particles";
  } // if MC
  
  //Brailsford 2017/10/16
  //Initially call the number of neutrinos to be stored the number of MCTruth objects.  This is not strictly true i.e. BNB + cosmic overlay but we will count the number of neutrinos later
  nMCNeutrinos = mclist.size();

  CreateData(); // tracker data is created with default constructor
  if (fSaveGenieInfo){
    fData->ResizeGenie(nGeniePrimaries);
    fData->ResizeMCNeutrino(nMCNeutrinos);
  }
  if (fSaveCryInfo)
    fData->ResizeCry(nCryPrimaries);
  if (fSaveGeantInfo)    
    fData->ResizeGEANT(nGEANTparticles);
  fData->ClearLocalData(); // don't bother clearing tracker data yet
  
//  const size_t Nplanes       = 3; // number of wire planes; pretty much constant...
  const size_t NTrackers = GetNTrackers(); // number of trackers passed into fTrackModuleLabel
  const size_t NHits     = hitlist.size(); // number of hits
  // make sure there is the data, the tree and everything;
  CreateTree();

  /// transfer the run and subrun data to the tree data object
//  fData->RunData = RunData;
  fData->SubRunData = SubRunData;

  fData->isdata = int(!isMC);
 
  /*****************************************************************************/
  //            Added by Rhiannon on 2nd October 2018 for MCP 0.9
  // If the user is requesting to save hierarchical information, the producers
  // must be from pandora and are therefore able to go via PFParticles to get
  // track and vertex associations and their corresponding hierarchy
  //
  // IF saving the hierarchy has been requested, check that the pfparticle
  // object exists, and that it has track and vertex associations 
  //  Fill maps with tracks and vertices as the key, and pfparticles as the object
  //  When saving track and vertex information, add additional query for 
  //  saving the hierarchy and save the information
  //
  /// Declare the track handle and vector for all the different tracking producers
  std::vector< art::Handle< std::vector<recob::Track> > >  trackListHandle(NTrackers);
  std::vector< art::Handle< std::vector<recob::Vertex> > > vtxListHandle(NTrackers);
  art::Handle< std::vector<recob::Shower> >                showerHandle;
  std::vector< std::vector<art::Ptr<recob::Track> > >      tracklist(NTrackers);
  std::vector< std::vector< art::Ptr<recob::Vertex> > >    vtxlist(NTrackers);
  std::vector< art::Ptr<recob::Shower> >                   shwlist;
  std::vector< int > NVertices;

  // Declare maps for track and vertex to pfparticle associations 
  //  this is the opposite way around to how they are produced
  typedef std::map<art::Ptr<recob::Track>,  art::Ptr<recob::PFParticle> > trkPfpMap;
  typedef std::map<art::Ptr<recob::Vertex>, art::Ptr<recob::PFParticle> > vtxPfpMap;
  typedef std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle> > shwPfpMap;
  typedef std::map<art::Ptr<recob::Track>,  art::Ptr<recob::PFParticle> >::const_iterator trkPfpMapIt;
  typedef std::map<art::Ptr<recob::Vertex>, art::Ptr<recob::PFParticle> >::const_iterator vtxPfpMapIt;
  typedef std::map<art::Ptr<recob::Shower>, art::Ptr<recob::PFParticle> >::const_iterator shwPfpMapIt;

  shwPfpMap showerPFParticleMap;
  lar_pandora::PFParticleVector pfplist;
  lar_pandora::PFParticleMap pfpmap;
  std::vector< trkPfpMap > trackerPFParticleMaps;
  std::vector< vtxPfpMap > verticesPFParticleMaps;

  art::Handle< std::vector<recob::PFParticle> > pfpHandle;
  if(evt.getByLabel(fPFParticleModuleLabel,pfpHandle)){
    if(pfpHandle->size()){
      art::fill_ptr_vector(pfplist, pfpHandle);
      lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfplist, pfpmap);
    }
  }
  trkPfpMap trackPFParticleMap;
  vtxPfpMap vertexPFParticleMap;
  // Loop over trackers and see if the hierarchy information needs saving
  for (unsigned int iTracker=0; iTracker < NTrackers; ++iTracker){
    verticesPFParticleMaps.push_back(vertexPFParticleMap);
    trackerPFParticleMaps.push_back(trackPFParticleMap);
  
    if(!fSaveHierarchyInfo[iTracker]) continue;

    // Get and check that the pfparticle handle is valid and exists
    if(pfpHandle->empty())
      mf::LogError("AnalysisTree:limits") << " Event has no PFParticle information ";

    // Vertex, shower and track associations
    art::FindManyP<recob::Vertex> fvtx(pfpHandle, evt, fVertexModuleLabel[iTracker]);
    art::FindManyP<recob::Track> fmtrk(pfpHandle, evt, fTrackModuleLabel[iTracker]);

    for(unsigned int i = 0; i < pfpHandle->size(); ++i) {
      const art::Ptr< recob::PFParticle > pfp( pfpHandle, i );

      // Get vertex association
      std::vector< art::Ptr<recob::Vertex> > vtxAssn = fvtx.at(pfp.key());

      // Make sure there is exactly 1 vertex associated to each pfparticle
      if(vtxAssn.size() > 1){
        throw cet::exception("AnalysisTree:limits") << "PFParticle has " << vtxAssn.size() << " associated vertices, should only have 1 or 0 ";
      }
      if(vtxAssn.size() == 0) continue;
      verticesPFParticleMaps.at(iTracker).emplace(vtxAssn.front(),pfp);
    }
    // Fill the vector of maps for each tracker
    for(unsigned int i = 0; i < pfpHandle->size(); ++i) {
      const art::Ptr< recob::PFParticle > pfp( pfpHandle, i );

      std::vector< art::Ptr<recob::Track> >  trkAssn = fmtrk.at(pfp.key());
      if(trkAssn.size()  > 1){
        throw cet::exception("AnalysisTree:limits") << "PFParticle has " << trkAssn.size() << " associated tracks, should only have 1 or 0 ";
      }
      if(trkAssn.size() == 0) continue;
      trackerPFParticleMaps.at(iTracker).emplace(trkAssn.front(),pfp);
    }
  }// end loop over trackers

  if(fSaveShowerHierarchyInfo){
    art::Handle< std::vector<recob::PFParticle> > pfpHandle;
    if(!evt.getByLabel(fPFParticleModuleLabel,pfpHandle)){
      if(!pfpHandle->size())
      mf::LogError("AnalysisTree:limits") << " Event has no PFParticle information ";
    }
    if(pfpHandle->size()){
      art::fill_ptr_vector(pfplist, pfpHandle);
      lar_pandora::LArPandoraHelper::BuildPFParticleMap(pfplist, pfpmap);
  
      art::FindManyP<recob::Shower> fshw(pfpHandle, evt, fShowerModuleLabel);
      for(unsigned int i = 0; i < pfpHandle->size(); ++i) {
        art::Ptr< recob::PFParticle > pfp( pfpHandle, i );

        // Get vertex association
        std::vector< art::Ptr<recob::Shower> > shwAssn = fshw.at(pfp->Self());

        // Make sure there is maximum 1 shower associated to each pfparticle
        if(shwAssn.size()  > 1){
          mf::LogError("AnalysisTree:limits") << "PFParticle has " << shwAssn.size() << " associated showers, should only have 1 or 0 ";
          continue;
        }
        if(shwAssn.size() == 0) continue;
        showerPFParticleMap.emplace(shwAssn[0],pfp);

      }
    }
  } // save maps for SaveShowerHierarchyInfo

  if (evt.getByLabel(fShowerModuleLabel,showerHandle))
    art::fill_ptr_vector(shwlist, showerHandle);
  const size_t NShowers = shwlist.size();
  
  for (unsigned int it = 0; it < NTrackers; ++it){
    if (evt.getByLabel(fTrackModuleLabel[it],trackListHandle[it]))
      art::fill_ptr_vector(tracklist[it], trackListHandle[it]);
    if (evt.getByLabel(fVertexModuleLabel[it],vtxListHandle[it]))
      art::fill_ptr_vector(vtxlist[it], vtxListHandle[it]);
  }

  art::Handle< std::vector<simb::MCFlux> > mcfluxListHandle;
  std::vector<art::Ptr<simb::MCFlux> > fluxlist;
  if (evt.getByLabel(fGenieGenModuleLabel,mcfluxListHandle))
    art::fill_ptr_vector(fluxlist, mcfluxListHandle);

  std::vector<const sim::AuxDetSimChannel*> fAuxDetSimChannels;
  if (fSaveAuxDetInfo && fSaveGeantInfo){
    evt.getView(fLArG4ModuleLabel, fAuxDetSimChannels);
  }

  std::vector<const sim::SimChannel*> fSimChannels;
  if (isMC && fSaveGeantInfo){
    evt.getView(fLArG4ModuleLabel, fSimChannels);
  }

  fData->run = evt.run();
  fData->subrun = evt.subRun();
  fData->event = evt.id().event();

  art::Timestamp ts = evt.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  fData->evttime = tts.AsDouble();

  //copied from MergeDataPaddles.cxx
  art::Handle< raw::BeamInfo > beam;
  if (evt.getByLabel("beam",beam)){
    fData->beamtime = (double)beam->get_t_ms();
    fData->beamtime/=1000.; //in second
  }

//  std::cout<<detprop->NumberTimeSamples()<<" "<<detprop->ReadOutWindowSize()<<std::endl;
//  std::cout<<geom->DetHalfHeight()*2<<" "<<geom->DetHalfWidth()*2<<" "<<geom->DetLength()<<std::endl;
//  std::cout<<geom->Nwires(0)<<" "<<geom->Nwires(1)<<" "<<geom->Nwires(2)<<std::endl;

  //hit information
  if (fSaveHitInfo){
    fData->no_hits = (int) NHits;
    if (NHits > kMaxHits) {
      // got this error? consider increasing kMaxHits
      // (or ask for a redesign using vectors)
      mf::LogError("AnalysisTree:limits") << "event has " << NHits
        << " hits, only kMaxHits=" << kMaxHits << " stored in tree";
    }
    for (size_t i = 0; i < NHits && i < kMaxHits ; ++i){//loop over hits
      fData->hit_channel[i] = hitlist[i]->Channel();
      fData->hit_tpc[i]     = hitlist[i]->WireID().TPC;
      fData->hit_plane[i]   = hitlist[i]->WireID().Plane;
      fData->hit_wire[i]    = hitlist[i]->WireID().Wire;
      fData->hit_peakT[i]   = hitlist[i]->PeakTime();
      fData->hit_charge[i]  = hitlist[i]->Integral();
      fData->hit_ph[i]  = hitlist[i]->PeakAmplitude();
      fData->hit_startT[i] = hitlist[i]->PeakTimeMinusRMS();
      fData->hit_endT[i] = hitlist[i]->PeakTimePlusRMS();
      /*
      for (unsigned int it=0; it<fTrackModuleLabel.size();++it){
        art::FindManyP<recob::Track> fmtk(hitListHandle,evt,fTrackModuleLabel[it]);
        if (fmtk.at(i).size()!=0){
          hit_trkid[it][i] = fmtk.at(i)[0]->ID();
        }
        else
          hit_trkid[it][i] = 0;
      }
      */

      if (!evt.isRealData()){
         fData -> hit_nelec[i] = 0;
         fData -> hit_energy[i] = 0;
         const sim::SimChannel* chan = 0;
         for(size_t sc = 0; sc < fSimChannels.size(); ++sc){
           if(fSimChannels[sc]->Channel() == hitlist[i]->Channel()) chan = fSimChannels[sc];
         }
         if (chan){
           for (auto const& tdcide: chan->TDCIDEMap()) {
             // loop over the vector of IDE objects.
             std::vector<sim::IDE> const& idevec = tdcide.second;
             for(auto const& ide: idevec) {
                fData -> hit_nelec[i] += ide.numElectrons;
                fData -> hit_energy[i] += ide.energy;
             } // all IDEs in the TDC tick
           } // all TDC ticks
         }
       }
    }

    if (evt.getByLabel(fHitsModuleLabel,hitListHandle)){
      //Find tracks associated with hits
      art::FindManyP<recob::Track> fmtk(hitListHandle,evt,fTrackModuleLabel[0]);
      for (size_t i = 0; i < NHits && i < kMaxHits ; ++i){//loop over hits
        if (fmtk.isValid()){
	  if (fmtk.at(i).size()!=0){
	    fData->hit_trkid[i] = fmtk.at(i)[0]->ID();
	  }
	  else
	    fData->hit_trkid[i] = -1;
        }
      }
    }
  }// end (fSaveHitInfo) 

  //vertex information
  if (fSaveVertexInfo){
    for (unsigned int iTracker=0; iTracker < NTrackers; ++iTracker){
      AnalysisTreeDataStruct::VertexDataStruct& VertexData = fData->GetVertexData(iTracker);
    
      size_t NVertices = vtxlist[iTracker].size();
      // allocate enough space for this number of tracks (but at least for one of them!)
      VertexData.SetMaxVertices(std::max(NVertices, (size_t) 1));
      VertexData.Clear(); // clear all the data
    
      VertexData.nvtx = static_cast<int>(NVertices);
    
      // now set the tree addresses to the newly allocated memory;
      // this creates the tree branches in case they are not there yet
      SetVerticesAddresses(iTracker);
      if (NVertices > VertexData.GetMaxVertices()) {
        // got this error? it might be a bug,
        // since we are supposed to have allocated enough space to fit all vertices
        mf::LogError("AnalysisTree:limits") << "event has " << NVertices
          << " " << fVertexModuleLabel[iTracker] << " vertices, only "
          << VertexData.GetMaxVertices() << " stored in tree";
      }

      for (size_t i = 0; i < NVertices && i < kMaxVertices ; ++i){//loop over vertices
        art::Ptr<recob::Vertex> pvertex(vtxListHandle[iTracker], i);
//      art::Ptr<recob::Vertex> pvertex = vtxlist[iTracker].at(i);
        Double_t xyz[3];
        pvertex->XYZ(xyz);
        for (size_t j = 0; j<3; ++j){
          VertexData.vtx[i][j] = xyz[j];
        }
        
        // If also saving the hierarchy info, set the primaryvtx boolean
        if(fSaveHierarchyInfo[iTracker]){
          vtxPfpMap vertexPFParticleMap = verticesPFParticleMaps.at(iTracker);
          vtxPfpMapIt it = vertexPFParticleMap.find(pvertex);
          // Check there is a map entry for this vertex
          
          if(it == vertexPFParticleMap.end()) continue;

          art::Ptr<recob::PFParticle> tempParticle = it->second;
          if(tempParticle->IsPrimary())
            VertexData.primaryvtx[i] = 1;
          else VertexData.primaryvtx[i] = 0;
        }// end (fSaveHierarchyInfo)
      }// end loop over vertices
    }// end loop over trackers
  }// end (fSaveVertexInfo)
  
  // shower information
  if (fSaveShowerInfo){
    AnalysisTreeDataStruct::ShowerDataStruct& ShowerData = fData->GetShowerData();
    if (NShowers > kMaxShowers){
      // got this error? consider increasing kMaxShowers
      // (or ask for a redesign using vectors)
      mf::LogError("AnalysisTree:limits") << "event has " << NShowers
        << " showers, only kMaxShowers=" << kMaxShowers << " stored in tree";
    }
    ShowerData.SetMaxShowers(std::max(NShowers, (size_t) 1));
    ShowerData.Clear();

    ShowerData.nshowers = (int) NShowers;
    
    SetShowerAddresses();

    for (size_t i = 0; i < NShowers && i < kMaxShowers ; ++i){//loop over showers
      const art::Ptr<recob::Shower> &shw = shwlist[i];
     
      ShowerData.shwId[i]        = shw->ID();
      ShowerData.shwbestplane[i] = shw->best_plane();
      ShowerData.shwlength[i]    = shw->Length();
      ShowerData.shwopenangle[i] = shw->OpenAngle();
      ShowerData.shwstartx[i]    = shw->ShowerStart()[0];
      ShowerData.shwstarty[i]    = shw->ShowerStart()[1];
      ShowerData.shwstartz[i]    = shw->ShowerStart()[2];
      ShowerData.shwdirx[i]      = shw->Direction()[0];
      ShowerData.shwdiry[i]      = shw->Direction()[1];
      ShowerData.shwdirz[i]      = shw->Direction()[2];
      
      // If also saving the hierarchy info, set the primaryvtx boolean
      if(fSaveShowerHierarchyInfo){
        shwPfpMapIt it;
        // Check there is a map entry for this vertex
        it = showerPFParticleMap.find(shw);
        if(it == showerPFParticleMap.end()) continue;

        art::Ptr<recob::PFParticle> tempParticle = it->second;
        ShowerData.shwisprimary[i]   = std::round(static_cast<Int_t>(lar_pandora::LArPandoraHelper::IsFinalState(pfpmap,tempParticle)));
        ShowerData.shwndaughters[i]  = tempParticle->NumDaughters();
        ShowerData.shwpfpid[i]       = tempParticle->Self();
        ShowerData.shwparentpfpid[i] = tempParticle->Parent();
      }// end (fSaveHierarchyInfo)
    }
  }// end (fSaveShowerInfo)

  //track information for multiple trackers
  if (fSaveTrackInfo){
    for (unsigned int iTracker=0; iTracker < NTrackers; ++iTracker){
      AnalysisTreeDataStruct::TrackDataStruct& TrackerData = fData->GetTrackerData(iTracker);
    
      size_t NTracks = tracklist[iTracker].size();
      // allocate enough space for this number of tracks (but at least for one of them!)
      TrackerData.SetMaxTracks(std::max(NTracks, (size_t) 1));
      TrackerData.Clear(); // clear all the data
    
      TrackerData.ntracks = (int) NTracks;
    
      // now set the tree addresses to the newly allocated memory;
      // this creates the tree branches in case they are not there yet
      SetTrackerAddresses(iTracker);
      if (NTracks > TrackerData.GetMaxTracks()) {
        // got this error? it might be a bug,
        // since we are supposed to have allocated enough space to fit all tracks
        mf::LogError("AnalysisTree:limits") << "event has " << NTracks
          << " " << fTrackModuleLabel[iTracker] << " tracks, only "
          << TrackerData.GetMaxTracks() << " stored in tree";
      }
      AnalysisTreeDataStruct::VertexDataStruct& VertexData = fData->GetVertexData(iTracker);
    
      size_t NVertices = vtxlist[iTracker].size();
    
      // now set the tree addresses to the newly allocated memory;
      // this creates the tree branches in case they are not there yet
      SetVerticesAddresses(iTracker);
      if (NVertices > VertexData.GetMaxVertices()) {
        // got this error? it might be a bug,
        // since we are supposed to have allocated enough space to fit all tracks
        mf::LogError("AnalysisTree:limits") << "event has " << NVertices
          << " " << fVertexModuleLabel[iTracker] << " vertices, only "
          << VertexData.GetMaxVertices() << " stored in tree";
      }
    
      //call the track momentum algorithm that gives you momentum based on track range
      trkf::TrackMomentumCalculator trkm;

      for(size_t iTrk=0; iTrk < NTracks; ++iTrk){//loop over tracks
        //Cosmic Tagger information
        if (fCosmicTaggerAssocLabel.size() > iTracker) {
          art::FindManyP<anab::CosmicTag> fmct(trackListHandle[iTracker],evt,fCosmicTaggerAssocLabel[iTracker]);
          if (fmct.isValid()){          
            TrackerData.trkncosmictags_tagger[iTrk]     = fmct.at(iTrk).size();
            if (fmct.at(iTrk).size()>0){
              if(fmct.at(iTrk).size()>1)
                std::cerr << "\n Warning : more than one cosmic tag per track in module! assigning the first tag to the track" << fCosmicTaggerAssocLabel[iTracker];
              TrackerData.trkcosmicscore_tagger[iTrk] = fmct.at(iTrk).at(0)->CosmicScore();
              TrackerData.trkcosmictype_tagger[iTrk] = fmct.at(iTrk).at(0)->CosmicType();
            }
          }
        } // if we have matching fCosmicTaggerAssocLabel

        //Flash match compatibility information
        if (fFlashMatchAssocLabel.size() > iTracker) {
          //Unlike CosmicTagger, Flash match doesn't assign a cosmic tag for every track. For those tracks, AnalysisTree initializes them with -9999 or -99999
          art::FindManyP<anab::CosmicTag> fmbfm(trackListHandle[iTracker],evt,fFlashMatchAssocLabel[iTracker]);
          if (fmbfm.isValid()){  
            TrackerData.trkncosmictags_flashmatch[iTrk] = fmbfm.at(iTrk).size();
            if (fmbfm.at(iTrk).size()>0){
              if(fmbfm.at(iTrk).size()>1) 
                std::cerr << "\n Warning : more than one cosmic tag per track in module! assigning the first tag to the track" << fFlashMatchAssocLabel[iTracker];
            TrackerData.trkcosmicscore_flashmatch[iTrk] = fmbfm.at(iTrk).at(0)->CosmicScore();
                TrackerData.trkcosmictype_flashmatch[iTrk] = fmbfm.at(iTrk).at(0)->CosmicType();
        //std::cout<<"\n"<<evt.event()<<"\t"<<iTrk<<"\t"<<fmbfm.at(iTrk).at(0)->CosmicScore()<<"\t"<<fmbfm.at(iTrk).at(0)->CosmicType();
            }
          }
        } // if we have matching fFlashMatchAssocLabel
        
        art::Ptr<recob::Track> ptrack(trackListHandle[iTracker], iTrk);
        const recob::Track& track = *ptrack;

        double tlen = 0.; 
        double mom = 0.;
        int TrackID = -1;
 
        int ntraj = track.NumberTrajectoryPoints();
        if (ntraj > 0) {
          const auto& pos       = track.Vertex();
          const auto& dir_start = track.VertexDirection();
          const auto& dir_end   = track.EndDirection();
          const auto& end       = track.End();

          tlen        = length(track);
          if(track.HasMomentum() > 0)
            mom = track.VertexMomentum();
          TrackID = track.ID();
          
          double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
          double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
          double dpos = bdist(pos);
          double dend = bdist(end);
        
          TrackerData.trkId[iTrk]                 = TrackID;
          TrackerData.trkstartx[iTrk]             = pos.X();
          TrackerData.trkstarty[iTrk]             = pos.Y();
          TrackerData.trkstartz[iTrk]             = pos.Z();
          TrackerData.trkstartd[iTrk]		          = dpos;
          TrackerData.trkendx[iTrk]		            = end.X();
          TrackerData.trkendy[iTrk]		            = end.Y();
          TrackerData.trkendz[iTrk]		            = end.Z();
          TrackerData.trkendd[iTrk]		            = dend;
          TrackerData.trktheta[iTrk]		          = dir_start.Theta();
          TrackerData.trkphi[iTrk]		            = dir_start.Phi();
          TrackerData.trkstartdcosx[iTrk]	        = dir_start.X();
          TrackerData.trkstartdcosy[iTrk]	        = dir_start.Y();
          TrackerData.trkstartdcosz[iTrk]	        = dir_start.Z();
          TrackerData.trkenddcosx[iTrk] 	        = dir_end.X();
          TrackerData.trkenddcosy[iTrk] 	        = dir_end.Y();
          TrackerData.trkenddcosz[iTrk] 	        = dir_end.Z();
          TrackerData.trkthetaxz[iTrk]  	        = theta_xz;
          TrackerData.trkthetayz[iTrk]  	        = theta_yz;
          TrackerData.trkmom[iTrk]		            = mom;
          TrackerData.trklen[iTrk]		            = tlen;
          TrackerData.trkmomrange[iTrk] 	        = trkm.GetTrackMomentum(tlen,13);
          TrackerData.trkmommschi2[iTrk]	        = trkm.GetMomentumMultiScatterChi2(ptrack);
          TrackerData.trkmommsllhd[iTrk]	        = trkm.GetMomentumMultiScatterLLHD(ptrack);
        } // if we have trajectory

        // find vertices associated with this track
       /* 
        art::FindMany<recob::Vertex> fmvtx(trackListHandle[iTracker], evt, fVertexModuleLabel[iTracker]);
        if(fmvtx.isValid()) {
          std::vector<const recob::Vertex*> verts = fmvtx.at(iTrk);
          // should have two at most
          for(size_t ivx = 0; ivx < verts.size(); ++ivx) {
            verts[ivx]->XYZ(xyz);
            // find the vertex in TrackerData to get the index
            short theVtx = -1;
            for(short jvx = 0; jvx < TrackerData.nvtx; ++jvx) {
              if(TrackerData.vtx[jvx][2] == xyz[2]) {
                theVtx = jvx;
                break;
              }
            } // jvx
            // decide if it should be assigned to the track Start or End.
            // A simple dz test should suffice
            if(fabs(xyz[2] - TrackerData.trkstartz[iTrk]) < 
               fabs(xyz[2] - TrackerData.trkendz[iTrk])) {
              TrackerData.trksvtxid[iTrk] = theVtx;
            } else {
              TrackerData.trkevtxid[iTrk] = theVtx;
            }
          } // vertices
        } // fmvtx.isValid()
        */
        Float_t minsdist = 10000;
        Float_t minedist = 10000;
        for (int ivx = 0; ivx < VertexData.nvtx && ivx < kMaxVertices; ++ivx){
          Float_t sdist = sqrt(pow(TrackerData.trkstartx[iTrk]-VertexData.vtx[ivx][0],2)+
                               pow(TrackerData.trkstarty[iTrk]-VertexData.vtx[ivx][1],2)+
                               pow(TrackerData.trkstartz[iTrk]-VertexData.vtx[ivx][2],2));
          Float_t edist = sqrt(pow(TrackerData.trkendx[iTrk]-VertexData.vtx[ivx][0],2)+
                               pow(TrackerData.trkendy[iTrk]-VertexData.vtx[ivx][1],2)+
                               pow(TrackerData.trkendz[iTrk]-VertexData.vtx[ivx][2],2));
          if (sdist<minsdist){
            minsdist = sdist;
            if (minsdist<10) TrackerData.trksvtxid[iTrk] = ivx;
          }
          if (edist<minedist){
            minedist = edist;
            if (minedist<10) TrackerData.trkevtxid[iTrk] = ivx;
          }
        }
        
        // find particle ID info
        art::FindMany<anab::ParticleID> fmpid(trackListHandle[iTracker], evt, fParticleIDModuleLabel[iTracker]);
        if(fmpid.isValid()) {
          
          std::vector<const anab::ParticleID*> pids = fmpid.at(iTrk);
          if(pids.size() > TrackerData.GetMaxPlanesPerTrack(iTrk)) {
            mf::LogError("AnalysisTree:limits")
              << "the " << fTrackModuleLabel[iTracker] << " track #" << iTrk
              << " has " << pids.size() 
              << " set of ParticleID variables. Only " << TrackerData.GetMaxPlanesPerTrack(iTrk) << " stored in the tree";
          }
          for (size_t ipid = 0; ipid < pids.size(); ++ipid){
            if (!pids[ipid]->PlaneID().isValid) continue;
            int planenum = pids[ipid]->PlaneID().Plane;
            if (planenum<0||planenum>2) continue;
            TrackerData.trkpidpdg[iTrk][planenum] = pids[ipid]->Pdg();
            TrackerData.trkpidchi[iTrk][planenum] = pids[ipid]->MinChi2();
            TrackerData.trkpidchipr[iTrk][planenum] = pids[ipid]->Chi2Proton();
            TrackerData.trkpidchika[iTrk][planenum] = pids[ipid]->Chi2Kaon();
            TrackerData.trkpidchipi[iTrk][planenum] = pids[ipid]->Chi2Pion();
            TrackerData.trkpidchimu[iTrk][planenum] = pids[ipid]->Chi2Muon();
            TrackerData.trkpidpida[iTrk][planenum] = pids[ipid]->PIDA();
          }
        } // fmpid.isValid()
      
        art::FindMany<anab::Calorimetry> fmcal(trackListHandle[iTracker], evt, fCalorimetryModuleLabel[iTracker]);
        if (fmcal.isValid()){
          std::vector<const anab::Calorimetry*> calos = fmcal.at(iTrk);
          if (calos.size() > TrackerData.GetMaxPlanesPerTrack(iTrk)) {
            // if you get this message, there is probably a bug somewhere since
            // the calorimetry planes should be 3.
            mf::LogError("AnalysisTree:limits")
              << "the " << fTrackModuleLabel[iTracker] << " track #" << iTrk
              << " has " << calos.size() << " planes for calorimetry , only "
              << TrackerData.GetMaxPlanesPerTrack(iTrk) << " stored in tree";
          }
          for (size_t ical = 0; ical<calos.size(); ++ical){
            if (!calos[ical]) continue;
            if (!calos[ical]->PlaneID().isValid) continue;
            int planenum = calos[ical]->PlaneID().Plane;
            if (planenum<0||planenum>2) continue;
            TrackerData.trkke[iTrk][planenum]    = calos[ical]->KineticEnergy();
            TrackerData.trkrange[iTrk][planenum] = calos[ical]->Range();
            //For now make the second argument as 13 for muons. 
            TrackerData.trkpitchc[iTrk][planenum]= calos[ical] -> TrkPitchC();
            const size_t NHits = calos[ical] -> dEdx().size();
            TrackerData.ntrkhits[iTrk][planenum] = (int) NHits;
            if (NHits > TrackerData.GetMaxHitsPerTrack(iTrk, planenum)) {
              // if you get this error, you'll have to increase kMaxTrackHits
              mf::LogError("AnalysisTree:limits")
                << "the " << fTrackModuleLabel[iTracker] << " track #" << iTrk
                << " has " << NHits << " hits on calorimetry plane #" << planenum
                <<", only "
                << TrackerData.GetMaxHitsPerTrack(iTrk, planenum) << " stored in tree";
            }
            if (!isCosmics){
              for(size_t iTrkHit = 0; iTrkHit < NHits && iTrkHit < TrackerData.GetMaxHitsPerTrack(iTrk, planenum); ++iTrkHit) {
                TrackerData.trkdedx[iTrk][planenum][iTrkHit]  = (calos[ical] -> dEdx())[iTrkHit];
                TrackerData.trkdqdx[iTrk][planenum][iTrkHit]  = (calos[ical] -> dQdx())[iTrkHit];
                TrackerData.trkresrg[iTrk][planenum][iTrkHit] = (calos[ical] -> ResidualRange())[iTrkHit];
                const auto& TrkPos = (calos[ical] -> XYZ())[iTrkHit];
                auto& TrkXYZ = TrackerData.trkxyz[iTrk][planenum][iTrkHit];
                TrkXYZ[0] = TrkPos.X();
                TrkXYZ[1] = TrkPos.Y();
                TrkXYZ[2] = TrkPos.Z();
              } // for track hits
            }
          } // for calorimetry info
          if(TrackerData.ntrkhits[iTrk][0] > TrackerData.ntrkhits[iTrk][1] && TrackerData.ntrkhits[iTrk][0] > TrackerData.ntrkhits[iTrk][2]) TrackerData.trkpidbestplane[iTrk] = 0;
          else if(TrackerData.ntrkhits[iTrk][1] > TrackerData.ntrkhits[iTrk][0] && TrackerData.ntrkhits[iTrk][1] > TrackerData.ntrkhits[iTrk][2]) TrackerData.trkpidbestplane[iTrk] = 1;
          else if(TrackerData.ntrkhits[iTrk][2] > TrackerData.ntrkhits[iTrk][0] && TrackerData.ntrkhits[iTrk][2] > TrackerData.ntrkhits[iTrk][1]) TrackerData.trkpidbestplane[iTrk] = 2;
          else if(TrackerData.ntrkhits[iTrk][2] == TrackerData.ntrkhits[iTrk][0] && TrackerData.ntrkhits[iTrk][2] > TrackerData.ntrkhits[iTrk][1]) TrackerData.trkpidbestplane[iTrk] = 2;
          else if(TrackerData.ntrkhits[iTrk][2] == TrackerData.ntrkhits[iTrk][1] && TrackerData.ntrkhits[iTrk][2] > TrackerData.ntrkhits[iTrk][0]) TrackerData.trkpidbestplane[iTrk] = 2;
          else if(TrackerData.ntrkhits[iTrk][1] == TrackerData.ntrkhits[iTrk][0] && TrackerData.ntrkhits[iTrk][1] > TrackerData.ntrkhits[iTrk][2]) TrackerData.trkpidbestplane[iTrk] = 0;
          else if(TrackerData.ntrkhits[iTrk][1] == TrackerData.ntrkhits[iTrk][0] && TrackerData.ntrkhits[iTrk][1] == TrackerData.ntrkhits[iTrk][2]) TrackerData.trkpidbestplane[iTrk] = 2;
        } // if has calorimetry info

        //track truth information
        if (isMC){
          //get the hits on each plane
          art::FindManyP<recob::Hit>      fmht(trackListHandle[iTracker], evt, fTrackModuleLabel[iTracker]);
          std::vector< art::Ptr<recob::Hit> > allHits = fmht.at(iTrk);
          std::vector< art::Ptr<recob::Hit> > hits[kNplanes];

          for(size_t ah = 0; ah < allHits.size(); ++ah){
            if (/* allHits[ah]->WireID().Plane >= 0 && */ // always true
              allHits[ah]->WireID().Plane <  3){
              hits[allHits[ah]->WireID().Plane].push_back(allHits[ah]);
            }
          }
          
          for (size_t ipl = 0; ipl < 3; ++ipl){
            TrackerData.trkidtruth_recoutils_totaltrueenergy[iTrk][ipl] = RecoUtils::TrueParticleIDFromTotalTrueEnergy(clockData, hits[ipl]);
            TrackerData.trkidtruth_recoutils_totalrecocharge[iTrk][ipl] = RecoUtils::TrueParticleIDFromTotalRecoCharge(clockData, hits[ipl]);
            TrackerData.trkidtruth_recoutils_totalrecohits[iTrk][ipl] = RecoUtils::TrueParticleIDFromTotalRecoHits(clockData, hits[ipl]);
            double maxe = 0;
            HitsPurity(clockData, hits[ipl],TrackerData.trkidtruth[iTrk][ipl],TrackerData.trkpurtruth[iTrk][ipl],maxe);
          //std::cout<<"\n"<<iTracker<<"\t"<<iTrk<<"\t"<<ipl<<"\t"<<trkidtruth[iTracker][iTrk][ipl]<<"\t"<<trkpurtruth[iTracker][iTrk][ipl]<<"\t"<<maxe;
            if (TrackerData.trkidtruth[iTrk][ipl]>0){
              const art::Ptr<simb::MCTruth> mc = pi_serv->TrackIdToMCTruth_P(TrackerData.trkidtruth[iTrk][ipl]);
              TrackerData.trkorigin[iTrk][ipl] = mc->Origin();
              const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(TrackerData.trkidtruth[iTrk][ipl]);
              double tote = 0;
              const std::vector<const sim::IDE*> vide(bt_serv->TrackIdToSimIDEs_Ps(TrackerData.trkidtruth[iTrk][ipl]));
              for (auto ide: vide) {
                 tote += ide->energy;
                 TrackerData.trksimIDEenergytruth[iTrk][ipl] = ide->energy;
                 TrackerData.trksimIDExtruth[iTrk][ipl] = ide->x;
                 TrackerData.trksimIDEytruth[iTrk][ipl] = ide->y;
                 TrackerData.trksimIDEztruth[iTrk][ipl] = ide->z;
              }
              TrackerData.trkpdgtruth[iTrk][ipl] = particle->PdgCode();
              TrackerData.trkefftruth[iTrk][ipl] = maxe/(tote/kNplanes); //tote include both induction and collection energies
            //std::cout<<"\n"<<trkpdgtruth[iTracker][iTrk][ipl]<<"\t"<<trkefftruth[iTracker][iTrk][ipl];
            }
          }
        }//end if (isMC)
        // If saving the hierarchy info
        if(fSaveHierarchyInfo[iTracker]){
          trkPfpMap trackPFParticleMap = trackerPFParticleMaps[iTracker];
          trkPfpMapIt it;
          // Check there is a map entry for this vertex
          it = trackPFParticleMap.find(ptrack);
          if(it == trackPFParticleMap.end()) continue;

          art::Ptr<recob::PFParticle> tempParticle = it->second;
          // Get information, find the neutrino and then call its daughters 
          // primary tracks
          TrackerData.trkisprimary[iTrk]   = std::round(static_cast<Int_t>(lar_pandora::LArPandoraHelper::IsFinalState(pfpmap,tempParticle)));
          TrackerData.trkndaughters[iTrk]  = tempParticle->NumDaughters();
          TrackerData.trkpfpid[iTrk]       = tempParticle->Self();
          TrackerData.trkparentpfpid[iTrk] = tempParticle->Parent();
        } // end save hierarchy info
      }//end loop over track
    }//end loop over track module labels
  }// end (fSaveTrackInfo) 
  

  //mc truth information
  if (isMC){
    if (fSaveCryInfo){ 
      //store cry (cosmic generator information) 
      fData->mcevts_truthcry = mclistcry.size();
      fData->cry_no_primaries = nCryPrimaries;
      //fData->cry_no_primaries;
      for(Int_t iPartc = 0; iPartc < mctruthcry->NParticles(); ++iPartc){
        const simb::MCParticle& partc(mctruthcry->GetParticle(iPartc));
        fData->cry_primaries_pdg[iPartc]=partc.PdgCode();
        fData->cry_Eng[iPartc]=partc.E();
        fData->cry_Px[iPartc]=partc.Px();
        fData->cry_Py[iPartc]=partc.Py();
        fData->cry_Pz[iPartc]=partc.Pz();
        fData->cry_P[iPartc]=partc.P();
        fData->cry_StartPointx[iPartc] = partc.Vx();
        fData->cry_StartPointy[iPartc] = partc.Vy();
        fData->cry_StartPointz[iPartc] = partc.Vz();	
        fData->cry_status_code[iPartc]=partc.StatusCode();
        fData->cry_mass[iPartc]=partc.Mass();
        fData->cry_trackID[iPartc]=partc.TrackId();
        fData->cry_ND[iPartc]=partc.NumberDaughters();
        fData->cry_mother[iPartc]=partc.Mother();
      } // for cry particles  
    }// end fSaveCryInfo   

    fData->mcevts_truth = mclist.size();
    //Brailsford 2017/10/16
    //Issue 17917
    //To keep a 1:1 between neutrinos and 'flux' we need the assns
    art::FindOneP<simb::MCFlux> fmFluxNeutrino(mctruthListHandle, evt, fGenieGenModuleLabel);
    // Get GTruth information for scattering code
    art::FindManyP< simb::GTruth > fmgt( mctruthListHandle, evt, fGenieGenModuleLabel );

    if (fData->mcevts_truth > 0){//at least one mc record
    if (fSaveGenieInfo){

      //Brailsford 2017/10/16 
      //Issue 17917
      //Loop over every truth in the spill rather than just looking at the first one.
      //Because MCTruth could be a neutrino OR something else (e.g. cosmics) we are going to have to count up how many neutrinos there are
      fData->mcevts_truth = 0;
      for (unsigned int i_mctruth = 0; i_mctruth < mclist.size(); i_mctruth++){
        //fetch an mctruth
        art::Ptr<simb::MCTruth> curr_mctruth = mclist[i_mctruth];
        //Check if it's a neutrino
        if (!curr_mctruth->NeutrinoSet()) continue;

        // Genie Truth association only for the neutrino
        std::vector< art::Ptr<simb::GTruth> > mcgtAssn = fmgt.at(i_mctruth);
        
        fData->nuScatterCode_truth[i_mctruth] = mcgtAssn[0]->fGscatter;

        fData->nuPDG_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().PdgCode();
        fData->ccnc_truth[i_mctruth] = curr_mctruth->GetNeutrino().CCNC();
        fData->mode_truth[i_mctruth] = curr_mctruth->GetNeutrino().Mode();
        fData->Q2_truth[i_mctruth] = curr_mctruth->GetNeutrino().QSqr();
        fData->W_truth[i_mctruth] = curr_mctruth->GetNeutrino().W();
        fData->hitnuc_truth[i_mctruth] = curr_mctruth->GetNeutrino().HitNuc();
        fData->enu_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().E();
        fData->nuvtxx_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().Vx();
        fData->nuvtxy_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().Vy();
        fData->nuvtxz_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().Vz();
        if (curr_mctruth->GetNeutrino().Nu().P()){
          fData->nu_dcosx_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().Px()/curr_mctruth->GetNeutrino().Nu().P();
          fData->nu_dcosy_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().Py()/curr_mctruth->GetNeutrino().Nu().P();
          fData->nu_dcosz_truth[i_mctruth] = curr_mctruth->GetNeutrino().Nu().Pz()/curr_mctruth->GetNeutrino().Nu().P();
        }
        fData->lep_mom_truth[i_mctruth] = curr_mctruth->GetNeutrino().Lepton().P();
        if (curr_mctruth->GetNeutrino().Lepton().P()){
          fData->lep_dcosx_truth[i_mctruth] = curr_mctruth->GetNeutrino().Lepton().Px()/curr_mctruth->GetNeutrino().Lepton().P();
          fData->lep_dcosy_truth[i_mctruth] = curr_mctruth->GetNeutrino().Lepton().Py()/curr_mctruth->GetNeutrino().Lepton().P();
          fData->lep_dcosz_truth[i_mctruth] = curr_mctruth->GetNeutrino().Lepton().Pz()/curr_mctruth->GetNeutrino().Lepton().P();
        }
        //Brailsford
        //2017/10/17
        //Issue 12918
        //Use the art::Ptr key as the neutrino's unique ID
        fData->nuID_truth[i_mctruth] = curr_mctruth.key();
        //We need to also store N 'flux' neutrinos per event so now check that the FindOneP is valid and, if so, use it!
        if (fmFluxNeutrino.isValid()){
          art::Ptr<simb::MCFlux> curr_mcflux = fmFluxNeutrino.at(i_mctruth);
          fData->tpx_flux[i_mctruth] = curr_mcflux->ftpx;
          fData->tpy_flux[i_mctruth] = curr_mcflux->ftpy;
          fData->tpz_flux[i_mctruth] = curr_mcflux->ftpz;
          fData->tptype_flux[i_mctruth] = curr_mcflux->ftptype;
        }

        //Let's increase the neutrino count
        fData->mcevts_truth++;
      }

      if (mctruth->NeutrinoSet()){
        //Brailsford 2017/10/16
        //Issue 17917
        //Bypass this section of filling code as there can now be multiple neutrinos per TTree::entry stored in the TTree
        /*
        fData->nuPDG_truth = mctruth->GetNeutrino().Nu().PdgCode();
        fData->ccnc_truth = mctruth->GetNeutrino().CCNC();
        fData->mode_truth = mctruth->GetNeutrino().Mode();
        fData->Q2_truth = mctruth->GetNeutrino().QSqr();
        fData->W_truth = mctruth->GetNeutrino().W();
        fData->hitnuc_truth = mctruth->GetNeutrino().HitNuc();
        fData->enu_truth = mctruth->GetNeutrino().Nu().E();
        fData->nuvtxx_truth = mctruth->GetNeutrino().Nu().Vx();
        fData->nuvtxy_truth = mctruth->GetNeutrino().Nu().Vy();
        fData->nuvtxz_truth = mctruth->GetNeutrino().Nu().Vz();
        if (mctruth->GetNeutrino().Nu().P()){
          fData->nu_dcosx_truth = mctruth->GetNeutrino().Nu().Px()/mctruth->GetNeutrino().Nu().P();
          fData->nu_dcosy_truth = mctruth->GetNeutrino().Nu().Py()/mctruth->GetNeutrino().Nu().P();
          fData->nu_dcosz_truth = mctruth->GetNeutrino().Nu().Pz()/mctruth->GetNeutrino().Nu().P();
        }
        fData->lep_mom_truth = mctruth->GetNeutrino().Lepton().P();
        if (mctruth->GetNeutrino().Lepton().P()){
          fData->lep_dcosx_truth = mctruth->GetNeutrino().Lepton().Px()/mctruth->GetNeutrino().Lepton().P();
          fData->lep_dcosy_truth = mctruth->GetNeutrino().Lepton().Py()/mctruth->GetNeutrino().Lepton().P();
          fData->lep_dcosz_truth = mctruth->GetNeutrino().Lepton().Pz()/mctruth->GetNeutrino().Lepton().P();
        }
        //flux information
        art::Ptr<simb::MCFlux>  mcflux = fluxlist[imc];
        fData->tpx_flux = mcflux->ftpx;
        fData->tpy_flux = mcflux->ftpy;
        fData->tpz_flux = mcflux->ftpz;
        fData->tptype_flux = mcflux->ftptype;
        */

        //genie particles information
        fData->genie_no_primaries = mctruth->NParticles();

        size_t StoreParticles = std::min((size_t) fData->genie_no_primaries, fData->GetMaxGeniePrimaries());
        if (fData->genie_no_primaries > (int) StoreParticles) {
          // got this error? it might be a bug,
          // since the structure should have enough room for everything
          mf::LogError("AnalysisTree:limits") << "event has "
            << fData->genie_no_primaries << " MC particles, only "
            << StoreParticles << " stored in tree";
        }
        for(size_t iPart = 0; iPart < StoreParticles; ++iPart){
          const simb::MCParticle& part(mctruth->GetParticle(iPart));
          fData->genie_primaries_pdg[iPart]=part.PdgCode();
          fData->genie_Eng[iPart]=part.E();
          fData->genie_Px[iPart]=part.Px();
          fData->genie_Py[iPart]=part.Py();
          fData->genie_Pz[iPart]=part.Pz();
          fData->genie_P[iPart]=part.P();
          fData->genie_status_code[iPart]=part.StatusCode();
          fData->genie_mass[iPart]=part.Mass();
          fData->genie_trackID[iPart]=part.TrackId();
          fData->genie_ND[iPart]=part.NumberDaughters();
          fData->genie_mother[iPart]=part.Mother();
        } // for particle
      } //if neutrino set
    }// end (fSaveGenieInfo)  

      //GEANT particles information
      if (fSaveGeantInfo){ 
        const sim::ParticleList& plist = pi_serv->ParticleList();
        
        std::string pri("primary");
        int primary=0;
        int active = 0;
        int geant_particle=0;
        sim::ParticleList::const_iterator itPart = plist.begin(),
          pend = plist.end(); // iterator to pairs (track id, particle)
        	  
        for(size_t iPart = 0; (iPart < plist.size()) && (itPart != pend); ++iPart){
          const simb::MCParticle* pPart = (itPart++)->second;
          if (!pPart) {
            throw art::Exception(art::errors::LogicError)
              << "GEANT particle #" << iPart << " returned a null pointer";
          }
          
          ++geant_particle;
          bool isPrimary = pPart->Process() == pri;
          if (isPrimary) ++primary;
          
          int TrackID = pPart->TrackId();

          TVector3 mcstart, mcend;
          double plen = length(*pPart, mcstart, mcend);
          bool isActive = plen != 0;
          if (plen) active++;

          if (iPart < fData->GetMaxGEANTparticles()) {
	   if (pPart->E()>fG4minE||isPrimary){
            fData->process_primary[iPart] = int(isPrimary);
            fData->processname[iPart]= pPart->Process();
            fData->Mother[iPart]=pPart->Mother();
            fData->TrackId[iPart]=TrackID;
            fData->pdg[iPart]=pPart->PdgCode();
            fData->status[iPart] = pPart->StatusCode();
            fData->Eng[iPart]=pPart->E();
	          fData->EndE[iPart]=pPart->EndE();
            fData->Mass[iPart]=pPart->Mass();
            fData->Px[iPart]=pPart->Px();
            fData->Py[iPart]=pPart->Py();
            fData->Pz[iPart]=pPart->Pz();
            fData->P[iPart]=pPart->Momentum().Vect().Mag();
            fData->StartPointx[iPart]=pPart->Vx();
            fData->StartPointy[iPart]=pPart->Vy();
            fData->StartPointz[iPart]=pPart->Vz();
            fData->StartT[iPart] = pPart->T();
            fData->EndPointx[iPart]=pPart->EndPosition()[0];
            fData->EndPointy[iPart]=pPart->EndPosition()[1];
            fData->EndPointz[iPart]=pPart->EndPosition()[2];
            fData->EndT[iPart] = pPart->EndT();
            fData->theta[iPart] = pPart->Momentum().Theta();
            fData->phi[iPart] = pPart->Momentum().Phi();
            fData->theta_xz[iPart] = std::atan2(pPart->Px(), pPart->Pz());
            fData->theta_yz[iPart] = std::atan2(pPart->Py(), pPart->Pz());
            fData->pathlen[iPart]  = plen;
            fData->NumberDaughters[iPart]=pPart->NumberDaughters();
            fData->inTPCActive[iPart] = int(isActive);
            if (isActive){	  
              fData->StartPointx_tpcAV[iPart] = mcstart.X();
              fData->StartPointy_tpcAV[iPart] = mcstart.Y();
              fData->StartPointz_tpcAV[iPart] = mcstart.Z();
              fData->EndPointx_tpcAV[iPart] = mcend.X();
              fData->EndPointy_tpcAV[iPart] = mcend.Y();
              fData->EndPointz_tpcAV[iPart] = mcend.Z();
            }		       
            //Brailsford
            //2017/10/17
            //Issue 17918
            //Get the mother neutrino of this particle and use the hosting art::Ptr key as the matched ID
            const art::Ptr<simb::MCTruth> matched_mctruth = pi_serv->ParticleToMCTruth_P(pPart);
            fData->MotherNuId[iPart] = matched_mctruth.key();
           } 
            //access auxiliary detector parameters
            if (fSaveAuxDetInfo) {
              unsigned short nAD = 0; // number of cells that particle hit
              
              // find deposit of this particle in each of the detector cells
              for (const sim::AuxDetSimChannel* c: fAuxDetSimChannels) {
        	
        	// find if this cell has a contribution (IDE) from this particle,
        	// and which one
        	const std::vector<sim::AuxDetIDE>& setOfIDEs = c->AuxDetIDEs();
        	// using a C++ "lambda" function here; this one:
        	// - sees only TrackID from the current scope
        	// - takes one parameter: the AuxDetIDE to be tested
        	// - returns if that IDE belongs to the track we are looking for
        	std::vector<sim::AuxDetIDE>::const_iterator iIDE
        	  = std::find_if(
        	    setOfIDEs.begin(), setOfIDEs.end(),
        	    [TrackID](const sim::AuxDetIDE& IDE){ return IDE.trackID == TrackID; }
        	  );
        	if (iIDE == setOfIDEs.end()) continue;
        	
        	// now iIDE points to the energy released by the track #i (TrackID)
        	
              // look for IDE with matching trackID
              // find trackIDs stored in setOfIDEs with the same trackID, but negative,
              // this is an untracked particle who's energy should be added as deposited by this original trackID
              float totalE = 0.; // total energy deposited around by the GEANT particle in this cell
              for(const auto& adtracks: setOfIDEs) {
                 if( fabs(adtracks.trackID) == TrackID )
                   totalE += adtracks.energyDeposited;
              } // for
              
              // fill the structure
              if (nAD < kMaxAuxDets) {
                fData->AuxDetID[iPart][nAD] = c->AuxDetID();
                fData->entryX[iPart][nAD]   = iIDE->entryX;
                fData->entryY[iPart][nAD]   = iIDE->entryY;
                fData->entryZ[iPart][nAD]   = iIDE->entryZ;
                fData->entryT[iPart][nAD]   = iIDE->entryT;
                fData->exitX[iPart][nAD]    = iIDE->exitX;
                fData->exitY[iPart][nAD]    = iIDE->exitY;
                fData->exitZ[iPart][nAD]    = iIDE->exitZ;
                fData->exitT[iPart][nAD]    = iIDE->exitT;
                fData->exitPx[iPart][nAD]   = iIDE->exitMomentumX;
                fData->exitPy[iPart][nAD]   = iIDE->exitMomentumY;
                fData->exitPz[iPart][nAD]   = iIDE->exitMomentumZ;
                fData->CombinedEnergyDep[iPart][nAD] = totalE;
              }
              ++nAD;
            } // for aux det sim channels
            fData->NAuxDets[iPart] = nAD; 
            
            if (nAD > kMaxAuxDets) {
              // got this error? consider increasing kMaxAuxDets
              mf::LogError("AnalysisTree:limits") << "particle #" << iPart
                << " touches " << nAD << " auxiliary detector cells, only "
                << kMaxAuxDets << " of them are saved in the tree";
            } // if too many detector cells
          } // if (fSaveAuxDetInfo) 
        }
        else if (iPart == fData->GetMaxGEANTparticles()) {
          // got this error? it might be a bug,
          // since the structure should have enough room for everything
          mf::LogError("AnalysisTree:limits") << "event has "
            << plist.size() << " MC particles, only "
            << fData->GetMaxGEANTparticles() << " will be stored in tree";
        }     
      } // for particles
      
      fData->geant_list_size_in_tpcAV = active;
      fData->no_primaries = primary;
      fData->geant_list_size = geant_particle;
      
      MF_LOG_DEBUG("AnalysisTree") << "Counted "
        << fData->geant_list_size << " GEANT particles ("
        << fData->geant_list_size_in_tpcAV << " in AV), "
        << fData->no_primaries << " primaries, "
        << fData->genie_no_primaries << " GENIE particles";
      
      FillWith(fData->MergedId, 0);

      // helper map track ID => index
      std::map<int, size_t> TrackIDtoIndex;
      const size_t nTrackIDs = fData->TrackId.size();
      for (size_t index = 0; index < nTrackIDs; ++index)
        TrackIDtoIndex.emplace(fData->TrackId[index], index);
      
      // for each particle, consider all the direct ancestors with the same
      // PDG ID, and mark them as belonging to the same "group"
      // (having the same MergedId)
      int currentMergedId = 1;
      for(size_t iPart = fData->geant_list_size; iPart-- > 0; ) {
        // if the particle already belongs to a group, don't bother
        if (fData->MergedId[iPart]) continue;
        // the particle starts its own group
        fData->MergedId[iPart] = currentMergedId;
        
        // look in the ancestry, one by one
        int currentMotherTrackId = fData->Mother[iPart];
        while(currentMotherTrackId > 0) {
          // find the mother (we have its track ID in currentMotherTrackId)
          std::map<int, size_t>::const_iterator iMother
            = TrackIDtoIndex.find(currentMotherTrackId);
          if (iMother == TrackIDtoIndex.end()) break; // no mother found
          size_t currentMotherIndex = iMother->second;
          // if the mother particle is of a different type,
          // don't bother with iPart ancestry any further
          if (fData->pdg[iPart] != fData->pdg[currentMotherIndex]) break;
          
          // group the "current mother" (actually, ancestor) with iPart
          fData->MergedId[currentMotherIndex] = currentMergedId;
          // then consider the grandmother (one generation earlier)
          currentMotherTrackId = fData->Mother[currentMotherIndex];
        } // while ancestry exists
        ++currentMergedId;
      } // for merging check
     } // if (fSaveGeantInfo) 
      
    }//if (mcevts_truth)
  }//if (isMC){
  //fData->taulife = detprop->ElectronLifetime();
  fTree->Fill();
  
  if (mf::isDebugEnabled()) {
    // use mf::LogDebug instead of MF_LOG_DEBUG because we reuse it in many lines;
    // thus, we protect this part of the code with the line above
    mf::LogDebug logStream("AnalysisTreeStructure");
    logStream
      << "Tree data structure contains:"
      << "\n - " << fData->no_hits << " hits (" << fData->GetMaxHits() << ")"
      << "\n - " << fData->genie_no_primaries << " genie primaries (" << fData->GetMaxGeniePrimaries() << ")"
      << "\n - " << fData->geant_list_size << " GEANT particles (" << fData->GetMaxGEANTparticles() << "), "
        << fData->no_primaries << " primaries"
      << "\n - " << fData->geant_list_size_in_tpcAV << " GEANT particles in AV "
      << "\n - " << ((int) fData->kNTracker) << " trackers:"
      ;
    
    size_t iTracker = 0;
    for (auto tracker = fData->TrackData.cbegin();
      tracker != fData->TrackData.cend(); ++tracker, ++iTracker
    ) {
      logStream
         << "\n -> " << tracker->ntracks << " " << fTrackModuleLabel[iTracker]
           << " tracks (" << tracker->GetMaxTracks() << ")"
         ;
      for (int iTrk = 0; iTrk < tracker->ntracks; ++iTrk) {
        logStream << "\n    [" << iTrk << "] "<< tracker->ntrkhits[iTrk][0];
        for (size_t ipl = 1; ipl < tracker->GetMaxPlanesPerTrack(iTrk); ++ipl)
          logStream << " + " << tracker->ntrkhits[iTrk][ipl];
        logStream << " hits (" << tracker->GetMaxHitsPerTrack(iTrk, 0);
        for (size_t ipl = 1; ipl < tracker->GetMaxPlanesPerTrack(iTrk); ++ipl)
          logStream << " + " << tracker->GetMaxHitsPerTrack(iTrk, ipl);
        logStream << ")";
      } // for tracks
    } // for trackers
  } // if logging enabled
  
  // if we don't use a permanent buffer (which can be huge),
  // delete the current buffer, and we'll create a new one on the next event
  if (!fUseBuffer) {
    MF_LOG_DEBUG("AnalysisTreeStructure") << "Freeing the tree data structure";
    DestroyData();
  }
} // sbnd::AnalysisTree::analyze()

void sbnd::AnalysisTree::HitsPurity(detinfo::DetectorClocksData const& clockData,
                                    std::vector< art::Ptr<recob::Hit> > const& hits, Int_t& trackid, Float_t& purity, double& maxe){

  trackid = -1;
  purity = -1;

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  std::map<int,double> trkide;

  for(size_t h = 0; h < hits.size(); ++h){

    art::Ptr<recob::Hit> hit = hits[h];
    std::vector<sim::IDE> ides;
    std::vector<sim::TrackIDE> eveIDs = bt_serv->HitToEveTrackIDEs(clockData, hit);

    for(size_t e = 0; e < eveIDs.size(); ++e){
      //std::cout<<h<<" "<<e<<" "<<eveIDs[e].trackID<<" "<<eveIDs[e].energy<<" "<<eveIDs[e].energyFrac<<std::endl;
      trkide[eveIDs[e].trackID] += eveIDs[e].energy;
    }
  }

  maxe = -1;
  double tote = 0;
  for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
    tote += ii->second;
    if ((ii->second)>maxe){
      maxe = ii->second;
      trackid = ii->first;
    }
  }

  //std::cout << "the total energy of this reco track is: " << tote << std::endl;
  

  if (tote>0){
    purity = maxe/tote;
  }
}

// Calculate distance to boundary.
double sbnd::AnalysisTree::bdist(const recob::Track::Point_t& pos)
{
  // Get geometry.
  art::ServiceHandle<geo::Geometry> geom;

  double d1 = pos.X();                             // Distance to right side (wires).
  double d2 = 2.*geom->DetHalfWidth() - pos.X();   // Distance to left side (cathode).
  double d3 = pos.Y() + geom->DetHalfHeight();     // Distance to bottom.
  double d4 = geom->DetHalfHeight() - pos.Y();     // Distance to top.
  double d5 = pos.Z();                             // Distance to front.
  double d6 = geom->DetLength() - pos.Z();           // Distance to back.

  double result = std::min(std::min(std::min(std::min(std::min(d1, d2), d3), d4), d5), d6);
  return result;
}

// Length of reconstructed track, trajectory by trajectory.
double sbnd::AnalysisTree::length(const recob::Track& track)
{
  return track.Length();
}

// Length of MC particle, trajectory by trajectory.
double sbnd::AnalysisTree::length(const simb::MCParticle& part, TVector3& start, TVector3& end)
{
  // Get geometry.
  art::ServiceHandle<geo::Geometry> geom;

  // Get active volume boundary.
  double xmin = -2.0 * geom->DetHalfWidth() - 1e-8;
  double xmax = 2.0 * geom->DetHalfWidth() + 1e-8;
  double ymin = -geom->DetHalfHeight() -1e-8;
  double ymax = geom->DetHalfHeight() + 1e-8;
  double zmin = 0. -1e-8;
  double zmax = geom->DetLength() + 1e-8;
  
  // Get number traj points
  int n = part.NumberTrajectoryPoints();
  if( n <= 1 ) return 0.;
  
  double  L	= 0.;
  bool	  first	= true; 

  // Loop over points (start with 2nd)
  for(int i = 1; i < n; ++i) {

    TVector3 p1(part.Vx(i),part.Vy(i),part.Vz(i));
    
    if(	  p1.X() >= xmin && p1.X() <= xmax
      &&  p1.Y() >= ymin && p1.Y() <= ymax
      &&  p1.Z() >= zmin && p1.Z() <= zmax ) {
    
      TVector3 p0(part.Vx(i-1),part.Vy(i-1),part.Vz(i-1));
      
      if(first)	start = p1; 
      else L += (p1-p0).Mag();
      first = false;
      end   = p1;
    }
  }

  return L;
}


namespace sbnd{

  DEFINE_ART_MODULE(AnalysisTree)

}

#endif

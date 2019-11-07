#ifndef OPT0FINDER_OPT0FINDERTYPES_H
#define OPT0FINDER_OPT0FINDERTYPES_H

#include <vector>
#include <numeric>
#include "OpT0FinderConstants.h"
#include <string>
namespace flashana {

  /// Index used to identify Flash_t/QPointCollection_t uniquely in an event
  typedef size_t ID_t;
  /// Invalid ID
  const ID_t kINVALID_ID = kINVALID_SIZE;

  /// Enumerator for different types of algorithm
  enum Algorithm_t {
    kTPCFilter,       ///< Algorithm type to filter out TPC objects from matching candidate list
    kFlashFilter,     ///< Algorithm type to filter out flash from matching candidate list
    kFlashMatch,      ///< Algorithm type to match flash hypothesis and reconstructed flash
    kMatchProhibit,   ///< Algorithm type to prohibit a match between a flash and a cluster
    kFlashHypothesis, ///< Algorithm type to make QCluster_t => Flash_t hypothesis
    kCustomAlgo,      ///< Algorithm type that does not play a role in the framework execution but inherits from BaseAlgorithm
    kAlgorithmTypeMax ///< enum flag for algorithm type count & invalid type
  };

  /// Struct to represent an optical flash
  struct Flash_t {
  public:

    std::vector<double> pe_v; ///< PE distribution over photo-detectors
    std::vector<double> pe_err_v; ///< PE value error
    double x,y,z;             ///< Flash position
    double x_err,y_err,z_err; ///< Flash position error
    double time;              ///< Flash timing, a candidate T0
    ID_t idx;                 ///< index from original larlite vector
    /// Default ctor assigns invalid values
    Flash_t() : pe_v() {
      x = y = z = kINVALID_DOUBLE;
      x_err = y_err = z_err = kINVALID_DOUBLE;
      time = kINVALID_DOUBLE;
      idx = kINVALID_ID;
    }
    /// Total PE calcualtion
    double TotalPE() const {
      double res=0.;
      for(auto const& v : pe_v) if(v>=0.) res+=v;
      return res;
    }
    /// Check validity
    bool Valid(size_t nopdet=0) const {
      return (nopdet ? (pe_v.size() == nopdet && pe_err_v.size() == nopdet) : (pe_v.size() == pe_err_v.size()));
    }
    //double TotalPE() const{ return std::accumulate(pe_v.begin(),pe_v.end(),0.0);}
  };

  /// Struct to represent an energy deposition point in 3D space
  struct QPoint_t{

    double x,y,z; ///< Spatial position in [cm]
    double q;     ///< Charge in an arbitrary unit
    /// Default ctor assigns invalid values
    QPoint_t()
      : x(kINVALID_DOUBLE)
      , y(kINVALID_DOUBLE)
      , z(kINVALID_DOUBLE)
      , q(kINVALID_DOUBLE)
      {}
    /// Alternative ctor
    QPoint_t(double xvalue,
             double yvalue,
             double zvalue,
             double qvalue)
      : x(xvalue)
      , y(yvalue)
      , z(zvalue)
      , q(qvalue)
      {}
  };

  /// Collection of charge deposition 3D point (cluster)
  class QCluster_t : public std::vector<QPoint_t>{
  public:
    ID_t idx;     ///< index from original larlite vector
    double time;  ///< assumed time w.r.t. trigger for reconstruction

    /// Default constructor
    QCluster_t() : idx(kINVALID_ID), time(0) {}
    ~QCluster_t() {}

    inline QCluster_t& operator+=(const QCluster_t& rhs) {
      this->reserve(rhs.size() + this->size());
      for(auto const& pt : rhs) this->push_back(pt);
      return (*this);
    }

    inline QCluster_t operator+(const QCluster_t& rhs) const {
      QCluster_t res((*this));
      res += rhs;
      return res;
    }

  };
  /// Collection of 3D point clusters (one use case is TPC object representation for track(s) and shower(s))
  typedef std::vector<flashana::QCluster_t> QClusterArray_t;
  /// Collection of Flash objects
  typedef std::vector<flashana::Flash_t> FlashArray_t;

  /// Index collection
  typedef std::vector<flashana::ID_t> IDArray_t;

  /// Flash-TPC match info
  struct FlashMatch_t {
    ID_t tpc_id;   ///< matched TPC object ID
    ID_t flash_id; ///< matched Flash ID
    double score;  ///< floating point representing the "goodness" (algorithm dependent)
    QPoint_t tpc_point; ///< estimated & matched 3D flash hypothesis point from TPC information
    QPoint_t tpc_point_err; ///< error on the estimated point
    std::vector<double> hypothesis;       ///< Hypothesis flash object
    /// Default ctor assigns invalid values
    FlashMatch_t() : hypothesis()
      { tpc_id = kINVALID_ID; flash_id = kINVALID_ID; score = -1; }
    /// Alternative ctor
    FlashMatch_t(const ID_t& tpc_id_value,
                 const ID_t& flash_id_value,
                 const double& score_value) : hypothesis()
      { tpc_id = tpc_id_value; flash_id = flash_id_value; score = score_value; }
#ifndef __CINT__ // hyde move from fucking CINT cuz it's fucked
    /// Alternative ctor
    FlashMatch_t(const ID_t& tpc_id_value,
                 const ID_t& flash_id_value,
                 const double& score_value,
                 std::vector<double>&& hypo) : hypothesis(std::move(hypo))
      { tpc_id = tpc_id_value; flash_id = flash_id_value; score = score_value; }
#endif
  };

  /// Enum to define MC source type (MCTrack or MCShower) for a given QCluster
  enum MCAncestor_t {
    kMCShowerAncestor,
    kMCTrackAncestor,
    kUnknownAncestor
  };

  /// Struct to represent the ancestor information for a specific interaction (QCluster)
  struct MCSource_t {
    int     index_id; ///< MCTrac/MCShower collection index ID of the ancestor
    double  g4_time;  ///< Interaction G4 time in micro-seconds
    double  energy_deposit;   ///< Deposited energy total
    MCAncestor_t source_type; ///< Ancestor source type
    MCSource_t()
      {
        index_id = -1;
        g4_time  = kINVALID_DOUBLE;
        source_type = kUnknownAncestor;
      }
  };

  namespace msg {
    /// Verbosity message level
    enum Level_t {
      kDEBUG,
      kINFO,
      kNORMAL,
      kWARNING,
      kERROR,
      kEXCEPTION,
      kMSG_TYPE_MAX
    };

    const std::string kStringPrefix[kMSG_TYPE_MAX] =
    {
      "\033[94m     [DEBUG]  \033[00m", ///< DEBUG message prefix
      "\033[92m      [INFO]  \033[00m", ///< INFO message prefix
      "\033[95m    [NORMAL]  \033[00m", ///< NORMAL message prefix
      "\033[93m   [WARNING]  \033[00m", ///< WARNING message prefix
      "\033[91m     [ERROR]  \033[00m", ///< ERROR message prefix
      "\033[5;1;33;41m [EXCEPTION]  \033[00m"  ///< CRITICAL message prefix
    };
    ///< Prefix of message
  }
}
#endif

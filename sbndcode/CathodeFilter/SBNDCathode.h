// SBNDCathode.h
// -----------------------------------------------------------------------------
// Data-driven representation of the SBND cathode (CPA) volume, extracted from
// cosmic-muon track endpoints in the SBND Gen2 run.  The cathode is stored as
// two surfaces on a regular (y, z) grid:
//
//     x_neg(y, z)  -  the -x face of the cathode
//     x_pos(y, z)  -  the +x face of the cathode
//
// A 3D point (x, y, z) [cm] is *inside* the cathode volume when
//
//     Y_MIN <= y <= Y_MAX,  Z_MIN <= z <= Z_MAX,
//     x_neg(y, z) <= x <= x_pos(y, z).
//
// The grids live in an ASCII file written by tools/export_volume.py.  NaN
// cells (no endpoint data in that bin) are treated as "outside" by default;
// use setNanPolicy() to change that.
//
// This header is header-only-safe and has no external deps beyond the C++17
// standard library, so it can be included from any LArSoft/art module.
// -----------------------------------------------------------------------------
#ifndef SBNDCODE_CATHODEFILTER_SBNDCATHODE_H
#define SBNDCODE_CATHODEFILTER_SBNDCATHODE_H

#include <string>
#include <vector>
#include <stdexcept>

#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

namespace sbnd {

class SBNDCathode {
public:
  enum class NanPolicy {
    Outside,   // NaN cell -> contains() returns false (default, safest)
    Inside,    // NaN cell -> treat as a thin slab centred on x = 0
    Throw      // NaN cell -> throw std::runtime_error in contains()
  };

  enum class Interp {
    Nearest,   // nearest-neighbour lookup (fastest)
    Bilinear   // bilinear interpolation in (y, z) (smoother, default)
  };

  // Load a cathode volume file (see tools/export_volume.py for format).
  // Throws std::runtime_error on parse errors.
  explicit SBNDCathode(const std::string& filename);

  // Test whether point (x, y, z) [cm] lies inside the cathode volume.
  // Any coordinate outside the (y, z) bounding box returns false immediately.
  bool contains(double x, double y, double z) const;

  // Get x_neg, x_pos [cm] at (y, z).  Returns NaN for NaN cells.
  double x_neg(double y, double z) const;
  double x_pos(double y, double z) const;

  // Local mesh midplane, 0.5 * (x_neg + x_pos).  NaN if either is NaN.
  double x_mesh(double y, double z) const;

  // Local thickness (x_pos - x_neg).  NaN if either is NaN.
  double thickness(double y, double z) const;

  // Nominal fiducial bounding box.
  double yMin() const { return fYMin; }
  double yMax() const { return fYMax; }
  double zMin() const { return fZMin; }
  double zMax() const { return fZMax; }

  // Grid shape (rows = NY, cols = NZ).
  int nY() const { return fNY; }
  int nZ() const { return fNZ; }

  // Configuration knobs.
  void setNanPolicy(NanPolicy p) { fNanPolicy = p; }
  void setInterp(Interp i)       { fInterp = i; }
  NanPolicy nanPolicy() const    { return fNanPolicy; }
  Interp    interp()    const    { return fInterp;    }

  // Summary string for logs.
  std::string describe() const;

private:
  int    fNY = 0, fNZ = 0;
  double fYMin = 0, fYMax = 0;
  double fZMin = 0, fZMax = 0;
  double fBinCm = 0;
  double fDy = 0, fDz = 0;

  std::vector<double> fXNeg;
  std::vector<double> fXPos;

  NanPolicy fNanPolicy = NanPolicy::Outside;
  Interp    fInterp    = Interp::Bilinear;

  double sampleNeg(double y, double z) const;
  double samplePos(double y, double z) const;
  double sample(const std::vector<double>& grid, double y, double z) const;
  bool   inBBox(double y, double z) const;
  static bool isnan_(double v);
};

} // namespace sbnd

#endif // SBNDCODE_CATHODEFILTER_SBNDCATHODE_H

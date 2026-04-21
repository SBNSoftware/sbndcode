// SBNDCathode.cxx
// Implementation of sbnd::SBNDCathode.  See header.
#include "sbndcode/CathodeFilter/SBNDCathode.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>

namespace sbnd {
namespace {

std::string nextNonComment(std::istream& is) {
  std::string line;
  while (std::getline(is, line)) {
    auto p = line.find_first_not_of(" \t\r\n");
    if (p == std::string::npos) continue;
    if (line[p] == '#') continue;
    return line;
  }
  return {};
}

double parseDouble(const std::string& s) {
  if (s == "nan" || s == "NaN" || s == "NAN")
    return std::numeric_limits<double>::quiet_NaN();
  return std::stod(s);
}

} // anon

bool SBNDCathode::isnan_(double v) {
  return std::isnan(v);
}

SBNDCathode::SBNDCathode(const std::string& filename) {
  std::ifstream f(filename);
  if (!f.is_open())
    throw std::runtime_error("SBNDCathode: cannot open file " + filename);

  auto expect = [&](const std::string& key, auto& out) {
    std::string line = nextNonComment(f);
    std::istringstream iss(line);
    std::string k; iss >> k;
    if (k != key)
      throw std::runtime_error(
        "SBNDCathode: expected key '" + key + "', got '" + k + "'");
    iss >> out;
  };

  int version = 0;
  expect("VERSION", version);
  if (version != 1)
    throw std::runtime_error(
      "SBNDCathode: unsupported VERSION " + std::to_string(version));
  expect("NY",    fNY);
  expect("NZ",    fNZ);
  expect("Y_MIN", fYMin);
  expect("Y_MAX", fYMax);
  expect("Z_MIN", fZMin);
  expect("Z_MAX", fZMax);
  expect("BIN_CM", fBinCm);

  fDy = (fYMax - fYMin) / fNY;
  fDz = (fZMax - fZMin) / fNZ;

  auto readGrid = [&](const std::string& tag, std::vector<double>& grid) {
    std::string line = nextNonComment(f);
    if (line.find(tag) == std::string::npos)
      throw std::runtime_error(
        "SBNDCathode: expected grid tag '" + tag + "', got '" + line + "'");
    grid.assign(fNY * fNZ, 0.0);
    for (int iy = 0; iy < fNY; ++iy) {
      std::string row;
      do { std::getline(f, row); } while (row.empty() || row[0] == '#');
      std::istringstream iss(row);
      for (int iz = 0; iz < fNZ; ++iz) {
        std::string tok;
        if (!(iss >> tok))
          throw std::runtime_error(
            "SBNDCathode: short row in " + tag +
            " at iy=" + std::to_string(iy));
        grid[iy * fNZ + iz] = parseDouble(tok);
      }
    }
  };

  readGrid("X_NEG", fXNeg);
  readGrid("X_POS", fXPos);
}

bool SBNDCathode::inBBox(double y, double z) const {
  return (y >= fYMin && y <= fYMax && z >= fZMin && z <= fZMax);
}

double SBNDCathode::sample(const std::vector<double>& grid,
                           double y, double z) const {
  double fy = (y - fYMin) / fDy - 0.5;
  double fz = (z - fZMin) / fDz - 0.5;

  if (fInterp == Interp::Nearest) {
    int iy = static_cast<int>(std::floor(fy + 0.5));
    int iz = static_cast<int>(std::floor(fz + 0.5));
    iy = std::clamp(iy, 0, fNY - 1);
    iz = std::clamp(iz, 0, fNZ - 1);
    return grid[iy * fNZ + iz];
  }

  int iy0 = static_cast<int>(std::floor(fy));
  int iz0 = static_cast<int>(std::floor(fz));
  int iy1 = iy0 + 1;
  int iz1 = iz0 + 1;
  double ty = fy - iy0;
  double tz = fz - iz0;
  iy0 = std::clamp(iy0, 0, fNY - 1);
  iy1 = std::clamp(iy1, 0, fNY - 1);
  iz0 = std::clamp(iz0, 0, fNZ - 1);
  iz1 = std::clamp(iz1, 0, fNZ - 1);

  const double v00 = grid[iy0 * fNZ + iz0];
  const double v01 = grid[iy0 * fNZ + iz1];
  const double v10 = grid[iy1 * fNZ + iz0];
  const double v11 = grid[iy1 * fNZ + iz1];

  auto pick = [&]() {
    double vs[4] = {v00, v01, v10, v11};
    int n = 0; double mean = 0.0;
    for (double v : vs) if (std::isfinite(v)) { mean += v; ++n; }
    return n ? mean / n : std::numeric_limits<double>::quiet_NaN();
  };
  if (isnan_(v00) || isnan_(v01) || isnan_(v10) || isnan_(v11)) return pick();

  return (1 - ty) * ((1 - tz) * v00 + tz * v01)
       +      ty  * ((1 - tz) * v10 + tz * v11);
}

double SBNDCathode::sampleNeg(double y, double z) const { return sample(fXNeg, y, z); }
double SBNDCathode::samplePos(double y, double z) const { return sample(fXPos, y, z); }

double SBNDCathode::x_neg(double y, double z) const {
  if (!inBBox(y, z)) return std::numeric_limits<double>::quiet_NaN();
  return sampleNeg(y, z);
}
double SBNDCathode::x_pos(double y, double z) const {
  if (!inBBox(y, z)) return std::numeric_limits<double>::quiet_NaN();
  return samplePos(y, z);
}
double SBNDCathode::x_mesh(double y, double z) const {
  const double n = x_neg(y, z);
  const double p = x_pos(y, z);
  if (isnan_(n) || isnan_(p)) return std::numeric_limits<double>::quiet_NaN();
  return 0.5 * (n + p);
}
double SBNDCathode::thickness(double y, double z) const {
  const double n = x_neg(y, z);
  const double p = x_pos(y, z);
  if (isnan_(n) || isnan_(p)) return std::numeric_limits<double>::quiet_NaN();
  return p - n;
}

bool SBNDCathode::contains(double x, double y, double z) const {
  if (!inBBox(y, z)) return false;

  const double xn = sampleNeg(y, z);
  const double xp = samplePos(y, z);

  if (isnan_(xn) || isnan_(xp)) {
    switch (fNanPolicy) {
      case NanPolicy::Outside: return false;
      case NanPolicy::Inside : {
        const double half = 0.5 * 3.2;
        return (x >= -half && x <= +half);
      }
      case NanPolicy::Throw  :
        throw std::runtime_error(
          "SBNDCathode::contains: no data at (y,z) = ("
          + std::to_string(y) + ", " + std::to_string(z) + ")");
    }
  }
  return (x >= xn && x <= xp);
}

std::string SBNDCathode::describe() const {
  std::ostringstream os;
  os << "SBNDCathode: " << fNY << " x " << fNZ
     << " grid  y in [" << fYMin << "," << fYMax << "] cm,"
     << "  z in [" << fZMin << "," << fZMax << "] cm,"
     << "  bin = " << fBinCm << " cm,"
     << "  interp = " << (fInterp == Interp::Bilinear ? "bilinear" : "nearest")
     << ",  NaN policy = "
     << (fNanPolicy == NanPolicy::Outside ? "outside"
       : fNanPolicy == NanPolicy::Inside  ? "inside" : "throw");
  return os.str();
}

} // namespace sbnd

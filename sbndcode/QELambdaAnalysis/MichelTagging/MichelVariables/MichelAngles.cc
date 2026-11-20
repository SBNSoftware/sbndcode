#include "MichelAngles.hh"


inline double getAngle(
    const double x,
    const double y,
    const double centerX,
    const double centerY,
    const double dirX,
    const double dirY)
{
    double vx = x - centerX;
    double vy = y - centerY;

    double norm = std::hypot(vx, vy);

    if (norm == 0) return -1;

    double dirNorm = std::hypot(dirX, dirY);

    double cosTheta = std::abs(vx * dirX + vy * dirY) / (norm * dirNorm);
    double angle = std::acos(cosTheta) * 180.0 / M_PI;

    if (angle > 90.0) {
        angle = 180.0 - angle;
    }

    return angle;
}

double quantile(const std::vector<double>& x, const double q)
{
    if (x.empty()) return -1;

    std::vector<double> sorted = x;
    std::sort(sorted.begin(), sorted.end());

    double pos = q * (sorted.size() - 1);

    size_t i = static_cast<size_t>(std::floor(pos));
    size_t j = static_cast<size_t>(std::ceil(pos));

    if (i == j)
        return sorted[i];

    double t = pos - i;

    return (1.0 - t) * sorted[i] + t * sorted[j]; 
}


double MichelAngles::AnglesQ80()
{

    for (const MuHit& hit: fmuHits) {

        if (hit.TpcID != fbragg.TpcID) continue;

        double dist = std::hypot(hit.ChannelID - fbragg.ChannelID, hit.PeakTime - fbragg.PeakTime);
        if (dist > fmaxRadius) continue;

        double angle = getAngle(hit.ChannelID, hit.PeakTime, fbragg.ChannelID, fbragg.PeakTime, fpca.dirX, fpca.dirY);
        if (angle >= 0) fangles.push_back(angle);
    }

    if (fangles.empty()) return -1;

    return quantile(fangles, 0.8);
}




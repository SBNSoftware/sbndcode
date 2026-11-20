#include "MichelRadial.hh"

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
    if (norm == 0) return 0.0;

    double dirNorm = std::hypot(dirX, dirY);

    double cosTheta = (vx * dirX + vy * dirY) / (norm * dirNorm);
    double angle = std::acos(std::clamp(cosTheta, -1.0, 1.0)) * 180.0 / M_PI;

    // Rango completo 0°-180°, sin plegar a 90°.
    return angle;
}



std::vector<MuHit> expandHit(const MuHit& hit, int nSubHits)
{
    std::vector<MuHit> outHits;

    std::vector<double> offsets(nSubHits);

    for (int i = 0; i < nSubHits; ++i) {
        double q = 0.1 + 0.8 * i / (nSubHits - 1);
        offsets[i] = ROOT::Math::normal_quantile(q, 1.0);
    }

    for (int i = 0; i < nSubHits; ++i) {
        MuHit h;

        h.PlaneID = hit.PlaneID;
        h.ChannelID = hit.ChannelID;
        h.PeakTime = hit.PeakTime + offsets[i] * hit.RMS;
        h.Integral = hit.Integral / nSubHits;

        outHits.push_back(h);
    }

    return outHits;
}

std::pair<std::vector<double>, std::vector<double>> MichelRadial::getChargeProfile(
    const std::vector<MuHit>& hits)
{
    std::vector<double> rVec;
    std::vector<double> chargeVec;

    double increment = (fParams.MaxRadius-0.1) / fParams.NRadialIntervals;
    double last_r = 0;
    double charge = 0;

    for (double r = 0.1; r < fParams.MaxRadius; r += increment)
    {   
        rVec.push_back(r);

        for (auto& hit : hits)
        {
            double dist = std::hypot(hit.ChannelID - fbragg.ChannelID, hit.PeakTime - fbragg.PeakTime);
            if ((dist > last_r) && (dist <= r)) {
                charge += hit.Integral;
            }
        }

        chargeVec.push_back(charge);
        last_r = r;
    }

    return {rVec, chargeVec};
}

std::vector<double> gradient(const std::vector<double>& x,
                             const std::vector<double>& y)
{
    size_t n = y.size();
    std::vector<double> g(n);

    if (n < 2)
        return g;

    g[0] = (y[1] - y[0]) / (x[1] - x[0]);
    g[n-1] = (y[n-1] - y[n-2]) / (x[n-1] - x[n-2]);

    for (size_t i = 1; i < n-1; ++i)
    {
        g[i] = (y[i+1] - y[i-1]) /
               (x[i+1] - x[i-1]);
    }

    return g;
}


std::vector<double> smoothGradient(const std::vector<double>& x,
                                   const std::vector<double>& y)
{
    std::vector<double> grad = gradient(x, y);

    TGraph g(grad.size(), x.data(), grad.data());

    TGraphSmooth smoother;
    TGraph* gs = smoother.SmoothKern(&g, "normal", 2.2);

    std::vector<double> gradSmooth(gs->GetN());

    double xx, yy;
    for (int i = 0; i < gs->GetN(); ++i)
    {
        gs->GetPoint(i, xx, yy);
        gradSmooth[i] = yy;
    }

    return gradSmooth;
}

struct DropIdxResult {
    int idx;
    int thrIdx;
    double cutoff;
    bool dropFound;
};

DropIdxResult getfdropIdx(const std::vector<double>& vec)
{   
    double max = *std::max_element(vec.begin(), vec.end());
    double min = std::min(100.0, max / 25.0);
    int n = vec.size();
    int idx = n - 1;
    int riseIdx = n;
    bool dropFound = false;

    // No charge observed
    if (max < 10){
        return {0, 0, 0.0, false};
    }

    for (int i = 0; i < n; ++i){
        if (vec[i] > min){
            riseIdx = i;
            break;
        }
    }

    if (riseIdx > n/8){
        return {0, 0, 0.0, false};
    }

    for (int i = riseIdx; i < static_cast<int>(vec.size()); ++i)
    {
        if (vec[i] < min)
        {
            idx = i;
            dropFound = true;
            break;
        }
    }

    return {static_cast<int>(idx), static_cast<int>(riseIdx), min, dropFound};
}


void MichelRadial::Compute() 
{

    std::vector<MuHit> selectedHits;
    double rDrop = -1;

    for (const MuHit& hit: fmuHits) {

        if (hit.TpcID != fbragg.TpcID) continue;

        double dist = std::hypot(hit.ChannelID - fbragg.ChannelID, hit.PeakTime - fbragg.PeakTime);
        if (dist > fParams.MaxRadius) continue;

        double angle = getAngle(hit.ChannelID, hit.PeakTime, fbragg.ChannelID, fbragg.PeakTime, fpca.dirX, fpca.dirY);
        if (angle <= fParams.AngleCut) continue;

        std::vector<MuHit> expandedHits = expandHit(hit, fParams.NSubHits);
        selectedHits.insert(selectedHits.end(), expandedHits.begin(), expandedHits.end());
    }

    auto radialInfo = getChargeProfile(selectedHits);
    fradii = radialInfo.first;
    fcharge = radialInfo.second;
    
    fdQdr = gradient(fradii, fcharge);
    fsmoothdQdr = smoothGradient(fradii, fcharge);

    DropIdxResult dropInfo = getfdropIdx(fsmoothdQdr);
    fdropIdx = dropInfo.idx;
    fdropThr = dropInfo.thrIdx;
    fmindQdr = dropInfo.cutoff;
    fdropFound = dropInfo.dropFound;

    if (fdropIdx == 0){
        fdropRadius = 0;
        return;
    }

    if (!fdropFound){
        fdropRadius = fradii.back();
        return;
    }

    if (fdropIdx > 0)
    {
        double r0 = fradii[fdropIdx];
        double dqdr0 = fsmoothdQdr[fdropIdx];
        double r1 = fradii[fdropIdx - 1];
        double dqdr1 = fsmoothdQdr[fdropIdx - 1];

        double denom = dqdr1 - dqdr0;

        if (std::abs(denom) > 1e-12)
        {
            rDrop = r0 + (fmindQdr - dqdr0) * (r1 - r0) / denom;
            rDrop = std::clamp(rDrop, fradii.front(), fradii.back());
        }
        else
        {
            rDrop = r0;
        }
    }

    fdropRadius = rDrop;

}
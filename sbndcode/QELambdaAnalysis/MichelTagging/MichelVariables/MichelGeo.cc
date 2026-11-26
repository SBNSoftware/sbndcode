#include "MichelGeo.hh"


std::vector<double> slidingMeans(
    const std::vector<MuCalo>& caloPoints, 
    int windowSize, 
    double maxdEdx, 
    double maxResRange
    )
{
    std::vector<double> means;
    int n = caloPoints.size();
    int w = std::min(windowSize, n);

    if (n == 0 || w == 0) return means;

    for (int i = 0; i <= n-w; i++){
        if (caloPoints[i+w].ResRange > maxResRange) break;

        double sum = 0.0;
        for (int j = i; j < i + w; j++){
            if (caloPoints[j].dEdx > maxdEdx) continue;

            sum += caloPoints[j].dEdx;
        }

        means.push_back(sum / w);
    }

    return means;
}


void MichelGeo::findBraggPeak()
{
    if (fmuCalos.empty()) {
        fbragg.Found = false;
        return;
    }

    std::vector<double> means = slidingMeans(fmuCalos, fParams.SlidingWindowSize, fParams.MaxDEDx, fParams.MaxResRange);

    if (means.empty()) {
        fbragg.Found = false;
        return;
    }

    auto it = std::max_element(means.begin(), means.end());
    double braggdEdx = *it;
    std::size_t braggIdx = std::distance(means.begin(), it);

    // Fill bragg information

    fbragg.Found = true;
    fbragg.Idx = braggIdx;
    fbragg.ChannelID = fmuCalos[braggIdx].ChannelID;
    fbragg.PeakTime = fmuCalos[braggIdx].PeakTime;
    fbragg.TpcID = fmuCalos[braggIdx].TpcID;
    fbragg.ResRange = fmuCalos[braggIdx].ResRange;
    fbragg.dEdx = braggdEdx;
    fbragg.X = fmuCalos[braggIdx].X;
    fbragg.Y = fmuCalos[braggIdx].Y;
    fbragg.Z = fmuCalos[braggIdx].Z;

    return;

}

bool firstPCADirection(
    const std::vector<double>& x, 
    const std::vector<double>& y,
    const double centerX,
    const double centerY,
    double& dirX,
    double& dirY)
{
    
    if (x.size() != y.size() || x.size() < 2) 
        return false;

    const size_t n = x.size();

    // Average x and y
    double meanX = std::accumulate(x.begin(), x.end(), 0.0) / n;
    double meanY = std::accumulate(y.begin(), y.end(), 0.0) / n;

    // Covariance matrix elements
    double covXX = 0.0, covXY = 0.0, covYY = 0.0;

    for (size_t i = 0; i < n; ++i) {
        double dx = x[i] - meanX;
        double dy = y[i] - meanY;
        covXX += dx * dx;
        covXY += dx * dy;
        covYY += dy * dy;
    }

    covXX /= (n-1);
    covXY /= (n-1);
    covYY /= (n-1);

    // Covariance matrix
    TMatrixDSym cov(2);
    cov(0,0) = covXX;
    cov(0,1) = covXY;
    cov(1,0) = covXY;
    cov(1,1) = covYY;

    // PCA
    TMatrixDSymEigen eigen(cov);

    TVectorD values = eigen.GetEigenValues();
    TMatrixD vectors = eigen.GetEigenVectors();

    int idx = (values[0] > values[1]) ? 0 : 1;

    dirX = vectors(0, idx);
    dirY = vectors(1, idx);

    double vx = meanX - centerX;
    double vy = meanY - centerY;

    if (dirX * vx + dirY * vy < 0) {
        dirX = -dirX;
        dirY = -dirY;
    }

    return true;

}


void MichelGeo::getMuonDirection() 
{

    if (!fbragg.Found) {
        fpcaResult.Valid = false;
        return;
    }

    std::vector<double> preBraggCh, preBraggT;

    for (int i=fbragg.Idx; i<static_cast<int>(fmuCalos.size()); i++){
        if (fmuCalos[i].ResRange > fParams.MaxResRange) break;
        if (fmuCalos[i].TpcID != fbragg.TpcID) continue;

        preBraggCh.push_back(fmuCalos[i].ChannelID);
        preBraggT.push_back(fmuCalos[i].PeakTime);
    }

    if (firstPCADirection(preBraggCh, preBraggT, fbragg.ChannelID, fbragg.PeakTime, fpcaResult.dirX, fpcaResult.dirY)) {
        fpcaResult.Valid = true;
    } else {
        fpcaResult.Valid = false;
        return;
    }
}


// MAIN FUNCTION

void MichelGeo::Compute()
{
    findBraggPeak();

    if (fbragg.Found) {
        getMuonDirection();
    } else {
        fpcaResult.Valid = false;
    }
}
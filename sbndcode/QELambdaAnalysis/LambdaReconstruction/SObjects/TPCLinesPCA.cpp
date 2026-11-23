////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesPCA.cpp
//
// \brief Definition of TPCLinesPCA
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCLinesPCA.h"

double mean(const std::vector<SHit>& data, double (SHit::*attribute)() const) {
    double sum = 0.0;
    for (const SHit& s : data) {
        sum += (s.*attribute)();
    }
    return sum / data.size();
}

struct LinearRegressionResult {
    double slope;
    double intercept;
};

LinearRegressionResult linearRegression(const std::vector<SHit>& data) {
    double sumX = 0.0;
    double sumY = 0.0;
    double sumXY = 0.0;
    double sumX2 = 0.0;

    for (const SHit& point : data) {
        sumX += point.X();
        sumY += point.Y();
        sumXY += point.X() * point.Y();
        sumX2 += point.X() * point.X();
    }

    const size_t n = data.size();
    double slope = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    double intercept = (sumY - slope * sumX) / n;

    return {slope, intercept};
}

double calculateR_squared(const std::vector<SHit>& data, const LinearRegressionResult& regression) {
    double totalSumOfSquares = 0.0;
    double residualSumOfSquares = 0.0;

    double yMean = 0.0;

    for (const SHit& point : data) {
        yMean += point.Y();
    }

    yMean /= data.size();

    for (const SHit& point : data) {
        double predictedY = regression.slope * point.X() + regression.intercept;
        totalSumOfSquares += (point.Y() - yMean) * (point.Y() - yMean);
        residualSumOfSquares += (point.Y() - predictedY) * (point.Y() - predictedY);
    }

    return 1.0 - (residualSumOfSquares / totalSumOfSquares);
}


double pearsonCorrelation(const std::vector<SHit>& data) {
    double meanX = 0.0;
    double meanY = 0.0;
    double sumXY = 0.0;
    double sumX2 = 0.0;
    double sumY2 = 0.0;

    for (const SHit& s : data) {
        meanX+=s.X();
        meanY+=s.Y();
    }
    meanX/=data.size();
    meanY/=data.size();

    for (const SHit& s : data) {
        double x = s.X();
        double y = s.Y();
        sumXY += (x - meanX) * (y - meanY);
        sumX2 += std::pow(x - meanX, 2);
        sumY2 += std::pow(y - meanY, 2);
    }

    return sumXY / (std::sqrt(sumX2) * std::sqrt(sumY2));
}


// Function to calculate the mean of a vector of SHit objects
SHit calculateMean(const std::vector<SHit>& data) {
    double sumX = 0.0;
    double sumY = 0.0;
    for (const SHit& point : data) {
        sumX += point.X();
        sumY += point.Y();
    }
    return SHit(-1, sumX / data.size(), sumY / data.size());
}

// Function to calculate SSR (Sum of Squares for Regression)
double calculateSSR(const std::vector<SHit>& data, double slope, double intercept) {
    SHit mean = calculateMean(data);
    double ssr = 0.0; // Regression sum of squares
    
    for (const SHit& point : data) {
        double predictedY = slope * point.X() + intercept;
        ssr += (predictedY - mean.Y()) * (predictedY - mean.Y());
    }
    
    return ssr;
}

// Function to calculate the closest distance of a point to the line
double calculateClosestDistance(const SHit& point, double slope, double intercept) {
    // Calculate the perpendicular distance between the point and the line
    double distance = std::abs(point.Y() - (slope * point.X() + intercept)) / std::sqrt(1 + slope * slope);
    return distance;
}

// Function to calculate Mean Squared Error (MSE) using closest distances
double calculateMSE(const std::vector<SHit>& data, double slope, double intercept) {
    double sum = 0.0;
    for (const SHit& point : data) {
        double closestDist = calculateClosestDistance(point, slope, intercept);
        sum += closestDist * closestDist;
    }
    return sum / data.size();
}

// Function to calculate Mean Absolute Error (MAE) for a linear regression model
double calculateMAE(const std::vector<SHit>& data, double slope, double intercept) {
    double sum = 0.0;
    for (const SHit& point : data) {
        double predictedY = slope * point.X() + intercept;
        sum += std::abs(point.Y() - predictedY);
    }
    return sum / data.size();
}



double calculateAverageWidth(const std::vector<SHit>& data) {
    double sum = 0.0;
    for (const SHit& point : data) {
        sum += point.Width();
    }
    return sum / data.size();
}


LineEquation TPCLinesPCA::PerformPCA2D(std::vector<SHit>& data) {

    // first center the data
    double meanX = 0.0, meanY = 0.0;
    for (SHit& h : data) {
        meanX += h.X();
        meanY += h.Y();
    }
    meanX /= data.size();
    meanY /= data.size();

    // get the covariance and fill centered vectors
    double covXX = 0.0, covYY = 0.0, covXY = 0.0;
    for (SHit& h : data) {
        double x = h.X()-meanX;
        double y = h.Y()-meanY;
        covXX += x * x;
        covYY += y * y;
        covXY += x * y;
    }
    covXX /= data.size();
    covYY /= data.size();
    covXY /= data.size();

    double eigenvalue1 = (covXX + covYY + std::sqrt((covXX - covYY) * (covXX - covYY) + 4 * covXY * covXY)) / 2.0;
    double eigenvectorX1 = eigenvalue1 - covYY;
    double eigenvectorY1 = covXY;

    /* double eigenvalue2 = (covXX + covYY - std::sqrt((covXX - covYY) * (covXX - covYY) + 4 * covXY * covXY)) / 2.0;
    double eigenvectorX2 = eigenvalue2 - covYY;
    double eigenvectorY2 = covXY;*/



    // Calculate slopes (m) and intercepts (b) for lines
    double slope1 = eigenvectorY1 / eigenvectorX1;
    double yIntercept1 = meanY - slope1 * meanX; // Since the line passes through the mean

    double mae = calculateMAE(data, slope1, yIntercept1);

    /*// Total variance (sum of covariances)
    double total_variance = covXX + covYY;
    double explained_variance = eigenvalue1 / total_variance;

    LinearRegressionResult LinRegResult = linearRegression(data);
    double RSquared = calculateR_squared(data, LinRegResult);
   
    double Pearson = pearsonCorrelation(data);
    double Residual = calculateSSR(data, slope1, yIntercept1);
    Residual = Residual/std::sqrt(Residual);
    double mae = calculateMAE(data, slope1, yIntercept1);
    double mse = calculateMSE(data, slope1, yIntercept1);
    
    std::cout << "\nLinear regression Slope: " << LinRegResult.slope << "Intercept: " << LinRegResult.intercept << std::endl;
    std::cout<<"Exaplined variance "<<explained_variance<<" Pearson: "<<Pearson<<" SSR:"<<Residual<<" RSq:"<<RSquared<<std::endl;
    std::cout << "MAE: " << mae << " MSE: " << mse << std::endl;*/

    double correlation = mae/calculateAverageWidth(data);

    /*std::cout << "Eigenvalues: " << eigenvalue1 << ", " << eigenvalue2 << std::endl;
    std::cout << "Eigenvectors: (" << eigenvectorX1 << ", " << eigenvectorY1 << "), (" << eigenvectorX2 << ", " << eigenvectorY2 << ")" << std::endl;
    double slope2 = eigenvectorY2 / eigenvectorX2;
    double yIntercept2 = 0.0; // Since the line passes through the origin (0, 0)*/
    /*std::cout << "Slope-Intercept Form of Lines:" << std::endl;
    std::cout << "Line 1: y = " << slope1 << "x + " << yIntercept1 << std::endl;*/

    return LineEquation(slope1, yIntercept1, correlation);
}


LineEquation TPCLinesPCA::PerformPCA2DThreshold(std::vector<SHit> &hitList, float th, bool reverse) {
    if(th > 1) th = 1;

    std::vector<SHit> points = hitList;

    // Order points by x
    if(reverse==true)
        std::sort(points.begin(), points.end(), [](SHit& a, SHit& b) {return a.X() > b.X();});
    else
        std::sort(points.begin(), points.end(), [](SHit& a, SHit& b) {return a.X() < b.X();});
    
    // Get first selected point
    points.resize(static_cast<size_t>(th * points.size()));

    return PerformPCA2D(points);
}

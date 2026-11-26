#include "MichelChi2.hh"

double getSumSquaredResiduals(const std::vector<double>& x,
                              const std::vector<double>& y,
                              TF1* func)
{
    if (x.size() != y.size())
        return -1.0;

    double sum = 0.0;

    for (size_t i = 0; i < x.size(); ++i)
    {
        double yPred = func->Eval(x[i]);
        double residual = y[i] - yPred;

        sum += residual * residual;
    }

    return sum;
}

double fitConstant(
    const std::vector<double>& x,
    const std::vector<double>& y)
{

    if (x.size() != y.size() || x.size() < 2) 
        return -1;

    TGraph gr(x.size(), x.data(), y.data());

    TF1 f("f", "[0]", x.front(), x.back());

    gr.Fit(&f, "Q");

    return getSumSquaredResiduals(x, y, &f);
}

double fitExp(
    const std::vector<double>& x,
    const std::vector<double>& y)
{

    if (x.size() != y.size() || x.size() < 3) 
        return -1;

    TGraph gr(x.size(), x.data(), y.data());

    TF1 f("f", "[0]*exp([1]*x)", x.front(), x.back());
    f.SetParameters(y.front(), -0.1);

    gr.Fit(&f, "Q");

    return getSumSquaredResiduals(x, y, &f);
}

double MichelChi2::LowerChi2Ratio()
{
    std::vector<double> x, y;

    for (const auto& calo : fmuCalos) {
        if (calo.ResRange >= fbragg.ResRange) break;
        if (calo.ResRange > fParams.MaxResRange) break;
        if (calo.dEdx >= fParams.MaxDEDx) continue;

        x.push_back(calo.ResRange);
        y.push_back(calo.dEdx);
    }

    if (x.size() < 3) return -1;

    double chi2Const = fitConstant(x, y);
    double chi2Exp = fitExp(x, y);

    return (chi2Exp > 0) ? chi2Exp / chi2Const : -1;
}

double MichelChi2::AllChi2Ratio()
{
    std::vector<double> x, y;

    for (const auto& calo : fmuCalos) {
        if (calo.ResRange >= fParams.MaxResRange) break;
        if (calo.dEdx > fParams.MaxDEDx) continue;

        x.push_back(calo.ResRange);
        y.push_back(calo.dEdx);
    }

    if (x.size() < 3) return -1;

    double chi2Const = fitConstant(x, y);
    double chi2Exp = fitExp(x, y);

    return (chi2Exp > 0) ? chi2Exp / chi2Const : -1;
}
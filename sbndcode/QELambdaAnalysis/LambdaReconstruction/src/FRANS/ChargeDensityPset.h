///////////////////////////////////////////////////////////////////////
/// File: ChargeDensityAlgConf.h
///
/// Created by Fran Nicolas, November 2022
////////////////////////////////////////////////////////////////////////

#ifndef SBND_CHARGEDENSITYPSET_H
#define SBND_CHARGEDENSITYPSET_H

struct FRANSPsetType {
  bool ApplyRawSmoothing;
  bool ApplySmoothing;
  bool ApplyCumulativeSmoothing;
  unsigned int NDriftPack;
  unsigned int NWirePack;
  float ExpoAvSmoothPar;
  int UnAvNeighbours;
  double CumulativeCut;
  int SlidingWindowN;
  int NSamplesBeginSlope;
  int MaxRadius;
  bool UseHitWidth;
  int Verbose;
  bool CalculateScore;
  bool UseAlpha;
  bool UseIota;
  std::string TMVAFilename;
  std::string OutputPath;

  // Constructor to initialize the members
  FRANSPsetType(bool rawSmoothing=false,
                bool smoothing=false,
                bool cumulativeSmoothing=false,
                unsigned int driftPack=4,
                unsigned int wirePack=1,
                float smoothParam=0.3,
                int unweightedNeighbors=1,
                double cut=0.8,
                int windowN=3,
                int samplesSlope=3,
                int radius=70,
                bool useHitWidth=true,
                int verbose=0,
                bool calculate=false,
                bool useAlpha=false,
                bool useIota=false,
                const std::string& filename="",
                const std::string& outputpath="" )
                  : ApplyRawSmoothing(rawSmoothing),
                    ApplySmoothing(smoothing),
                    ApplyCumulativeSmoothing(cumulativeSmoothing),
                    NDriftPack(driftPack),
                    NWirePack(wirePack),
                    ExpoAvSmoothPar(smoothParam),
                    UnAvNeighbours(unweightedNeighbors),
                    CumulativeCut(cut),
                    SlidingWindowN(windowN),
                    NSamplesBeginSlope(samplesSlope),
                    MaxRadius(radius),
                    UseHitWidth(useHitWidth),
                    Verbose(verbose),
                    CalculateScore(calculate),
                    UseAlpha(useAlpha),
                    UseIota(useIota),
                    TMVAFilename(filename),
                    OutputPath(outputpath)
  {}

  void Print() const {
    std::cout << "ApplyRawSmoothing: " << ApplyRawSmoothing << std::endl;
    std::cout << "ApplySmoothing: " << ApplySmoothing << std::endl;
    std::cout << "ApplyCumulativeSmoothing: " << ApplyCumulativeSmoothing << std::endl;
    std::cout << "NDriftPack: " << NDriftPack << std::endl;
    std::cout << "NWirePack: " << NWirePack << std::endl;
    std::cout << "ExpoAvSmoothPar: " << ExpoAvSmoothPar << std::endl;
    std::cout << "UnAvNeighbours: " << UnAvNeighbours << std::endl;
    std::cout << "CumulativeCut: " << CumulativeCut << std::endl;
    std::cout << "SlidingWindowN: " << SlidingWindowN << std::endl;
    std::cout << "NSamplesBeginSlope: " << NSamplesBeginSlope << std::endl;
    std::cout << "MaxRadius: " << MaxRadius << std::endl;
    std::cout << "UseHitWidth: " << UseHitWidth << std::endl;
    std::cout << "Verbose: " << Verbose << std::endl;
    std::cout << "CalculateScore: " << CalculateScore << std::endl;
    std::cout << "UseAlpha: " << UseAlpha << std::endl;
    std::cout << "UseIota: " << UseIota << std::endl;
    std::cout << "TMVAFilename: " << TMVAFilename << std::endl;
    std::cout << "OutputPath: " << OutputPath << std::endl;
  }
};



#endif

////////////////////////////////////////////////////////////////////////////
//
// \file TPCLineParameters.h
//
// \brief Definition of TPCLinesParameters
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_PARSERCONFIG_PARAMETERS_H
#define TPC_PARSERCONFIG_PARAMETERS_H

#if LAMBDAANA_LARSOFT == 1
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/TPCLinesAlgo/TPCLinesParameters.h"
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/FRANS/ChargeDensityPset.h"
#else
#include "ChargeDensityPset.h"
#include "TPCLinesParameters.h"
#endif

#include <boost/property_tree/info_parser.hpp>
#include <filesystem>
#include <fstream>
#include <map>

namespace fs = std::filesystem;

std::string FindFile(const std::string& fileName, const char* envVarName="LAMBDAANA_SEARCHPATH") {
    
    const char* envVarValue = std::getenv(envVarName);

    std::string fullPath="";

    if (envVarValue!=nullptr) {
        
        std::istringstream paths(envVarValue);
        std::string path;
        while (std::getline(paths, path, ':')) { 
            fs::path directoryPath(path);

            if (fs::is_directory(directoryPath)) {
                for (const auto& entry : fs::directory_iterator(directoryPath)) {
                    
                    if (fs::is_regular_file(entry) && entry.path().filename() == fileName) {
                        fullPath = entry.path().string(); // Return the full path of the file.
                    }
                }
            }
        }
    }
    else{
        std::cerr << "PATH environmental variable not found." << std::endl;
        return "";
    }

    return fullPath;
}


boost::property_tree::ptree GetPropertyTreeFromFileName(std::string filename, std::string blockName){
    
    // Define the section names
    const std::string begin_prolog = "BEGIN_PROLOG";
    const std::string end_prolog = "END_PROLOG";

    // Create a property tree
    boost::property_tree::ptree pt;

    // Open the configuration file
    std::ifstream config_file(filename);
    if (config_file.is_open()) {

        std::string line;
        bool inside_prolog = false;
        bool inside_block = false;

        while (std::getline(config_file, line)) {
            if (line == begin_prolog){
                inside_prolog = true;
            }
            else if (line == end_prolog) {
                inside_prolog = false;
            }
            else if (line ==  blockName) {
                inside_block = true;
            }
            else if (inside_prolog && inside_block) {
                
                if (!line.empty() && line[0] != '#') {  // Ignore empty lines and comments
                    size_t separator_pos = line.find(":");
                    if (separator_pos != std::string::npos) {
                        std::string key = line.substr(0, separator_pos);
                        std::string value = line.substr(separator_pos + 1);
                        
                        // Remove spaces
                        key.erase(std::remove_if(key.begin(), key.end(), ::isspace), key.end());
                        value.erase(std::remove_if(value.begin(), value.end(), ::isspace), value.end());

                        // Add the key-value pair to the property tree
                        pt.put(key, value);
                    }
                }
                
            }
        }

        config_file.close();

    }



    return pt;
}


FRANSPsetType ReadFRANSPset(std::string filename, std::string blockName){

    // Create the parameter set
    FRANSPsetType fPsetFRANS;
    
    // Create the boost property and read all the parameter in the given file
    boost::property_tree::ptree pt = GetPropertyTreeFromFileName(filename, blockName);

    // Fill the parameter set
    fPsetFRANS.ApplyRawSmoothing = pt.get<bool>("ApplyRawSmoothing");
    fPsetFRANS.ApplySmoothing = pt.get<bool>("ApplySmoothing");
    fPsetFRANS.ApplyCumulativeSmoothing = pt.get<bool>("ApplyCumulativeSmoothing");
    fPsetFRANS.NDriftPack = pt.get<unsigned int>("NDriftPack");
    fPsetFRANS.NWirePack = pt.get<unsigned int>("NWirePack");
    fPsetFRANS.ExpoAvSmoothPar = pt.get<double>("ExpoAvSmoothPar");
    fPsetFRANS.UnAvNeighbours = pt.get<int>("UnAvNeighbours");
    fPsetFRANS.CumulativeCut = pt.get<double>("CumulativeCut");
    fPsetFRANS.SlidingWindowN = pt.get<int>("SlidingWindowN");
    fPsetFRANS.NSamplesBeginSlope = pt.get<int>("NSamplesBeginSlope");
    fPsetFRANS.MaxRadius = pt.get<int>("MaxRadius");
    fPsetFRANS.UseHitWidth = pt.get<bool>("UseHitWidth");
    fPsetFRANS.Verbose = pt.get<int>("Verbose");
    fPsetFRANS.CalculateScore = pt.get<bool>("CalculateScore");
    fPsetFRANS.UseAlpha = pt.get<bool>("UseAlpha");
    fPsetFRANS.UseIota = pt.get<bool>("UseIota");
    fPsetFRANS.TMVAFilename = pt.get<std::string>("TMVAFilename");
    fPsetFRANS.TMVAFilename = fPsetFRANS.TMVAFilename.substr(1, fPsetFRANS.TMVAFilename.size()-2);
    fPsetFRANS.OutputPath = pt.get<std::string>("OutputPath");
    
    return fPsetFRANS;
}


TrackFinderAlgorithmPsetType ReadTrackFinderAlgorithmPset(std::string filename, std::string blockName){

    // Create the parameter set
    TrackFinderAlgorithmPsetType fPsetTrackFinder;
    
    // Create the boost property and read all the parameter in the given file
    boost::property_tree::ptree pt = GetPropertyTreeFromFileName(filename, blockName);

    // Fill the parameter set
    fPsetTrackFinder.MaxDTube = pt.get<int>("MaxDTube");
    fPsetTrackFinder.MaxDCluster = pt.get<double>("MaxDCluster");
    fPsetTrackFinder.SingleWireMode = pt.get<bool>("SingleWireMode");
    fPsetTrackFinder.MinClusterHits = pt.get<int>("MinClusterHits");
    fPsetTrackFinder.DCleaning = pt.get<double>("DCleaning");
    fPsetTrackFinder.ClusterCompletenessCut = pt.get<double>("ClusterCompletenessCut");
    fPsetTrackFinder.ClusterAngleCut = pt.get<double>("ClusterAngleCut");
    fPsetTrackFinder.CaptureMissingHits = pt.get<bool>("CaptureMissingHits");
    fPsetTrackFinder.MinTrackHits = pt.get<int>("MinTrackHits");
    fPsetTrackFinder.HitDensityThreshold = pt.get<float>("HitDensityThreshold");
    fPsetTrackFinder.UseCompactness = pt.get<bool>("UseCompactness");
    fPsetTrackFinder.ConnectednessTol = pt.get<double>("ConnectednessTol");
    fPsetTrackFinder.ConnectednessWidthTol = pt.get<double>("ConnectednessWidthTol");
    fPsetTrackFinder.CompactnessTol = pt.get<double>("CompactnessTol");
    fPsetTrackFinder.Verbose = pt.get<int>("Verbose");
    
    
    return fPsetTrackFinder;
}


VertexFinderAlgorithmPsetType ReadVertexFinderAlgorithmPset(std::string filename, std::string blockName){

    // Create the parameter set
    VertexFinderAlgorithmPsetType fPsetVertexFinder;
    
    // Create the boost property and read all the parameter in the given file
    boost::property_tree::ptree pt = GetPropertyTreeFromFileName(filename, blockName);

    // Fill the parameter set
    fPsetVertexFinder.MaxDistToEdge = pt.get<double>("MaxDistToEdge");
    fPsetVertexFinder.RefineVertexIntersection = pt.get<bool>("RefineVertexIntersection");
    fPsetVertexFinder.SlidingWindowN = pt.get<int>("SlidingWindowN");
    fPsetVertexFinder.UseEdgesDiscard = pt.get<bool>("UseEdgesDiscard");
    fPsetVertexFinder.MaxTrackFractionInMain = pt.get<float>("MaxTrackFractionInMain");
    fPsetVertexFinder.DecideMainTrack = pt.get<bool>("DecideMainTrack");
    fPsetVertexFinder.AddCollinearLines = pt.get<bool>("AddCollinearLines");
    fPsetVertexFinder.VertexDistanceROI = pt.get<float>("VertexDistanceROI");
    fPsetVertexFinder.MinWires = pt.get<int>("MinWires");
    fPsetVertexFinder.AngleTolerance = pt.get<float>("AngleTolerance");
    fPsetVertexFinder.TriangleInequalityTol = pt.get<float>("TriangleInequalityTol");
    fPsetVertexFinder.VertexHitsTol = pt.get<float>("VertexHitsTol");
    fPsetVertexFinder.VertexHitsMinHits = pt.get<int>("VertexHitsMinHits");
    fPsetVertexFinder.VertexCompactnessTol = pt.get<float>("VertexCompactnessTol");
    fPsetVertexFinder.MinTrackOccupancy = pt.get<float>("MinTrackOccupancy");
    fPsetVertexFinder.MinTrackGoodness = pt.get<float>("MinTrackGoodness");
    fPsetVertexFinder.MakeCalorimetry = pt.get<bool>("MakeCalorimetry");
    fPsetVertexFinder.Verbose = pt.get<int>("Verbose");
    
    
    return fPsetVertexFinder;
}

HoughAlgorithmPsetType ReadHoughAlgorithmPset(std::string filename, std::string blockName){

    // Create the parameter set
    HoughAlgorithmPsetType fPsetHough;
    
    // Create the boost property and read all the parameter in the given file
    boost::property_tree::ptree pt = GetPropertyTreeFromFileName(filename, blockName);

    // Fill the parameter set
    fPsetHough.MaxRadiusLineHypothesis = pt.get<double>("MaxRadiusLineHypothesis");
    fPsetHough.ThetaRes = pt.get<double>("ThetaRes");
    fPsetHough.MaxDistanceTube = pt.get<double>("MaxDistanceTube");
    fPsetHough.MinHoughHits = pt.get<int>("MinHoughHits");
    fPsetHough.Verbose = pt.get<int>("Verbose");
    
    
    return fPsetHough;
}

TPCLinesAlgoPsetType ReadTPCLinesAlgoPset(std::string filename, std::string blockName){

    // Create the parameter set
    TPCLinesAlgoPsetType fPsetTPCLines;
    
    // Create the boost property and read all the parameter in the given file
    boost::property_tree::ptree pt = GetPropertyTreeFromFileName(filename, blockName);

    // Fill the parameter set
    fPsetTPCLines.MaxRadius = pt.get<double>("MaxRadius");
    fPsetTPCLines.DriftConversion = pt.get<double>("DriftConversion");
    fPsetTPCLines.MaxMeanOccupation = pt.get<double>("MaxMeanOccupation");
    fPsetTPCLines.MaxHoughTracks = pt.get<int>("MaxHoughTracks");
    fPsetTPCLines.MinTrackHits = pt.get<int>("MinTrackHits");
    fPsetTPCLines.RemoveIsolatedHits = pt.get<bool>("RemoveIsolatedHits");
    fPsetTPCLines.MaxNeighbourDistance = pt.get<double>("MaxNeighbourDistance");
    fPsetTPCLines.MinNeighboursHits = pt.get<int>("MinNeighboursHits");
    fPsetTPCLines.MinTrackGoodness = pt.get<float>("MinTrackGoodness");
    fPsetTPCLines.CustomKinkPoint = pt.get<bool>("CustomKinkPoint");
    fPsetTPCLines.ParallelTracksMaxD = pt.get<int>("ParallelTracksMaxD");
    fPsetTPCLines.VertexAlgorithm = pt.get<int>("VertexAlgorithm");
    fPsetTPCLines.View = pt.get<int>("View");
    fPsetTPCLines.OutputPath = pt.get<std::string>("OutputPath");
    fPsetTPCLines.Verbose = pt.get<int>("Verbose");
    fPsetTPCLines.DebugMode = pt.get<int>("DebugMode");
    
    
    return fPsetTPCLines;
}

CaloAlgorithmPsetType ReadCaloAlgorithmPset(std::string filename, std::string blockName){

    // Create the parameter set
    CaloAlgorithmPsetType fPsetCalo;

    // Create the boost property and read all the parameter in the given file
    boost::property_tree::ptree pt = GetPropertyTreeFromFileName(filename, blockName);

    // Fill the parameter set
    fPsetCalo.MaxTrackHits = pt.get<int>("MaxTrackHits");
    fPsetCalo.MakeFit = pt.get<bool>("MakeFit");
    fPsetCalo.UseHitRMSInFit = pt.get<bool>("UseHitRMSInFit");
    fPsetCalo.HitWidthTolInFit = pt.get<double>("HitWidthTolInFit");
    fPsetCalo.Verbose = pt.get<int>("Verbose");

    return fPsetCalo;
}

#endif // TPC_LINES_PARAMETERS_H
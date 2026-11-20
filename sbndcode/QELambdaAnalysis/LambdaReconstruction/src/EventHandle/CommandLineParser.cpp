#include "CommandLineParser.h"


CommandLineParser::CommandLineParser(int argc, char* argv[]) {
    
    // Initialize default values
    Debug = 0;
    ConfPsetPath = "";
    DebugMode = -1;
    n = 1e6;
    nskip = -1;
    event = -1;
    sr = -1;
    file_name = "";
    directory_path = ".";
    ext = ".root";
    vertexOption = 0;
    treeName = "ana";
    plotFRANS = false;
    view = 2;
    threeViews = false;
    
    // Loop through command-line arguments
    for (int i = 1; i < argc; ++i) {
        std::string argument = argv[i];

        // Check for specific flags or labels with single hyphen prefix
        if (argument == "-d") {
            if (i + 1 < argc) {
                Debug = std::stoi(argv[i + 1]);
                i++; // Skip the next argument since it has been processed
            }
        } else if (argument == "-c") {
            if (i + 1 < argc) {
                ConfPsetPath = argv[i + 1];
                i++;
            }
        } else if (argument == "-m") {
            if (i + 1 < argc) {
                DebugMode = std::stoi(argv[i + 1]);
                i++;
            }
        } else if (argument == "-n") {
            if (i + 1 < argc) {
                n = std::stoi(argv[i + 1]);
                i++;
            }
        } else if (argument == "-nskip") {
            if (i + 1 < argc) {
                nskip = std::stoi(argv[i + 1]);
                i++;
            }
        } else if (argument == "-e") {
            if (i + 1 < argc) {
                event = std::stoi(argv[i + 1]);
                i++;
            }
        } else if (argument == "-sr") {
            if (i + 1 < argc) {
                sr = std::stoi(argv[i + 1]);
                i++;
            }
        } else if (argument == "-s") {
            if (i + 1 < argc) {
                file_name = argv[i + 1];
                i++;
            }
        } else if (argument == "-dir") {
            if (i + 1 < argc) {
                directory_path = argv[i + 1];
                i++;
            }
        } else if (argument == "-ext") {
            if (i + 1 < argc) {
                ext = argv[i + 1];
                i++;
            }
        } else if (argument == "-vx") {
            if (i + 1 < argc) {
                vertexOption = std::stoi(argv[i + 1]);
                i++;
            }
        } else if (argument == "-t") {
            if (i + 1 < argc) {
                treeName = argv[i + 1];
                i++;
            }
        } else if (argument == "-frans") {
            if (i + 1 < argc) {
                plotFRANS = true;
                i++;
            }
        }
        else if (argument == "-v") {
            if (i + 1 < argc) {
                view = std::stoi(argv[i + 1]);
                i++;
            }
        }
        else if (argument == "-3v") {
            if (i + 1 < argc) {
                threeViews = true;
                i++;
            }
        }



    }
}
#include <iostream>
#include <string>
#include <vector>

class FRANSAnaCommandLineParser {

    public:
        FRANSAnaCommandLineParser(int argc, char* argv[]);

        std::string getFileName() const { return fFileName; }

    private:
        std::string fFileName;

};



FRANSAnaCommandLineParser::FRANSAnaCommandLineParser(int argc, char* argv[]) {
    
    // Initialize default values
    fFileName = "";
    
    // Loop through command-line arguments
    for (int i = 1; i < argc; ++i) {
        std::string argument = argv[i];

        // Check for specific flags or labels with single hyphen prefix
        if (argument == "-s") {
            if (i + 1 < argc) {
                fFileName = argv[i + 1];
                i++;
            }
        }
    }


}
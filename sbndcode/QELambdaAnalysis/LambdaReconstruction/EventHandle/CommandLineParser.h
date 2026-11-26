#include <iostream>
#include <string>
#include <vector>

class CommandLineParser {

    public:
        CommandLineParser(int argc, char* argv[]);

        int getDebug() const { return Debug; }
        std::string getPsetPath() const { return ConfPsetPath; }
        int getDebugMode() const { return DebugMode; }
        int getN() const { return n; }
        int getNskip() const { return nskip; }
        int getEvent() const { return event; }
        int getSr() const { return sr; }
        std::string getFileName() const { return file_name; }
        std::string getDirectoryPath() const { return directory_path; }
        std::string getExtension() const { return ext; }
        int getVertexOption() const { return vertexOption; }
        std::string getTreeName() const { return treeName; }
        bool getPlotFRANS() const { return plotFRANS; }
        int getView() const { return view; }
        bool getThreeViews() const { return threeViews; }

    private:
        int Debug;
        std::string ConfPsetPath;
        int DebugMode;
        int n;
        int nskip;
        int event;
        int sr;
        std::string file_name;
        std::string directory_path;
        std::string ext;
        int vertexOption;
        std::string treeName;
        bool plotFRANS;
        int view;
        bool threeViews;

};

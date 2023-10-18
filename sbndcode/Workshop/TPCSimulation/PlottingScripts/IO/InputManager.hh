#ifndef INPUTMANAGER_HH
#define INPUTMANAGER_HH

#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <iostream>
#include <string>

class InputManager{
protected:
  TFile*      f_Input;
  std::string treename;
  std::string filename;

public:
  TTree* t_Input;

  InputManager(){};
  virtual ~InputManager(){};

  int GetEntries() const {
    return t_Input->GetEntries();
  };

  virtual void GetEntry (const int i=0) {
    t_Input->GetEntry(i);
  };

  std::string GetInputFile()                       const { return filename; };
  std::string GetInputTree()                       const { return treename; };
  void        SetInputFile(const std::string s="")       { filename=s; };
  void        SetInputTree(const std::string s="")       { treename=s; };

  virtual void LoadTree() {};
protected:
  void Initialise() {
    if (filename == "") {
      std::cerr << "Need to provide an input file" << std::endl;
      exit(1);
    }

    f_Input = new TFile(filename.c_str(), "READ");
    if (!f_Input->IsOpen()) {
      std::cerr << "The file " << filename.c_str() << " does not exist." << std::endl;
      exit(1);
    }

    if (f_Input->Get(treename.c_str())) {
      t_Input = (TTree*)f_Input->Get(treename.c_str());
    } else {
      std::cerr << "The requested tree (" << treename.c_str() << ") does not exist in file " << file\
name << std::endl;
      exit(1);
    }
  }
};

class FullGeoAnaInputManager: public InputManager {

public:
  int Run   ;
  int SubRun;
  int Event ;

  std::vector<int>         * True_Bck_PDG          ;
  std::vector<int>         * True_Bck_Mother       ;
  std::vector<int>         * True_Bck_EndProcess   ;
  std::vector<int>         * True_Bck_ID           ;
  std::vector<double>      * True_Bck_VertX        ;
  std::vector<double>      * True_Bck_VertY        ;
  std::vector<double>      * True_Bck_VertZ        ;
  std::vector<double>      * True_Bck_Time         ;
  std::vector<double>      * True_Bck_Energy       ;
  std::vector<double>      * True_Bck_EndE         ;
  std::vector<double>      * True_Bck_EndX         ;
  std::vector<double>      * True_Bck_EndY         ;
  std::vector<double>      * True_Bck_EndZ         ;
  std::vector<double>      * True_Bck_EndT         ;
  std::vector<std::string> * True_Bck_StartMaterial;
  std::vector<std::string> * True_Bck_EndMaterial  ;

  std::vector<int>   * Hit_View;
  std::vector<int>   * Hit_Chan;
  std::vector<int>   * Hit_TPC ;
  std::vector<int>   * Hit_Size;
  std::vector<float> * Hit_Time;

  ~FullGeoAnaInputManager() {
    delete True_Bck_PDG                 ; True_Bck_PDG                  = NULL;
    delete True_Bck_Mother              ; True_Bck_Mother               = NULL;
    delete True_Bck_EndProcess          ; True_Bck_EndProcess           = NULL;
    delete True_Bck_ID                  ; True_Bck_ID                   = NULL;
    delete True_Bck_VertX               ; True_Bck_VertX                = NULL;
    delete True_Bck_VertY               ; True_Bck_VertY                = NULL;
    delete True_Bck_VertZ               ; True_Bck_VertZ                = NULL;
    delete True_Bck_Time                ; True_Bck_Time                 = NULL;
    delete True_Bck_Energy              ; True_Bck_Energy               = NULL;
    delete True_Bck_EndE                ; True_Bck_EndE                 = NULL;
    delete True_Bck_EndX                ; True_Bck_EndX                 = NULL;
    delete True_Bck_EndY                ; True_Bck_EndY                 = NULL;
    delete True_Bck_EndZ                ; True_Bck_EndZ                 = NULL;
    delete True_Bck_EndT                ; True_Bck_EndT                 = NULL;
    delete True_Bck_StartMaterial       ; True_Bck_StartMaterial        = NULL;
    delete True_Bck_EndMaterial         ; True_Bck_EndMaterial          = NULL;

    delete Hit_View              ; Hit_View = NULL;
    delete Hit_Chan              ; Hit_Chan = NULL;
    delete Hit_TPC               ; Hit_TPC  = NULL;
    delete Hit_Size              ; Hit_Size = NULL;
    delete Hit_Time              ; Hit_Time = NULL;

    if (f_Input) f_Input->Close();
  };

  FullGeoAnaInputManager():
    Run   (0),
    SubRun(0),
    Event (0),

    True_Bck_PDG                (NULL),
    True_Bck_Mother             (NULL),
    True_Bck_EndProcess         (NULL),
    True_Bck_ID                 (NULL),
    True_Bck_VertX              (NULL),
    True_Bck_VertY              (NULL),
    True_Bck_VertZ              (NULL),
    True_Bck_Time               (NULL),
    True_Bck_Energy             (NULL),
    True_Bck_EndE               (NULL),
    True_Bck_EndX               (NULL),
    True_Bck_EndY               (NULL),
    True_Bck_EndZ               (NULL),
    True_Bck_EndT               (NULL),
    True_Bck_StartMaterial      (NULL),
    True_Bck_EndMaterial        (NULL),

    Hit_View                    (NULL),
    Hit_Chan                    (NULL),
    Hit_TPC                     (NULL),
    Hit_Size                    (NULL),
    Hit_Time                    (NULL) {

    f_Input = NULL;
    t_Input = NULL;
    filename = "";
    treename = "";
  };

  void LoadTree() {
    Initialise();
    t_Input->SetBranchAddress("Run"   , &Run   );
    t_Input->SetBranchAddress("SubRun", &SubRun);
    t_Input->SetBranchAddress("Event" , &Event );

    if (t_Input->GetListOfBranches()->FindObject("True_Bck_PDG")) {
      t_Input->SetBranchAddress("True_Bck_PDG"                  , &True_Bck_PDG                 );
      t_Input->SetBranchAddress("True_Bck_Mother"               , &True_Bck_Mother              );
      t_Input->SetBranchAddress("True_Bck_EndProcess"           , &True_Bck_EndProcess          );
      t_Input->SetBranchAddress("True_Bck_ID"                   , &True_Bck_ID                  );
      t_Input->SetBranchAddress("True_Bck_VertX"                , &True_Bck_VertX               );
      t_Input->SetBranchAddress("True_Bck_VertY"                , &True_Bck_VertY               );
      t_Input->SetBranchAddress("True_Bck_VertZ"                , &True_Bck_VertZ               );
      t_Input->SetBranchAddress("True_Bck_Time"                 , &True_Bck_Time                );
      t_Input->SetBranchAddress("True_Bck_Energy"               , &True_Bck_Energy              );
      t_Input->SetBranchAddress("True_Bck_EndE"                 , &True_Bck_EndE                );
      t_Input->SetBranchAddress("True_Bck_EndX"                 , &True_Bck_EndX                );
      t_Input->SetBranchAddress("True_Bck_EndY"                 , &True_Bck_EndY                );
      t_Input->SetBranchAddress("True_Bck_EndZ"                 , &True_Bck_EndZ                );
      t_Input->SetBranchAddress("True_Bck_EndT"                 , &True_Bck_EndT                );
      t_Input->SetBranchAddress("True_Bck_StartMaterial"        , &True_Bck_StartMaterial       );
      t_Input->SetBranchAddress("True_Bck_EndMaterial"          , &True_Bck_EndMaterial         );

      t_Input->SetBranchAddress("Hit_View", &Hit_View);
      t_Input->SetBranchAddress("Hit_Chan", &Hit_Chan);
      t_Input->SetBranchAddress("Hit_TPC" , &Hit_TPC );
      t_Input->SetBranchAddress("Hit_Time", &Hit_Time);
      t_Input->SetBranchAddress("Hit_Size", &Hit_Size);

    }
  };
};

#endif



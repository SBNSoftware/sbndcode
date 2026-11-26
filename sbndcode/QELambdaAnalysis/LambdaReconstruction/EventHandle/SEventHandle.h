////////////////////////////////////////////////////////////////////////////
//
// \file TPCLineParameters.h
//
// \brief Definition of TPCLinesParameters
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_EVENTHANDLE_PARAMETERS_H
#define TPC_EVENTHANDLE_PARAMETERS_H

#include <vector>

// ROOT includes
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TSystemDirectory.h"

#if LAMBDAANA_LARSOFT == 1
#include "sbndcode/QELambdaAnalysis/LambdaReconstruction/SObjects/TPCSimpleEvents.h"
#else
#include "TPCSimpleEvents.h"
#endif


std::vector<TPad*> buildpadcanvas(int nx, int ny);

std::vector<TString> GetInputFileList(std::string file_name,  std::string ext, std::string dirPath=".");

class SEventId {
    public:
        SEventId(int r=0, int sr=0, int e=0):
            fRun(r),
            fSubRun(sr),
            fEvent(e),
            fEventLabel(std::to_string(r)+"_"+std::to_string(sr)+"_"+std::to_string(e))
            {};


        SEventId(const SEventId& other) : SEventId(other.fRun, other.fSubRun, other.fEvent) {}

        std::string Label(){return fEventLabel;};
        int Run(){return fRun;};
        int SubRun(){return fSubRun;};
        int Event(){return fEvent;};

        friend std::ostream& operator<<(std::ostream& os, const SEventId& statusPrinter) {
        os << " ***** R=" << statusPrinter.fRun << " SR=" << statusPrinter.fSubRun << " E=" << statusPrinter.fEvent << " *****" << std::endl;
        return os;
    }

    private:
        int fRun, fSubRun, fEvent;
        std::string fEventLabel;
};

class SEventSelection {
    public:
        SEventSelection(SEventId ev):
            fEvent(ev),
            fNSlices(0),
            fNSelected(0),
            fNNotSelected(0),
            fNSkipped(0)
            {};

        SEventSelection(){
            fEvent = SEventId();
            fNSlices=0;
            fNSelected=0;
            fNNotSelected=0;
            fNSkipped=0;
        }
        
        SEventSelection(const SEventSelection& other) : SEventSelection(other.fEvent) {}

        SEventId EventId(){return fEvent;};

        void AddSelected(){
            fNSlices++;
            fNSelected++;
        };

        void AddSkipped(){
            fNSlices++;
            fNSkipped++;
        }

        void AddNotSelected(){
            fNSlices++;
            fNNotSelected++;
        }

        int NSlices(){return fNSlices;};
        int NSelected(){return fNSelected;};
        int NNotSelected(){return fNNotSelected;};
        int NSkipped(){return fNSkipped;};

    private:
        SEventId fEvent;
        int fNSlices;
        int fNSelected;
        int fNNotSelected;
        int fNSkipped;
};


class EfficiencyCalculator {
public:
    EfficiencyCalculator();

    void AddEvent(SEventId ev){
        if(fEventList.find(ev.Label()) == fEventList.end()){
            SEventSelection sEv(ev);
            fEventList[ev.Label()] = SEventSelection(sEv);
        }
    }

    void UpdateSkipped(SEventId ev){
        AddEvent(ev);
        fEventList[ev.Label()].AddSkipped();
    }

    void UpdateSelected(SEventId ev){
        AddEvent(ev);
        fEventList[ev.Label()].AddSelected();
    }

    void UpdateNotSelected(SEventId ev){
        AddEvent(ev);
        fEventList[ev.Label()].AddNotSelected();

    }

    // return NProcessed events
    int NProcessedEvents(){return nProcessedEvents;};
    // return NSelected events
    int NSelectedEvents(){return nEventsSelected;};
    // return NEvents
    int NEvents(){return nEvents;};


    void UpdateHistograms(SEvent recoEv);

    void UpdateValues();

    friend std::ostream& operator<<(std::ostream& os, EfficiencyCalculator& statusPrinter) {
        statusPrinter.UpdateValues();
        os << " ********** Final status-" << std::endl;
        os << " ... NTotalEvents=" << statusPrinter.nEvents << " NSkipped=" << statusPrinter.nEventsSkipped << std::endl;
        os << " ... NProcessed=" << statusPrinter.nProcessedEvents << " NSelected=" << statusPrinter.nEventsSelected;
        os << " efficiency=" << static_cast<double>(statusPrinter.nEventsSelected) / statusPrinter.nProcessedEvents;
        os << " efficiency all=" << static_cast<double>(statusPrinter.nEventsSelected) / statusPrinter.nEvents << std::endl;
        return os;
    }

    void DrawHistograms(TCanvas *c);


private:
    std::map<std::string, SEventSelection> fEventList;
    int nEvents;
    int nEventsSkipped;
    int nProcessedEvents;
    int nEventsSelected;

    TH1F *hNOrigins;
    TH1F *hNOriginsMult1;
    TH1F *hNOriginsMult2;
    TH1F *hNOriginsMultGt3;
    TH2F *hProf2DMatrix;
    TH1F *hHitDensity;
};


#endif // TPC_LINES_PARAMETERS_H
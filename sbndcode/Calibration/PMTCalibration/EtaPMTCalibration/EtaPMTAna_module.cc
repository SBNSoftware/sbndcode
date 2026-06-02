#include "sbndcode/Calibration/PMTCalibration/EtaPMTCalibration/EtaPMTAna_module.hh"
#include "larcorealg/Geometry/OpDetGeo.h"

// -------- Constructor --------
opdet::EtaPMTAna::EtaPMTAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fOpFlashesModuleLabel( p.get<std::vector<std::string>>("OpFlashesModuleLabel") ),
  fCRTTrackModuleLabel( p.get<std::string>("CRTTrackModuleLabel") ),
  fMaxTimeDiff( p.get<double>("MaxTimeDiff") ),
  fMinTrackLength( p.get<double>("MinTrackLength") ),
  fMakePDSGeoTree( p.get<bool>("MakePDSGeoTree") ),
  fVerbosity( p.get<int>("Verbosity") ), 
  fAVmin_tpc0( p.get<std::vector<double>>("AVminTPC0") ),
  fAVmax_tpc0( p.get<std::vector<double>>("AVmaxTPC0") ),
  fAVmin_tpc1( p.get<std::vector<double>>("AVminTPC1") ),
  fAVmax_tpc1( p.get<std::vector<double>>("AVmaxTPC1") ),
  fIsochronousWidth ( p.get<double>("IsochronousWidth") )
{
}

// -------- Begi job function --------
void opdet::EtaPMTAna::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("OpAnaTree", "Optical Analyzer Tree");

  fTree->Branch("eventID", &_eventID, "eventID/i");
  fTree->Branch("runID", &_runID, "runID/i");
  fTree->Branch("subrunID", &_subrunID, "subrunID/i");

  fTree->Branch("flash_time","std::vector<double>", &_flash_time);
  fTree->Branch("flash_total_pe", "std::vector<double>", &_flash_total_pe);
  fTree->Branch("flash_tpc", "std::vector<int>", &_flash_tpc);
  fTree->Branch("flash_pmt_ratio", "std::vector<double>", &_flash_pmt_ratio);

  fTree->Branch("tr_midpoint_x", "std::vector<double>", &_tr_midpoint_x);
  fTree->Branch("tr_length", "std::vector<double>", &_tr_length);


  if(fMakePDSGeoTree){
    fPDSMapTree = tfs->make<TTree>("PDSMapTree", "PDS Map Tree");
    fPDSMapTree->Branch("OpDetID", "std::vector<int>", &_opDetID);
    fPDSMapTree->Branch("OpDetX", "std::vector<double>", &_opDetX);
    fPDSMapTree->Branch("OpDetY", "std::vector<double>", &_opDetY);
    fPDSMapTree->Branch("OpDetZ", "std::vector<double>", &_opDetZ);
    fPDSMapTree->Branch("OpDetType", "std::vector<int>", &_opDetType);

    if(fVerbosity>0)
      std::cout << " -- Dumping SBND PDS mapping -- \n";
    for(unsigned int opch=0; opch<fWireReadout.NOpChannels(); opch++){
      auto pdCenter = fWireReadout.OpDetGeoFromOpChannel(opch).GetCenter();
      _opDetID.push_back(opch);
      _opDetX.push_back(pdCenter.X());
      _opDetY.push_back(pdCenter.Y());
      _opDetZ.push_back(pdCenter.Z());
      if(fPDSMap.pdType(opch)=="pmt_coated") _opDetType.push_back(kPMTCoated);
      else if(fPDSMap.pdType(opch)=="pmt_uncoated") _opDetType.push_back(kPMTUncoated);
      else if(fPDSMap.pdType(opch)=="xarapuca_vuv") _opDetType.push_back(kXARAPUCAVUV);
      else if(fPDSMap.pdType(opch)=="xarapuca_vis") _opDetType.push_back(kXARAPUCAVIS);
      else _opDetType.push_back(kPDUnknown);
      if(fVerbosity>0)
        std::cout << opch << " " << pdCenter.X() << " " << pdCenter.Y() << " " << pdCenter.Z() << " " << fPDSMap.pdType(opch) << std::endl;
    }
    fPDSMapTree->Fill();
  }
  for(size_t oc=0; oc<fPDSMap.size(); oc++){
    fPDSBoxIDs.insert( fPDSMap.pdBox(oc) );
  }
}


// -------- Main function --------
void opdet::EtaPMTAna::analyze(art::Event const& e)
{
  _flash_pmt_ratio.clear();
  _flash_tpc.clear();
  _flash_total_pe.clear();
  _flash_time.clear();

  _tr_midpoint_x.clear();
  _tr_length.clear();


  // --- Event General Info
  _eventID = e.id().event();
  _runID = e.id().run();
  _subrunID = e.id().subRun();
  if(fVerbosity>0) std::cout << " -- Running SBNDPDSAnalyzer -- \n Run=" << _runID << " Subrun=" << _subrunID << " Event=" << _eventID << std::endl;

    //Read opflash
    art::Handle< std::vector<recob::OpFlash> > opflashListHandle;
    //Read crttracks
    art::Handle<std::vector<sbnd::crt::CRTTrack>> CRTTrackHandle;
    e.getByLabel(fCRTTrackModuleLabel, CRTTrackHandle);
    std::vector<art::Ptr<sbnd::crt::CRTTrack>> CRTTrackVec;
    art::fill_ptr_vector(CRTTrackVec, CRTTrackHandle);

    // Read OpFlash Time
    for (size_t s = 0; s < fOpFlashesModuleLabel.size(); s++) {
      e.getByLabel(fOpFlashesModuleLabel[s], opflashListHandle);
      if(!opflashListHandle.isValid()){
        std::cout << "OpFlash with label " << fOpFlashesModuleLabel[s] << " not found..." << std::endl;
        throw std::exception();
      }
      if(fVerbosity>0)
        std::cout << "Saving OpFlashes from " << fOpFlashesModuleLabel[s] << std::endl;
      for (unsigned int i = 0; i < opflashListHandle->size(); ++i) {
        // Get OpFlash
        art::Ptr<recob::OpFlash> FlashPtr(opflashListHandle, i);
        recob::OpFlash Flash = *FlashPtr;
        std::vector<size_t> associatedTracks;
        FillAssociatedTracks(Flash, CRTTrackVec);
      }
    }
    fTree->Fill();
}



void opdet::EtaPMTAna::FillAssociatedTracks(const recob::OpFlash& flash, const std::vector<art::Ptr<sbnd::crt::CRTTrack>>& crtTracks) {

  size_t _tpc = flash.XCenter() < 0 ? 0 : 1; // Determine TPC based on flash X position

  for (size_t i = 0; i < crtTracks.size(); ++i) {
        const auto& track = crtTracks[i];
        double time_diff = std::abs(track->Ts0() / 1000.0 - flash.AbsTime());
        if(time_diff > fMaxTimeDiff) continue; // Skip if time difference is greater than min time coincidence requirement
        TVector3 EntryPoint;
        TVector3 ExitPoint;
        CRTTrackCrossesAV(_tpc, *track, EntryPoint, ExitPoint);
        double _track_delta_x = abs(ExitPoint.X() - EntryPoint.X());
        if(_track_delta_x>fIsochronousWidth) continue; // Skip if track is not approximately parallel to x
        double _track_x = (ExitPoint.X() + EntryPoint.X()) / 2.0;
        double _track_length_in_av = (ExitPoint - EntryPoint).Mag();
        if(_track_length_in_av < fMinTrackLength) continue; // Skip if track length is less than the minimum
        double _pmt_ratio = GetPMTRatio(flash.PEs());
        _flash_pmt_ratio.push_back(_pmt_ratio);
        _flash_tpc.push_back(_tpc);
        _flash_total_pe.push_back(flash.TotalPE());
        _flash_time.push_back(flash.AbsTime());
        _tr_midpoint_x.push_back(_track_x);
        _tr_length.push_back(_track_length_in_av);
    }
}

// -------- Function to check if a CRTTrack crosses the AV of TPC(0/1) and return then entry/exit point --------
bool opdet::EtaPMTAna::CRTTrackCrossesAV(const int tpc, const sbnd::crt::CRTTrack &crttrack, TVector3& EntryPoint, TVector3& ExitPoint){

  TVector3 AV_min, AV_max;
  // Define the AV limits for TPC0 and TPC1
  if (tpc == 0) {
      AV_min = TVector3(fAVmin_tpc0[0], fAVmin_tpc0[1], fAVmin_tpc0[2]);
      AV_max = TVector3(fAVmax_tpc0[0], fAVmax_tpc0[1], fAVmax_tpc0[2]);
  }
  else if (tpc == 1) {
      AV_min = TVector3(fAVmin_tpc1[0], fAVmin_tpc1[1], fAVmin_tpc1[2]);
      AV_max = TVector3(fAVmax_tpc1[0], fAVmax_tpc1[1], fAVmax_tpc1[2]);
  }
  else std::runtime_error("Invalid TPC number. Only 0 and 1 are valid.");

  const geo::Point_t start = crttrack.Start();
  const geo::Point_t end = crttrack.End();

  TVector3 p1 = {start.X(), start.Y(), start.Z()};
  TVector3 p2 = {end.X(), end.Y(), end.Z()};

  TVector3 d = p2 - p1;
  for (int i = 0; i < 3; ++i) {
      if (d[i] == 0.0)
          d[i] = 1e-10;
  }

  std::array<double, 3> t_min, t_max, t1, t2;
  for (int i = 0; i < 3; ++i) {
      t_min[i] = (AV_min[i] - p1[i]) / d[i];
      t_max[i] = (AV_max[i] - p1[i]) / d[i];
      t1[i] = std::min(t_min[i], t_max[i]);
      t2[i] = std::max(t_min[i], t_max[i]);
  }

  double t_entry = *std::max_element(t1.begin(), t1.end());
  double t_exit  = *std::min_element(t2.begin(), t2.end());

  if (t_entry <= t_exit && t_exit >= 0.0 && t_entry <= 1.0) {
    TVector3 entry = p1 + d * t_entry;
    TVector3 exit = p1 + d * t_exit;
    EntryPoint = entry;
    ExitPoint = exit;
    return true;
  }
  else
  {
    EntryPoint={-9999., -9999., -9999.};
    ExitPoint={-9999., -9999., -9999.};
    return false;
  }
}



  double opdet::EtaPMTAna::GetPMTRatio(std::vector<double> PE_v){

    std::map<int, double> BoxMap_PECoated;
    std::map<int, double> BoxMap_PEUncoated;
    std::map<int, int> BoxMap_NCoatedCh;
    std::map<int, int> BoxMap_NUncoatedCh;

    for(size_t oc=0; oc<PE_v.size(); oc++){
      // skip 0 PE channels
      if(PE_v[oc]==0) continue;

      std::string pd_type=fPDSMap.pdType(oc);
      // exclude xarapucas by now
      if(pd_type=="xarapuca_vuv" || pd_type=="xarapuca_vis") continue;

      // get PDS box
      int box_id = fPDSMap.pdBox(oc);

      // we store the pe in each box per PMT flavour
      // and the number of "triggered" PMTs
      if(pd_type=="pmt_coated") {
        BoxMap_PECoated[box_id]+=PE_v[oc];
        BoxMap_NCoatedCh[box_id]+=1;
      }
      else if(pd_type=="pmt_uncoated") {
        BoxMap_PEUncoated[box_id]+=PE_v[oc];
        BoxMap_NUncoatedCh[box_id]+=1;
      }
    }

    double pmtratio=-1.;
    // compute PMTRatio metric

    double TotalPE=0;
    double RatioPerBoxWeight=0;
    for(size_t boxID=0; boxID<fPDSBoxIDs.size(); boxID++){
      double RatioPerBox;
      double PECoated=0, PEUncoated=0;
      //we need the uncoated PMT in each window and at least one coated
      if( BoxMap_NUncoatedCh[boxID]==1 && BoxMap_NCoatedCh[boxID]>=1){
        PECoated+= BoxMap_PECoated[boxID];
        PEUncoated+=BoxMap_PEUncoated[boxID];
        TotalPE += PECoated + PEUncoated;
        RatioPerBox = (PEUncoated/PECoated)*BoxMap_NCoatedCh[boxID];
        RatioPerBoxWeight += RatioPerBox * (PECoated + PEUncoated);
      }
    }
    if(TotalPE!=0) pmtratio = RatioPerBoxWeight/TotalPE;
  
    return pmtratio;

  }
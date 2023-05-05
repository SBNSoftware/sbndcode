////////////////////////////////////////////////////////////////////////
// Class:       EventDisplayAndDump
// Plugin Type: analyzer
// File:        EventDisplayAndDump_module.cc
// Author:      Henry Lay (h.lay@lancaster.ac.uk)
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"

#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TCanvas.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TSystem.h"

typedef std::map<int, std::vector<art::Ptr<recob::Hit>>> HitMap;

constexpr double wireAngleU = 1.04719758034;
constexpr double wireAngleV = -1.04719758034;
constexpr double wireAngleW = 0;

const double cosU = TMath::Cos(wireAngleU);
const double sinU = TMath::Sin(wireAngleU);
const double cosV = TMath::Cos(wireAngleV);
const double sinV = TMath::Sin(wireAngleV);
const double cosW = TMath::Cos(wireAngleW);
const double sinW = TMath::Sin(wireAngleW);

namespace sbnd {
  class EventDisplayAndDump;
}

class sbnd::EventDisplayAndDump : public art::EDAnalyzer {
public:
  explicit EventDisplayAndDump(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  EventDisplayAndDump(EventDisplayAndDump const&) = delete;
  EventDisplayAndDump(EventDisplayAndDump&&) = delete;
  EventDisplayAndDump& operator=(EventDisplayAndDump const&) = delete;
  EventDisplayAndDump& operator=(EventDisplayAndDump&&) = delete;

  // Required functions.
  void analyze(const art::Event &e) override;

  bool SignalEvent(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles);

  void DumpTrue(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles, std::ofstream &log);

  void DumpReco(const art::Event &e, const std::vector<art::Ptr<recob::Slice>> &sliceVec, std::ofstream &log);

  std::vector<HitMap> GetTrueMaps(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &hitVec);

  std::vector<HitMap> GetRecoMaps(const art::Event &e, const std::vector<art::Ptr<recob::Slice>> &sliceVec);

  void DrawTrueFull(const art::Event &e, const std::vector<HitMap> &maps);

  void DrawTrueSep(const art::Event &e, const std::vector<HitMap> &maps);

  void DrawTrueView(const HitMap &hitMap0, const HitMap &hitMap1, const TString &name, const bool draw_tpc_edges = true);

  void DrawRecoFull(const art::Event &e, const std::vector<HitMap> &maps);

  void DrawRecoSep(const art::Event &e, const std::vector<HitMap> &maps);

  void DrawRecoView(const HitMap &hitMap0, const HitMap &hitMap1, const TString &name, const bool draw_tpc_edges = true);

  bool VolumeCheck(const TVector3 &pos, const double &walls = 0., const double &cath = 0., const double &front = 0., const double &back = 0.);

  double YZtoU(const double y, const double z);

  double YZtoV(const double y, const double z);

  double YZtoW(const double y, const double z);

private:

  art::ServiceHandle<cheat::ParticleInventoryService> particleInv;
  art::ServiceHandle<cheat::BackTrackerService>       backTracker;

  std::string fMCParticleModuleLabel, fHitModuleLabel, fSliceModuleLabel,
    fPFPModuleLabel, fTrackModuleLabel, fShowerModuleLabel;

  bool fAVOnly, fPrintTrue, fLogTrue, fPrintReco, fLogReco, fDrawTrueFull, fDrawTrueSep, fDrawRecoFull, fDrawRecoSep,
    fSavePDF, fSavePNG;

  std::string fSaveDir, fEventSaveDir;

  int  fRun;
  int  fSubRun;
  int  fEvent;
  bool fSignal;

  std::set<int> fTrackIds;

  std::vector<int> fColours;
  std::map<int, int> fTrueColourMap, fRecoColourMap;
  int fTrueColourCounter, fRecoColourCounter;

  std::map<int, TString> fPDGString;
};

sbnd::EventDisplayAndDump::EventDisplayAndDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
  , fMCParticleModuleLabel(p.get<std::string>("MCParticleModuleLabel", "largeant"))
  , fHitModuleLabel(p.get<std::string>("HitModuleLabel", "gaushit"))
  , fSliceModuleLabel(p.get<std::string>("SliceModuleLabel", "pandoraSCE"))
  , fPFPModuleLabel(p.get<std::string>("PFPModuleLabel", "pandoraSCE"))
  , fTrackModuleLabel(p.get<std::string>("TrackModuleLabel", "pandoraSCETrack"))
  , fShowerModuleLabel(p.get<std::string>("ShowerModuleLabel", "pandoraSCEShower"))
  , fAVOnly(p.get<bool>("AVOnly", true))
  , fPrintTrue(p.get<bool>("PrintTrue", true))
  , fLogTrue(p.get<bool>("LogTrue", true))
  , fPrintReco(p.get<bool>("PrintReco", true))
  , fLogReco(p.get<bool>("LogReco", true))
  , fDrawTrueFull(p.get<bool>("DrawTrueFull", true))
  , fDrawTrueSep(p.get<bool>("DrawTrueSep", true))
  , fDrawRecoFull(p.get<bool>("DrawRecoFull", true))
  , fDrawRecoSep(p.get<bool>("DrawRecoSep", true))
  , fSavePDF(p.get<bool>("SavePDF", true))
  , fSavePNG(p.get<bool>("SavePNG", true))
  , fSaveDir(p.get<std::string>("SaveDir", "/sbnd/data/users/hlay/ncpizero/2022A/plots/spring2023/evds"))
  {
    fEventSaveDir = fSaveDir;

    fColours = {kRed+2, kBlue+2, kGreen+2, kMagenta+2, kCyan+2, kYellow+2, kAzure+1, kSpring+9, kPink+9, kTeal+9, kOrange+2, kViolet-8};

    fPDGString[13]    = "#mu^{-}";
    fPDGString[-13]   = "#mu^{+}";
    fPDGString[11]    = "e^{-}";
    fPDGString[-11]   = "e^{+}";
    fPDGString[211]   = "#pi^{+}";
    fPDGString[-211]  = "#pi^{-}";
    fPDGString[2212]  = "p";
    fPDGString[22]    = "#gamma";
    fPDGString[22111] = "#gamma (#pi^{0})";

    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineWidth(4);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetStatColor(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFont(62);
    gStyle->SetLegendTextSize(0.09);

    // gStyle->SetPaperSize(20,26);
    // gStyle->SetCanvasDefH(1000);
    // gStyle->SetCanvasDefW(1400);
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetPadRightMargin(0.02);
    gStyle->SetPadBottomMargin(0.1);
    gStyle->SetPadLeftMargin(0.02);

    gStyle->SetTextFont(62);
    gStyle->SetTextSize(0.07);
    gStyle->SetLabelFont(62,"x");
    gStyle->SetLabelFont(62,"y");
    gStyle->SetLabelFont(62,"z");
    gStyle->SetLabelSize(0,"x");
    gStyle->SetTitleSize(0,"x");
    gStyle->SetLabelSize(0,"y");
    gStyle->SetTitleSize(0,"y");
    gStyle->SetLabelSize(0,"z");
    gStyle->SetTitleSize(0,"z");
    gStyle->SetLabelFont(62,"t");
    gStyle->SetTitleFont(62,"x");
    gStyle->SetTitleFont(62,"y");
    gStyle->SetTitleFont(62,"z");
    gStyle->SetTitleFont(62,"t");
    gStyle->SetTitleFillColor(kWhite);
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleAlign(23);
    gStyle->SetTitleOffset(1,"y");
    gStyle->SetTitleOffset(1,"y");
    gStyle->SetTitleOffset(1,"z");
    gStyle->SetTitleFontSize(0.07);
    gStyle->SetTitleFont(62,"pad");
    gStyle->SetTitleBorderSize(0);

    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerSize(.25);
    gStyle->SetHistLineWidth(0);
    gStyle->SetLineStyleString(2,"[12 12]");

    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    gStyle->SetNdivisions(5000, "X");
    gStyle->SetNdivisions(1000, "Y");

    gStyle->SetPalette(kRainBow);
  }

void sbnd::EventDisplayAndDump::analyze(const art::Event &e)
{
  fRun    = e.id().run();
  fSubRun = e.id().subRun();
  fEvent  =  e.id().event();

  fEventSaveDir = Form("%s/r%i_s%i_e%i", fSaveDir.c_str(), fRun, fSubRun, fEvent);
  if(fLogTrue || fDrawTrueFull || fDrawTrueSep || fDrawRecoFull || fDrawRecoSep)
    gSystem->Exec(Form("mkdir -p %s", fEventSaveDir.c_str()));

  std::ofstream log(Form("%s/evd_r%i_s%i_e%i_log.txt", fEventSaveDir.c_str(), fRun, fSubRun, fEvent), std::ios_base::out | std::ios_base::app);

  if(fPrintTrue) std::cout << "This is event " << fRun << "-" << fSubRun << "-" << fEvent << std::endl;
  if(fLogTrue) log << "This is event " << fRun << "-" << fSubRun << "-" << fEvent << std::endl;

  // Get MCTruths
  std::vector<art::Handle<std::vector<simb::MCTruth>>> MCTruthHandles = e.getMany<std::vector<simb::MCTruth>>();
  fSignal = SignalEvent(e, MCTruthHandles);

  if(fSignal)
    {
      if(fPrintTrue) std::cout << "Is signal event" << std::endl;
      if(fLogTrue) log << "Is signal event" << std::endl;
    }
  else
    {
      if(fPrintTrue) std::cout << "Is not signal event" << std::endl;
      if(fLogTrue) log << "Is not signal event" << std::endl;
    }

  if(fPrintTrue || fLogTrue)
    DumpTrue(e, MCTruthHandles, log);
  
  // Get Hits
  art::Handle<std::vector<recob::Hit>> hitHandle;
  e.getByLabel(fHitModuleLabel, hitHandle);

  std::vector<art::Ptr<recob::Hit>> hitVec;
  art::fill_ptr_vector(hitVec, hitHandle);

  if(fDrawTrueFull || fDrawTrueSep)
    {
      std::vector<HitMap> maps = GetTrueMaps(e, hitVec);

      if(fDrawTrueFull)
        DrawTrueFull(e, maps);

      if(fDrawTrueSep)
        DrawTrueSep(e, maps);
    }

  // Get Slices
  art::Handle<std::vector<recob::Slice>> sliceHandle;
  e.getByLabel(fSliceModuleLabel, sliceHandle);

  std::vector<art::Ptr<recob::Slice>> sliceVec;
  art::fill_ptr_vector(sliceVec, sliceHandle);

  if(fPrintReco || fLogReco)
    DumpReco(e, sliceVec, log);

  if(fDrawTrueFull || fDrawTrueSep)
    {
      std::vector<HitMap> maps = GetRecoMaps(e, sliceVec);

      if(fDrawRecoFull)
        DrawRecoFull(e, maps);

      if(fDrawRecoSep)
        DrawRecoSep(e, maps);
    }
}

bool sbnd::EventDisplayAndDump::SignalEvent(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles)
{
  std::vector<std::pair<bool, double>> events;

  for(auto const& MCTruthHandle : MCTruthHandles)
    {
      std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
      art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

      art::FindManyP<simb::MCParticle> MCTruthToMCParticles(MCTruthHandle, e, fMCParticleModuleLabel);

      for(auto const& mct : MCTruthVec)
        {
          if(mct->Origin() != 1)
            continue;

          const simb::MCNeutrino mcn = mct->GetNeutrino();
          const simb::MCParticle nu  = mcn.Nu();

          const bool nc = mcn.CCNC() == 1;
          const bool fv = VolumeCheck(nu.Position().Vect(), 20., 5., 10., 50.);

          unsigned pizeros = 0;

          for(int i = 0; i < mct->NParticles(); ++i)
            {
              const auto mcp = mct->GetParticle(i);

              if(mcp.PdgCode() == 111 && mcp.StatusCode() != 1)
                ++pizeros;
            }

          const bool pizero = pizeros > 0;

          const std::vector<art::Ptr<simb::MCParticle>> MCParticleVec = MCTruthToMCParticles.at(mct.key());
          double total_en = 0.;

          for(auto const& mcp : MCParticleVec)
            {
              std::vector<const sim::IDE*> ides = backTracker->TrackIdToSimIDEs_Ps(mcp->TrackId());

              for(auto const& ide : ides)
                total_en += ide->energy / 1000.;
            }

          events.push_back({nc && fv && pizero, total_en});
        }
    }

  if(events.size() == 0)
    return false;

  std::sort(events.begin(), events.end(),
            [](const auto &a, const auto &b)
            { return a.second > b.second; });

  return events.at(0).first;
}

void sbnd::EventDisplayAndDump::DumpTrue(const art::Event &e, const std::vector<art::Handle<std::vector<simb::MCTruth>>> &MCTruthHandles, std::ofstream &log)
{
  for(auto const& MCTruthHandle : MCTruthHandles)
    {
      std::vector<art::Ptr<simb::MCTruth>> MCTruthVec;
      art::fill_ptr_vector(MCTruthVec, MCTruthHandle);

      art::FindManyP<simb::MCParticle> MCTruthToMCParticles(MCTruthHandle, e, fMCParticleModuleLabel);

      for(auto const& mct : MCTruthVec)
        {
          if(mct->Origin() != 1)
            continue;

          const simb::MCNeutrino mcn = mct->GetNeutrino();
          const simb::MCParticle nu  = mcn.Nu();

          if(fAVOnly && !VolumeCheck(nu.Position().Vect()))
            continue;

          const bool nc = mcn.CCNC() == 1;
	  
          if(fPrintTrue) std::cout << '\n'
				   << "===== Neutrino Event =====\n"
				   << "PDG: " << nu.PdgCode() << '\n'
				   << "CCNC: " << (nc ? "NC" : "CC") << '\n'
				   << "E: " << nu.E() << '\n'
				   << "Vtx: (" << nu.Vx() << ", " << nu.Vy() << ", " << nu.Vz() << ") cm\n";

          if(fLogTrue) log << '\n'
			   << "===== Neutrino Event =====\n"
			   << "PDG: " << nu.PdgCode() << '\n'
			   << "CCNC: " << (nc ? "NC" : "CC") << '\n'
			   << "E: " << nu.E() << '\n'
			   << "Vtx: (" << nu.Vx() << ", " << nu.Vy() << ", " << nu.Vz() << ") cm\n";

          const std::vector<art::Ptr<simb::MCParticle>> MCParticleVec = MCTruthToMCParticles.at(mct.key());

          int i = 0;

          for(auto const& mcp : MCParticleVec)
            {
              if(mcp->StatusCode() != 1)
                continue;

              fTrackIds.insert(mcp->TrackId());

              if(mcp->Mother() != nu.TrackId() + 10000000)
                continue;

              if(fPrintTrue) std::cout << "\tParticle " << i << '\n'
				       << "\t\tPDG: " << mcp->PdgCode() << '\n'
				       << "\t\tStatusCode: " << mcp->StatusCode() << '\n'
				       << "\t\tTrackID: " << mcp->TrackId() << '\n'
				       << "\t\tMother: " << mcp->Mother() << '\n'
				       << "\t\tNumDaughters: " << mcp->NumberDaughters() << '\n'
				       << "\t\tLength: " << (mcp->EndPosition().Vect() - mcp->Position().Vect()).Mag() << '\n'
				       << "\t\tE: " << mcp->E() << '\n'
				       << "\t\tp: (" << mcp->Px() << ", " << mcp->Py() << ", " << mcp->Pz() << ") cm\n";

              if(fLogTrue) log << "\tParticle " << i << '\n'
			       << "\t\tPDG: " << mcp->PdgCode() << '\n'
			       << "\t\tStatusCode: " << mcp->StatusCode() << '\n'
			       << "\t\tTrackID: " << mcp->TrackId() << '\n'
			       << "\t\tMother: " << mcp->Mother() << '\n'
			       << "\t\tNumDaughters: " << mcp->NumberDaughters() << '\n'
			       << "\t\tLength: " << (mcp->EndPosition().Vect() - mcp->Position().Vect()).Mag() << '\n'
			       << "\t\tE: " << mcp->E() << '\n'
			       << "\t\tp: (" << mcp->Px() << ", " << mcp->Py() << ", " << mcp->Pz() << ") cm\n";

              if(mcp->PdgCode() == 111 && mcp->NumberDaughters() == 2)
                {
                  const simb::MCParticle* gamma1 = particleInv->TrackIdToParticle_P(mcp->Daughter(0));
                  const simb::MCParticle* gamma2 = particleInv->TrackIdToParticle_P(mcp->Daughter(1));

                  const TVector3 gamma1_mom = gamma1->Momentum().Vect();
                  const TVector3 gamma2_mom = gamma2->Momentum().Vect();
                  const double open_angle   = TMath::RadToDeg() * gamma1_mom.Angle(gamma2_mom);

                  if(fPrintTrue) std::cout << "\t\t\tPDG: " << gamma1->PdgCode() << '\n'
					   << "\t\t\tStatusCode: " << gamma1->StatusCode() << '\n'
					   << "\t\t\tTrackID: " << gamma1->TrackId() << '\n'
					   << "\t\t\tE: " << gamma1->E() << '\n'
					   << "\t\t\tp: (" << gamma1->Px() << ", " << gamma1->Py() << ", " << gamma1->Pz() << ") cm\n"
					   << "\t\t\tLength: " << (gamma1->EndPosition().Vect() - gamma1->Position().Vect()).Mag() << '\n'
					   << "\t\t\tPDG: " << gamma2->PdgCode() << '\n'
					   << "\t\t\tStatusCode: " << gamma2->StatusCode() << '\n'
					   << "\t\t\tTrackID: " << gamma2->TrackId() << '\n'
					   << "\t\t\tE: " << gamma2->E() << '\n'
					   << "\t\t\tp: (" << gamma2->Px() << ", " << gamma2->Py() << ", " << gamma2->Pz() << ") cm\n"
					   << "\t\t\tLength: " << (gamma2->EndPosition().Vect() - gamma2->Position().Vect()).Mag() << '\n'
					   << "\t\tOpen Angle: " << open_angle << '\n';

                  if(fLogTrue) log << "\t\t\tPDG: " << gamma1->PdgCode() << '\n'
				   << "\t\t\tStatusCode: " << gamma1->StatusCode() << '\n'
				   << "\t\t\tTrackID: " << gamma1->TrackId() << '\n'
				   << "\t\t\tE: " << gamma1->E() << '\n'
				   << "\t\t\tp: (" << gamma1->Px() << ", " << gamma1->Py() << ", " << gamma1->Pz() << ") cm\n"
				   << "\t\t\tLength: " << (gamma1->EndPosition().Vect() - gamma1->Position().Vect()).Mag() << '\n'
				   << "\t\t\tPDG: " << gamma2->PdgCode() << '\n'
				   << "\t\t\tStatusCode: " << gamma2->StatusCode() << '\n'
				   << "\t\t\tTrackID: " << gamma2->TrackId() << '\n'
				   << "\t\t\tE: " << gamma2->E() << '\n'
				   << "\t\t\tp: (" << gamma2->Px() << ", " << gamma2->Py() << ", " << gamma2->Pz() << ") cm\n"
				   << "\t\t\tLength: " << (gamma2->EndPosition().Vect() - gamma2->Position().Vect()).Mag() << '\n'
				   << "\t\tOpen Angle: " << open_angle << '\n';
                }
              ++i;
            }
          if(fPrintTrue) std::cout << std::endl;
	  if(fLogTrue) log << std::endl;
        }
    }

  if(fPrintTrue) std::cout << std::endl;
}

void sbnd::EventDisplayAndDump::DumpReco(const art::Event &e, const std::vector<art::Ptr<recob::Slice>> &sliceVec, std::ofstream &log)
{
  art::Handle<std::vector<recob::Slice>> sliceHandle;
  e.getByLabel(fSliceModuleLabel, sliceHandle);

  art::FindManyP<recob::PFParticle> slicesToPFPs(sliceHandle, e, fPFPModuleLabel);

  for(auto const &slice : sliceVec)
    {
      const std::vector<art::Ptr<recob::PFParticle>> pfps = slicesToPFPs.at(slice.key());

      int pdg = 0;
      for(auto const& pfp : pfps)
        {
          if(pfp->IsPrimary())
            {
              pdg = pfp->PdgCode();
              break;
            }
        }

      if(pdg == 0)
        std::cout << "No primary found" << std::endl;
      else if(pdg == 11 || pdg == 13)
        continue;

      if(fPrintReco) std::cout << '\n'
			       << "===== Reconstructed Neutrino Slice =====\n"
			       << "PDG: " << pdg << '\n'
			       << "nPFPs: " << pfps.size() << '\n';

      if(fLogReco) log << '\n'
		       << "===== Reconstructed Neutrino Slice =====\n"
		       << "PDG: " << pdg << '\n'
		       << "nPFPs: " << pfps.size() << '\n';

      for(auto const& pfp : pfps)
        {
          const int self = pfp->Self();
          std::vector<art::Ptr<recob::Hit>> hits;

	  if(fPrintReco) std::cout << "\tPFP: " << self << '\n'
				   << "\t\tPDG: " << pfp->PdgCode() << '\n';

	  if(fLogReco) log << "\tPFP: " << self << '\n'
			   << "\t\tPDG: " << pfp->PdgCode() << '\n';
	}
    }
}

std::vector<HitMap> sbnd::EventDisplayAndDump::GetTrueMaps(const art::Event &e, const std::vector<art::Ptr<recob::Hit>> &hitVec)
{
  const detinfo::DetectorClocksData clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);

  std::vector<HitMap> maps;
  maps.resize(6);

  for(auto const &hit : hitVec)
    {
      const int view    = hit->View();
      const int tpc     = hit->WireID().asTPCID().deepestIndex();
      const int trackid = TruthMatchUtils::TrueParticleID(clockData,hit,true);

      if(tpc == 0)
        {
          if(view == geo::kU)
            maps[0][trackid].push_back(hit);
          else if(view == geo::kV)
            maps[1][trackid].push_back(hit);
          else if(view == geo::kW)
            maps[2][trackid].push_back(hit);
          else
            std::cout << "Hit with view: " << view << std::endl;
        }
      else if(tpc == 1)
        {
          if(view == geo::kU)
            maps[3][trackid].push_back(hit);
          else if(view == geo::kV)
            maps[4][trackid].push_back(hit);
          else if(view == geo::kW)
            maps[5][trackid].push_back(hit);
          else
            std::cout << "Hit with view: " << view << std::endl;
        }
      else
        std::cout << "TPC " << tpc << "?????" << std::endl;
    }

  return maps;
}

std::vector<HitMap> sbnd::EventDisplayAndDump::GetRecoMaps(const art::Event &e, const std::vector<art::Ptr<recob::Slice>> &sliceVec)
{
  art::Handle<std::vector<recob::Slice>> sliceHandle;
  e.getByLabel(fSliceModuleLabel, sliceHandle);

  art::Handle<std::vector<recob::PFParticle>> PFPHandle;
  e.getByLabel(fPFPModuleLabel, PFPHandle);

  art::Handle<std::vector<recob::Track>> trackHandle;
  e.getByLabel(fTrackModuleLabel, trackHandle);

  art::Handle<std::vector<recob::Shower>> showerHandle;
  e.getByLabel(fShowerModuleLabel, showerHandle);

  art::FindManyP<recob::PFParticle> slicesToPFPs(sliceHandle, e, fPFPModuleLabel);
  art::FindManyP<recob::Track> pfpsToTracks(PFPHandle, e, fTrackModuleLabel);
  art::FindManyP<recob::Shower> pfpsToShowers(PFPHandle, e, fShowerModuleLabel);
  art::FindManyP<recob::Hit> tracksToHits(trackHandle, e, fTrackModuleLabel);
  art::FindManyP<recob::Hit> showersToHits(showerHandle, e, fShowerModuleLabel);

  std::vector<HitMap> maps;
  maps.resize(6);

  for(auto const &slice : sliceVec)
    {
      const std::vector<art::Ptr<recob::PFParticle>> pfps = slicesToPFPs.at(slice.key());

      int pdg = 0;
      for(auto const& pfp : pfps)
        {
          if(pfp->IsPrimary())
            {
              pdg = pfp->PdgCode();
              break;
            }
        }

      if(pdg == 0)
        std::cout << "No primary found" << std::endl;
      else if(pdg == 11 || pdg == 13)
        continue;

      for(auto const& pfp : pfps)
        {
          const int self = pfp->Self();
          std::vector<art::Ptr<recob::Hit>> hits;

          if(pfp->PdgCode() == 13)
            {
              const std::vector<art::Ptr<recob::Track>> tracks = pfpsToTracks.at(pfp.key());
              hits = tracksToHits.at(tracks[0].key());
            }
          else if(pfp->PdgCode() == 11)
            {
              const std::vector<art::Ptr<recob::Shower>> showers = pfpsToShowers.at(pfp.key());
              hits = showersToHits.at(showers[0].key());
            }

          for(auto const &hit : hits)
            {
              const int view = hit->View();
              const int tpc  = hit->WireID().asTPCID().deepestIndex();

              if(tpc == 0)
                {
                  if(view == geo::kU)
                    maps[0][self].push_back(hit);
                  else if(view == geo::kV)
                    maps[1][self].push_back(hit);
                  else if(view == geo::kW)
                    maps[2][self].push_back(hit);
                  else
                    std::cout << "Hit with view: " << view << std::endl;
                }
              else if(tpc == 1)
                {
                  if(view == geo::kU)
                    maps[3][self].push_back(hit);
                  else if(view == geo::kV)
                    maps[4][self].push_back(hit);
                  else if(view == geo::kW)
                    maps[5][self].push_back(hit);
                  else
                    std::cout << "Hit with view: " << view << std::endl;
                }
              else
                std::cout << "TPC " << tpc << "?????" << std::endl;
            }
        }
    }

  return maps;
}

void sbnd::EventDisplayAndDump::DrawTrueFull(const art::Event &e, const std::vector<HitMap> &maps)
{
  fTrueColourMap.clear();
  fTrueColourCounter = 0;

  TCanvas *c = new TCanvas("c", "", 2100, 4200);
  c->cd();
  c->Divide(1,3);

  c->cd(1);
  DrawTrueView(maps[0], maps[3], "U view", false);

  c->cd(2);
  DrawTrueView(maps[1], maps[4], "V view", false);

  c->cd(3);
  DrawTrueView(maps[2], maps[5], "W view", true);

  if(fSavePDF) c->SaveAs(Form("%s/evd_r%i_s%i_e%i_true.pdf", fEventSaveDir.c_str(), fRun, fSubRun, fEvent));
  if(fSavePNG) c->SaveAs(Form("%s/evd_r%i_s%i_e%i_true.png", fEventSaveDir.c_str(), fRun, fSubRun, fEvent));

  delete c;
}

void sbnd::EventDisplayAndDump::DrawTrueSep(const art::Event &e, const std::vector<HitMap> &maps)
{
  fTrueColourMap.clear();
  fTrueColourCounter = 0;

  TCanvas *cU = new TCanvas("cU", "", 2100, 1400);
  cU->cd();
  DrawTrueView(maps[0], maps[3], "U view", false);

  if(fSavePDF) cU->SaveAs(Form("%s/evd_r%i_s%i_e%i_u_view_true.pdf", fEventSaveDir.c_str(), fRun, fSubRun, fEvent));
  if(fSavePNG) cU->SaveAs(Form("%s/evd_r%i_s%i_e%i_u_view_true.png", fEventSaveDir.c_str(), fRun, fSubRun, fEvent));

  delete cU;

  TCanvas *cV = new TCanvas("cV", "", 2100, 1400);
  cV->cd();
  DrawTrueView(maps[1], maps[4], "V view", false);

  if(fSavePDF) cV->SaveAs(Form("%s/evd_r%i_s%i_e%i_v_view_true.pdf", fEventSaveDir.c_str(), fRun, fSubRun, fEvent));
  if(fSavePNG) cV->SaveAs(Form("%s/evd_r%i_s%i_e%i_v_view_true.png", fEventSaveDir.c_str(), fRun, fSubRun, fEvent));

  delete cV;

  TCanvas *cW = new TCanvas("cW", "", 2100, 1400);
  cW->cd();
  DrawTrueView(maps[2], maps[5], "W view", true);

  if(fSavePDF) cW->SaveAs(Form("%s/evd_r%i_s%i_e%i_w_view_true.pdf", fEventSaveDir.c_str(), fRun, fSubRun, fEvent));
  if(fSavePNG) cW->SaveAs(Form("%s/evd_r%i_s%i_e%i_w_view_true.png", fEventSaveDir.c_str(), fRun, fSubRun, fEvent));

  delete cW;
}

void sbnd::EventDisplayAndDump::DrawRecoFull(const art::Event &e, const std::vector<HitMap> &maps)
{
  fRecoColourMap.clear();
  fRecoColourCounter = 0;

  TCanvas *c = new TCanvas("c", "", 2100, 4200);
  c->cd();
  c->Divide(1,3);

  c->cd(1);
  DrawRecoView(maps[0], maps[3], "U view", false);

  c->cd(2);
  DrawRecoView(maps[1], maps[4], "V view", false);

  c->cd(3);
  DrawRecoView(maps[2], maps[5], "W view", true);

  if(fSavePDF) c->SaveAs(Form("%s/evd_r%i_s%i_e%i_reco.pdf", fEventSaveDir.c_str(), fRun, fSubRun, fEvent));
  if(fSavePNG) c->SaveAs(Form("%s/evd_r%i_s%i_e%i_reco.png", fEventSaveDir.c_str(), fRun, fSubRun, fEvent));

  delete c;
}

void sbnd::EventDisplayAndDump::DrawRecoSep(const art::Event &e, const std::vector<HitMap> &maps)
{
  fRecoColourMap.clear();
  fRecoColourCounter = 0;

  TCanvas *cU = new TCanvas("cU", "", 2100, 1400);
  cU->cd();
  DrawRecoView(maps[0], maps[3], "U view", false);

  if(fSavePDF) cU->SaveAs(Form("%s/evd_r%i_s%i_e%i_u_view_reco.pdf", fEventSaveDir.c_str(), fRun, fSubRun, fEvent));
  if(fSavePNG) cU->SaveAs(Form("%s/evd_r%i_s%i_e%i_u_view_reco.png", fEventSaveDir.c_str(), fRun, fSubRun, fEvent));

  delete cU;

  TCanvas *cV = new TCanvas("cV", "", 2100, 1400);
  cV->cd();
  DrawRecoView(maps[1], maps[4], "V view", false);

  if(fSavePDF) cV->SaveAs(Form("%s/evd_r%i_s%i_e%i_v_view_reco.pdf", fEventSaveDir.c_str(), fRun, fSubRun, fEvent));
  if(fSavePNG) cV->SaveAs(Form("%s/evd_r%i_s%i_e%i_v_view_reco.png", fEventSaveDir.c_str(), fRun, fSubRun, fEvent));

  delete cV;

  TCanvas *cW = new TCanvas("cW", "", 2100, 1400);
  cW->cd();
  DrawRecoView(maps[2], maps[5], "W view", true);

  if(fSavePDF) cW->SaveAs(Form("%s/evd_r%i_s%i_e%i_w_view_reco.pdf", fEventSaveDir.c_str(), fRun, fSubRun, fEvent));
  if(fSavePNG) cW->SaveAs(Form("%s/evd_r%i_s%i_e%i_w_view_reco.png", fEventSaveDir.c_str(), fRun, fSubRun, fEvent));

  delete cW;
}

void sbnd::EventDisplayAndDump::DrawTrueView(const HitMap &hitMap0, const HitMap &hitMap1, const TString &name, const bool draw_tpc_edges)
{
  TPaveText *t = new TPaveText(.9, .95, .98, 1, "NB NDC");
  t->SetTextAlign(32);
  t->SetTextSize(0.085);
  t->SetTextColor(kBlack);
  t->SetFillStyle(4000);
  t->AddText(name);

  TPaveText *t0 = new TPaveText(.03, .92, .2, .95, "NB NDC");
  t0->SetTextAlign(12);
  t0->SetTextSize(0.06);
  t0->SetTextColor(kBlack);
  t0->SetFillStyle(4000);
  t0->AddText("TPC 0");

  TPaveText *t1 = new TPaveText(.03, .82, .2, .85, "NB NDC");
  t1->SetTextAlign(12);
  t1->SetTextSize(0.06);
  t1->SetTextColor(kBlack);
  t1->SetFillStyle(4000);
  t1->AddText("TPC 1");

  TPad *top    = new TPad("top", "", 0.005, 0.505, 0.995, 0.995);
  TPad *bottom = new TPad("bottom", "", 0.005, 0.005, 0.995, 0.495);
  top->SetBottomMargin(0.002);
  bottom->SetTopMargin(0.002);
  top->Draw();
  bottom->Draw();

  TLegend *leg = new TLegend(0.1, 0, .9, .05);

  top->cd();

  const double n_wires = name == "W view" ? 1664 : 1984;

  TH1F* base_hist = new TH1F("base_hist", "", 10, 0.5 - 0.5*(1984 - n_wires), 1984.5 - 0.5*(1984 - n_wires));
  base_hist->SetMaximum(3400.0);
  base_hist->Draw();
  t->Draw();
  t1->Draw();

  while(hitMap0.size() + hitMap1.size() > fColours.size())
    fColours.insert(fColours.end(), fColours.begin(), fColours.end());

  std::set<int> tpc1_used;

  for(auto const& [trackid, hits] : hitMap1)
    {
      if(hits.size() < 5)
	continue;

      const simb::MCParticle* part = particleInv->TrackIdToParticle_P(trackid);
      int pdg = part == NULL ? 2 : part->PdgCode();
      if(part != NULL && fTrackIds.count(part->TrackId()) == 0)
        pdg = 2;

      if(pdg == 22)
	{
	  const simb::MCParticle *mother = particleInv->TrackIdToParticle_P(part->Mother());
	  if(mother != NULL && mother->PdgCode() == 111)
	    pdg = 22111;
	}

      TGraph *g = new TGraph();

      for(auto const& hit : hits)
        g->SetPoint(g->GetN(), hit->WireID().deepestIndex(), 3400-hit->PeakTime());

      if(pdg == 2)
	g->SetMarkerColor(kGray);
      else
	{
	  const int colour = fTrueColourMap.count(trackid) == 0 ? fColours[fTrueColourCounter] : fTrueColourMap[trackid];
	  g->SetMarkerColor(colour);
	  TLegendEntry *le = leg->AddEntry(g, fPDGString[pdg], "");
	  le->SetTextColor(colour);
	  if(fTrueColourMap.count(trackid) == 0)
	    ++fTrueColourCounter;
	  fTrueColourMap[trackid] = colour;
	}

      g->Draw("Psame");
      tpc1_used.insert(trackid);
    }

  if(draw_tpc_edges)
    {
      TLine *l1 = new TLine(0, 0, 0, 3400);
      l1->SetLineWidth(3);
      l1->SetLineColor(kBlack);
      l1->Draw();

      TLine *l2 = new TLine(n_wires, 0, n_wires, 3400);
      l2->SetLineWidth(3);
      l2->SetLineColor(kBlack);
      l2->Draw();
    }

  bottom->cd();
  base_hist->Draw();
  t0->Draw();

  for(auto const& [trackid, hits] : hitMap0)
    {
      if(hits.size() < 5)
	continue;

      const simb::MCParticle* part = particleInv->TrackIdToParticle_P(trackid);
      int pdg = part == NULL ? 2 : part->PdgCode();
      if(part != NULL && fTrackIds.count(part->TrackId()) == 0)
        pdg = 2;

      if(pdg == 22)
	{
	  const simb::MCParticle *mother = particleInv->TrackIdToParticle_P(part->Mother());
	  if(mother != NULL && mother->PdgCode() == 111)
	    pdg = 22111;
	}

      TGraph *g = new TGraph();

      for(auto const& hit : hits)
        g->SetPoint(g->GetN(), hit->WireID().deepestIndex(), hit->PeakTime());

      if(pdg == 2)
	g->SetMarkerColor(kGray);
      else
	{
	  const int colour = fTrueColourMap.count(trackid) == 0 ? fColours[fTrueColourCounter] : fTrueColourMap[trackid];
	  g->SetMarkerColor(colour);

	  if(tpc1_used.count(trackid) == 0)
	    {
	      TLegendEntry *le = leg->AddEntry(g, fPDGString[pdg], "");
	      le->SetTextColor(colour);
	    }
	  if(fTrueColourMap.count(trackid) == 0)
	    ++fTrueColourCounter;
	    
	  fTrueColourMap[trackid] = colour;
	}

      g->Draw("Psame");
    }

  if(draw_tpc_edges)
    {
      TLine *l1 = new TLine(0, 0, 0, 3400);
      l1->SetLineWidth(3);
      l1->SetLineColor(kBlack);
      l1->Draw();

      TLine *l2 = new TLine(n_wires, 0, n_wires, 3400);
      l2->SetLineWidth(3);
      l2->SetLineColor(kBlack);
      l2->Draw();
    }

  leg->SetNColumns(fTrueColourCounter);
  leg->Draw();
}

void sbnd::EventDisplayAndDump::DrawRecoView(const HitMap &hitMap0, const HitMap &hitMap1, const TString &name, const bool draw_tpc_edges)
{
  TPaveText *t = new TPaveText(.9, .95, .98, 1, "NB NDC");
  t->SetTextAlign(32);
  t->SetTextSize(0.085);
  t->SetTextColor(kBlack);
  t->SetFillStyle(4000);
  t->AddText(name);

  TPaveText *t0 = new TPaveText(.03, .92, .2, .95, "NB NDC");
  t0->SetTextAlign(12);
  t0->SetTextSize(0.06);
  t0->SetTextColor(kBlack);
  t0->SetFillStyle(4000);
  t0->AddText("TPC 0");

  TPaveText *t1 = new TPaveText(.03, .82, .2, .85, "NB NDC");
  t1->SetTextAlign(12);
  t1->SetTextSize(0.06);
  t1->SetTextColor(kBlack);
  t1->SetFillStyle(4000);
  t1->AddText("TPC 1");

  TPad *top    = new TPad("top", "", 0.005, 0.505, 0.995, 0.995);
  TPad *bottom = new TPad("bottom", "", 0.005, 0.005, 0.995, 0.495);
  top->SetBottomMargin(0.002);
  bottom->SetTopMargin(0.002);
  top->Draw();
  bottom->Draw();

  top->cd();

  const double n_wires = name == "W view" ? 1664 : 1984;

  TH1F* base_hist = new TH1F("base_hist", "", 10, 0.5 - 0.5*(1984 - n_wires), 1984.5 - 0.5*(1984 - n_wires));
  base_hist->SetMaximum(3400.0);
  base_hist->Draw();
  t->Draw();
  t1->Draw();

  while(hitMap0.size() + hitMap1.size() > fColours.size())
    fColours.insert(fColours.end(), fColours.begin(), fColours.end());

  for(auto const& [pfpid, hits] : hitMap1)
    {
      TGraph *g = new TGraph();

      for(auto const& hit : hits)
        g->SetPoint(g->GetN(), hit->WireID().deepestIndex(), 3400-hit->PeakTime());

      const int colour = fRecoColourMap.count(pfpid) == 0 ? fColours[fRecoColourCounter] : fRecoColourMap[pfpid];
      g->SetMarkerColor(colour);
      if(fRecoColourMap.count(pfpid) == 0)
	++fRecoColourCounter;
      fRecoColourMap[pfpid] = colour;

      g->Draw("Psame");
    }

  if(draw_tpc_edges)
    {
      TLine *l1 = new TLine(0, 0, 0, 3400);
      l1->SetLineWidth(3);
      l1->SetLineColor(kBlack);
      l1->Draw();

      TLine *l2 = new TLine(n_wires, 0, n_wires, 3400);
      l2->SetLineWidth(3);
      l2->SetLineColor(kBlack);
      l2->Draw();
    }

  bottom->cd();
  base_hist->Draw();
  t0->Draw();

  for(auto const& [pfpid, hits] : hitMap0)
    {
      TGraph *g = new TGraph();

      for(auto const& hit : hits)
        g->SetPoint(g->GetN(), hit->WireID().deepestIndex(), hit->PeakTime());

      const int colour = fRecoColourMap.count(pfpid) == 0 ? fColours[fRecoColourCounter] : fRecoColourMap[pfpid];
      g->SetMarkerColor(colour);
      if(fRecoColourMap.count(pfpid) == 0)
	++fRecoColourCounter;
      fRecoColourMap[pfpid] = colour;

      g->Draw("Psame");
    }

  if(draw_tpc_edges)
    {
      TLine *l1 = new TLine(0, 0, 0, 3400);
      l1->SetLineWidth(3);
      l1->SetLineColor(kBlack);
      l1->Draw();

      TLine *l2 = new TLine(n_wires, 0, n_wires, 3400);
      l2->SetLineWidth(3);
      l2->SetLineColor(kBlack);
      l2->Draw();
    }
}

bool sbnd::EventDisplayAndDump::VolumeCheck(const TVector3 &pos, const double &walls, const double &cath, const double &front, const double &back)
{
  const bool xedges = pos.X() < (200. - walls) && pos.X() > (-200. + walls);
  const bool yedges = pos.Y() < (200. - walls) && pos.Y() > (-200. + walls);
  const bool zedges = pos.Z() < (500. - back)  && pos.Z() > (0. + front);
  const bool caths  = pos.X() > cath || pos.X() < -cath;

  return xedges && yedges && zedges && caths;
}

double sbnd::EventDisplayAndDump::YZtoU(const double y, const double z)
{
  return z * cosU - y * sinU;
}

double sbnd::EventDisplayAndDump::YZtoV(const double y, const double z)
{
  return z * cosV - y * sinV;
}

double sbnd::EventDisplayAndDump::YZtoW(const double y, const double z)
{
  return z * cosW - y * sinW;
}

DEFINE_ART_MODULE(sbnd::EventDisplayAndDump)

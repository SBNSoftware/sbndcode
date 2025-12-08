#include <string>
#include <vector>

#include "TF1.h"
#include "TGraph.h"
#include "TH1D.h"

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "sbncode/CAFMaker/FillTrue.h"
#include "sbncode/WireMod/Utility/WireModUtility.hh"

namespace wiremod {


bool ide_in_tpc(const sim::IDE* ide_ptr) {
    // avoid larsoft exception if IDE is outside TPC
    return (ide_ptr->x > -199 && ide_ptr->x < 199)
        && (ide_ptr->y > -199 && ide_ptr->y < 199)
        && (ide_ptr->z > 1 && ide_ptr->z < 499);
}

// version of this function in sbncode folds the angle and assumes degrees, and
// our convention is different for the plane angles -- larsoft returns pi/6 for
// plane angles but splines are done using pi/3 we're also using radians until
// the final spline lookup
double ThetaXW(double dxdr, double dydr, double dzdr, int tpc, int ip, bool fold=false) {
    // TODO use the version in WireMod utility library
    double planeAngle = 0;
    if (tpc == 0 && ip == 0) {
        dzdr *= -1;
        planeAngle = util::pi() / 3;
    }
    else if (tpc == 0 && ip == 1) {
        dzdr *= 1;
        planeAngle = util::pi() / 3;
    }
    else if (tpc == 1 && ip == 0) {
        dzdr *= 1;
        planeAngle = util::pi() / 3;
    }
    else if (tpc == 1 && ip == 1) {
        dzdr *= -1;
        planeAngle = util::pi() / 3;
    }
    double sinPlaneAngle = std::sin(planeAngle);
    double cosPlaneAngle = std::cos(planeAngle);

    double dzdrPlaneRel = dzdr * cosPlaneAngle + dydr * sinPlaneAngle;

    double theta = std::atan(dxdr / dzdrPlaneRel);
    if (!fold) {
        return theta;
    }
    return (std::abs(theta) > 0.5 * util::pi()) ? util::pi() - std::abs(theta) : std::abs(theta);
}



class WireModifier : public art::EDProducer
{
public:
    struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> HitLabel {
            Name("HitLabel"), Comment("Label for hits")
        };
        fhicl::Atom<art::InputTag> WireLabel {
            Name("WireLabel"), Comment("Label for wires")
        };
        fhicl::Atom<art::InputTag> SimChannelLabel {
            Name("SimChannelLabel"), Comment("Label for sim channels")
        };
        fhicl::Atom<art::InputTag> MCPartAssnLabel {
            Name("MCPartAssnLabel"), Comment("Label for MC particle,hit association")
        };
        fhicl::Atom<art::InputTag> WireAssnLabel {
            Name("WireAssnLabel"), Comment("Label for wire,hit association")
        };
        fhicl::Atom<bool> ApplyYZScale {
            Name("ApplyYZScale"), Comment("Apply scaling based on Y and Z")
        };
        fhicl::Atom<bool> ApplyXXZAngleScale {
            Name("ApplyXXZAngleScale"), Comment("Apply scaling based on X and theta_XW")
        };
        fhicl::Atom<std::string> SplineFileXTXW_Q {
            Name("SplineFileXTXW_Q"), Comment("ROOT file containing TGraph2D for ADC scaling in each plane by X and Theta_XW")
        };
        fhicl::Atom<std::string> SplineFileXTXW_W {
            Name("SplineFileXTXW_W"), Comment("ROOT file containing TGraph2D for width scaling in each plane by X and Theta_XW")
        };
        fhicl::Atom<std::string> SplineFileYZ_Q {
            Name("SplineFileYZ_Q"), Comment("ROOT file containing TGraph2D for ADC scaling in each plane by Y and Z")
        };
        fhicl::Atom<std::string> SplineFileYZ_W {
            Name("SplineFileYZ_W"), Comment("ROOT file containing TGraph2D for width scaling in each plane by Y and Z")
        };
        fhicl::Atom<float> MaxThetaXW {
            Name("MaxThetaXW"), Comment("Maximum value of a hit's Theta-XW to apply WireMod")
        };
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit WireModifier(Parameters const& config);

    void produce(art::Event& evt) override;

protected:

private:
    const geo::GeometryCore* fGeometry = lar::providerFrom<geo::Geometry>();
    const geo::WireReadoutGeom* fWireReadout = &(art::ServiceHandle<geo::WireReadout const>()->Get());

    std::string fSplineFileXTXW_Q;
    std::string fSplineFileXTXW_W;
    std::string fSplineFileYZ_Q;
    std::string fSplineFileYZ_W;
    float fMaxThetaXW;
    art::InputTag fHitLabel; 
    art::InputTag fWireLabel; 
    art::InputTag fSimChannelLabel; 
    art::InputTag fMCPartAssnLabel; 
    art::InputTag fWireAssnLabel; 
    bool applyYZScale;
    bool applyXXZAngleScale;

    std::vector<TGraph2D*> splines_x_txw_q;
    std::vector<TGraph2D*> splines_x_txw_w;
    std::vector<TGraph2D*> splines_y_z_q;
    std::vector<TGraph2D*> splines_y_z_w;

}; // end WireModifier class

WireModifier::WireModifier(Parameters const& config) :
    EDProducer{config}, 
    fSplineFileXTXW_Q(config().SplineFileXTXW_Q()),
    fSplineFileXTXW_W(config().SplineFileXTXW_W()),
    fSplineFileYZ_Q(config().SplineFileYZ_Q()),
    fSplineFileYZ_W(config().SplineFileYZ_W()),
    fMaxThetaXW(config().MaxThetaXW()),
    fHitLabel(config().HitLabel()),
    fWireLabel(config().WireLabel()),
    fSimChannelLabel(config().SimChannelLabel()),
    fMCPartAssnLabel(config().MCPartAssnLabel()),
    fWireAssnLabel(config().WireAssnLabel()),
    applyYZScale(config().ApplyYZScale()), 
    applyXXZAngleScale(config().ApplyXXZAngleScale())
{

    cet::search_path sp("FW_SEARCH_PATH");

    auto load_splines = [&sp](const std::string& fname, std::vector<TGraph2D*>& spl) {
        std::string fname_fullpath;
        if (!sp.find_file(fname, fname_fullpath)) {
            // try WireMod subdirectory within FW_SEARCH_PATH
            if (!sp.find_file("WireMod/" + fname, fname_fullpath)) {
                throw std::runtime_error("Could not find spline file with name " + fname + " in FW_SEARCH_PATH");
            }
        }
        TFile* f = TFile::Open(fname_fullpath.c_str());
        MF_LOG_INFO("WireModifier") <<
            "Loading splines from " << fname_fullpath.c_str();
        for (UInt_t i = 0; i < 6; i++) {
            std::string spline_name(Form("splines_%d", i));
            TGraph2D* spline = static_cast<TGraph2D*>(f->Get(spline_name.c_str()));
            if (!spline) {
                throw std::runtime_error("Could not read spline with name " + spline_name + " from file " + fname_fullpath);
            }
            spl.push_back(static_cast<TGraph2D*>(spline->Clone()));
        }
        f->Close();
    };
    load_splines(fSplineFileXTXW_Q, splines_x_txw_q);
    load_splines(fSplineFileXTXW_W, splines_x_txw_w);
    load_splines(fSplineFileYZ_Q, splines_y_z_q);
    load_splines(fSplineFileYZ_W, splines_y_z_w);

    produces<std::vector<recob::Wire>>();

    if (fMaxThetaXW < 0. || fMaxThetaXW > 90.) {
        throw std::runtime_error("MaxThetaXW value should be between 0 and 90 (not " + std::to_string(fMaxThetaXW) + ")");
    }
}


void WireModifier::produce(art::Event& evt) {
    art::ServiceHandle<art::TFileService> tfs;

    const auto det_prop =
        art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt);

    // TODO can we refactor this into the constructor?
    sys::WireModUtility wmUtil(fGeometry, fWireReadout, det_prop);
    wmUtil.applyYZScale = applyYZScale;
    wmUtil.applyXXZAngleScale = applyXXZAngleScale;

    art::Handle<std::vector<recob::Wire>> wire_handle;
    evt.getByLabel(fWireLabel, wire_handle);
    auto const& wire_vec(*wire_handle);

    art::Handle<std::vector<recob::Hit>> hit_handle;
    evt.getByLabel(fHitLabel, hit_handle);
    auto const& hit_vec(*hit_handle);

    wmUtil.FillROIMatchedHitMap(hit_vec, wire_vec);
    MF_LOG_INFO("WireModifier")
        << "Got Hit Map.";

    // ------
    const auto det_clocks =
        art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

    art::ServiceHandle<cheat::BackTrackerService> bt_serv;
    const cheat::BackTrackerService* bt = bt_serv.get();

    // sim channels
    const auto simchannel_handle =
        evt.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelLabel);
    std::vector<art::Ptr<sim::SimChannel>> simch_ptr_vec;
    art::fill_ptr_vector(simch_ptr_vec, simchannel_handle);

    // hits
    const auto hit_ptr_handle =
        evt.getValidHandle<std::vector<recob::Hit>>(fHitLabel);
    std::vector<art::Ptr<recob::Hit>> hit_ptr_vec;
    art::fill_ptr_vector(hit_ptr_vec, hit_ptr_handle);

    // TODO this changes interface between versions
    // std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE*>>> id_to_ide_map =
    std::map<int, std::vector<sbn::ReadoutIDE>> id_to_ide_map =
        sbn::PrepSimChannels(simch_ptr_vec, *fWireReadout);
        // TODO this changes namespace between versions
        // caf::PrepSimChannels(simch_ptr_vec, *fWireReadout);

    std::map<int, std::vector<art::Ptr<recob::Hit>>> id_to_truehit_map =
        // TODO this changes namespace between versions
        // caf::PrepTrueHits(hit_ptr_vec, det_clocks, *bt);
        sbn::PrepTrueHits(hit_ptr_vec, det_clocks, *bt);


    art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mc_part_for_hits(
            hit_handle, evt, fMCPartAssnLabel);
    // ------

    // loop over wires
    std::unique_ptr<std::vector<recob::Wire>> new_wires(new std::vector<recob::Wire>());
    for (size_t i_w = 0; i_w < wire_vec.size(); i_w++) {
        const auto& wire = wire_vec.at(i_w);

        recob::Wire::RegionsOfInterest_t new_rois;
        new_rois.resize(wire.SignalROI().size());

        // loop over ROIs
        // std::cout << "\nROI Loop\n";
        size_t nranges = wire.SignalROI().get_ranges().size();
        for (size_t i_r = 0; i_r < nranges; i_r++) {
            const auto& roi_range = wire.SignalROI().get_ranges()[i_r];

            // std::cout << "-";
            auto roi_properties = wmUtil.CalcROIProperties(wire, i_r);
            MF_LOG_INFO("WireModifier")
                << "    ROI Properties:" << '\n'
                << "                    key:     (" << roi_properties.key.first << ", " << roi_properties.key.second << ")" << '\n'
                << "                    view:    " << roi_properties.view << '\n'
                << "                    begin:   " << roi_properties.begin << '\n'
                << "                    end:     " << roi_properties.end << '\n'
                << "                    total_q: " << roi_properties.total_q << '\n'
                << "                    center:  " << roi_properties.center << '\n'
                << "                    sigma:   " << roi_properties.sigma;

            // find hits within the ROI
            sys::WireModUtility::ROI_Key_t roi_key(wire.Channel(), i_r);
            auto it_hit_map = wmUtil.ROIMatchedHitMap.find(roi_key);
            if (it_hit_map == wmUtil.ROIMatchedHitMap.end()) {
                MF_LOG_INFO("WireModifier") << "No hits within ROI.";
                // TODO: Handle this case, in case a new hit would appear
                // for now: just keep the original ROI
                new_rois.add_range(roi_properties.begin, roi_range.data());
                continue;
            }

            std::vector<const recob::Hit*> roi_hit_vec;
            for (auto i_h : it_hit_map->second) {
                roi_hit_vec.push_back(&hit_vec[i_h]);
            }
            auto sub_roi_properties = wmUtil.CalcSubROIProperties(roi_properties, roi_hit_vec);

            std::map<sys::WireModUtility::SubROI_Key_t, sys::WireModUtility::ScaleValues_t> scale_map;
            for (size_t i_sr = 0; i_sr != sub_roi_properties.size(); i_sr++) {
                // calculate truth & perform spline interpolation
                // first set some defaults
                sys::WireModUtility::ScaleValues_t scale_vals;
                scale_vals.r_Q = 1.;
                scale_vals.r_sigma = 1.;
                auto subroi_prop = sub_roi_properties.at(i_sr);
                auto key = subroi_prop.key;
                scale_map[key] = scale_vals;

                auto hit = roi_hit_vec.at(i_sr);
                art::Ptr<recob::Hit> hit_ptr;
                for (const auto& hp : hit_ptr_vec) {
                    if (hp->Channel() != hit->Channel() || hp->PeakTime() != hit->PeakTime()) continue;
                    hit_ptr = hp;
                    break;
                }

                // hit associated with particle?
                auto mc_parts = mc_part_for_hits.at(hit_ptr.key());
                if (mc_parts.size() == 0) continue;

                // particle in IDE map?
                const auto particle = mc_parts.at(0);
                if (!id_to_ide_map.count(particle->TrackId())) continue;

                // TODO particle direction taken as the true hit direction
                // this will be replaced in future productions when we have
                // access to SimEnergyDeposits
                const auto mcp_dir = particle->Momentum().Vect().Unit();

                // TODO this interface changes with version
                // std::vector<std::pair<geo::WireID, const sim::IDE*>> id_to_ide =
                const std::vector<sbn::ReadoutIDE> id_to_ide =
                    id_to_ide_map.at(particle->TrackId());

                double factor_q = 1.0;
                double factor_w = 1.0;

                // defined here just printing below
                double factor_yz_q = 1.0;
                double factor_yz_w = 1.0;
                double factor_xtxw_q = 1.0;
                double factor_xtxw_w = 1.0;

                // TODO for SBND SPRING PRODUCTION ONLY using sim::IDEs
                // We will use WireMod Utility functions for most of this once
                // we have access to SimEnergyDeposits
                for (const auto& id_ide : id_to_ide) {
                    const auto id = id_ide.wire;
                    if (wire.Channel() != fWireReadout->PlaneWireToChannel(id)) continue;
                    const auto ide_ptr = id_ide.ide;
                    if (!wiremod::ide_in_tpc(ide_ptr)) continue;

                    // get the plane & corresponding spline
                    const geo::TPCGeo& tpc_geom = fGeometry->PositionToTPC({ ide_ptr->x, ide_ptr->y, ide_ptr->z });
                    const auto plane = fWireReadout->Plane(tpc_geom.ID(), wire.View());
                    unsigned int ip = plane.ID().Plane;

                    size_t plane_idx = ip + 3 * tpc_geom.ID().TPC;

                    // spline is data/MC ratio, we want to scale MC to match data
                    if (wmUtil.applyYZScale) { 
                        TGraph2D* spline_yz_q = splines_y_z_q.at(plane_idx);
                        TGraph2D* spline_yz_w = splines_y_z_w.at(plane_idx);
                        factor_yz_q = spline_yz_q->Interpolate(ide_ptr->y, ide_ptr->z);
                        factor_yz_w = spline_yz_w->Interpolate(ide_ptr->y, ide_ptr->z);
                        factor_q *= factor_yz_q;
                        factor_w *= factor_yz_w;
                    }

                    Double_t txw = std::numeric_limits<double>::quiet_NaN();
                    if (wmUtil.applyXXZAngleScale) { 
                        // x vs theta_xw angle
                        txw = (180.0 / util::pi()) * wiremod::ThetaXW(mcp_dir.X(), mcp_dir.Y(), mcp_dir.Z(), tpc_geom.ID().TPC, ip);
                        if (std::abs(txw) < fMaxThetaXW) {
                            TGraph2D* spline_xtxw_q = splines_x_txw_q.at(plane_idx);
                            TGraph2D* spline_xtxw_w = splines_x_txw_w.at(plane_idx);
                            factor_xtxw_q = spline_xtxw_q->Interpolate(ide_ptr->x, txw);
                            factor_xtxw_w = spline_xtxw_w->Interpolate(ide_ptr->x, txw);
                            factor_q *= factor_xtxw_q;
                            factor_w *= factor_xtxw_w;
                        }
                    }
                
                    MF_LOG_INFO("WireModifier")
                        << "IDE info\n"
                        << " - x: " << ide_ptr->x << "\n"
                        << " - y: " << ide_ptr->y << "\n"
                        << " - z: " << ide_ptr->z << "\n"
                        << " - dirx: " << mcp_dir.X() << "\n"
                        << " - diry: " << mcp_dir.Y() << "\n"
                        << " - dirz: " << mcp_dir.Z() << "\n"
                        << " - PlaneAng: " << plane.ThetaZ() << "\n"
                        << " - txw: " << txw << "\n"
                        // << " - tyz: " << tyz << "\n"
                        << " - TPC: " << tpc_geom.ID() << "\n"
                        << " - Plane: " << ip << "\n"
                        << " - spline index: " << plane_idx << "\n"
                        << " - YZ Spline (q, w): " << factor_yz_q << " " << factor_yz_w << "\n"
                        << " - XThetaXW Spline (q, w): " << factor_xtxw_q << " " << factor_xtxw_w << "\n"
                        << " - Total Spline (q, w): " << factor_q << " " << factor_w << "\n";
                    break;
                }

                scale_vals.r_Q *= factor_q;
                scale_vals.r_sigma *= factor_w;
                scale_map[key] = scale_vals;
            }

            std::vector<float> modified_data(roi_range.data());
            wmUtil.ModifyROI(modified_data, roi_properties, sub_roi_properties, scale_map);
            new_rois.add_range(roi_properties.begin, modified_data);
        }
        new_wires->emplace_back(new_rois, wire.Channel(), wire.View());
    }

    evt.put(std::move(new_wires));
}

DEFINE_ART_MODULE(WireModifier)
}


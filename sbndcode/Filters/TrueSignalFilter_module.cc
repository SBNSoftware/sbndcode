/*
 * Filter module for common true signal definition options
 * Each filter block can check
 * - True nu flavors
 * - CC/NC
 * - Modes
 * - Fiducial vertex
 * - Final state primary PDG
 * - Final state primary KE thresholds
 * - Disallowed PDG list
 * - Exact counts? Exact match to primary list, or allow others (except those listed in exclude list or below threshold)
 *
 * Multiple filters blocks per event are supported. Event is kept if it passes
 * any filter block
 */

#include <iostream>
#include <algorithm>
#include <optional>
#include <array>

#include "TGeoManager.h"

#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Tuple.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "sbndcode/RecoUtils/RecoUtils.h"

#include "messagefacility/MessageLogger/MessageLogger.h"


namespace filt{

// All optional. Empty block config will pass all events
struct FilterBlockConfig {
    fhicl::OptionalSequence<int> NuPDGs {
        fhicl::Name("NuPDGs"),
        fhicl::Comment("Allowed PDG codes for the neutrino")
    };

    fhicl::OptionalSequence<int> TargetPDGs {
        fhicl::Name("TargetPDGs"),
        fhicl::Comment("Allowed PDG codes for the target")
    };

    fhicl::OptionalSequence<float, 2> WRange {
        fhicl::Name("WRange"),
        fhicl::Comment("Range of allowed invariant hadronic mass (W) as (min, max). Use -1 for min or max to indicate one-sided bound.")
    };

    fhicl::OptionalAtom<bool> InTPC {
        fhicl::Name("InTPC"),
        fhicl::Comment("Require interaction vertex in the TPC")
    };

    fhicl::OptionalAtom<bool> IsCC {
        fhicl::Name("IsCC"),
        fhicl::Comment("If true, only CC events are accepted. If false, only NC events are accepted. If ommitted, no requirement")
    };

    fhicl::OptionalAtom<bool> KeepBadStatus {
        fhicl::Name("KeepBadStatus"),
        fhicl::Comment("If true, keep particles with a bad status code. This is needed for the production of unstable primaries like the eta meson")
    };

    fhicl::OptionalSequence<int> Modes {
        fhicl::Name("Modes"),
        fhicl::Comment("List of allowed interaction modes")
    };

    fhicl::OptionalSequence<int> RequiredPDGs {
        fhicl::Name("RequiredPDGs"),
        fhicl::Comment("List of PDG codes that must appear in the event as primaries. May repeat PDG codes to require multiple particles of the same type")
    };

    fhicl::OptionalAtom<bool> ExactCounts {
        fhicl::Name("ExactCounts"),
        fhicl::Comment("If true, event must contain the exact count of particles listed in the required PDG list with no extra primaries, except those in the ignored list or below KE threshold")
    };
    
    fhicl::OptionalSequence<int> IgnoredPDGs {
        fhicl::Name("IgnoredPDGs"),
        fhicl::Comment("List of PDG codes which are not considered when counting primaries, even if they are above KE thresholds")
    };

    fhicl::OptionalSequence<int> DisallowedPDGs {
        fhicl::Name("DisallowedPDGs"),
        fhicl::Comment("List of PDG codes which must not be present in the list of primaries if they are above KE thresholds")
    };

    fhicl::OptionalSequence<fhicl::Tuple<int, float>> KEThresholds {
        fhicl::Name("KEThresholds"),
        fhicl::Comment("Minimum KE required for particles to count towards the list of primaries. Format is [PDG code, threshold (MeV)]")
    };
};


struct FilterBlock {
    std::optional<std::vector<int>> nu_pdgs;
    std::optional<std::vector<int>> target_pdgs;
    std::optional<std::array<float, 2>> wrange;
    std::optional<bool> in_tpc;
    std::optional<bool> iscc;
    std::optional<bool> keep_bad_status;
    std::optional<std::vector<int>> modes;
    std::optional<std::vector<int>> required_pdgs;
    std::optional<std::vector<int>> ignored_pdgs;
    std::optional<std::vector<int>> disallowed_pdgs;
    std::optional<std::map<int, float>> ke_thresholds;
    std::optional<bool> exact_counts;

    // computed from options
    // (pdg, nrequired)
    std::map<int, int> required_counts;
    float wmin = -1.0;
    float wmax = -1.0;
};


class TrueSignalFilter : public art::EDFilter {
public:
    struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag> GENIELabel {
            Name("GENIELabel"), Comment("Label for MC truth (GENIE)")
        };
        fhicl::Sequence<fhicl::Table<FilterBlockConfig>> Filters {
             Name("Filters"),
             Comment("List of filter-condition blocks; event passes if it matches any block")
         };
    };
    using Parameters = art::EDFilter::Table<Config>;
    explicit TrueSignalFilter(const Parameters& config);

    virtual bool filter(art::Event& e) override;

protected:
    bool PassBlock(const art::Ptr<simb::MCTruth>&, const FilterBlock&, mf::LogDebug&) const;
    void PrintBlock(const FilterBlock&) const;

private:
    art::InputTag fGenieModuleLabel;
    std::vector<FilterBlock> fFilterBlocks;
};


TrueSignalFilter::TrueSignalFilter(const Parameters& pset) :
    EDFilter{pset}, fGenieModuleLabel{pset().GENIELabel()}
{
    // loop over filter blocks
    auto const& cfg_blocks = pset().Filters();
    fFilterBlocks.reserve(cfg_blocks.size());

    for (std::size_t i = 0; i < cfg_blocks.size(); ++i) {
        auto const& cfg = cfg_blocks.at(i);
        FilterBlock block;

        std::vector<int> nu_pdgs;
        if (cfg.NuPDGs(nu_pdgs)) {
            block.nu_pdgs = nu_pdgs;
        }

        std::vector<int> target_pdgs;
        if (cfg.TargetPDGs(target_pdgs)) {
            block.target_pdgs = target_pdgs;
        }

        std::array<float, 2> wrange;
        if (cfg.WRange(wrange)) {
            block.wrange = wrange;
            block.wmin = wrange.at(0);
            block.wmax = wrange.at(1);

            if (block.wmin > block.wmax && block.wmax > 0) {
                throw cet::exception("TrueSignalFilter")
                    << "In Filters[" << i << "]: WRange ("
                    << block.wmin << ", " << block.wmax << ") invalid\n";
            }

        }

        bool in_tpc;
        if (cfg.InTPC(in_tpc)) {
            block.in_tpc = in_tpc;
        }

        bool iscc;
        if (cfg.IsCC(iscc)) {
            block.iscc = iscc;
        }
        bool keep_bad_status;
        if (cfg.KeepBadStatus(keep_bad_status)) {
            block.keep_bad_status = keep_bad_status;
        }

        std::vector<int> modes;
        if (cfg.Modes(modes)) {
            block.modes = modes;
        }

        std::vector<int> required_pdgs;
        if (cfg.RequiredPDGs(required_pdgs)) {
            block.required_pdgs = required_pdgs;

            // Count multiplicities in the requirements.
            for (int pdg : required_pdgs) {
                block.required_counts[pdg]++;
            }
        }

        std::vector<int> disallowed_pdgs;
        if (cfg.DisallowedPDGs(disallowed_pdgs)) {
            block.disallowed_pdgs = disallowed_pdgs;
        }

        std::vector<std::tuple<int, float>> ke_thresholds;
        if (cfg.KEThresholds(ke_thresholds)) {
            std::map<int, float> ke_threshold_map;
            for (const auto& pdg_threshold : ke_thresholds) {
                int pdg = std::get<0>(pdg_threshold);
                float threshold = std::get<1>(pdg_threshold);
                auto [it, inserted] = ke_threshold_map.emplace(pdg, threshold);

                // error on duplicates
                if (!inserted) {
                    throw cet::exception("TrueSignalFilter")
                        << "In Filters[" << i << "]: KEThreshold duplicate "
                        << "entry for PDG " << pdg << ".\n";
                }
            }
            block.ke_thresholds = ke_threshold_map;
        }

        bool exact_counts;
        if (cfg.ExactCounts(exact_counts)) {
            block.exact_counts = exact_counts;
        }

        PrintBlock(block);
        fFilterBlocks.push_back(std::move(block));
    }
}


bool TrueSignalFilter::PassBlock(const art::Ptr<simb::MCTruth>& mc, const FilterBlock& block, mf::LogDebug& debug_log) const {
    const std::string kSpace = "    ";
    debug_log << kSpace << "Checking MCTruth is neutrino... ";
    if (mc->Origin() != simb::kBeamNeutrino) return false;
    debug_log << "PASS\n";

    const auto& nu = mc->GetNeutrino();

    if (block.nu_pdgs) {
        int nu_pdg = nu.Nu().PdgCode();
        debug_log << kSpace << "Checking MCTruth has appropriate neutrino PDG (Nu PDG=" << nu_pdg << ")... ";
        if (std::find(block.nu_pdgs->begin(), block.nu_pdgs->end(), nu_pdg)
                == block.nu_pdgs->end()) return false;
        debug_log << "PASS\n";
    }

    if (block.target_pdgs) {
        int target_pdg = nu.Target();
        debug_log << kSpace << "Checking MCTruth has appropriate target PDG (Target PDG=" << target_pdg << ")... ";
        if (std::find(block.target_pdgs->begin(), block.target_pdgs->end(), target_pdg)
                == block.target_pdgs->end()) return false;
        debug_log << "PASS\n";
    }

    if (block.wrange) {
        float w = nu.W();
        debug_log << kSpace << "Checking MCTruth has W within range (W=" << w << ")... ";
        if (block.wmin >= 0. && w < block.wmin) return false;
        if (block.wmax >= 0. && w > block.wmax) return false;
        debug_log << "PASS\n";
    }

    if (block.modes) {
        int mode = nu.Mode();
        debug_log << kSpace << "Checking MCTruth has appropriate mode (mode=" << mode << ")... ";
        if (std::find(block.modes->begin(), block.modes->end(), mode)
                == block.modes->end()) return false;
        debug_log << "PASS\n";
    }

    if (block.iscc) {
        debug_log << kSpace << "Checking MCTruth CC/NC (CCNC=" << nu.CCNC() << ")... ";
        if (block.iscc.value() && (nu.CCNC() != simb::kCC)) return false;
        if (!block.iscc.value() && (nu.CCNC() != simb::kNC)) return false;
        debug_log << "PASS\n";
    }

    if (block.in_tpc) {
        TVector3 vtx{ nu.Nu().Vx(), nu.Nu().Vy(), nu.Nu().Vz() };
        debug_log << kSpace << "Checking MCTruth is inside the TPC (x=" << vtx(0) << ", y=" << vtx(1) << ", z=" << vtx(2) << ")... ";
        // technically user could request vertices outside the tpc?
        if (RecoUtils::IsInsideTPC(vtx, 0) != block.in_tpc.value()) return false;
        debug_log << "PASS\n";
    }

    //get a vector of final state particles for the next few checks
    std::vector<int> primaries;
    bool has_thresholds = block.ke_thresholds.has_value();
    bool has_ignored = block.ignored_pdgs.has_value();
    for (int i = 0; i < mc->NParticles(); ++i) {
        const auto& part(mc->GetParticle(i));
        debug_log << kSpace << "MC particle list " << i << " " << part.PdgCode() << " STATUS=" << part.StatusCode() << ", E=" << part.E() << "\n";
        // only consider primaries
        if (part.StatusCode() != 1 && !block.keep_bad_status) continue;

        // don't count particles in the ignored list
        int pdg = part.PdgCode();
        if (has_ignored) {
            if (std::find(block.ignored_pdgs->begin(), block.ignored_pdgs->end(), pdg)
                    != block.ignored_pdgs->end()) continue;
        }

        // don't count particles below threshold
        if (has_thresholds) {
            if (block.ke_thresholds->find(pdg) != block.ke_thresholds->end()) {
                float ke = 1.0e3 * (part.E() - part.Mass());
                if (ke < block.ke_thresholds->at(pdg)) continue;
            }
        }
        primaries.push_back(part.PdgCode());
    }

    if (block.required_pdgs) {
        debug_log << kSpace << "Checking list of primaries [";
        // check that primaries contains all the required PDGs
        // count multiplicities in the event
        std::map<int,int> counts;
        for (int pdg : primaries) {
            counts[pdg]++;
        }

        for (const auto& [k, v] : counts) {
            float threshold = 0.;
            if (has_thresholds) {
                if (block.ke_thresholds->find(k) != block.ke_thresholds->end()) {
                    threshold = block.ke_thresholds->at(k);
                }
            }
            debug_log << "\n" << kSpace << kSpace << "PDG=" << k
                << ", count above " << threshold << " MeV threshold=" << v;
        }
        debug_log << "\n" << kSpace << "]... ";

        // Check that for every required PDG, the event has enough.
        // if exact_counts: check equality
        bool exact_counts = block.exact_counts.value_or(false);
        for (const auto& [required_pdg, nrequired] : block.required_counts) {
            auto it = counts.find(required_pdg);
            if (it == counts.end() || it->second < nrequired) return false;
            if (exact_counts && it->second != nrequired) return false;
        }

        // check if there are extra particles not in the required list if exact_counts is set
        // note: counts already excludes below threshold particles & ignored particles
        if (exact_counts) {
            for (const auto& [primary_pdg, count] : counts) {
                // extra particle not in the list
                if (block.required_counts.find(primary_pdg) == block.required_counts.end()) return false;
            }
        }

        debug_log << "PASS\n";
    }

    if (block.disallowed_pdgs) {
        debug_log << kSpace << "Checking no primaries are disallowed... ";
        // check that the primaries list doesn't contain anything disallowed
        // note: primaries still are only counted if they are above threshold
        for (int disallowed_pdg : *block.disallowed_pdgs) {
            if (std::find(primaries.begin(), primaries.end(), disallowed_pdg)
                    != primaries.end()) return false;
        }
        debug_log << "PASS\n";
    }

    return true;
}


bool TrueSignalFilter::filter(art::Event & e) {
    // pass debug log reference to PassBlock so that we can print messages to
    // the same line without adding too many timestamps to the output logs, and
    // we also only have to write "FAIL" once in the code.
    // We get a logger message once for every event
    mf::LogDebug debug_log{"TrueSignalFilter"};

    auto mclists = e.getMany<std::vector<simb::MCTruth>>();
    for (size_t i = 0; i != mclists.size(); i++) {
        const auto& mclist(mclists.at(i));
        for (size_t j = 0; j != mclist->size(); j++) {
            art::Ptr<simb::MCTruth> mc(mclist, j);
            debug_log << "Filtering MCTruth " << j << "...\n";
            for (size_t ib = 0; ib != fFilterBlocks.size(); ib++) {
                const auto& block = fFilterBlocks.at(ib);
                if (PassBlock(mc, block, debug_log)) {
                    // mctruth passes at least one filter, so we are done
                    debug_log << "MCTruth " << j << " passed filter " << ib << "!\n";
                    return true;
                }
                else {
                    debug_log << "FAIL\n";
                    debug_log << "MCTruth " << j << " did not pass filter " << ib << "\n";
                }
            }
        }
    }

    // did not pass any filters
    return false;
}


void TrueSignalFilter::PrintBlock(const FilterBlock& block) const {
    const std::string kNotSet("<not set>");

    // print helpers
    mf::LogInfo log{"TrueSignalFilter"};
    log << "Filter configuration:\n";
    auto print_int_vec = [&](const std::string& name, std::optional<std::vector<int>> const& val) {
        log << " - " << name << ": ";
        if (!val) {
            log << kNotSet << "\n";
        }
        else {
            log << "{ ";
            const auto& vec = val.value();
            for (size_t i = 0; i != vec.size(); i++) {
                log << vec.at(i);
                if (i < vec.size() - 1) log << ", ";
            }
            log << " }\n";
        }
    };
    auto print_bool = [&](const std::string& name, std::optional<bool> const& val) {
        log << " - " << name << ": ";
        if (!val) {
            log << kNotSet << "\n";
        }
        else {
            log << (val.value() ? "True" : "False") << "\n";
        }
    };
    auto print_ke_thresholds = [&](const std::string& name, const std::optional<std::map<int, float>>& val) {
        log << " - " << name << ": ";
        if (!val) {
            log << kNotSet << "\n";
        }
        else {
            log << "\n";
            for (const auto& [k, v] : val.value()) {
                log << "    - " << k << ": " << v << " MeV\n";
            }
        }
    };

    print_int_vec("NuPDGs", block.nu_pdgs);
    print_int_vec("TargetPDGs", block.target_pdgs);
    print_bool("InTPC", block.in_tpc);
    print_bool("IsCC", block.iscc);
    print_int_vec("Modes", block.modes);
    print_int_vec("RequiredPDGs", block.required_pdgs);
    print_int_vec("IgnoredPDGs", block.ignored_pdgs);
    print_int_vec("DisallowedPDGs", block.disallowed_pdgs);
    print_bool("ExactCounts", block.exact_counts);
    print_ke_thresholds("KEThresholds", block.ke_thresholds);
    // TODO WRange
}

DEFINE_ART_MODULE(TrueSignalFilter)

}

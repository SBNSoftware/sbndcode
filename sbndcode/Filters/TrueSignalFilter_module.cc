/*
 * Filter module for common true signal definition options
 * - True nu flavors
 * - CC/NC
 * - Modes
 * - Fiducial vertex
 * - Final state primary PDG
 * - Final state primary KE
 * - Exclude PDG list
 * - Exclusive? Exact match to primary list, or allow others (except those listed in exclude list)
 *
 * Can pass multiple filter options with fcl parameters
 * Example: Pass CC events with no final state proton OR NC events with final
 * state proton (not sure why you would want to do this, but it's possible!)
 * {
 *   nu: 14
 *   cc: true
 *   NoFinalStatePDG: [ 2112 ] 
 * }
 * {
 *   nu: 14
 *   cc: False
 *   FinalStatePDG: [ 2112 ] 
 * }
 */

#include <iostream>
#include <algorithm>
#include <optional>

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
        fhicl::Comment("PDG codes of the neutrino")
    };

    fhicl::OptionalAtom<bool> InTPC {
        fhicl::Name("InTPC"),
        fhicl::Comment("Require interaction vertex in the TPC")
    };


    fhicl::OptionalAtom<bool> IsCC {
        fhicl::Name("IsCC"),
        fhicl::Comment("If true, only CC events are accepted. If false, only NC events are accepted. If ommitted, no requirement")
    };

    fhicl::OptionalSequence<int> Modes {
        fhicl::Name("Modes"),
        fhicl::Comment("List of interaction modes")
    };

    fhicl::OptionalSequence<int> RequiredPDGs {
        fhicl::Name("RequiredPDGs"),
        fhicl::Comment("List of PDG codes that must appear in the event as primaries. May repeat PDG codes for multiple particles of the same type")
    };

    fhicl::OptionalSequence<fhicl::Tuple<int, float>> KEThresholds {
        fhicl::Name("KEThresholds"),
        fhicl::Comment("Minimum KE for particles to count towards either required or disallowed particle lists. Format is [PDG code, threshold (MeV)]")
    };

    fhicl::OptionalAtom<bool> Exclusive {
        fhicl::Name("Exclusive"),
        fhicl::Comment("If true, event must contain the exact required PDGs with no extra primary particles")
    };

    fhicl::OptionalSequence<int> DisallowedPDGs {
        fhicl::Name("DisallowedPDGs"),
        fhicl::Comment("List of PDG codes which must not be present in the primaries (not used if \"Exclusive\" option is true)")
    };
};


struct FilterBlock { 
    std::optional<std::vector<int>> nu_pdgs;
    std::optional<bool> in_tpc;
    std::optional<bool> iscc; 
    std::optional<std::vector<int>> modes;
    std::optional<std::vector<int>> required_pdgs;
    std::optional<std::vector<int>> disallowed_pdgs;
    std::optional<std::map<int, float>> ke_thresholds;
    std::optional<bool> exclusive;
    
    // computed from options
    // (pdg, nrequired)
    std::map<int, int> required_counts;
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
    bool PassBlock(const art::Ptr<simb::MCTruth>, const FilterBlock&) const;
    void PrintBlock(const FilterBlock&) const;

private:
    art::InputTag fGenieModuleLabel; 
    std::vector<FilterBlock> fFilterBlocks;
};


TrueSignalFilter::TrueSignalFilter(const Parameters& pset) :
    EDFilter{pset}, fGenieModuleLabel(pset().GENIELabel())
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

        bool in_tpc; 
        if (cfg.InTPC(in_tpc)) {
            block.in_tpc = in_tpc;
        }

        bool iscc; 
        if (cfg.IsCC(iscc)) {
            block.iscc = iscc;
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
            for (auto& pdg_threshold : ke_thresholds) {
                std::map<int, float> ke_thresholds;
                int pdg = std::get<0>(pdg_threshold);
                float threshold = std::get<1>(pdg_threshold);
                auto [it, inserted] = ke_thresholds.emplace(pdg, threshold);

                // error on duplicates
                if (!inserted) {
                    throw cet::exception("TrueSignalFilter")
                        << "In Filters[" << i << "]: KEThreshold duplicate "
                        << "entry for PDG " << pdg << ".\n";
                }

                block.ke_thresholds = ke_thresholds;
            }
        }

        bool exclusive; 
        if (cfg.Exclusive(exclusive)) {
            block.exclusive = exclusive;
        }

        PrintBlock(block);
        fFilterBlocks.push_back(std::move(block));
    }
}


bool TrueSignalFilter::PassBlock(const art::Ptr<simb::MCTruth> mc, const FilterBlock& block) const {
    if (mc->Origin() != simb::kBeamNeutrino) return false;

    const auto& nu = mc->GetNeutrino();

    if (block.nu_pdgs) {
        if (std::find(block.nu_pdgs->begin(), block.nu_pdgs->end(), nu.Nu().PdgCode())
                == block.nu_pdgs->end()) return false;
    }

    if (block.modes) {
        if (std::find(block.modes->begin(), block.modes->end(), nu.Mode())
                == block.modes->end()) return false;
    }

    if (block.in_tpc) {
        // technically user could request vertices outside the tpc?
        TVector3 vtx{ nu.Nu().Vx(), nu.Nu().Vy(), nu.Nu().Vz() };
        if (RecoUtils::IsInsideTPC(vtx, 0) != block.in_tpc.value()) return false;
    }

    if (block.iscc) {
        if (block.iscc.value() && (nu.CCNC() != simb::kCC)) return false;
        if (!block.iscc.value() && (nu.CCNC() != simb::kNC)) return false;
    }

    //get a vector of final state particles for the next few checks
    std::vector<int> primaries;
    for (int i = 0; i < mc->NParticles(); ++i) {
        const auto& part(mc->GetParticle(i));
        // only consider primaries
        if (part.StatusCode() != 1) continue;
        if (part.Mother() > 0) continue;

        // check against threshold map. Primaries not counted if below threshold
        if (block.ke_thresholds) {
            int pdg = part.PdgCode();
            if (block.ke_thresholds->find(pdg) != block.ke_thresholds->end()) {
                float ke = part.E() - part.Mass();
                if (ke < block.ke_thresholds->at(pdg)) continue;
            }
        }
        primaries.push_back(part.PdgCode());
    }

    if (block.required_pdgs) {
        // check that primaries contains all the required PDGs
        // count multiplicities in the event
        std::map<int,int> counts;
        for (int pdg : primaries) {
            counts[pdg]++;
        }

        // Check that for every required PDG, the event has enough.
        // if exclusive: check equality
        bool exclusive = block.exclusive.value_or(false);
        for (auto const& [required_pdg, nrequired] : block.required_counts) {
            auto it = counts.find(required_pdg);
            if (it == counts.end() || it->second < nrequired) return false;
            if (exclusive && it->second != nrequired) return false;
        }

        // TODO check if there are extra particles not in the required list if exclusive
    }

    if (block.disallowed_pdgs) {
        // check that the primaries list doesn't contain anything disallowed
        // note: primaries still are only counted if they are above threshold
        for (int disallowed_pdg : *block.disallowed_pdgs) {
            if (std::find(primaries.begin(), primaries.end(), disallowed_pdg)
                    != primaries.end()) return false;
        }
    }

    return true;
}


bool TrueSignalFilter::filter(art::Event & e) {
    auto mclists = e.getMany<std::vector<simb::MCTruth>>();
    for (size_t i = 0; i != mclists.size(); i++) {
        const auto& mclist(mclists.at(i));
        for (size_t j = 0; j != mclist->size(); j++) {
            art::Ptr<simb::MCTruth> mc(mclist, j);
            for (auto const& block : fFilterBlocks) {
                if (PassBlock(mc, block)) {
                    // mctruth passes at least one filter, so we are done
                    return true;
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
            for (int i : val.value()) {
                log << i << ", ";
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

    print_int_vec("NuPDGs", block.nu_pdgs);
    print_bool("InTPC", block.in_tpc);
    print_bool("IsCC", block.iscc);
    print_int_vec("Modes", block.modes);
    print_int_vec("RequiredPDGs", block.required_pdgs);
    print_int_vec("DisallowedPDGs", block.disallowed_pdgs);
    print_bool("Exclusive", block.exclusive);
}

DEFINE_ART_MODULE(TrueSignalFilter)

}

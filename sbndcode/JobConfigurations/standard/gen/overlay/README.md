# GENIE+Cosmics Event Generation

There are several fhicl configurations for producing GENIE events. The are divided in several sub-directories `other_flux_config/flux_config*/`, depending on the type of flux files that are used as input. The flux configurations are described here: https://sbnsoftware.github.io/sbndcode_wiki/The_SBND_flux_files.html.

If you are uncertain, use the fhicl files in this directory.

| FHiCL Name | Neutrinos | Cosmics | Comments |
| ------------- | ------------- | ------------- | ------------- |
| `prodgenie_corsika_cmc_nu_spill_cryo_sbnd.fcl` | GENIE, BNB flux, cryostat only | CORSIKA CMC model |  |
| `prodgenie_corsika_cmc_nu_spill_tpc_sbnd.fcl` | GENIE, BNB flux, TPC only | CORSIKA CMC model |  |
| `prodgenie_corsika_p_nu_spill_cryo_sbnd.fcl` | GENIE, BNB flux, cryostat only | CORSIKA proton model |  |
| `prodgenie_corsika_p_nu_spill_tpc_sbnd.fcl` | GENIE, BNB flux, TPC only | CORSIKA proton model |  |
| `prodgenie_corsika_p_rockbox_sbnd.fcl` | GENIE, BNB flux | CORSIKA proton model | Dirt generation with rockbox model |
| `prodgenie_corsika_p_rockbox_intrnue_sbnd.fcl` | GENIE, nue only | CORSIKA proton model | Dirt generation with rockbox model |
| `prodgenie_corsika_p_rockbox_sce_sbnd.fcl` | GENIE, nue only | CORSIKA proton model | Dirt generation with rockbox model, includes Space Charge Effect |
| `prodgenie_corsika_p_rockbox_sce_keep_corsika_trajectories_sbnd.fcl` | GENIE, BNB flux | CORSIKA proton model | Dirt generation with rockbox model, keeps all particle trajectories for GENIE and CORSIKA, includes Space Charge Effect |
| `prodgenie_corsika_p_rockbox_sce_no_shower_rollup_sbnd.fcl` | GENIE, BNB flux | CORSIKA proton model | Dirt generation with rockbox model, doesn't drop secondary e+/e-/gamma in showers, includes Space Charge Effect |
| - | - | - | - |
| `prodgibuu_corsika_p_dirtpropagation_sbnd.fcl` | GiBUU, BNB flux | CORSIKA proton model | To be used with GiBUU workflow |


# GENIE+Cosmics Event Generation

There are several fhicl configurations for producing GENIE events. The are divided in several sub-directories `other_flux_config/flux_config*/`, depending on the type of flux files that are used as input. The flux configurations are described here: https://sbnsoftware.github.io/sbndcode_wiki/The_SBND_flux_files.html.

If you are uncertain, use the fhicl files in this directory.

| FHiCL Name | Generator | Nu Flavor | Nu Generation Volume | Comments |
| ------------- | ------------- | ------------- | ------------- |
| `prodgenie_corsika_cmc_nu_spill_cryo_sbnd.fcl` | GENIE + CORSIKA CMC | BNB flux | Cryostat | CORSIKA CMC model |
| `prodgenie_corsika_cmc_nu_spill_tpc_sbnd.fcl` | GENIE + CORSIKA CMC | BNB flux | TPC | CORSIKA CMC model |
| `prodgenie_corsika_p_nu_spill_cryo_sbnd.fcl` | GENIE + CORSIKA proton | BNB flux | Cryostat | CORSIKA proton model |
| `prodgenie_corsika_p_nu_spill_tpc_sbnd.fcl` | GENIE + CORSIKA proton | BNB flux | TPC | CORSIKA proton model |
| `prodgenie_corsika_p_rockbox_sbnd.fcl` | GENIE + CORSIKA proton | BNB flux | Dirt | Dirt generation with rockbox model |
| `prodgenie_corsika_p_rockbox_intrnue_sbnd.fcl` | GENIE + CORSIKA proton | nue | Dirt | Dirt generation with rockbox model |
| `prodgenie_corsika_p_rockbox_sce_sbnd.fcl` | GENIE + CORSIKA proton | nue | Dirt | Dirt generation with rockbox model, includes Space Charge Effect |
| `prodgenie_corsika_p_rockbox_sce_keep_corsika_trajectories_sbnd.fcl` | GENIE + CORSIKA proton | BNB flux | Dirt | Dirt generation with rockbox model, keeps all particle trajectories for GENIE and CORSIKA, includes Space Charge Effect |
| `prodgenie_corsika_p_rockbox_sce_no_shower_rollup_sbnd.fcl` | GENIE + CORSIKA proton | BNB flux | Dirt | Dirt generation with rockbox model, doesn't drop secondary e+/e-/gamma in showers, includes Space Charge Effect |
| - | - | - | - |
| `prodgibuu_corsika_p_dirtpropagation_sbnd.fcl` | GiBUU + CORSIKA proton | BNB flux | - | To be used with GiBUU workflow |


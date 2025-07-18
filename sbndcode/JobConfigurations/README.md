# SBND Core FHiCL Files

- 25th February 2025 (Marco Del Tutto)
`standard_g4_[rockbox_,intime_]sbnd.fcl` has been promoted back to an up-to-date fcl and has rejoined the standard workflow. SCE services are enabled by default.
`standard_detsim_sbnd.fcl` has been promoted back to an up-to-date fcl and has rejoined the standard workflow. SCE services are enabled by default.
- 4th October 2024 (Dom Brailsford)
`standard_reco2_sbnd.fcl` has been promoted back to an up-to-date fcl and has rejoined the standard workflow.  SCE services are enabled by default
- 15th July 2024 (Dom Brailsford)
`standard_reco1_sbnd.fcl` has been promoted back to an up-to-date fcl so can be used as part of any standard workflow.  The below suggested workflows have been updated to include this information.


The intended operation of the fcl workflows is that the `standard-*` fcls run the standard workflow. The WireCell 2D TPC simulation/signal processing workflow has now been implemented into the `standard-*` fcls, and the 1D simulation is now deprecated. The *standard* workflow also includes the dropping of some heavy data products.

At the time of writing, the core workflow (for BNB + Dirt + Cosmics) is the following:

- `prodgenie_corsika_proton_rockbox_sbnd.fcl`
- `standard_g4_rockbox_sbnd.fcl` (does not include TPC electron drift simulation)
- `standard_detsim_sbnd.fcl` (includes TPC drift simulation, TPC electronics simulation, and signal processing)
- `standard_reco1_sbnd.fcl`
- `standard_reco2_sbnd.fcl`

The intime workflow is as follows:

- `prodcorsika_proton_intime_sbnd.fcl`
- `standard_g4_intime_sbnd.fcl`
- `standard_detsim_sbnd.fcl`
- `standard_reco1_sbnd.fcl`
- `standard_reco2_sbnd.fcl`

For single generator workflows (like intrinsic neutrino samples) the workflow is the following:

- `<your-gen>.fcl`
- `standard_g4_sbnd.fcl`
- `standard_detsim_sbnd.fcl`
- `standard_reco1_sbnd.fcl`
- `standard_reco2_sbnd.fcl`

This may well change over the coming months, and this README should be updated to reflect this.

** ALERT **

Due to changes implemented in sbndcode PRs #408 and #409 the 1D simulation fcls will not work out of the box, they will need editing by experts!

Henry Lay, Lynn Tung, Bear Carlson - Feb 2024

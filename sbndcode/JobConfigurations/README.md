# SBND Core FHiCL Files

The intended operation of the fcl workflows is that the `standard-*` fcls run the standard workflow. This isn't currently true and should be acknowledged here. The WireCell 2D TPC simulation/signal processing workflow has now been implemented into the `standard-*` fcls, and the 1D simulation is now deprecated. The main deviation of the *core* workflow (described below) from the `standard-*` fcls is the inclusion of the space charge simulation. The *core* workflow also includes the dropping of some heavy data products (hence the `lite` suffix).

At the time of writing, the core workflow (for BNB + Dirt + Cosmics) is the following:

- `prodoverlay_corsika_cosmics_proton_genie_rockbox_sce.fcl`
- `g4_sce_dirt_filter_lite.fcl` (does not include TPC electron drift simulation)
- `detsim_sce_overlay_lite.fcl` (includes TPC drift simulation, TPC electronics simulation, and signal processing)
- `reco1_sce_lite.fcl`
- `reco2_sce.fcl`

For single generator workflows (like intrinsic neutrino samples) the workflow is the follow:

- `<your-gen>.fcl`
- `g4_sce_lite.fcl`
- `detsim_sce_lite.fcl`
- `reco1_sce_no_overlay_lite.fcl`
- `reco2_sce.fcl`

This may well change over the coming months, and this README should be updated to reflect this.

** ALERT **

Due to changes implemented in sbndcode PRs #408 and #409 the 1D simulation fcls will not work out of the box, they will need editing by experts!

Henry Lay, Lynn Tung, Bear Carlson - Feb 2024

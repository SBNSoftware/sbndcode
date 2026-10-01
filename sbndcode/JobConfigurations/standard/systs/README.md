# DENT detector variation

DENT variation: data-driven cathode (cathode filter), data-driven distorted
E-field / displacement map (`SBND_DataMap_v4.root`) and measured drift
velocity (0.15532 cm/us).

## Workflow

1. CV:

   ```
   prodgenie_corsika_proton_rockbox0p1_sbnd.fcl
   standard_g4_rockbox_sbnd.fcl
   standard_detsim_sbnd.fcl
   standard_reco1_sbnd.fcl
   standard_reco2_sbnd.fcl
   cafmakerjob_sbnd_sce_systtools_and_fluxwgt_and_g4rw.fcl
   ```

2. DENT, starting from the gen-stage output of step 1:

   ```
   standard_g4_rockbox_sbnd_cathodefilt.fcl
   standard_detsim_sbnd_cathodefilt.fcl
   standard_reco1_sbnd.fcl
   standard_reco2_sbnd.fcl
   cafmakerjob_sbnd_sce_systtools_and_fluxwgt_and_g4rw.fcl
   ```

DENT is run separately from the gen-stage output, even though the
SimEnergyDeposits are saved, because the WireCell simulation is
non-deterministic. The DENT variation therefore includes the WireCell
randomness, to avoid a statistics shortage.

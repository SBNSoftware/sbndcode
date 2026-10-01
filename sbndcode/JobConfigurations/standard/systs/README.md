# DENT detector systematic (cathode filter)

FHiCL wrappers for producing the DENT detector-variation MC samples. Compared
with the CV workflow, a DENT sample changes three things:

1. **Data-driven cathode.** `CathodeFilterSBND` (`sbndcode/CathodeFilter`)
   removes `sim::SimEnergyDeposit`s that lie inside the measured 3D cathode
   volume (`sbnd_cathode_v8.txt`). The filter runs in the g4 stage between
   `largeant` and `ionandscint`, so it acts on un-shifted deposit positions,
   before SCE is applied.
2. **Data-driven distorted E-field / displacement map:**
   `SCEoffsets/SBND_DataMap_v4.root` (from sbnd_data v01_43_00).
3. **Measured drift velocity, 0.15532 cm/us:** `DriftVelFudgeFactor`
   0.99373 in g4 and detsim, and WireCell `driftSpeed` 1.5532 in detsim.

The stock geometry (`sbnd_v02_06`) is used. Its nominal CPA lies entirely
inside the data-driven cathode volume.

## Production workflow

| Stage  | fcl                                              |
|--------|--------------------------------------------------|
| gen    | `prodgenie_corsika_proton_rockbox0p1_sbnd.fcl`   |
| g4     | `standard_g4_rockbox_sbnd_cathodefilt.fcl`       |
| detsim | `standard_detsim_sbnd_cathodefilt.fcl`           |
| reco1  | `standard_reco1_sbnd.fcl`                        |
| reco2  | `standard_reco2_sbnd.fcl`                        |
| caf    | `cafmakerjob_sbnd_sce_systtools_and_fluxwgt.fcl` |

Only g4 and detsim differ from CV; gen, reco1, reco2 and caf are the standard
fcls. For example:

```
lar -c prodgenie_corsika_proton_rockbox0p1_sbnd.fcl -n 10 -o gen.root
lar -c standard_g4_rockbox_sbnd_cathodefilt.fcl -s gen.root -o g4.root
lar -c standard_detsim_sbnd_cathodefilt.fcl     -s g4.root  -o detsim.root
lar -c standard_reco1_sbnd.fcl                  -s detsim.root -o reco1.root
lar -c standard_reco2_sbnd.fcl                  -s reco1.root  -o reco2.root
lar -c cafmakerjob_sbnd_sce_systtools_and_fluxwgt.fcl -s reco2.root
```

The gen stage of `prodgenie_corsika_proton_rockbox0p1_sbnd.fcl` runs Geant4
with a dirt filter, so fewer events reach g4 than were generated.

### The chain must start from gen

Do not scrub existing CV reco1 files (`scrub_g4_wcls_detsim_reco1.fcl`) and
re-run these wrappers on them:

- In the rockbox workflow Geant4 runs at gen, and the cathode filter reads its
  `largeant:LArG4DetectorServicevolTPCActive` deposits. The CV detsim output
  drops all `largeant` SimEnergyDeposits (`detsim_drops.fcl`), so they are not
  in the reco1 files.
- A scrubbed file still remembers the `G4`, `DetSim` and `Reco1` process
  names, so these stages cannot reuse them.

## Files

| File | Purpose |
|------|---------|
| `geometry_sbnd_cathodefilt.fcl` | Shared overrides included by the wrappers: drift velocity and SCE map. No geometry change despite the name, which is kept from the v10_06 version. |
| `standard_g4_rockbox_sbnd_cathodefilt.fcl` | g4 for rockbox gen. Inserts `cathodefiltered` before `ionandscint` and drops the unfiltered `largeant` TPC-active deposits from the output. |
| `standard_g4_sbnd_cathodefilt.fcl` | Same as above for generators that do not run `largeant` at gen, such as `singlegen_anode_muons.fcl`. |
| `standard_detsim_sbnd_cathodefilt.fcl` | Standard detsim with the measured drift velocity. |
| `standard_reco1_sbnd_cathodefilt.fcl`, `standard_reco2_sbnd_cathodefilt.fcl` | Pass-throughs to the standard reco fcls. |
| `prodgenie_corsika_proton_rockbox_sbnd_cathodefilt.fcl` | Pass-through to the standard gen fcl, kept for the v10_06 fcl names. |
| `standard_cathodefilter_sbnd.fcl` | Deprecated no-op stage, kept for the v10_06 fcl names. |
| `singlegen_anode_muons.fcl` | Particle gun of muons crossing both anodes, for studies of the cathode model. |

## Validation

Ran with sbndcode v10_14_02_0602 and sbnd_data v01_43_00, 2 events through
the full workflow above: every stage exits 0, and the CAF files have 2
`recTree` entries. The filter removes about 1–1.4% of TPC-active deposits.

Ported from SBNSoftware/sbndcode#926 (production/v10_06_00).

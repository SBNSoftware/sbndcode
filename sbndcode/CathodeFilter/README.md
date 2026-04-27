# sbndcode/CathodeFilter

Data-driven representation of the SBND cathode (CPA) volume, reconstructed
from cosmic-muon track endpoints, and an art::EDProducer that uses it to
drop `sim::SimEnergyDeposit`s that lie inside the cathode volume.

## Files
- `SBNDCathode.h / .cxx`  — C++17 class that loads the ASCII volume file and
  exposes `contains(x, y, z)`, `x_neg`, `x_pos`, `thickness`, etc.
- `CathodeSimEnergyDepositFilter_module.cc` — art producer; reads a
  `sim::SimEnergyDeposit` collection and writes a filtered collection
  that excludes every deposit whose MidPoint (and optionally Start/End)
  lies inside the cathode volume.
- `fcl/cathode_simenergydeposit_filter.fcl` — default PROLOG.
- `fcl/run_cathode_filter.fcl` — drop-in job that wires the producer
  around any SBND simulation input file.
- `data/sbnd_cathode_v6.txt` — 160 × 200 grid of (x_neg, x_pos) in cm,
  installed to `$FW_SEARCH_PATH/sbndcode/CathodeFilter/sbnd_cathode_v6.txt`.

## Build
Append `add_subdirectory(CathodeFilter)` to `sbndcode/sbndcode/CMakeLists.txt`.
Then rebuild sbndcode:

```
mrb z
mrbsetenv
mrb i -j4
```

## Run
```
lar -c run_cathode_filter.fcl -s <input.root> -n 10
```

## Parameters (FHiCL)
| name | type | default | meaning |
| --- | --- | --- | --- |
| `SimEnergyDepositLabel` | `art::InputTag` | `"ionandscint"` | input sim::SimEnergyDeposit collection |
| `CathodeVolumeFile`     | `std::string`   | *required*       | fw-search-path-resolved ASCII volume file |
| `UseStepEndpoints`      | `bool`          | `false`          | also reject deposits whose Start or End lies inside the cathode (in addition to the MidPoint test) |
| `NanPolicy`             | `std::string`   | `"outside"`      | `outside` / `inside` / `throw` at (y,z) bins with no endpoint data |
| `Interp`                | `std::string`   | `"bilinear"`     | `bilinear` / `nearest` |
| `Verbose`               | `bool`          | `false`          | log kept/dropped counts per event |

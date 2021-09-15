# Crossing muon tracks

These modules will run over samples, create crossing muon tracks products, filter over events, and finally create a hitdumper tree with branches containing 3D endpoints, t0 (for applicable tracks), and trajcetory info.

## Run command

To run the job, you only have to use one command:

`lar -c run_muontrack.fcl <sample.root>`

Output: `muon_hitdumper.root`

## Relevant info

These modules use a new data product that is located within sbnobj. Make sure that your local version of sbnobj contains the Commissioning directory in the SBND subdirectory, and has the files for the MuonTrack data product.

Relevant files (mostly located within subdirectory `/sbndcode/Commissioning/...`)

- **MuonTrackProducer_module.cc**: produces MuonTrack objects
- **MuonTrackFilter_module.cc**: filters for events that contain MuonTrack objects
- **muontrackmodule.fcl**: contains parameters for MuonTrackProducer + MuonTrackFilter
  - by default, the producer will create *all* types of muon tracks. To create only anode-cathode crossing muon tracks for example, set `KeepMuonTypes: [0]`
- **run_muontrack.fcl**: fcl for running the job, also overwrites default for the hitdumper tree

**HitDumper_module.cc** and **hitdumpermodule.fcl** have both been changed to include parameters and branches for the muon tracks. The hitdumpertree now includes branches that contain track information such as t0, 3d endpoints, trajectory, tpc, and track type. These branches are labeled: `muontrk_<var>`.

## Caveats

- t0 information is only available for tracks with a unique t0. Top-bottom,  up-downstream, and uncategorized tracks do not have unique t0's, and will be assigned a default value of -999.
- similarly, 3D endpoint information is not accurate for tracks of type 5 (uncategorized). The endpoint information in the tree for these types of tracks are the wire intersection endpoints on the collection/induction planes.

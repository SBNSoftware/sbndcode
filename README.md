# Anode-cathode crossing muon tracks 
Relevant files (all located within subdirectory `/sbndcode/Commissioning/...`)

- `ACProducer_module.cc`: produces information for AC tracks contained within the CRTTrack data product.
- `ACFilter_module.cc`: Filters for events that contain AC tracks 
- `run_acproducer.fcl`: contains relevant parameters for Hough Transform. Path contains: fasthit, ACProducer, and ACFilter 
    - The defaults are appropriate ONLY for performing an *anode-cathode crossing* muon selection cut, not any other crossing tracks. 

`HitDumper_module.cc` and `hitdumpermodule.fcl` have both been changed to include parameters and branches for the AC tracks. The hitdumpertree includes branches that contain AC information such as t0, 3d endpoints, and trajectory. 

NOTE: A new data product is under development for AC tracks, but the information is currently stored in a CRTTrack data product. The endpoints are intuitively stored in variables `x1_pos` and `z2_pos` etc, but t0 (in us) for the AC muon is stored in variable `ts1_ns` and theta_xz and theta_yz are stored in `thetaxy` and `phizy` respectively. If you are looking only at the hitdumper output, this shouldn't interfere with anything.

# Congifuration files
- run_decohitfinder.fcl: runs deconvolution+ophit finder + opflash finder (input is artroot file with opdaq)
- opdeconvolution_sbnd.fcl: configuration file for SBNDOpDeconvolution_module (produces deconvolved raw::OpDetWaveform)
- flashfinder_deco_sbnd.fcl: contains OpHit and OpFlass algorithm configurations to be used with the deconvolved waveforms

# Directories
- Alg: deconvolution algorithms and modules
- reco_deco_sbnd: conatins refactored reco1 and reco2 fhicl files including deconvolution for PDS

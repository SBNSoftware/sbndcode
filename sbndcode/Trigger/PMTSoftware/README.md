# PMT Software Trigger

This directory contains the producer, filter, and analyzer modules for SBND's PMT Software Trigger (also known as PMT Beam metrics). 

- `PMTMetricProducer` is nearly identical to the version of the producer that runs in the DAQ. The DAQ version can be found in: `sbndaq-artdaq/ArtModules/SBND/SoftwareTrigger/PMTMetricProducer_module.cc`. The input is CAEN1730 Fragments. 
- `PMTMCMetricProducer` runs on the detector simulation. The input is `raw::OpDetWaveform`s from the PMTs. 
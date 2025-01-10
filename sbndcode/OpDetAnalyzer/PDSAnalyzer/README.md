# Overview

`SBNDPDSAnalyzer` is an `art::analyzer` that saves useful information for PDS MC/reco studies at different stages of the simulation:
- Gen: `MCTruth`
- G4: `MCParticles` and deposited energy
- Detsim: `OpDetWaveforms`
- Reco1: deconvolved waveforms (`OpDetWaveforms`), `OpHits` and `OpFlashes`
- Reco2: flash matchers scores (`SimpleFlash` and `OpT0Finder`) and `CRUMBS` score

Three TTrees are created:
- `PDSMapTree`: Creates a tree with the optical channel IDs, position (X, Y and Z in cm) and the photon detector type (0: CoatedPMT, 1: UncoatedPMT, 2: VUV XARAPUCA and 3: VIS XARAPUCA).
- `OpAnaTree`: TTree with the analysis variables (described below). Each entry in the TTree corresponds to an event.
- `OpAnaPerTrackTree`: TTree with the number of SimPhotons produced by each MC particle. Each entry in the TTree corresponds to one MC particle in one event. Requires SimPhotons (not SimPhotonsLite).

Configuration parameters: 
- `SaveMCTruth`, `SaveMCParticles`, ...: booleans to save different data products at different stages of the simulation/reconstruction chain.
- `Verbosity`: shows some printous while running the analyzer.
- `MakePerTrackTree`: option to make a tree with the number of `SimPhotons` per `MCParticle` (`OpAnaPerTrackTree`).
- `MakePDSGeoTree`: option to dump the PDS mapping in a TTree (`PDSMapTree`)
- `UseSimPhotonsLite`: use `SimPhotonsLite` or `SimPhotons` (default is true, i.e. use `SimPhotonsLite`)
- `KeepPDGCode`: specify what PDG codes will be stored in the TTree from the MC particles list. If no PDGs are specified, it saves all the MC particles (default is [], so it saves everything).
- `MCTruthOrigin`: specify the `MCTruth` origins that will be saved. Default is [1], i.e. neutrino generated events.
- `MCTruthPDG`: specify the PDG codes of the candidate vertex particles for the `MCTruth`. Defaults is [12, -12, 14, -14], i.e. BNB neutrinos.
- `PDTypes`: vector specifying the PD types for which it will save the `SimPhotons` (default is ["pmt_coated", "pmt_uncoated"])
- `G4BufferBoxX`, `G4BufferBoxY`, `G4BufferBoxZ`: only store MC particle trajectories inside the boundaries defined by the G4BufferBox variables.
- `G4BeamWindow`: only store MC particles and trajectories in the time window defined by this parameter. Default is [-10000, 12000] #ns.


# Variables in the TTree

### Variables at the particle generation stage
  
  | Branch name                 |    Type                                     |  Description |             
  |-----------------------------|---------------------------------------------|--------------|  
  | nuvX, nuvY, nuvZ | std::vector\<double\> | True neutrino interaction vertex (in cm) |
  | nuvT | std::vector\<double\> | True neutrino interaction time (in ns) |
  | nuvE | std::vector\<double\> | True neutrino energy (in GeV) |

  
  They store information regarging all the primary neutrino interactions in a given readout window. Loop over the std::vector to get the energy, interaction vertex and interaction time corresponding to the $i_{th}$ neutrino interaction.

### Variables at the LArG4 stage (particle propagation+ionization+scintillation)
  
  #### Variables regarding the MC particles propagated by G4
  | Branch name                 |    Type                                     |  Description |             
  |-----------------------------|---------------------------------------------|--------------|  
  | E | std::vector\<double\> | Primary energy of each MC particle (in GeV)|
  | process | std::vector\<std::string\> | Primary process of MC each particle|
  | trackID | std::vector\<int\> |  MC particle ID|
  | motherID | std::vector\<int\> | ID of the mother MC particle|
  | PDGCode | std::vector\<int\> | PDG code |
  | InTimeCosmics | bool | Returns true if there is a cosmic interaction during the BNB spill |
  
  Example: Consider a 1207 MeV beam $\nu_\mu$. It undergoes a charged current interaction, leading to a $\mu^-$ (E=574 MeV) and a proton (E=1543 MeV). The first entry in the previous vectors will have the following values:
  E[0]=0.574, process[0]=primary, trackID[0]=1, motherID[0]=0, PDGCode[0]=13,  E[1]=1.543, PDGCode[1]=2212... If e.g. the muon decays to a Michel electron, there will be an entry with PDGCode=11, motherID=1 and process=Decay.

  #### Variables regarding the deposited energy
  
  | Branch name                 |    Type                                     |  Description |             
  |-----------------------------|---------------------------------------------|--------------|  
  | energydep | std::vector\<std::vector<double\>\> | Energy deposition (in MeV) at each G4 tracking step. It's saved for each MC particle|
  | energydepX, energydepY, energydepZ | std::vector\<std::vector<double\>\> | Location (in cm) of each energy deposition. |
  | dEpromX, dEpromY, dEpromZ | std::vector<double\> | Average X, Y, Z (in cm) location of the energy depositions. It's saved for the two TPCs (vector size is 2)|
  | dEspreadX, dEspreadY, dEspreadZ | std::vector<double\> | X, Y, Z standard deviation of the energy depositions. It's saved for the two TPCs (vector size is 2)|
  | dElowedges, dEmaxedges | std::vector\<std::vector<double\>\> | (X, Y, Z) coordinates of the lowest (max) energy deposition. It's saved for the two TPCs (vector size is 2)|
  
  Following the previous example, to read the energy deposition values and their locations induced by the primary proton you can take energydep[1], energydepX[1]...
  
  #### Variables regarding the scintillation photons

  
  | Branch name                 |    Type                                     |  Description |             
  |-----------------------------|---------------------------------------------|--------------|  
  | SimPhotonsperOpChVUV | std::vector\<double\> | Number of true photons at each PD (VUV wavelength)|
  | SimPhotonsperOpChVIS | std::vector\<double\> | Number of true photons at each PD (visible wavelength)|
  | SimPhotonsLiteVUV | std::vector\<std::vector<double\>\> | Photon arrival times at G4 stage (VUV). In ns.|  
  | SimPhotonsLiteVIS | std::vector\<std::vector<double\>\> | Photon arrival times at G4 stage (VIS). In ns.|  
  | NPhotons* variables| double | Integrated number of photons in the events per PD type. |  
  
  SBND has 312 PDs, hence the size of the SimPhotonsperOpChVUV(VIS) is 312. You can obtain the number of VUV photons reaching the coated PMT with ID 144 by taking SimPhotonsperOpChVUV[144]. The size of the SimPhotonsLiteVUV(VIS) is also 312. Each vector in the 'vector of vecrtors' contains the times (in ns) in which each photon gets to the given PD. Imagine 567 VUV photons reach the coated PMT with ID 144. The size of the vector SimPhotonsLiteVUV[144] will be 567.
  
  
### Variables at the digitization stage
| Branch name                 |    Type                                     |  Description |             
|-----------------------------|---------------------------------------------|--------------|  
| SignalsDigi | std::vector\<std::vector<double\>\> | Digitized signals (ADC values) |  
| StampTime | std::vector<double\> | Start time of each digitized waveform (in $\mu$s) |  
| OpChDigi | std::vector<int\> | Associated PD ID |  

The PMT/XARAPUCA output signals (including electronic response) are stored in the previous vectors. We only save the regions of the waveforms going above a certain thereshold (region of interest or ROIs). The ADC values of each identifeid ROI correspond to an entry in the SignalsDigi "vector of vectors". Note that we may have more than one ROI per PD, so the size of the SignalsDigi branch will be in general different than the number of PDs (312). To get the start time and the channel corresponding to the $i_{th}$ ROI get the StampTime and OpChDigi with index $i_{th}$.

### Variables at the reconstruction stage

#### Deconvolution
| Branch name                 |    Type                                     |  Description |             
|-----------------------------|---------------------------------------------|--------------|  
| SignalsDeco | std::vector\<std::vector\<double\>\> | Deconvolved signals |  
| StampTimeDeco | std::vector\<double\> | Start time of each digitized waveform (in $\mu s$) |  
| OpChDeco | std::vector\<int\> | Associated PD ID | 

Same scheme followed, but storing the deconvolved signals instead of the raw waveforms.

#### Pulse finder (a.k.a. OpHits)

It dumps all the reconstructed OpHits in the event.

| Branch name                 |    Type                                     |  Description |             
|-----------------------------|---------------------------------------------|--------------|  
| nophits | int | Total number of reconstructed OpHits |  
| ophit_opch | std::vector<int\> | Optical channel corresponding to the reconstructed OpHit |  
| ophit_peakT | std::vector<double\> | Waveform bin in which the OpHit gets the maximum value (in $\mu s$)  |  
| ophit_startT | std::vector<double\>| Start of the OpHit (in $\mu s$)  |  
| ophit_peakT | std::vector<double\> | OpHit rise time, relative to the StartTime (in $\mu s$)  |  
| ophit_width | std::vector<double\> | With of the OpHit (in $\mu s$)  |  
| ophit_amplitude | std::vector<double\> | Amplitude of the OpHit (in ADC units)  |  
| ophit_area | std::vector<double\> | Area of the OpHit (in $\mu s$ x ADC units)  |  
| ophit_pe | std::vector<double\> | Reconstructed number of PE  |  


#### Clustering among different PDs (a.k.a. OpFlash)

It dumps all the reconstructed OpFlash in the event and the associated OpHits associated to each OpFlash. Each entry in the following vectors correspond to one OpFlash.

| Branch name                 |    Type                                     |  Description |             
|-----------------------------|---------------------------------------------|--------------|  
| nopflash | int | Total number of reconstructed OpFlash objects |  
| flash_time| std::vector<double\> | t0 of the reconstructed OpFlashes |  
| flash_total_pe | std::vector<double\> | Integrated (all optical channels) number of photoelectrons in each OpFlash |  
| flash_pe_v |  std::vector\<std::vector<double\>\> | Vector containing the reconstructed number of photoelectron in each optical channel for each OpFlash |  
| flash_x, flash_y, flash_z | std::vector<double\> | X, Y, Z position of the reconstructed OpFlash |  
| flash_ophit_* | std::vector\<std::vector<double\>\> | Save the attributes of the OpHits associated to each OpFlash  |  

  #### Legend
  - PD: Photon Detector
  - PE: Photoelectron
  - BNB: Booster Neutrino Beam

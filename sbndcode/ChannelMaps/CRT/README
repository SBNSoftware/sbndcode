# SBND CRT Channel Maps

The SBND CRT channel maps are stored in text files that contain 3 numbers:

GDML ID | FEB MAC5 | Inversion

These numbers are:
- _GDML ID_ the number used to describe the active scintillator area in the simulated geometry
- _FEB MAC5_ the hardware address of the readout board (reported in every data packet)
- _Inversion_ a boolean describing whether the 32 channels are ordered in the same direction reality and simulation (0) or opposite (1)


During commissioning a number of boards have been replaced for a variety of reasons. This means the channel map has changed in a time/run dependent fashion. The channel maps here each represent a distinct period in this history.

- v1 before run 13511 (28th May 2024)
- v2 before run 15947 (9th August 2024)
- v3 before run 16955 (26th September 2024)
- v4 before run 17192 (11th October 2024)
- v5 Current - run 17192 onwards

Useful links to understand more detail of the SBND CRT Hardware-Simulation mapping.

- [CRT Master Database](https://docs.google.com/spreadsheets/d/1ReXP3Q2DuU-mO_vaQXZWj1h8NVdt6r13sBigYDvihXo/edit?usp=sharing) - a database containing the various numbering schemes and timing delays associated with the CRT hardware.
- [SBN DocDB 36401](https://sbn-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=36401) - maps of the original hardware installations - before any FEB swaps.
- [SBN DocDB 34844](https://sbn-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=34844) - visualisations showing the locations of each of the simulated modules by their GDML ID.


All of these channel maps use the number scheme associated with `sbnd_v02_03.gdml`. If using `sbnd_v02_02.gdml` or earlier these will not work.
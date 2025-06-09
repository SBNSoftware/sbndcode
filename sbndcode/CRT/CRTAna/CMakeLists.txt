art_make(
         MODULE_LIBRARIES
         ROOT::Tree
         ROOT::Graf3d
         artdaq_core::artdaq-core_Data
         larsim::Utils
         sbndaq_artdaq_core::sbndaq-artdaq-core_Obj_SBND
         sbnobj::SBND_CRT
         sbnobj::SBND_Timing
         sbndcode_GeoWrappers
         sbndcode_CRT_CRTBackTracker
)

install_fhicl()

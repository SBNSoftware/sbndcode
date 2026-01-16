art_make(
         MODULE_LIBRARIES
         ROOT::Tree
         ROOT::Graf3d
         art_root_io::TFileService_service
         larsim::Utils
         lardata::DetectorClocksService
         sbndaq_artdaq_core::sbndaq-artdaq-core_Obj_SBND
         sbnobj::SBND_CRT
         sbnobj::SBND_Timing
         sbndcode_CRT_CRTBackTracker
)

install_fhicl()

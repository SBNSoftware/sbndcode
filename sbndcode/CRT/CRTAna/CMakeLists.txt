art_make(
         MODULE_LIBRARIES
         ROOT::Tree
         art_root_io::TFileService_service
         larsim::Utils
         lardata::DetectorClocksService
         sbnobj::SBND_CRT
         sbnobj::SBND_Timing
         sbndcode_CRT_CRTBackTracker
)

install_fhicl()

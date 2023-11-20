art_make_library(
                 SOURCE SecondShowerFinderAlg.cxx
                 LIBRARIES
                 ROOT::Graf3d
                 art::Framework_Principal
                 lardata::DetectorPropertiesService
                 larcore::Geometry_Geometry_service
                 lardataobj::RecoBase
)

install_fhicl()

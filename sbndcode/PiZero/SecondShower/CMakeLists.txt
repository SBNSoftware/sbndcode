art_make_library(
                 SOURCE SecondShowerFinderAlg.cxx
                 LIBRARIES
                 art::Framework_Principal
                 lardataobj::RecoBase
)

install_fhicl()

cet_build_plugin(DriftEstimatorPMTRatio art::tool SOURCE DriftEstimatorPMTRatio_tool.cc LIBRARIES sbndcode_OpDetReco_OpFlash_FlashFinder )

cet_build_plugin(FlashGeoBarycenter art::tool SOURCE FlashGeoBarycenter_tool.cc LIBRARIES sbndcode_OpDetReco_OpFlash_FlashFinder )
cet_build_plugin(FlashGeoThreshold art::tool SOURCE FlashGeoThreshold_tool.cc LIBRARIES sbndcode_OpDetReco_OpFlash_FlashFinder )

cet_build_plugin(FlashT0FirstHit art::tool SOURCE FlashT0FirstHit_tool.cc LIBRARIES sbndcode_OpDetReco_OpFlash_FlashFinder )
cet_build_plugin(FlashT0SelectedChannels art::tool SOURCE FlashT0SelectedChannels_tool.cc LIBRARIES sbndcode_OpDetReco_OpFlash_FlashFinder )


install_headers()
install_fhicl()
install_source()

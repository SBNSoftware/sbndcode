file(GLOB json_files *.json)
install_fw(LIST ${json_files} SUBDIRNAME NuIntNNOpReco/3-inference-larsoft-module/tf/v0127_trained_w_160k_resnet18_fall_prod)
install(DIRECTORY saved_model DESTINATION scripts/NuIntNNOpReco/3-inference-larsoft-module/tf/v0127_trained_w_160k_resnet18_fall_prod)

file(GLOB json_files *.json)
install_fw(LIST ${json_files} SUBDIRNAME NuIntNNOpReco/3-inference-larsoft-module/tf/v0901_trained_w_165k_resnet18)
install(DIRECTORY saved_model DESTINATION scripts/NuIntNNOpReco/3-inference-larsoft-module/tf/v0901_trained_w_165k_resnet18)
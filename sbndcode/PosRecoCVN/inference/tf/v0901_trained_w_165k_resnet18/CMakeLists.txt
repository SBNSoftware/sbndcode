file(GLOB json_files *.json)
install_fw(LIST ${json_files} SUBDIRNAME PosRecoCVN/inference/tf/v0901_trained_w_165k_resnet18)
install(DIRECTORY saved_model DESTINATION scripts/PosRecoCVN/inference/tf/v0901_trained_w_165k_resnet18)
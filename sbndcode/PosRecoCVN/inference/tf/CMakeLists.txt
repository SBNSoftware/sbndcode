cet_make_library(SOURCE
  tf_graph.cc
  LIBRARIES PRIVATE
  TensorFlow::cc
  TensorFlow::framework
)

add_subdirectory(v0901_trained_w_165k_resnet18)
add_subdirectory(v0127_trained_w_160k_resnet18_fall_prod)
install_headers()
install_source()

cet_make_library(SOURCE
  tf_graph.cc
  LIBRARIES PRIVATE
  TensorFlow::cc
  TensorFlow::framework
)

install_headers()
install_source()

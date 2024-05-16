file(GLOB lib_srcs *.cxx)

cet_make_exec(NAME flatten_caf
              SOURCE flatten_caf.cc
              LIBRARIES sbnanaobj::StandardRecordFlat
                        ROOT::Core)

target_include_directories( flatten_caf PRIVATE $ENV{SRPROXY_INC} )

install_headers()
install_source()
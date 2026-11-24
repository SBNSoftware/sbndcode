# EventGeneratorListAssembler.xml holding only the "Default+SingleKaon" list
# (the AR25_20i Default list and DIS-CC-SINGLEK). Installed in a subdirectory so that it shadows
# GENIE's copy only for jobs that put this directory on GXMLPATH.
install_fw(LIST EventGeneratorListAssembler.xml SUBDIRNAME GENIE/SingleKaon)

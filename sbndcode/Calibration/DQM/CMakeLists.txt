cet_build_plugin(FilterChoppy art::module 
        lardataobj::RecoBase
        ROOT::Core
        ROOT::Tree
        art::Framework_Core
        art::Framework_Principal
        art::Utilities
        art_root_io::tfile_support
        art_root_io::TFileService_service
        canvas::canvas
        fhiclcpp::fhiclcpp
        cetlib::cetlib
        cetlib_except::cetlib_except
        messagefacility::MF_MessageLogger
        lardataobj::RawData
        larcore::Geometry_Geometry_service
        fhiclcpp::fhiclcpp
        cetlib::cetlib
        CLHEP::Random
        ROOT::Geom
        ROOT::XMLIO
        ROOT::Gdml
        ROOT::FFTW
)

install_headers()
install_source()
install_fhicl()

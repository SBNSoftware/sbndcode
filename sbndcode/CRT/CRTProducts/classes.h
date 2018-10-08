#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"
#include "sbndcode/CRT/CRTProducts/CRTData.hh"
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include <vector>
#include <utility>

template class art::Wrapper<sbnd::crt::CRTData>;
template class std::vector<sbnd::crt::CRTData>;
template class art::Wrapper<std::vector<sbnd::crt::CRTData> >;

template class art::Wrapper<sbnd::crt::CRTHit>;
template class std::vector<sbnd::crt::CRTHit>;
template class art::Wrapper<std::vector<sbnd::crt::CRTHit> >;

template class std::vector< std::pair<int,float> >;
template class std::map< unsigned char, std::vector< std::pair<int,float> > >;


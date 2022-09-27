#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/Source.h"
#include "ImportCRTData.h"

namespace crt {
    typedef art::Source<ImportCRTData> ImportCRTDataSource;
}

DEFINE_ART_INPUT_SOURCE(crt::ImportCRTDataSource)


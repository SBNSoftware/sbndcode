/**
 * @file   sbndcode/Decoders/DecoderTools/Dumpers/classes.h
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 * @date   Feb 7, 2023
 * 
 * Enables dictionary definitions for:
 * 
 * * `sbn::OpDetWaveformMeta`
 *   (and its associations with `raw::OpDetWaveform`)
 * 
 * See also `sbndcode/Decoders/DecoderTools/Dumpers/classes_def.xml`.
 */

// SBND libraries
#include "sbndcode/Decoders/DecoderTools/Dumpers/OpDetWaveformMeta.h"

// LArSoft libraries
#include "lardataobj/RawData/OpDetWaveform.h"

// framework libraries
#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"

// C++ libraries
#include <ostream>

namespace {

  sbn::OpDetWaveformMeta opdetwaveformmeta;

} // local namespace

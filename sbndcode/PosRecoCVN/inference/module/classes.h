#include "sbndcode/PosRecoCVN/inference/module/PixelMapVars.h"
#include "canvas/Persistency/Common/Wrapper.h"

// This file is needed to tell ROOT about our custom data structures
// so they can be written to and read from ROOT files properly

// Explicit instantiation of wrapper for our custom types
template class art::Wrapper<PixelMapVars>;
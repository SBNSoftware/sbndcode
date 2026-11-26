#include "MichelBlips.hh"



int MichelBlips::NBlips()
{
    int nBlips = 0;

    float muLength = std::hypot(fbragg.X - frecoVtxX, fbragg.Y - frecoVtxY, fbragg.Z - frecoVtxZ);

    float maxDist = std::min(muLength, fmaxBlipDist);
  
    for (auto & blip : fmuBlips)
    {   
        if (!blip.IsValid) continue;
        if (blip.TpcID != fbragg.TpcID) continue;
        if(blip.MaxWireSpan > fmaxWireSpan) continue;

        double vtxDist = std::hypot(blip.X - frecoVtxX, blip.Y - frecoVtxY, blip.Z - frecoVtxZ);
        double braggDist = std::hypot(blip.X - fbragg.X, blip.Y - fbragg.Y, blip.Z - fbragg.Z);

        
        if (vtxDist < braggDist) continue; // Only consider blips closer to bragg than to interaction vertex
        if (braggDist > maxDist) continue; // Reject blips too far away

        nBlips++;
    }

    return nBlips;

}

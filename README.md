# sbndcode

## How to run and what changed:

Run like:

`lar -c prodgenie_nu_spill_tpc_lowEen_sbnd.fcl -n 10`

**Workflow:**

```bash
#Made new directory and set up LArSoft newDev
mkdir && cd /exp/sbnd/app/users/lucykots/dk2nu_to_marley/
                                                                                                                                                  
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh                                                                                                       
mrb newDev -v v10_15_00 -q prof:e26 #Fall production tag                                                                                                                 
source localProducts_larsoft_v09_91_02_01_prof_e26/setup                                                                                                                 
cd srcs                                                                                                                                                                  
mrb g sbndcode                                                                                                                                                           
cd $MRB_SOURCE/sbndcode                                                                                                                                                  
git checkout feature/lucykots_dk2nutomarley                                                                                                                                       
mrbsetenv                                                                                                                                                                
mrb i -j4                                                                                                                                                                
mrbslp  
```

```bash
#Added this to genie_sbnd.fcl
#At L298
#/exp/sbnd/app/users/lucykots/dk2nu_to_marley/srcs/sbndcode/sbndcode/LArSoftConfigurations/gen/genie_sbnd.fcl

#                                                                                                                                                                        
# Booster Neutrino Beam, FHC mode, low energy uDAR & piDAR enabled, G4BNB v1.1.1, A                                                                                      
#                                                                                                                                                                        
sbnd_flux_g4bnb_1_1_1_lowE_enabled: {                                                                                                                                    
  FluxType:         "dk2nu"                                                                                                                                              
  DetectorLocation: "SBND"                                                                                                                                               
  FluxCopyMethod:   "DIRECT"                                                                                                                                             
  FluxSearchPaths: "/pnfs/sbnd/persistent/users/lucykots/G4BNB_muonDAR/1000n_5000p_pidar_udar/production_muDAR_piDAR_BooNE_50m_I174000A/all_NuBeam_files_fhc/"           
  FluxFiles: [ "NuBeam_production_muDAR_piDAR_BooNE_50m_I174000A_*.dk2nu.root" ]                                                                                         
}

#And added this to L
#
# Booster Neutrino Beam, low energy uDAR & piDAR enabled
#
sbnd_flux_bnb_lowE_en: @local::sbnd_flux_g4bnb_1_1_1_lowE_enabled
```

```cpp
#Copied Marco's fcl prodgenie_nu_spill_tpc_beamdump_sbnd.fcl and made my own here:
/exp/sbnd/app/users/lucykots/dk2nu_to_marley/srcs/sbndcode/sbndcode/JobConfigurations/standard/gen/genie/prodgenie_nu_spill_tpc_lowEen_sbnd.fcl

#Changed one line:
# Flux in from dk2nu files with low energy uDAR and piDAR enabled                                                                                                        
#include "set_flux_lowE_enabled.fcl"
```

```cpp
#Also copied Marco's fcl  that I included above and made my own here:
/exp/sbnd/app/users/lucykots/dk2nu_to_marley/srcs/sbndcode/sbndcode/LArSoftConfigurations/gen/set_flux_lowE_enabled.fcl

#Changed this line:
physics.producers.generator: {                                                                                                                                           
     @table::physics.producers.generator                                                                                                                                 
     @table::sbnd_flux_bnb_lowE_en                                                                                                                                       
}

```

I rebuilt and ran fcl command:
